# Contém as funções algorítmicas que não existem em structures.jl/utils.jl:
#   - populate_basins
#   - map_cells
#   - search_for_attractor!
#   - join_results!
#
# Structs, funções utilitárias e dinâmicas NÃO estão aqui (seria uma repetição, já que vêm de
# structures.jl, utils.jl e do próprio tests_milestone5.jl).

function populate_basins(bp :: BasinProblem) :: BasinResult
    total_cells = prod(bp.region.elements)

    p = Progress(total_cells)
    if bp.threads == 1
        results = map_cells(bp, 1:total_cells, p)
        return results
    else
        cell_ranges = Iterators.partition(1:total_cells, ceil(Int64, total_cells / bp.threads))
        tasks = map(cr -> Threads.@spawn(map_cells(bp, cr, p)), cell_ranges)
        results = map(t -> fetch(t), tasks)
        return foldl(join_results!, results)
    end
end

function map_cells(bp :: BasinProblem, cell_range, p) :: BasinResult
    grid = zeros(Int64, Tuple(x for x in bp.region.elements))

    attractors = []
    attractors_found = 1

    T = zeros(Int64, bp.maximum_cycles + 1)
    A = zeros(Int64, 20)

    u = zeros(Float64, length(bp.region.elements))

    ode_problem = ODEProblem(
        bp.f,
        zero(u),
        (0.0, bp.period * bp.maximum_cycles),
        bp.params,
    )
    integrator = init(
        ode_problem;
        dense=false,
        save_everystep=false,
        save_start=false,
        maxiters=1e10,
    )

    for c in cell_range
        if (grid[c] != 0)
            next!(p)
            continue
        end

        store_cell_center!(u, c, bp.region)
        set_integrator!(integrator; u, t=0.0)
        
        attractor = search_for_attractor!(;
            bp,
            grid,
            integrator,
            next_attr_num=attractors_found,
            T_prealloc=T,
            A_prealloc=A,
            u=u
        )

        if !isnothing(attractor)
            attractors_found += 1
            push!(attractors, attractor)
        end

        next!(p)
    end

    BasinResult(grid, attractors, bp.region)
end

function search_for_attractor!(; bp :: BasinProblem, grid, integrator, next_attr_num::Int64, 
                                T_prealloc::Vector{Int64}, A_prealloc::Vector{Int64}, u :: Vector{Float64})

    T = T_prealloc
    A = A_prealloc
    u0 = u

    T_len = A_len = 0
    cycles = 0
    extended_cycles = 0
    
    T[T_len += 1] = get_cell_number(u0, bp.region)

    step!(integrator, bp.period * bp.transient_cycles, true)
    adjust_cyclic(u, bp.region.range, bp.region.is_cyclic)
    
    set_integrator!(integrator; u, t=0.0)

    attr_number = -1
    attr_type = none

    found_new_attractor = false

    while cycles < bp.maximum_cycles
        if !is_inside_range(integrator.u, bp.region.range)
            if (!is_inside_range(integrator.u, bp.region.extended_range)
                || cycles > bp.maximum_cycles)
                attr_number = -1
                break
            end
            extended_cycles += 1
        else
            extended_cycles = 0
            cell = get_cell_number(integrator.u, bp.region)

            if grid[cell] != 0
                attr_number = grid[cell]
                break
            end

            T[T_len += 1] = cell

            occurrences = count(c -> c == cell, @view T[1:(T_len - 1)])
            if occurrences == bp.attractor_tolerance
                if length(A) < A_len + 1
                    push!(A, cell)
                else
                    A[A_len+=1] = cell
                end
            elseif occurrences > bp.attractor_tolerance
                found_new_attractor = true
                attr_number = next_attr_num
                attr_type = regular
                break
            end
        end

        step!(integrator, bp.period, true)
        adjust_cyclic(u, bp.region.range, bp.region.is_cyclic)
        cycles += 1
    end

    if cycles == bp.maximum_cycles
        found_new_attractor = true
        attr_type = quasi_periodic
        attr_number = next_attr_num
        A = @view(T[ceil(Int64, T_len/2):T_len])
        A_len = ceil(Int64, T_len/2)
    end

    foreach(c -> grid[c] = attr_number, @view T[1:T_len])

    if found_new_attractor
        return make_attractor_from_cells(@view(A[1:A_len]), type=attr_type, number=attr_number, region=bp.region)
    else
        return nothing
    end
end

function join_results!(resultA::BasinResult, resultB::BasinResult)
    convert_dict = Dict{eltype(resultA.cells),eltype(resultA.cells)}(-1 => -1)

    next_attr_num = maximum(map(a -> a.number, resultA.attractors))

    for attr in resultB.attractors
        other_attr = findfirst(a -> Set(a.points) == Set(attr.points),
            resultA.attractors)

        if isnothing(other_attr)
            convert_dict[attr.number] = (next_attr_num += 1)
            new_attr = @set attr.number = next_attr_num
            push!(resultA.attractors, new_attr)
        else
            convert_dict[attr.number] = resultA.attractors[other_attr].number
        end
    end

    for i in 1:length(resultA.cells)
        if iszero(resultA.cells[i]) && haskey(convert_dict, resultB.cells[i])
            resultA.cells[i] = convert_dict[resultB.cells[i]]
        end
    end

    return resultA
end