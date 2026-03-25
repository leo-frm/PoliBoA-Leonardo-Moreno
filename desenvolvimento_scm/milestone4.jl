function process_cell_range(
    cell_range,
    bp::BasinProblem
) :: SimpleCellMap

    scm_local = SimpleCellMap(bp.region)

    ode_problem = ODEProblem(
        bp.f,
        zeros(Float64, length(bp.region.elements)),
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

    for cell in cell_range
        scm_local.target[cell] = compute_single_cell_mapping(cell, integrator, bp)
        scm_local.computed_cells += 1
    end

    return scm_local
end

"""
    join_scmaps!(scm_a, scm_b) -> scm_a

Une dois SimpleCellMaps copiando os valores não-zero de scm_b para scm_a.
Lógica:
- scm_a.target[i] == 0 → célula ainda não processada em scm_a
- scm_b.target[i] != 0 → célula foi processada em scm_b
"""
function join_scmaps!(scm_a::SimpleCellMap, scm_b::SimpleCellMap) :: SimpleCellMap
    for cell in eachindex(scm_a.target)
        if scm_a.target[cell] == 0 && scm_b.target[cell] != 0
            scm_a.target[cell] = scm_b.target[cell]
            scm_a.computed_cells += 1
        end
    end
    return scm_a
end

function build_scmap_parallel(bp::BasinProblem)
    total_cells = prod(bp.region.elements)

    if bp.threads == 1
        ode_problem = ODEProblem(
            bp.f,
            zeros(Float64, length(bp.region.elements)),
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

        scm = SimpleCellMap(bp.region)

        for cell in 1:total_cells
            scm.target[cell] = compute_single_cell_mapping(cell, integrator, bp)
            scm.computed_cells += 1
        end

        return scm

    else
        cell_ranges = Iterators.partition(
            1:total_cells,
            ceil(Int64, total_cells / bp.threads)
        )

        tasks = map(cr -> Threads.@spawn(process_cell_range(cr, bp)), cell_ranges)
        results = map(fetch, tasks)

        return foldl(join_scmaps!, results)
    end
end

function compute_single_cell_mapping(
    cell_id::Int64,
    integrator,
    bp::BasinProblem
) :: Int64
    u = zeros(Float64, length(bp.region.elements))
    store_cell_center!(u, cell_id, bp.region)

    set_integrator!(integrator; u, t=0.0)
    step!(integrator, bp.period * bp.transient_cycles, true)

    set_integrator!(integrator; u=integrator.u, t=0.0)
    step!(integrator, bp.period, true)

    adjust_cyclic(integrator.u, bp.region.range, bp.region.is_cyclic)

    # Saiu da região estendida → divergente
    if !is_inside_range(integrator.u, bp.region.extended_range)
        return -1
    end

    # Está na região estendida mas fora do grid principal →
    # reintegra até maximum_extended_cycles tentando retornar (idêntico à M1)
    if !is_inside_range(integrator.u, bp.region.range)
        for _ in 1:bp.maximum_extended_cycles
            set_integrator!(integrator; u=integrator.u, t=0.0)
            step!(integrator, bp.period, true)
            adjust_cyclic(integrator.u, bp.region.range, bp.region.is_cyclic)

            if !is_inside_range(integrator.u, bp.region.extended_range)
                return -1
            end

            if is_inside_range(integrator.u, bp.region.range)
                return get_cell_number(integrator.u, bp.region)
            end
        end
        return -1
    end

    return get_cell_number(integrator.u, bp.region)
end