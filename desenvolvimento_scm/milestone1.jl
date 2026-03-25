# Constrói o mapa célula‑a‑célula integrando um ponto de cada célula
function build_scmap(bp :: BasinProblem)
    scm = SimpleCellMap(bp.region)
    total_cells = prod(bp.region.elements)

    # Define o problema de EDO para integração de cada centro de célula
    ode_problem = ODEProblem(
        bp.f,
        zeros(Float64, length(bp.region.elements)),
        (0.0, bp.period * bp.maximum_cycles),
        bp.params,
    )

    # Cria integrador
    integrator = init(
        ode_problem;
        dense=false,
        save_everystep=false,
        save_start=false,
        maxiters=1e10,
    )
    
    # Integra cada célula uma vez e registra seu destino
    for cell in 1:total_cells
        scm.target[cell] = compute_single_cell_mapping(cell, integrator, bp)
        scm.computed_cells += 1 # computed_cells? paralelizavel no futuro -> passa um range
    end

    return scm
end

# Computa o destino de uma única célula via integração da EDO
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

    # Se o ponto final saiu da região estendida:  divergência
    if !is_inside_range(integrator.u, bp.region.extended_range)
        return -1
    end

    # Se está dentro da estendida mas fora da região principal: 
    # reintegra até maximum_extended_cycles tentando retornar ao grid
    if !is_inside_range(integrator.u, bp.region.range)
        for ext_cycle in 1:bp.maximum_extended_cycles
            set_integrator!(integrator; u=integrator.u, t=0.0)
            step!(integrator, bp.period, true)
            adjust_cyclic(integrator.u, bp.region.range, bp.region.is_cyclic)

            # Saiu da região estendida:  divergente
            if !is_inside_range(integrator.u, bp.region.extended_range)
                return -1
            end

            # Voltou ao grid principal:  sucesso, retorna célula destino
            if is_inside_range(integrator.u, bp.region.range)
                return get_cell_number(integrator.u, bp.region)
            end
        end
        # Esgotou tentativas na região estendida: divergente
        return -1
    end
    
    # Caso geral: determina a célula destino
    target_cell = get_cell_number(integrator.u, bp.region)
    return target_cell
end