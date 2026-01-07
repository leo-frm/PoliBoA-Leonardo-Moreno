using Accessors
using DifferentialEquations
using Plots
using PlotThemes
using Plots.PlotMeasures
using DelimitedFiles
using LaTeXStrings
using Printf
using Colors
using Distributed
using ProgressMeter

function process_cell_range(
    cell_range,
    bp::BasinProblem
) :: SimpleCellMap
    # Cada thread cria SimpleCellMap
    scm_local = SimpleCellMap(bp.region)
    
    # Cada thread cria integrador
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
    
    # Processa apenas as células do range desta thread
    for cell in cell_range
        scm_local.target[cell] = compute_single_cell_mapping(cell, integrator, bp)
        scm_local.computed_cells += 1
    end
    
    # Retorna mapa particular da thread
    return scm_local
end

function join_scmaps!(scm_a::SimpleCellMap, scm_b::SimpleCellMap) :: SimpleCellMap
    """
    Une dois SimpleCellMaps copiando os valores não-zero de scm_b para scm_a.
    
    Lógica:
    - scm_a.target[i] == 0 → célula ainda não processada em scm_a
    - scm_b.target[i] != 0 → célula foi processada em scm_b
    - Copia o valor de scm_b para scm_a
    """
    
    total_cells = length(scm_a.target)
    
    for cell in 1:total_cells
        # Se scm_a não processou esta célula, mas scm_b sim...
        if scm_a.target[cell] == 0 && scm_b.target[cell] != 0
            scm_a.target[cell] = scm_b.target[cell]
            scm_a.computed_cells += 1
        end
    end
    
    return scm_a
end

function build_scmap_parallel(bp :: BasinProblem)
    scm = SimpleCellMap(bp.region)
    total_cells = prod(bp.region.elements)

    if bp.threads == 1
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

    # Multiple threads
    else
        # Divide the ranges for the especified number of threads
        cell_ranges = Iterators.partition(
            1:total_cells, 
            ceil(Int64, total_cells / bp.threads)
        )
        # Creates the tasks to be used (based on the number of threads)
        tasks = map(cr -> Threads.@spawn(process_cell_range_ccm_style(cr, bp)), cell_ranges)
        # Fetch the results of each thread
        results = map(t -> fetch(t), tasks)
        # results = [scm₁, scm₂, scm₃, scm₄] 
        # Join the results basins up to one result
        scm = foldl(join_scmaps!, results) 
        return scm
    end
end