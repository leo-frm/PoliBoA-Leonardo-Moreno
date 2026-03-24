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

    # Se o ponto final saiu da região estendida → divergência
    if !is_inside_range(integrator.u, bp.region.extended_range)
        return -1
    end

    # Se está dentro da estendida mas fora da região principal → fora da bacia (Na verdade, integra por um novo periodo de ciclos para tomar a decisão se ainda esta fora da bacia ou se realmente divergiu)
    if !is_inside_range(integrator.u, bp.region.range)
        return -2
    end
    
    # Caso geral: determina a célula destino
    target_cell = get_cell_number(integrator.u, bp.region)
    return target_cell
end

## Atencao paralelizavel:
# function populate_basins(bp :: BasinProblem) :: BasinResult #Coordena a paralelização da análise de bacias
#     # Total cells in the basin
#     total_cells = prod(bp.region.elements)

#     p = Progress(total_cells) #Inicializa a barra de progresso
#     # Single thread
#     if bp.threads == 1
#         # Map the cells for a single thread - no split necessary
#         results = map_cells(bp, 1:total_cells, p)
#         # Return the results - no join necessary
#         return results
#     # Multiple threads
#     else
#         # Divide the ranges for the especified number of threads
#         cell_ranges = Iterators.partition(1:total_cells, ceil(Int64, total_cells / bp.threads))
#         # Creates the tasks to be used (based on the number of threads)
#         tasks = map(cr -> Threads.@spawn(map_cells(bp, cr, p)), cell_ranges)
#         # Fetch the results of each thread
#         results = map(t -> fetch(t), tasks)
#         # Join the results basins up to one result
#         return foldl(join_results!, results)
#     end
# end