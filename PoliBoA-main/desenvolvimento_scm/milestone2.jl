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

# celula 1 -> 32 : 1 esta 0? esta em 0, 1 recebe 1, 32 esta em 0? recebe 1 também.
# capturar retornos na própria trajetória
# definir os pontos do atrator

function find_attractors_from_scm(scm::SimpleCellMap, bp::BasinProblem) :: BasinResult
    grid = zeros(Int64, Tuple(bp.region.elements))
    attractors = Attractor[]
    next_attractor_id = 1
    total_cells = prod(bp.region.elements)
    
    # Para cada célula ainda não classificada, segue o grafo do mapa de células
    for start_cell in 1:total_cells
        if grid[start_cell] != 0
            continue
        end

        trajectory = Int64[]
        current = start_cell
        max_steps = bp.maximum_cycles
        # trajectory guarda o caminho percorrido; current é a célula atual seguida
        
        # 
        for step in 1:max_steps
            push!(trajectory, current)
            # Cada passo adiciona a célula atual à trajetória; isso constrói o path

            if grid[current] != 0
                attractor_id = grid[current]

                foreach(c -> grid[c] = attractor_id, trajectory)
                break
            end
            # Se a célula atual já tem destino conhecido, toda a trajetória tem esse mesmo destino

        end
        next_cell = scm.target[current]

        if next_cell < 0
            foreach(c -> grid[c] = -1, trajectory)
            break
        end
        # Caso o mapeamento retorne -1, marca-se a trajetória inteira como divergente

        cycle_start = findfirst(==(next_cell), trajectory)
        if !isnothing(cycle_start)  
            attractor_cells = trajectory[cycle_start:end]
            attractor_type = length(attractor_cells) < 20 ? regular : quasi_periodic
            # Se next_cell já apareceu na trajetória, encontra-se um ciclo (atrator)

            push!(attractors, make_attractor_from_cells(
                attractor_cells,
                number = next_attractor_id,
                type = attractor_type,
                region = bp.region
            ))
            # Cria-se o objeto Attractor correspondente ao ciclo
        
            foreach(c -> grid[c] = next_attractor_id, trajectory)
            next_attractor_id += 1
            break
        end
        # Marcam-se todas as células visitadas como pertencendo ao novo atrator encontrado

        current = next_cell
    end

    return BasinResult(grid, attractors, bp.region)
end

#function trace_trajectory_in_map(u, region::mapSimpleCell) :: Int64 -> Embutida em find_attractors_from_scm por enquanto

#function detect_cycle_in_trajectory(trajectory::Vector{Int64}) :: Union{Int64, Nothing} -> Embutida em find_attractors_from_scm por enquanto

#function classify_attractor_type(u, region::mapSimpleCell) :: Int64 -> Embutida em find_attractors_from_scm por enquanto

#function mark_basin_cells(trajectory, attractor_id, grid) -> Embutida em find_attractors_from_scm por enquanto