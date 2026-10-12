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
        attractor_candidate_cells = Int64[]
        current = start_cell
        max_steps = bp.maximum_cycles
        found = false
        # trajectory guarda o caminho percorrido; current é a célula atual seguida
        # attractor_candidate_cells acumula células que atingiram attractor_tolerance visitas
        
        for step in 1:max_steps
            push!(trajectory, current)

            # Se a célula atual já tem destino conhecido, toda a trajetória herda esse destino
            if grid[current] != 0
                attractor_id = grid[current]
                for c in trajectory
                    grid[c] = attractor_id
                end
                found = true
                break
            end

            next_cell = scm.target[current]

            # Caso o mapeamento retorne negativo, marca-se a trajetória inteira como divergente
            if next_cell < 0
                foreach(c -> grid[c] = -1, trajectory)
                found = true
                break
            end

            # Conta quantas vezes next_cell já apareceu na trajetória
            occurrences = count(==(next_cell), trajectory)

            if occurrences == bp.attractor_tolerance
                # Célula visitada attractor_tolerance vezes → candidata a parte do atrator
                push!(attractor_candidate_cells, next_cell)
            elseif occurrences > bp.attractor_tolerance
                # Atrator confirmado: as células candidatas formam o atrator
                attractor_type = length(attractor_candidate_cells) < 20 ? regular : quasi_periodic

                push!(attractors, make_attractor_from_cells(
                    attractor_candidate_cells,
                    number = next_attractor_id,
                    type = attractor_type,
                    region = bp.region
                ))

                foreach(c -> grid[c] = next_attractor_id, trajectory)
                next_attractor_id += 1
                found = true
                break
            end

            current = next_cell
        end

        # Trajetória esgotou max_steps sem convergir
        if !found && length(trajectory) > 0 && grid[trajectory[1]] == 0
            half = ceil(Int64, length(trajectory) / 2)
            fallback_cells = trajectory[half:end]

            foreach(c -> grid[c] = next_attractor_id, trajectory)
            push!(attractors, make_attractor_from_cells(
                fallback_cells,
                number = next_attractor_id,
                type = quasi_periodic,
                region = bp.region
            ))
            next_attractor_id += 1
        end
    end

    # Análise de agrupamento (Hsu §8.2): une as soluções periódicas que representam o mesmo atrator
    attractors = merge_clustered_attractors!(grid, attractors, bp)

    return BasinResult(grid, attractors, bp.region)
end

# Análise de agrupamento das soluções periódicas, parte do método SCM (Hsu §8.2):
# as soluções periódicas encontradas tendem a se agrupar em torno de cada atrator, e a
# bacia do atrator é a união das bacias do agrupamento. Quando o ponto de um atrator fica
# muito perto da fronteira entre duas células, o centro da célula vizinha também
# pode ser mapeado em si mesmo, e o mapa passa a ter duas soluções periódicas para
# um único atrator; a bacia do atrator é a união das bacias delas. Para reconhecê-las,
# integra-se a trajetória contínua a partir de cada atrator, como no PoliBoA original:
# se ela termina na bacia de outro atrator, os dois são o mesmo e as bacias são unidas.
# Atratores distintos não são unidos, por mais próximos que estejam: a trajetória
# contínua de cada um permanece nele.
function merge_clustered_attractors!(grid, attractors::Vector{Attractor}, bp::BasinProblem) :: Vector{Attractor}
    # Atrator em cuja bacia termina a trajetória contínua de cada atrator
    target = collect(1:length(attractors))
    for attr in attractors
        cell = continuous_trajectory_cell(attr.points[1], bp)
        if cell > 0 && grid[cell] > 0
            target[attr.number] = grid[cell]
        end
    end

    representative = cluster_representatives(target)

    # Renumera os atratores que sobram em ordem crescente e reatribui as bacias
    new_number = Dict(old => k for (k, old) in enumerate(sort(unique(representative))))
    for c in eachindex(grid)
        if grid[c] > 0
            grid[c] = new_number[representative[grid[c]]]
        end
    end

    merged = Attractor[]
    for attr in attractors
        if representative[attr.number] == attr.number
            attr.number = new_number[attr.number]
            push!(merged, attr)
        end
    end

    return merged
end

#function trace_trajectory_in_map(u, region::mapSimpleCell) :: Int64 -> Embutida em find_attractors_from_scm por enquanto

#function detect_cycle_in_trajectory(trajectory::Vector{Int64}) :: Union{Int64, Nothing} -> Embutida em find_attractors_from_scm por enquanto

#function classify_attractor_type(u, region::mapSimpleCell) :: Int64 -> Embutida em find_attractors_from_scm por enquanto

#function mark_basin_cells(trajectory, attractor_id, grid) -> Embutida em find_attractors_from_scm por enquanto