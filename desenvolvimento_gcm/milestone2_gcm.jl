using SparseArrays

# Detecção de atratores pelo mapa probabilístico, segundo o algoritmo de Hsu
# (Cap. 11). O mapa é tratado como uma cadeia de Markov: cada célula tem várias
# imagens simultâneas, com probabilidades. Atratores passam a ser grupos
# persistentes — conjuntos fechados de células que se comunicam entre si — e
# cada célula transiente é atribuída ao grupo em que tem maior probabilidade de
# ser absorvida. Retorna BasinResult, no mesmo formato do SCM.
function find_attractors_from_gcm(gcm::GeneralizedCellMap, gp::GcmProblem) :: BasinResult
    analysis = analyze_gcm(gcm, gp)
    return basin_result_from_analysis(analysis, gp)
end

# Executa as etapas do Cap. 11 em sequência e devolve a análise completa,
# com as informações que o BasinResult não comporta (domicílios, absorção)
function analyze_gcm(gcm::GeneralizedCellMap, gp::GcmProblem) :: GcmAnalysis
    N = prod(gp.base.region.elements)

    # §11.2: quem mapeia para cada célula (necessário para marchar para trás)
    preimages = build_preimage_array(gcm)

    # §11.3: análise preliminar pelo SCM dominante, cujas células periódicas
    # servem de semente para a busca dos grupos persistentes (Teorema 11.3.1)
    periodic_solutions = analyze_dominant_scm(gcm)

    # §11.4 e §11.5: células absorventes e demais grupos persistentes
    groups, group_of = identify_persistent_groups(gcm, preimages, periodic_solutions)

    # §8.2: grupos persistentes que representam o mesmo atrator
    groups, group_of = merge_clustered_groups(groups, group_of, gp.base)

    # §11.6: período de cada grupo persistente
    periods = [group_period(gcm, members) for members in groups]

    # §11.10: domicílios de cada célula, marchando para trás a partir de cada grupo
    reach = determine_domiciles(preimages, groups, N)

    # §11.9: probabilidade de absorção de cada célula em cada grupo
    absorption = assign_basin_probabilities(gcm, group_of, reach)

    return GcmAnalysis(
        length(periodic_solutions),
        groups,
        periods,
        group_of,
        vec(sum(reach, dims=2)),
        absorption,
    )
end

# Constrói o BasinResult a partir da análise. Cada célula recebe o grupo de
# maior probabilidade de absorção. O grupo 1 é o sumidouro (-1, divergência);
# os grupos 2..K viram os atratores 1..K-1.
function basin_result_from_analysis(analysis::GcmAnalysis, gp::GcmProblem) :: BasinResult
    bp = gp.base
    N = prod(bp.region.elements)
    grid = zeros(Int64, Tuple(bp.region.elements))

    for cell in 1:N
        g = argmax(view(analysis.absorption, cell, :))
        grid[cell] = g == 1 ? -1 : g - 1
    end

    attractors = Attractor[]
    for g in 2:length(analysis.groups)
        push!(attractors, make_attractor_from_cells(
            analysis.groups[g],
            number = g - 1,
            type = classify_markov_group_type(analysis.groups[g], analysis.periods[g]),
            region = bp.region
        ))
    end

    return BasinResult(grid, attractors, bp.region)
end

# Devolve a célula-imagem de maior probabilidade de transição. Define o SCM
# dominante de Hsu (§11.3): um SCM selecionado a partir do próprio GCM, e por
# isso sempre compatível com ele.
function dominant_image(gcm::GeneralizedCellMap, cell::Int64) :: Int64
    images = images_of(gcm, cell)
    counts = counts_of(gcm, cell)
    return Int64(images[argmax(counts)])
end

# Imagem dominante já convertida em nó do grafo; o sumidouro mapeia para si mesmo
function dominant_node(gcm::GeneralizedCellMap, z::Int64, N::Int64) :: Int64
    return z == N + 1 ? N + 1 : gcm_node(dominant_image(gcm, z), N)
end

# Monta o arranjo de pré-imagens (Hsu §11.2). Primeiro conta quantas vezes cada
# nó aparece como imagem, para dimensionar os offsets; depois preenche.
function build_preimage_array(gcm::GeneralizedCellMap) :: PreImageArray
    N = length(gcm.cell_range)

    n_pre = zeros(Int32, N + 1)
    for z in 1:N, img in images_of(gcm, z)
        n_pre[gcm_node(img, N)] += 1
    end

    offsets = Vector{Int32}(undef, N + 2)
    offsets[1] = 1
    for n in 1:(N + 1)
        offsets[n + 1] = offsets[n] + n_pre[n]
    end

    cells = Vector{Int32}(undef, offsets[end] - 1)
    next_free = offsets[1:(N + 1)]
    for z in 1:N, img in images_of(gcm, z)
        n = gcm_node(img, N)
        cells[next_free[n]] = z
        next_free[n] += 1
    end

    return PreImageArray(offsets, cells)
end

# Análise global do SCM dominante pelo algoritmo de Hsu §8.2: segue-se o ponteiro
# único de cada célula até encontrar uma célula já processada (a sequência herda
# a solução dela) ou uma célula da própria sequência (nova solução periódica).
# Por ser um mapa determinístico sobre um conjunto finito, todo caminho termina
# num ciclo exato. Devolve as células periódicas de cada solução (conjunto E_pc).
function analyze_dominant_scm(gcm::GeneralizedCellMap) :: Vector{Vector{Int64}}
    N = length(gcm.cell_range)
    scm_group = zeros(Int64, N + 1) # Gr(z): 0 virgem, -1 em processamento, k>0 solução periódica k
    periodic_solutions = Vector{Vector{Int64}}()
    sequence = Int64[]

    for start_cell in 1:(N + 1)
        if scm_group[start_cell] != 0
            continue
        end

        empty!(sequence)
        current = start_cell
        while scm_group[current] == 0
            scm_group[current] = -1
            push!(sequence, current)
            current = dominant_node(gcm, current, N)
        end

        if scm_group[current] == -1
            # Célula repetida na sequência atual: as células a partir dela formam o ciclo
            first_repeat = findfirst(==(current), sequence)
            push!(periodic_solutions, sequence[first_repeat:end])
            solution_id = length(periodic_solutions)
        else
            # Sequência chegou a uma célula já processada: herda a solução dela
            solution_id = scm_group[current]
        end

        for c in sequence
            scm_group[c] = solution_id
        end
    end

    return periodic_solutions
end

# Busca dos grupos persistentes (Hsu §11.4 e §11.5). Rótulos de G(z), como na
# eq. (11.5.1): 0 célula ainda não classificada, -2 transiente, g>0 membro do
# grupo g. Ao final, os transientes são normalizados para 0.
#
# Pelo Teorema 11.3.1, todo grupo persistente contém ao menos uma solução
# periódica do SCM dominante; por isso basta partir de uma célula de cada
# solução, em vez de varrer as N células.
#
# Para cada semente usam-se as propriedades (ii) a (iv) de §11.5.1:
#   - expande-se o conjunto de todas as células alcançáveis a partir da semente,
#     seguindo TODAS as imagens (não apenas a dominante);
#   - se esse conjunto toca uma célula já processada (transiente ou de outro
#     grupo), a semente é transiente;
#   - se não toca, o conjunto é fechado; ele é um grupo persistente se todas as
#     suas células conseguem voltar à semente (busca para trás pelas pré-imagens).
function identify_persistent_groups(
    gcm::GeneralizedCellMap,
    preimages::PreImageArray,
    periodic_solutions::Vector{Vector{Int64}}
)
    N = length(gcm.cell_range)
    group_of = zeros(Int64, N + 1)
    groups = Vector{Vector{Int64}}()

    # §11.4: o sumidouro é absorvente e constitui o primeiro grupo persistente
    push!(groups, [N + 1])
    group_of[N + 1] = 1

    in_closure = falses(N + 1)
    returns_to_seed = falses(N + 1)
    closure = Int64[]
    backward = Int64[]

    for solution in periodic_solutions
        seed = first(solution)
        if group_of[seed] != 0
            continue
        end

        # §11.4: célula absorvente, I(z) = 1 e C(z,1) = z
        if n_images(gcm, seed) == 1 && gcm_node(first(images_of(gcm, seed)), N) == seed
            push!(groups, [seed])
            group_of[seed] = length(groups)
            continue
        end

        # Expansão para frente: todas as células alcançáveis a partir da semente
        empty!(closure)
        push!(closure, seed)
        in_closure[seed] = true
        leads_to_processed = false
        head = 1

        while head <= length(closure) && !leads_to_processed
            z = closure[head]
            head += 1
            for img in images_of(gcm, z)
                n = gcm_node(img, N)
                if group_of[n] != 0
                    # Propriedades (iii) e (iv): alcança célula já classificada → semente transiente
                    leads_to_processed = true
                    break
                end
                if !in_closure[n]
                    in_closure[n] = true
                    push!(closure, n)
                end
            end
        end

        # Conjunto fechado: verifica se todas as células retornam à semente
        is_persistent = false
        if !leads_to_processed
            empty!(backward)
            push!(backward, seed)
            returns_to_seed[seed] = true
            head = 1

            while head <= length(backward)
                z = backward[head]
                head += 1
                for p in preimages_of(preimages, z)
                    pz = Int64(p)
                    if in_closure[pz] && !returns_to_seed[pz]
                        returns_to_seed[pz] = true
                        push!(backward, pz)
                    end
                end
            end

            is_persistent = length(backward) == length(closure)
            for z in backward
                returns_to_seed[z] = false
            end
        end

        if is_persistent
            # Propriedade (ii): o conjunto alcançável a partir de uma célula persistente é o próprio grupo
            push!(groups, sort(closure))
            for c in closure
                group_of[c] = length(groups)
            end
        else
            # Teorema 11.3.3: todas as células periódicas desta solução são transientes
            for c in solution
                group_of[c] = -2
            end
        end

        for c in closure
            in_closure[c] = false
        end
    end

    # Toda célula fora dos grupos é transiente
    group_of[group_of .< 0] .= 0

    return groups, group_of
end

# Agrupamento de grupos persistentes (Hsu §8.2), pelo mesmo critério do SCM
# (merge_clustered_attractors!): integra-se a trajetória contínua a partir de cada
# grupo e, se ela termina numa célula de outro grupo, os dois representam o mesmo
# atrator e são unidos. Com m = 1 o GCM é o próprio SCM e herda dele os grupos
# duplicados; com m maior as amostras costumam ligar as células e nada muda.
# O sumidouro (grupo 1) nunca é unido.
function merge_clustered_groups(groups::Vector{Vector{Int64}}, group_of::Vector{Int64}, bp::BasinProblem)
    u = zeros(Float64, length(bp.region.elements))

    # Grupo em que termina a trajetória contínua de cada grupo
    target = collect(1:length(groups))
    for g in 2:length(groups)
        store_cell_center!(u, first(groups[g]), bp.region)
        cell = continuous_trajectory_cell(u, bp)
        if cell > 0 && group_of[cell] > 1
            target[g] = group_of[cell]
        end
    end

    representative = cluster_representatives(target)

    # Grupos que sobram, em ordem crescente; o sumidouro continua sendo o grupo 1
    merged_groups = [sort(vcat([groups[g] for g in eachindex(groups) if representative[g] == r]...))
                     for r in sort(unique(representative))]

    merged_group_of = zeros(Int64, length(group_of))
    for (g, members) in enumerate(merged_groups), c in members
        merged_group_of[c] = g
    end

    return merged_groups, merged_group_of
end

# Período de um grupo persistente (Hsu §11.6). Um grupo de período d se divide em
# d subgrupos visitados ciclicamente. Numera-se cada célula pelo nível em que é
# alcançada a partir da primeira, e o período é o máximo divisor comum das
# diferenças (nível de origem + 1 - nível de destino) em todas as transições
# internas. Uma célula que tem a si mesma como imagem torna o grupo acíclico (d = 1).
function group_period(gcm::GeneralizedCellMap, members::Vector{Int64}) :: Int64
    N = length(gcm.cell_range)
    if members == [N + 1]
        return 1
    end

    member_set = Set(members)
    level = Dict{Int64, Int64}(first(members) => 0)
    queue = [first(members)]
    head = 1
    d = 0

    while head <= length(queue)
        z = queue[head]
        head += 1
        for img in images_of(gcm, z)
            n = gcm_node(img, N)
            if !(n in member_set)
                continue
            end
            if !haskey(level, n)
                level[n] = level[z] + 1
                push!(queue, n)
            else
                d = gcd(d, abs(level[z] + 1 - level[n]))
            end
        end
    end

    return d == 0 ? 1 : d
end

# Classificação do grupo persistente em termos do AttractorType do SCM.
# Pela §10.6, um grupo periódico ou um grupo acíclico pequeno corresponde a uma
# solução periódica estável; um grupo acíclico grande corresponde a um atrator
# não periódico. Mantém-se o mesmo corte de 20 células usado no SCM.
function classify_markov_group_type(members::Vector{Int64}, period::Int64) :: AttractorType
    return (period > 1 || length(members) < 20) ? regular : quasi_periodic
end

# Domicílios de cada célula (Hsu §11.10). Um domicílio de uma célula é também
# domicílio de todas as suas pré-imagens; então, partindo dos membros de cada
# grupo e marchando para trás pelo arranjo de pré-imagens, encontram-se todas
# as células que têm aquele grupo como domicílio. reach[z, g] = true se o grupo
# g é domicílio da célula z. Dm(z) é a soma da linha z.
function determine_domiciles(preimages::PreImageArray, groups::Vector{Vector{Int64}}, N::Int64) :: BitMatrix
    reach = falses(N + 1, length(groups))
    queue = Int64[]

    for (g, members) in enumerate(groups)
        empty!(queue)
        for m in members
            reach[m, g] = true
            push!(queue, m)
        end

        head = 1
        while head <= length(queue)
            z = queue[head]
            head += 1
            for p in preimages_of(preimages, z)
                pz = Int64(p)
                if !reach[pz, g]
                    reach[pz, g] = true
                    push!(queue, pz)
                end
            end
        end
    end

    return reach
end

# Probabilidades de absorção (Hsu §11.9 e §10.4). Para células persistentes e
# transientes de domicílio único a absorção é trivial: probabilidade 1 no único
# grupo possível. Só as células de domicílio múltiplo exigem cálculo: para elas
# resolve-se o sistema linear esparso (I - Q) X = B, em que Q contém as
# transições entre células de domicílio múltiplo e B as probabilidades de, num
# passo, chegar a células cuja absorção já é conhecida (Hsu eq. 11.9.1 e §11.11.2).
function assign_basin_probabilities(
    gcm::GeneralizedCellMap,
    group_of::Vector{Int64},
    reach::BitMatrix
) :: Matrix{Float64}

    N = length(gcm.cell_range)
    K = size(reach, 2)
    absorption = zeros(Float64, N + 1, K)
    domicile_count = vec(sum(reach, dims=2))

    # Absorções triviais
    multiple = Int64[]
    for z in 1:(N + 1)
        if group_of[z] > 0
            absorption[z, group_of[z]] = 1.0
        elseif domicile_count[z] == 1
            absorption[z, findfirst(reach[z, :])] = 1.0
        elseif domicile_count[z] > 1
            push!(multiple, z)
        end
    end

    if isempty(multiple)
        return absorption
    end

    # Índice local de cada célula de domicílio múltiplo no sistema linear
    M = length(multiple)
    index = zeros(Int64, N + 1)
    for (k, z) in enumerate(multiple)
        index[z] = k
    end

    rows = Int64[]
    cols = Int64[]
    vals = Float64[]
    B = zeros(Float64, M, K)

    for (k, z) in enumerate(multiple)
        # Diagonal da identidade
        push!(rows, k)
        push!(cols, k)
        push!(vals, 1.0)

        images = images_of(gcm, z)
        counts = counts_of(gcm, z)
        for i in eachindex(images)
            n = gcm_node(images[i], N)
            p = counts[i] / gcm.samples_per_cell
            if index[n] > 0
                # Transição para outra célula de domicílio múltiplo: entra em -Q
                push!(rows, k)
                push!(cols, index[n])
                push!(vals, -p)
            else
                # Transição para célula de absorção conhecida: entra em B
                for g in 1:K
                    B[k, g] += p * absorption[n, g]
                end
            end
        end
    end

    X = sparse(rows, cols, vals, M, M) \ B

    for (k, z) in enumerate(multiple)
        absorption[z, :] = X[k, :]
    end

    return absorption
end