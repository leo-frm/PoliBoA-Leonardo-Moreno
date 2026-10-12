# Constrói o mapa probabilístico integrando m amostras sorteadas em cada célula
function build_gcmap(gp :: GcmProblem; verbose :: Bool = true)
    bp = gp.base
    total_cells = prod(bp.region.elements)
    M = samples_per_cell(gp)

    gcm = GeneralizedCellMap(bp.region, M)
    integrator = make_gcm_integrator(bp)

    if verbose
        @printf(">>> GCM: grid %s, subgrid %s, m=%d (%s, %d integrações)\n",
            join(bp.region.elements, "x"), join(gp.subgrid, "x"), M, gp.sampling, total_cells * M)
    end

    step_report = max(1, div(total_cells, 10))

    # A ordem crescente é obrigatória: o formato CSR exige que as listas de
    # imagens sejam concatenadas na mesma ordem em que os offsets são registrados
    for cell in 1:total_cells
        images, counts = compute_cell_transition_probabilities(cell, integrator, gp)

        append!(gcm.images, images)
        append!(gcm.counts, counts)
        gcm.offsets[cell + 1] = Int32(length(gcm.images) + 1)
        gcm.computed[cell] = true

        if verbose && cell % step_report == 0
            @printf("    %3d%% (%d/%d células)\n",
                round(Int, 100 * cell / total_cells), cell, total_cells)
        end
    end

    return gcm
end

# Computa as transições de uma única célula: integra as m subcélulas sorteadas
# no seu interior e agrupa os destinos. Devolve as células-imagem distintas em
# ordem crescente e o número de amostras que chegou a cada uma, de modo que
# sum(counts) == m e as probabilidades counts/m somam 1 (Hsu, eq. 11.1.1).
# O destino -1 é a célula sumidouro, que absorve as amostras divergentes: tratá-la
# como imagem legítima é o que permite a uma célula ter destino misto, situação
# que o SCM não consegue representar.
function compute_cell_transition_probabilities(
    cell_id :: Int64,
    integrator,
    gp :: GcmProblem
) :: Tuple{Vector{Int32}, Vector{Int32}}

    bp = gp.base

    u = zeros(Float64, length(bp.region.elements))

    images = Int32[]
    counts = Int32[]

    for s in select_samples(gp, cell_id)
        store_cell_sample!(u, cell_id, bp.region, s, gp.subgrid)
        target = map_point(u, integrator, bp)

        # Busca linear: o número de imagens distintas é pequeno (tipicamente 1 a 3),
        # de modo que a varredura é mais rápida que um dicionário
        idx = findfirst(==(Int32(target)), images)
        if idx === nothing
            push!(images, Int32(target))
            push!(counts, Int32(1))
        else
            counts[idx] += Int32(1)
        end
    end

    # Ordena por número da célula-imagem, tornando o mapa independente da ordem
    # em que as amostras foram integradas
    perm = sortperm(images)

    return images[perm], counts[perm]
end