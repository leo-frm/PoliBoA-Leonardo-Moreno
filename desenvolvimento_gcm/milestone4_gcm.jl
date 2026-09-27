function process_cell_range_gcm(
    cell_range,
    gp::GcmProblem
) :: GeneralizedCellMap

    bp = gp.base
    gcm_local = GeneralizedCellMap(
        bp.region,
        samples_per_cell(gp);
        cell_range = first(cell_range):last(cell_range),
    )

    integrator = make_gcm_integrator(bp)

    # Mesma gravação do build_gcmap, mas com o índice local da faixa: a k-ésima
    # célula da faixa ocupa a posição k dos offsets deste mapa parcial
    for cell in cell_range
        images, counts = compute_cell_transition_probabilities(cell, integrator, gp)
        k = local_index(gcm_local, cell)

        append!(gcm_local.images, images)
        append!(gcm_local.counts, counts)
        gcm_local.offsets[k + 1] = Int32(length(gcm_local.images) + 1)
        gcm_local.computed[k] = true
    end

    return gcm_local
end

"""
    join_gcmaps!(gcm_a, gcm_b) -> gcm_a

Une dois GeneralizedCellMaps parciais concatenando gcm_b ao final de gcm_a.
Lógica:
- as faixas precisam ser contíguas e em ordem crescente (o CSR concatena as
  listas de imagens na mesma ordem das células);
- os offsets de gcm_b são deslocados pelo número de imagens já presentes em gcm_a.

Diferente do join_scmaps!, que preenchia posições vazias em qualquer ordem, aqui
a ordem de união importa. O foldl em build_gcmap_parallel preserva essa ordem.
"""
function join_gcmaps!(gcm_a::GeneralizedCellMap, gcm_b::GeneralizedCellMap) :: GeneralizedCellMap
    @assert last(gcm_a.cell_range) + 1 == first(gcm_b.cell_range) "As faixas precisam ser contíguas e em ordem crescente"

    shift = Int32(length(gcm_a.images))

    append!(gcm_a.offsets, gcm_b.offsets[2:end] .+ shift)
    append!(gcm_a.images, gcm_b.images)
    append!(gcm_a.counts, gcm_b.counts)
    append!(gcm_a.computed, gcm_b.computed)
    gcm_a.cell_range = first(gcm_a.cell_range):last(gcm_b.cell_range)

    return gcm_a
end

function build_gcmap_parallel(gp::GcmProblem)
    bp = gp.base
    total_cells = prod(bp.region.elements)

    if bp.threads == 1
        integrator = make_gcm_integrator(bp)

        gcm = GeneralizedCellMap(bp.region, samples_per_cell(gp))

        for cell in 1:total_cells
            images, counts = compute_cell_transition_probabilities(cell, integrator, gp)

            append!(gcm.images, images)
            append!(gcm.counts, counts)
            gcm.offsets[cell + 1] = Int32(length(gcm.images) + 1)
            gcm.computed[cell] = true
        end

        return gcm

    else
        cell_ranges = Iterators.partition(
            1:total_cells,
            ceil(Int64, total_cells / bp.threads)
        )

        tasks = map(cr -> Threads.@spawn(process_cell_range_gcm(cr, gp)), cell_ranges)
        results = map(fetch, tasks)

        return foldl(join_gcmaps!, results)
    end
end
