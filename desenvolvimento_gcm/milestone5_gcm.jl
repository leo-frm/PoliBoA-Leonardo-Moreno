# Funções de validação para o método GCM:
#   - benchmark_method_gcm:        mede tempo médio, desvio padrão, bytes alocados e nº de integrações
#   - memory_scaling_gcm:          mede como a memória escala com o tamanho do grid
#   - thread_scaling_gcm:          mede como a performance escala com o número de threads
#   - compare_basins_by_domicile:  separa a concordância entre células de domicílio único e múltiplo
#   - refinement_reference_gcm:    calcula a referência fina do estudo de refinamento (SCM refinado)
#   - refinement_study_gcm:        compara GCMs grossos com a referência fina (custo vs precisão)
#
# Reaproveitadas sem alteração de desenvolvimento_scm/milestone5.jl:
#   BenchmarkResult, compare_basins, compare_benchmarks

using Printf
using Statistics

"""
    benchmark_method_gcm(run_fn; trials=5, warmup=1) -> BenchmarkResult

Executa `run_fn()` (uma função sem argumentos que retorna um `GeneralizedCellMap`)
`warmup` vezes para pré-compilar e `trials` vezes para medir.

O número de integrações do GCM é o número de células vezes o número de amostras
por célula (m), e não apenas o número de células como no SCM.
"""
function benchmark_method_gcm(run_fn; trials::Int=5, warmup::Int=1)
    # Warmup (pré-compilação)
    for _ in 1:warmup
        GC.gc()
        run_fn()
    end

    times = Float64[]
    last_bytes = Int64(0)
    last_integrations = Int64(0)
    last_grid_size = ""

    for i in 1:trials
        GC.gc()
        stats = @timed run_fn()

        push!(times, stats.time)
        last_bytes = stats.bytes

        # Extrai nº de integrações do resultado
        result = stats.value
        if result isa GeneralizedCellMap
            last_integrations = computed_cells(result) * result.samples_per_cell
            last_grid_size = join(result.region.elements, "x")
        elseif result isa SimpleCellMap
            last_integrations = result.computed_cells
            last_grid_size = join(result.range.elements, "x")
        elseif result isa BasinResult
            last_integrations = prod(result.region.elements)
            last_grid_size = join(result.region.elements, "x")
        end

        @printf("  trial %d/%d: %.3f s | %.1f MB\n", i, trials, stats.time, stats.bytes / 1e6)
    end

    return BenchmarkResult(
        mean(times),
        std(times),
        last_bytes,
        last_integrations,
        last_grid_size,
        trials,
    )
end

"""
    memory_scaling_gcm(build_fn, find_fn, gp_list; method_name)

Mede como a memória alocada escala com o tamanho do grid.
`gp_list` é um vetor de GcmProblem com resoluções crescentes.

Além das colunas do SCM, reporta Ī (número médio de imagens por célula) e a
memória ocupada pelo próprio mapa: se Ī se mantém estável quando o grid cresce,
o custo do mapa é O(N · Ī), linear em N.
"""
function memory_scaling_gcm(
    build_fn    :: Function,
    find_fn     :: Function,
    gp_list     :: Vector{GcmProblem};
    method_name :: String = "GCM"
)
    w = 90
    println("\n" * "=" ^ w)
    @printf("  ESCALABILIDADE DE MEMÓRIA  ·  %s\n", method_name)
    println("=" ^ w)
    @printf("  %-16s %10s %12s %12s %14s %8s %12s\n",
        "Grid", "Células", "Memória(MB)", "MB/kCell", "Integrações", "Ī", "Mapa(MB)")
    println("-" ^ w)

    for gp in gp_list
        GC.gc()
        total_cells = prod(gp.base.region.elements)
        grid_str    = join(gp.base.region.elements, "×")

        local gcm, result
        alloc_bytes = @allocated begin
            gcm    = build_fn(gp)
            result = find_fn(gcm, gp)
        end

        mb_total     = alloc_bytes / 1024^2
        mb_per_kcell = mb_total / (total_cells / 1_000.0)
        mean_images  = length(gcm.images) / computed_cells(gcm)

        @printf("  %-16s %10d %12.1f %12.3f %14d %8.3f %12.3f\n",
            grid_str, total_cells, mb_total, mb_per_kcell,
            computed_cells(gcm) * gcm.samples_per_cell, mean_images,
            memory_footprint(gcm) / 1024^2)
        GC.gc()
    end

    println("=" ^ w * "\n")
end

"""
    thread_scaling_gcm(build_fn, find_fn, make_gp_fn; thread_counts, trials)

Mede o tempo de execução do GCM para diferentes quantidades de threads
e calcula o speedup relativo à execução single-thread.

`make_gp_fn(threads)` deve retornar um `GcmProblem` configurado para
o número de threads dado.

Retorna vetor de NamedTuples com (threads, mean_time, speedup).
"""
function thread_scaling_gcm(
    build_fn    :: Function,
    find_fn     :: Function,
    make_gp_fn;
    thread_counts = [1, 2, 4],
    trials::Int = 3
)
    results = NamedTuple{(:threads, :mean_time, :speedup), Tuple{Int, Float64, Float64}}[]
    base_time = 0.0

    println("=" ^ 60)
    println("Análise de Escalabilidade por Threads (GCM)")
    println("=" ^ 60)
    @printf("  %-10s %12s %10s\n", "Threads", "Tempo médio", "Speedup")
    println("-" ^ 35)

    for nt in thread_counts
        gp = make_gp_fn(nt)
        bench = benchmark_method_gcm(trials=trials, warmup=1) do
            gcm = build_fn(gp)
            find_fn(gcm, gp)
            return gcm
        end

        if nt == first(thread_counts)
            base_time = bench.mean_time_s
        end

        speedup = base_time / bench.mean_time_s
        push!(results, (threads=nt, mean_time=bench.mean_time_s, speedup=speedup))
        @printf("  %-10d %10.4f s %10.2fx\n", nt, bench.mean_time_s, speedup)
    end

    println("=" ^ 60)
    return results
end

"""
    basin_id_mapping(result_a, result_b) -> Dict{Int64, Int64}

Mapeamento entre os identificadores de atratores de dois resultados, pela mesma
regra de compare_basins: -1 casa com -1, e os demais pares são atribuídos de
forma gulosa bijetiva por sobreposição decrescente.
"""
function basin_id_mapping(result_a::BasinResult, result_b::BasinResult) :: Dict{Int64, Int64}
    cells_a = result_a.cells
    cells_b = result_b.cells

    pair_counts = Dict{Tuple{Int64, Int64}, Int64}()
    for i in eachindex(cells_a)
        key = (cells_a[i], cells_b[i])
        pair_counts[key] = get(pair_counts, key, 0) + 1
    end

    mapping = Dict{Int64, Int64}(-1 => -1)

    candidates = [(a_id, b_id, count) for ((a_id, b_id), count) in pair_counts
                  if a_id != -1 && b_id != -1]
    sort!(candidates, by = x -> -x[3])

    used_a = Set{Int64}()
    used_b = Set{Int64}()
    for (a_id, b_id, _) in candidates
        if a_id ∉ used_a && b_id ∉ used_b
            mapping[a_id] = b_id
            push!(used_a, a_id)
            push!(used_b, b_id)
        end
    end

    return mapping
end

"""
    compare_basins_by_domicile(result_gcm, result_ref, analysis; label_ref)

Compara as bacias do GCM com as de um método de referência (SCM ou PoliBoA
original) separando as células pelo número de domicílios (Hsu §10.4.4, §11.10):

- domicílio único: células que só podem terminar num atrator (interior das bacias);
- domicílio múltiplo: células que podem terminar em mais de um (fronteira entre bacias).

O mapeamento de atratores é o mesmo de compare_basins, calculado sobre o grid
inteiro. A hipótese testada (Hsu §10.6) é que as discordâncias entre métodos se
concentram nas células de domicílio múltiplo.

Retorna NamedTuple com (overall, single, multiple, frac_multiple).
"""
function compare_basins_by_domicile(
    result_gcm :: BasinResult,
    result_ref :: BasinResult,
    analysis   :: GcmAnalysis;
    label_ref  :: String = "Referência"
)
    cells_gcm = vec(result_gcm.cells)
    cells_ref = vec(result_ref.cells)
    N = length(cells_gcm)

    mapping = basin_id_mapping(result_gcm, result_ref)
    agrees = [get(mapping, cells_gcm[i], nothing) == cells_ref[i] for i in 1:N]

    is_multiple = view(analysis.domicile_count, 1:N) .> 1
    n_multiple = count(is_multiple)
    n_single = N - n_multiple

    overall = count(agrees) / N
    single = n_single == 0 ? 1.0 : count(agrees .& .!is_multiple) / n_single
    multiple = n_multiple == 0 ? 1.0 : count(agrees .& is_multiple) / n_multiple

    w = 60
    println("-" ^ w)
    @printf("  GCM vs %s — concordância por domicílio\n", label_ref)
    println("-" ^ w)
    @printf("  %-22s %10s %12s\n", "", "Células", "Concordância")
    @printf("  %-22s %10d %11.2f%%\n", "Grid inteiro", N, 100 * overall)
    @printf("  %-22s %10d %11.2f%%\n", "Domicílio único", n_single, 100 * single)
    @printf("  %-22s %10d %11.2f%%\n", "Domicílio múltiplo", n_multiple, 100 * multiple)
    println("-" ^ w)

    return (
        overall       = overall,
        single        = single,
        multiple      = multiple,
        frac_multiple = n_multiple / N,
    )
end

# Converte o índice de uma célula fina no índice da célula grossa que a contém.
# Exige grids alinhados: fine_elements[i] = coarse_elements[i] * subgrid[i].
# Mesma numeração de get_cell_number: a primeira dimensão varia mais rápido.
function fine_to_coarse_cell(
    f :: Int64,
    fine_elements :: Vector{Int64},
    subgrid :: Vector{Int64},
    coarse_elements :: Vector{Int64}
) :: Int64
    rest = f - 1
    c = 0
    stride = 1
    for i in eachindex(fine_elements)
        idx = rest % fine_elements[i]
        rest = div(rest, fine_elements[i])
        c += div(idx, subgrid[i]) * stride
        stride *= coarse_elements[i]
    end
    return c + 1
end

"""
    refinement_reference_gcm(make_bp_fn, fine_elements) -> NamedTuple

Referência do estudo de refinamento: SCM refinado no grid `fine_elements`,
calculado como um GCM com uma única amostra no centro de cada célula
(subgrid 1, m = 1). Assim a referência e os GCMs grossos usam exatamente o
mesmo integrador. É calculada uma vez por sistema e reaproveitada em todas as
linhas do estudo.

`make_bp_fn(elements)` deve retornar um `BasinProblem` com o grid `elements`.
Retorna NamedTuple com (fine_elements, region, result, time_s).
"""
function refinement_reference_gcm(make_bp_fn, fine_elements :: Vector{Int64})
    Nf = prod(fine_elements)

    @printf("  Referência: SCM refinado %s, calculado como GCM com m = 1 (%d integrações)...\n", join(fine_elements, "×"), Nf)
    gp_ref = GcmProblem(make_bp_fn(fine_elements), ones(Int64, length(fine_elements)), 1)
    local gcm_ref, result_ref
    time_ref = @elapsed begin
        gcm_ref = build_gcmap_parallel(gp_ref)
        result_ref = find_attractors_from_gcm(gcm_ref, gp_ref)
    end

    ref_cells = vec(result_ref.cells)
    @printf("  Referência: %d atratores, %.1f s\n", length(result_ref.attractors), time_ref)
    for id in sort(unique(ref_cells))
        @printf("    bacia %3d: %.4f da área\n", id, count(==(id), ref_cells) / Nf)
    end

    return (fine_elements = fine_elements, region = gp_ref.base.region, result = result_ref, time_s = time_ref)
end

"""
    refinement_study_gcm(make_bp_fn, reference, coarse_list, m_list; sampling, seed, label, full_sampling)

Estudo de custo vs precisão proposto pelo orientador, contra a referência fina
calculada por refinement_reference_gcm.

Para cada grid grosso de `coarse_list` (alinhado ao fino) e cada m de `m_list`,
constrói o GCM sorteando m das subcélulas de cada célula. Como as subcélulas são
as células finas, todo ponto integrado pelo GCM grosso também foi integrado pela
referência: as diferenças vêm do método, não da integração. Com
`full_sampling = true`, acrescenta para cada grid grosso uma linha com m = K
(todas as subcélulas), que mostra o teto da concordância ponto a ponto.

Métricas, por linha:
  - atratores encontrados (vs referência);
  - integrações e custo relativo à referência;
  - erro máximo na fração de área de cada bacia, com a área estimada pela
    probabilidade de absorção (p.p.);
  - erro médio por célula: distância de variação total entre a probabilidade de
    absorção e a fração real das subcélulas que a referência manda a cada atrator;
  - concordância ponto a ponto: cada célula fina recebe o atrator mais provável
    da sua célula grossa e é comparada com a referência (compare_basins). Como uma
    célula grossa recebe um único rótulo, essa métrica tem um teto abaixo de 100%.

`make_bp_fn(elements)` deve retornar um `BasinProblem` com o grid `elements`.
Retorna NamedTuple com a referência e um vetor de linhas.
"""
function refinement_study_gcm(
    make_bp_fn,
    reference :: NamedTuple,
    coarse_list :: Vector{Vector{Int64}},
    m_list :: Vector{Int64};
    sampling :: Symbol = :random,
    seed :: Int64 = 1,
    label :: String = "",
    full_sampling :: Bool = false,
)
    fine_elements = reference.fine_elements
    Nf = prod(fine_elements)
    result_ref = reference.result

    ref_cells = vec(result_ref.cells)
    ref_ids = sort(unique(ref_cells))
    n_attr_ref = length(result_ref.attractors)

    w = 118
    println("\n" * "=" ^ w)
    @printf("  ESTUDO DE REFINAMENTO  ·  %s  ·  sorteio %s\n", label, sampling)
    println("=" ^ w)
    @printf("  %-10s %-8s %5s %12s %8s %10s %14s %12s %14s %9s\n",
        "Grid", "Subgrid", "m", "Integrações", "Custo", "Atratores", "Erro fração", "Erro/célula", "Ponto a ponto", "Tempo")
    println("-" ^ w)

    rows = []

    for coarse in coarse_list
        @assert all(fine_elements .% coarse .== 0) "O grid fino precisa ser múltiplo do grid grosso em cada dimensão"
        subgrid = fine_elements .÷ coarse
        K = prod(subgrid)
        Nc = prod(coarse)
        coarse_of = [fine_to_coarse_cell(f, fine_elements, subgrid, coarse) for f in 1:Nf]

        for m in (full_sampling ? vcat(m_list, K) : m_list)
            mm = min(m, K)
            gp = GcmProblem(make_bp_fn(coarse), subgrid, mm; seed=seed, sampling=sampling)

            local gcm, analysis
            elapsed = @elapsed begin
                gcm = build_gcmap_parallel(gp)
                analysis = analyze_gcm(gcm, gp)
            end
            result = basin_result_from_analysis(analysis, gp)
            coarse_cells = vec(result.cells)

            # Projeção do resultado grosso no grid fino
            projected = BasinResult(
                reshape([coarse_cells[c] for c in coarse_of], Tuple(fine_elements)),
                result.attractors,
                reference.region,
            )
            point_match = compare_basins(projected, result_ref)
            mapping = basin_id_mapping(projected, result_ref)

            # Rótulo de cada grupo persistente no grid de referência
            n_groups = length(analysis.groups)
            group_label = [g == 1 ? -1 : get(mapping, g - 1, 10_000 + g) for g in 1:n_groups]

            labels = sort(unique(vcat(ref_ids, group_label)))
            col = Dict(l => i for (i, l) in enumerate(labels))

            # Fração real (pela referência) e estimada (pela absorção) em cada célula grossa
            true_frac = zeros(Float64, Nc, length(labels))
            for f in 1:Nf
                true_frac[coarse_of[f], col[ref_cells[f]]] += 1.0 / K
            end
            est_frac = zeros(Float64, Nc, length(labels))
            for c in 1:Nc, g in 1:n_groups
                est_frac[c, col[group_label[g]]] += analysis.absorption[c, g]
            end

            cell_error = mean(0.5 .* vec(sum(abs.(true_frac .- est_frac), dims=2)))
            fraction_error = maximum(abs.(vec(mean(true_frac, dims=1)) .- vec(mean(est_frac, dims=1))))

            n_integrations = Nc * mm
            @printf("  %-10s %-8s %5d %12d %7.1f%% %6d de %d %11.2f p.p. %12.4f %13.2f%% %8.1f s\n",
                join(coarse, "×"), join(subgrid, "×"), mm, n_integrations, 100 * n_integrations / Nf,
                length(result.attractors), n_attr_ref, 100 * fraction_error, cell_error,
                100 * point_match, elapsed)

            push!(rows, (
                coarse = coarse,
                subgrid = subgrid,
                m = mm,
                n_integrations = n_integrations,
                cost = n_integrations / Nf,
                n_attractors = length(result.attractors),
                fraction_error = fraction_error,
                cell_error = cell_error,
                point_match = point_match,
                time_s = elapsed,
            ))
        end
    end

    println("=" ^ w)
    println("  Erro fração: maior diferença na área de uma bacia, com a área estimada pela probabilidade de absorção.")
    println("  Erro/célula: diferença média entre a probabilidade de absorção e a fração real de cada célula (0 = igual, 1 = oposto).")
    println("=" ^ w * "\n")

    return (n_attractors_ref = n_attr_ref, time_ref = reference.time_s, rows = rows)
end