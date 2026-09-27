# Funções de validação para o método GCM:
#   - benchmark_method_gcm:        mede tempo médio, desvio padrão, bytes alocados e nº de integrações
#   - memory_scaling_gcm:          mede como a memória escala com o tamanho do grid
#   - thread_scaling_gcm:          mede como a performance escala com o número de threads
#   - compare_basins_by_domicile:  separa a concordância entre células de domicílio único e múltiplo
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
por célula (S^d), e não apenas o número de células como no SCM.
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
