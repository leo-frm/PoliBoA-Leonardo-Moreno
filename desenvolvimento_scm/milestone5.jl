# Funções de validação para o método SCM:
#   - benchmark_method:    mede tempo médio, desvio padrão, bytes alocados e nº de integrações
#   - compare_basins:      compara dois BasinResult e retorna a fração de células coincidentes
#   - compare_benchmarks:  compara dois resultados de benchmark (tempo, bytes, integrações)
#   - memory_scaling:      mede como a memória escala com o tamanho do grid
#   - thread_scaling:      mede como a performance escala com o número de threads

using Printf
using Statistics

struct BenchmarkResult
    mean_time_s     :: Float64      # Tempo médio em segundos
    std_time_s      :: Float64      # Desvio padrão do tempo em segundos
    bytes_allocated :: Int64        # Bytes alocados (da última execução medida)
    n_integrations  :: Int64        # Número de integrações realizadas (= computed_cells do SCM)
    grid_size       :: String       # Descrição do grid (ex: "200x222")
    trials          :: Int64        # Número de execuções realizadas
end

"""
    benchmark_method(run_fn; trials=5, warmup=1) -> BenchmarkResult

Executa `run_fn()` (uma função sem argumentos que retorna um `SimpleCellMap`)
`warmup` vezes para pré-compilar e `trials` vezes para medir.

Retorna: tempo médio, desvio padrão, bytes alocados e nº de integrações.
"""
function benchmark_method(run_fn; trials::Int=5, warmup::Int=1)
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
        if result isa SimpleCellMap
            last_integrations = result.computed_cells
            last_grid_size = join(result.range.elements, "x")
        elseif result isa BasinResult
            last_integrations = prod(result.region.elements)
            last_grid_size = join(result.region.elements, "x")
        elseif result isa Tuple
            # Suporta retornos tipo (scm, basin_result)
            for r in result
                if r isa SimpleCellMap
                    last_integrations = r.computed_cells
                    last_grid_size = join(r.range.elements, "x")
                    break
                end
            end
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

function Base.show(io::IO, b::BenchmarkResult)
    println(io, "BenchmarkResult (grid $(b.grid_size), $(b.trials) trials):")
    @printf(io, "  Tempo médio:    %.4f ± %.4f s\n", b.mean_time_s, b.std_time_s)
    @printf(io, "  Bytes alocados: %.2f MB\n", b.bytes_allocated / 1e6)
    @printf(io, "  Integrações:    %d\n", b.n_integrations)
end

"""
    compare_basins(result_a::BasinResult, result_b::BasinResult) -> Float64

Compara célula a célula os dois resultados e retorna a fração de coincidência
(`match_fraction`), ignorando a numeração absoluta dos atratores: o que importa
é se duas células que pertencem ao mesmo atrator em `a` também pertencem ao mesmo
atrator em `b`.
"""
function compare_basins(result_a::BasinResult, result_b::BasinResult) :: Float64
    cells_a = result_a.cells
    cells_b = result_b.cells

    @assert length(cells_a) == length(cells_b) "Os dois resultados precisam ter o mesmo número de células"

    total = length(cells_a)
    if total == 0
        return 1.0
    end

    # Constrói mapeamento: atrator_a -> atrator_b mais frequente
    # Para cada par (a_id, b_id), conta quantas células compartilham
    pair_counts = Dict{Tuple{Int64, Int64}, Int64}()
    for i in eachindex(cells_a)
        key = (cells_a[i], cells_b[i])
        pair_counts[key] = get(pair_counts, key, 0) + 1
    end

    # Para cada atrator de A, encontra o melhor match em B
    best_match = Dict{Int64, Tuple{Int64, Int64}}() # a_id -> (b_id, count)
    for ((a_id, b_id), count) in pair_counts
        if !haskey(best_match, a_id) || count > best_match[a_id][2]
            best_match[a_id] = (b_id, count)
        end
    end

    # Conta células que batem com o melhor mapeamento
    mapping = Dict(a_id => bm[1] for (a_id, bm) in best_match)
    matches = count(i -> get(mapping, cells_a[i], -1) == cells_b[i], eachindex(cells_a))

    return matches / total
end

"""
    compare_benchmarks(bench_a::BenchmarkResult, bench_b::BenchmarkResult)

Compara dois BenchmarkResults e imprime: razão de tempo, razão de bytes e
razão de integrações. Retorna NamedTuple com os valores.
"""
function compare_benchmarks(bench_a::BenchmarkResult, bench_b::BenchmarkResult)
    time_ratio  = bench_a.mean_time_s / bench_b.mean_time_s
    bytes_ratio = bench_a.bytes_allocated / bench_b.bytes_allocated
    integ_ratio = bench_a.n_integrations / bench_b.n_integrations

    println("=" ^ 60)
    println("Comparação de Benchmarks")
    println("=" ^ 60)
    @printf("  %-25s %12s %12s %10s\n", "", "Método A", "Método B", "Razão A/B")
    println("-" ^ 60)
    @printf("  %-25s %10.4f s %10.4f s %10.2fx\n",
        "Tempo médio", bench_a.mean_time_s, bench_b.mean_time_s, time_ratio)
    @printf("  %-25s %9.2f MB %9.2f MB %10.2fx\n",
        "Bytes alocados", bench_a.bytes_allocated / 1e6, bench_b.bytes_allocated / 1e6, bytes_ratio)
    @printf("  %-25s %12d %12d %10.2fx\n",
        "Integrações", bench_a.n_integrations, bench_b.n_integrations, integ_ratio)
    println("=" ^ 60)

    return (
        time_ratio  = time_ratio,
        bytes_ratio = bytes_ratio,
        integ_ratio = integ_ratio,
    )
end

"""
    memory_scaling(build_fn, find_fn, bp_list; method_name)

Mede como a memória alocada escala com o tamanho do grid.
`bp_list` é um vetor de BasinProblem com resoluções crescentes.

Métrica derivada `MB/kCell` permite comparar a eficiência de memória
independentemente da resolução.
"""
function memory_scaling(
    build_fn    :: Function,
    find_fn     :: Function,
    bp_list     :: Vector{<:BasinProblem};
    method_name :: String = "SCM"
)
    w = 72
    println("\n" * "=" ^ w)
    @printf("  ESCALABILIDADE DE MEMÓRIA  ·  %s\n", method_name)
    println("=" ^ w)
    @printf("  %-16s %10s %12s %12s %14s\n",
        "Grid", "Células", "Memória(MB)", "MB/kCell", "Integrações")
    println("-" ^ w)

    for bp in bp_list
        GC.gc()
        total_cells = prod(bp.region.elements)
        grid_str    = join(bp.region.elements, "×")

        local scm, result
        alloc_bytes = @allocated begin
            scm    = build_fn(bp)
            result = find_fn(scm, bp)
        end

        mb_total     = alloc_bytes / 1024^2
        mb_per_kcell = mb_total / (total_cells / 1_000.0)

        @printf("  %-16s %10d %12.1f %12.3f %14d\n",
            grid_str, total_cells, mb_total, mb_per_kcell, scm.computed_cells)
        GC.gc()
    end

    println("=" ^ w * "\n")
end

"""
    thread_scaling(make_bp_fn, thread_counts; divs=200, trials=3)

Mede o tempo de execução do SCM para diferentes quantidades de threads
e calcula o speedup relativo à execução single-thread.

`make_bp_fn(threads)` deve retornar um `BasinProblem` configurado para
o número de threads dado.

Retorna vetor de NamedTuples com (threads, mean_time, speedup).
"""
function thread_scaling(
    build_fn    :: Function,
    find_fn     :: Function,
    make_bp_fn;
    thread_counts = [1, 2, 4],
    trials::Int = 3
)
    results = NamedTuple{(:threads, :mean_time, :speedup), Tuple{Int, Float64, Float64}}[]
    base_time = 0.0

    println("=" ^ 60)
    println("Análise de Escalabilidade por Threads")
    println("=" ^ 60)
    @printf("  %-10s %12s %10s\n", "Threads", "Tempo médio", "Speedup")
    println("-" ^ 35)

    for nt in thread_counts
        bp = make_bp_fn(nt)
        bench = benchmark_method(trials=trials, warmup=1) do
            scm = build_fn(bp)
            find_fn(scm, bp)
            return scm
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