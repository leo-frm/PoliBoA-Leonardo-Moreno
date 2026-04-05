#   1. Roda benchmark_method no SCM (helmholtz-duffing, 200 divs)
#   2. Roda benchmark_method no código original (populate_basins)
#   3. Compara os dois benchmarks (compare_benchmarks)
#   4. Compara as bacias célula a célula (compare_basins → match_fraction)
#   5. Escalabilidade de memória com tamanho do grid (memory_scaling)
#   6. (Opcional) Análise de escalabilidade por threads (thread_scaling)

using Accessors
using DifferentialEquations
using ProgressMeter
using Printf

include("structures.jl")
include("utils.jl")
include("milestone2.jl")
include("milestone4.jl")
include("milestone5.jl")
include("poliboa_original.jl")

#Teste
# Helmholtz-Duffing
function helmholtz_duffing!(du, u, (a,b,c,d,e,ω), t)
    du[1] = u[2]
    du[2] = -a*u[2] - b*u[1] - c*u[1]^2 - d*u[1]^3 + e*sin(ω*t)
end

# Spring Column
function spring_column!(dx, x, (c, p, alpha, q, q1, omega), t)
    dx[1] = x[2]
    dx[2] = -c * x[2] + p * sin(x[1]) - (1 - 1 / sqrt(1 + alpha * sin(x[1])) - (q + q1 * sin(omega * t))) * cos(x[1])
end

function make_helmduff_bp(divs; e=0.077, threads=1)
    region = BasinRegion(
        [[-1.2, 1.5], [-1.5, 1.5]],
        [divs, floor(Int64, (3.0/2.7) * divs)],
        [false, false],
        [[-2.0, 2.5], [-2.0, 2.0]],
    )

    return BasinProblem(
        helmholtz_duffing!,
        [0.1, -1.2, -0.3, 2.0, e, 1.17],
        2*pi/1.17,
        10,
        1000,
        20,
        80,
        region,
        threads,
    )
end

function make_spring_column_bp(divs; omega=0.0, threads=1)
    region = BasinRegion(
        [[-3.14, 3.14], [-1.3, 1.3]],
        [divs, divs],
        [true, false],
        [[-3.14, 3.14], [-2.0, 2.0]],
    )

    return BasinProblem(
        spring_column!,
        [0.01, 0.05, 0.8, 0.01, 0.05, omega],
        7.853981634,
        100,
        1000,
        20,
        80,
        region,
        threads,
    )
end

"""Roda o pipeline SCM completo (build + find_attractors) e retorna o SCM."""
function run_scm(bp::BasinProblem)
    scm = build_scmap_parallel(bp)
    find_attractors_from_scm(scm, bp)
    return scm
end

# =============================================================================
# TESTE 1: Benchmark do método SCM
# =============================================================================

function test_benchmark_scm(; divs=200, trials=5, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 1: Benchmark do SCM — $label (divs=$divs)")
    println("=" ^ 60)

    bp = make_bp(divs)
    bench = benchmark_method(trials=trials) do
        run_scm(bp)
    end

    println(bench)
    return bench
end

# =============================================================================
# TESTE 2: Comparação de bacias SCM vs código original
# =============================================================================

function test_compare_basins(; divs=200, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 2: Comparação de bacias SCM vs Original — $label (divs=$divs)")
    println("=" ^ 60)

    bp = make_bp(divs)

    println("  Rodando SCM...")
    scm = build_scmap_parallel(bp)
    result_scm = find_attractors_from_scm(scm, bp)

    println("  Rodando método original...")
    result_og = populate_basins(bp)

    match = compare_basins(result_scm, result_og)
    @printf("  → match_fraction = %.4f (%.1f%% de coincidência)\n", match, match * 100)

    return match
end

# =============================================================================
# TESTE 3: Comparação de benchmarks SCM vs Original
# =============================================================================

function test_compare_benchmarks(; divs=200, trials=3, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 3: Comparação de Benchmarks SCM vs Original — $label (divs=$divs)")
    println("=" ^ 60)

    bp = make_bp(divs)

    println("\n  >> Benchmarking SCM:")
    bench_scm = benchmark_method(trials=trials) do
        run_scm(bp)
    end

    println("\n  >> Benchmarking Original (populate_basins):")
    bench_og = benchmark_method(trials=trials) do
        populate_basins(bp)
    end

    ratios = compare_benchmarks(bench_scm, bench_og)
    return ratios
end

# =============================================================================
# TESTE 4: Escalabilidade de memória
# =============================================================================

function test_memory_scaling(; divs_range=[100, 200, 400], make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 4: Escalabilidade de Memória — $label (divs=$divs_range)")
    println("=" ^ 60)

    bp_list = [make_bp(d) for d in divs_range]

    memory_scaling(
        build_scmap_parallel,
        find_attractors_from_scm,
        bp_list;
        method_name = "SCM — $label",
    )
end

# =============================================================================
# TESTE 5: Escalabilidade por threads
# =============================================================================

function test_thread_scaling(; divs=200, trials=3, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
    max_t = Threads.nthreads()
    println("\n" * "=" ^ 60)
    println("TESTE 5: Escalabilidade por Threads — $label (divs=$divs, max_threads=$max_t)")
    println("=" ^ 60)

    if max_t == 1
        println("  ⚠ Julia iniciada com 1 thread. Use `julia -t N` para testar scaling.")
        println("  Pulando teste.")
        return nothing
    end

    counts = [1]
    t = 2
    while t <= max_t
        push!(counts, t)
        t *= 2
    end
    if last(counts) != max_t
        push!(counts, max_t)
    end

    results = thread_scaling(
        build_scmap_parallel,
        find_attractors_from_scm,
        nt -> make_bp(divs, threads=nt);
        thread_counts = counts,
        trials = trials,
    )

    return results
end

# =============================================================================
# EXECUÇÃO — HELMHOLTZ-DUFFING
# =============================================================================

println("\n" * "▓" ^ 60)
println("▓  HELMHOLTZ-DUFFING")
println("▓" ^ 60)

bench_hd  = test_benchmark_scm(divs=200, trials=3, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
match_hd  = test_compare_basins(divs=200, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
ratios_hd = test_compare_benchmarks(divs=200, trials=3, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
test_memory_scaling(divs_range=[100, 200, 400], make_bp=make_helmduff_bp, label="Helmholtz-Duffing")
test_thread_scaling(divs=200, trials=3, make_bp=make_helmduff_bp, label="Helmholtz-Duffing")

# =============================================================================
# EXECUÇÃO — SPRING COLUMN
# =============================================================================

println("\n" * "▓" ^ 60)
println("▓  SPRING COLUMN")
println("▓" ^ 60)

bench_sc  = test_benchmark_scm(divs=200, trials=3, make_bp=make_spring_column_bp, label="Spring Column")
match_sc  = test_compare_basins(divs=200, make_bp=make_spring_column_bp, label="Spring Column")
ratios_sc = test_compare_benchmarks(divs=200, trials=3, make_bp=make_spring_column_bp, label="Spring Column")
test_memory_scaling(divs_range=[100, 200, 400], make_bp=make_spring_column_bp, label="Spring Column")
test_thread_scaling(divs=200, trials=3, make_bp=make_spring_column_bp, label="Spring Column")

# =============================================================================
# OUTPUT
# =============================================================================
println("\n" * "=" ^ 60)
println("Resumo Final")
println("=" ^ 60)
@printf("  [HD] match_fraction:  %.4f | tempo médio: %.4f s\n", match_hd, bench_hd.mean_time_s)
@printf("  [SC] match_fraction:  %.4f | tempo médio: %.4f s\n", match_sc, bench_sc.mean_time_s)
println("=" ^ 60)
println("Done.")