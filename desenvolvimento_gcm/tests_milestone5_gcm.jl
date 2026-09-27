#   1. Roda benchmark_method_gcm no GCM (helmholtz-duffing, 200 divs)
#   2. Compara as bacias GCM vs original e GCM vs SCM (compare_basins → match_fraction),
#      separando a concordância por domicílio (compare_basins_by_domicile)
#   3. Compara os benchmarks GCM vs original e GCM vs SCM (compare_benchmarks)
#   4. Escalabilidade de memória com tamanho do grid (memory_scaling_gcm)
#   5. (Opcional) Análise de escalabilidade por threads (thread_scaling_gcm)

using Accessors
using DifferentialEquations
using ProgressMeter
using Printf

include("../desenvolvimento_scm/structures.jl")
include("../desenvolvimento_scm/utils.jl")
include("../desenvolvimento_scm/milestone2.jl")
include("../desenvolvimento_scm/milestone4.jl")
include("../desenvolvimento_scm/milestone5.jl")
include("../desenvolvimento_scm/poliboa_original.jl")
include("structures_gcm.jl")
include("utils_gcm.jl")
include("milestone1_gcm.jl")
include("milestone2_gcm.jl")
include("milestone4_gcm.jl")
include("milestone5_gcm.jl")

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

make_helmduff_gp(divs; S=3, threads=1) = GcmProblem(make_helmduff_bp(divs, threads=threads), S)
make_spring_column_gp(divs; S=3, threads=1) = GcmProblem(make_spring_column_bp(divs, threads=threads), S)

"""Roda o pipeline GCM completo (build + find_attractors) e retorna o GCM."""
function run_gcm(gp::GcmProblem)
    gcm = build_gcmap_parallel(gp)
    find_attractors_from_gcm(gcm, gp)
    return gcm
end

"""Roda o pipeline SCM completo (build + find_attractors) e retorna o SCM."""
function run_scm(bp::BasinProblem)
    scm = build_scmap_parallel(bp)
    find_attractors_from_scm(scm, bp)
    return scm
end

# =============================================================================
# TESTE 1: Benchmark do método GCM
# =============================================================================

function test_benchmark_gcm(; divs=200, S=3, trials=5, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 1: Benchmark do GCM — $label (divs=$divs, S=$S)")
    println("=" ^ 60)

    gp = make_gp(divs, S=S)
    bench = benchmark_method_gcm(trials=trials) do
        run_gcm(gp)
    end

    println(bench)
    return bench
end

# =============================================================================
# TESTE 2: Comparação de bacias GCM vs código original e GCM vs SCM
# =============================================================================

function test_compare_basins(; divs=200, S=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 2: Comparação de bacias GCM vs Original e vs SCM — $label (divs=$divs, S=$S)")
    println("=" ^ 60)

    gp = make_gp(divs, S=S)
    bp = gp.base

    println("  Rodando GCM...")
    gcm = build_gcmap_parallel(gp)
    analysis = analyze_gcm(gcm, gp)
    result_gcm = basin_result_from_analysis(analysis, gp)
    report_gcm_analysis(analysis)

    println("  Rodando SCM...")
    scm = build_scmap_parallel(bp)
    result_scm = find_attractors_from_scm(scm, bp)

    println("  Rodando método original...")
    result_og = populate_basins(bp)

    match_og = compare_basins(result_gcm, result_og)
    @printf("  → GCM vs Original: match_fraction = %.4f (%.1f%% de coincidência)\n", match_og, match_og * 100)

    match_scm = compare_basins(result_gcm, result_scm)
    @printf("  → GCM vs SCM:      match_fraction = %.4f (%.1f%% de coincidência)\n", match_scm, match_scm * 100)

    dom_og = compare_basins_by_domicile(result_gcm, result_og, analysis, label_ref="Original")
    dom_scm = compare_basins_by_domicile(result_gcm, result_scm, analysis, label_ref="SCM")

    return (
        match_og  = match_og,
        match_scm = match_scm,
        dom_og    = dom_og,
        dom_scm   = dom_scm,
    )
end

# =============================================================================
# TESTE 3: Comparação de benchmarks GCM vs Original e GCM vs SCM
# =============================================================================

function test_compare_benchmarks(; divs=200, S=3, trials=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 3: Comparação de Benchmarks GCM vs Original e vs SCM — $label (divs=$divs, S=$S)")
    println("=" ^ 60)

    gp = make_gp(divs, S=S)
    bp = gp.base

    println("\n  >> Benchmarking GCM:")
    bench_gcm = benchmark_method_gcm(trials=trials) do
        run_gcm(gp)
    end

    println("\n  >> Benchmarking SCM:")
    bench_scm = benchmark_method(trials=trials) do
        run_scm(bp)
    end

    println("\n  >> Benchmarking Original (populate_basins):")
    bench_og = benchmark_method(trials=trials) do
        populate_basins(bp)
    end

    println("\n  GCM (A) vs Original (B):")
    ratios_og = compare_benchmarks(bench_gcm, bench_og)

    println("\n  GCM (A) vs SCM (B):")
    ratios_scm = compare_benchmarks(bench_gcm, bench_scm)

    return (
        ratios_og  = ratios_og,
        ratios_scm = ratios_scm,
    )
end

# =============================================================================
# TESTE 4: Escalabilidade de memória
# =============================================================================

function test_memory_scaling(; divs_range=[100, 200, 400], S=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
    println("\n" * "=" ^ 60)
    println("TESTE 4: Escalabilidade de Memória — $label (divs=$divs_range, S=$S)")
    println("=" ^ 60)

    gp_list = [make_gp(d, S=S) for d in divs_range]

    memory_scaling_gcm(
        build_gcmap_parallel,
        find_attractors_from_gcm,
        gp_list;
        method_name = "GCM — $label",
    )
end

# =============================================================================
# TESTE 5: Escalabilidade por threads
# =============================================================================

function test_thread_scaling(; divs=200, S=3, trials=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
    max_t = Threads.nthreads()
    println("\n" * "=" ^ 60)
    println("TESTE 5: Escalabilidade por Threads — $label (divs=$divs, S=$S, max_threads=$max_t)")
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

    results = thread_scaling_gcm(
        build_gcmap_parallel,
        find_attractors_from_gcm,
        nt -> make_gp(divs, S=S, threads=nt);
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

bench_hd  = test_benchmark_gcm(divs=200, S=3, trials=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
match_hd  = test_compare_basins(divs=200, S=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
ratios_hd = test_compare_benchmarks(divs=200, S=3, trials=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
test_memory_scaling(divs_range=[100, 200, 400], S=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")
test_thread_scaling(divs=200, S=3, trials=3, make_gp=make_helmduff_gp, label="Helmholtz-Duffing")

# =============================================================================
# EXECUÇÃO — SPRING COLUMN
# =============================================================================

println("\n" * "▓" ^ 60)
println("▓  SPRING COLUMN")
println("▓" ^ 60)

bench_sc  = test_benchmark_gcm(divs=200, S=3, trials=3, make_gp=make_spring_column_gp, label="Spring Column")
match_sc  = test_compare_basins(divs=200, S=3, make_gp=make_spring_column_gp, label="Spring Column")
ratios_sc = test_compare_benchmarks(divs=200, S=3, trials=3, make_gp=make_spring_column_gp, label="Spring Column")
test_memory_scaling(divs_range=[100, 200, 400], S=3, make_gp=make_spring_column_gp, label="Spring Column")
test_thread_scaling(divs=200, S=3, trials=3, make_gp=make_spring_column_gp, label="Spring Column")

# =============================================================================
# OUTPUT
# =============================================================================
println("\n" * "=" ^ 60)
println("Resumo Final")
println("=" ^ 60)
@printf("  [HD] match vs Original: %.4f | vs SCM: %.4f | tempo médio: %.4f s\n",
    match_hd.match_og, match_hd.match_scm, bench_hd.mean_time_s)
@printf("  [HD] concordância no domicílio único vs Original: %.4f | fronteira: %.2f%% das células\n",
    match_hd.dom_og.single, 100 * match_hd.dom_og.frac_multiple)
@printf("  [SC] match vs Original: %.4f | vs SCM: %.4f | tempo médio: %.4f s\n",
    match_sc.match_og, match_sc.match_scm, bench_sc.mean_time_s)
@printf("  [SC] concordância no domicílio único vs Original: %.4f | fronteira: %.2f%% das células\n",
    match_sc.dom_og.single, 100 * match_sc.dom_og.frac_multiple)
println("=" ^ 60)
println("Done.")
