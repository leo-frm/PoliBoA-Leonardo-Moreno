# =============================================================================
# tests_milestone5.jl — Testes de Validação (Milestone 5)
# =============================================================================
#
# COMO RODAR:
#   cd desenvolvimento_scm/
#   julia tests_milestone5.jl
#
# PRÉ-REQUISITOS:
#   - Julia com pacotes: DifferentialEquations, Plots, Accessors, etc.
#   - Arquivos no mesmo diretório: structures.jl, utils.jl, milestone2.jl, milestone4.jl, milestone5.jl
#   - Para testes de thread_scaling: iniciar Julia com threads, ex: julia -t 4 tests_milestone5.jl
#
# O QUE ESTE SCRIPT FAZ:
#   1. Roda benchmark_method no SCM (helmholtz-duffing, 200 divs)
#   2. Roda benchmark_method no código original (populate_basins)
#   3. Compara os dois benchmarks (compare_benchmarks)
#   4. Compara as bacias célula a célula (compare_basins → match_fraction)
#   5. Gera perfil de alocação (allocation_profile)
#   6. (Opcional) Análise de escalabilidade por threads (thread_scaling)
#
# =============================================================================

using Accessors
using DifferentialEquations
using Printf

include("structures.jl")
include("utils.jl")
include("milestone2.jl")
include("milestone4.jl")
include("milestone5.jl")

# =============================================================================
# Dinâmicas (mesmo setup dos outros testes)
# =============================================================================

function helmholtz_duffing!(du, u, (a,b,c,d,e,ω), t)
    du[1] = u[2]
    du[2] = -a*u[2] - b*u[1] - c*u[1]^2 - d*u[1]^3 + e*sin(ω*t)
end

# =============================================================================
# Helpers para criar BasinProblems reutilizáveis
# =============================================================================

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

# =============================================================================
# Funções de execução que serão passadas ao benchmark_method
# =============================================================================

"""Roda o pipeline SCM completo (build + find_attractors) e retorna o SCM."""
function run_scm(bp::BasinProblem)
    scm = build_scmap_parallel(bp)
    find_attractors_from_scm(scm, bp)
    return scm
end

# =============================================================================
# TESTE 1: Benchmark do método SCM
# =============================================================================

function test_benchmark_scm(; divs=200, trials=5)
    println("\n" * "=" ^ 60)
    println("TESTE 1: Benchmark do SCM (divs=$divs)")
    println("=" ^ 60)

    bp = make_helmduff_bp(divs)
    bench = benchmark_method(trials=trials) do
        run_scm(bp)
    end

    println(bench)
    return bench
end

# =============================================================================
# TESTE 2: Comparação de bacias SCM vs código original
# =============================================================================

function test_compare_basins(; divs=200)
    println("\n" * "=" ^ 60)
    println("TESTE 2: Comparação de bacias SCM vs Original (divs=$divs)")
    println("=" ^ 60)

    bp = make_helmduff_bp(divs)

    # Resultado SCM
    println("  Rodando SCM...")
    scm = build_scmap_parallel(bp)
    result_scm = find_attractors_from_scm(scm, bp)

    # Resultado Original (populate_basins)
    println("  Rodando método original...")
    result_og = populate_basins(bp)

    # Compara
    match = compare_basins(result_scm, result_og)
    @printf("  → match_fraction = %.4f (%.1f%% de coincidência)\n", match, match * 100)

    return match
end

# =============================================================================
# TESTE 3: Comparação de benchmarks SCM vs Original
# =============================================================================

function test_compare_benchmarks(; divs=200, trials=3)
    println("\n" * "=" ^ 60)
    println("TESTE 3: Comparação de Benchmarks SCM vs Original (divs=$divs)")
    println("=" ^ 60)

    bp = make_helmduff_bp(divs)

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
# TESTE 4: Perfil de alocação
# =============================================================================

function test_allocation_profile(; divs=200)
    println("\n" * "=" ^ 60)
    println("TESTE 4: Perfil de Alocação (divs=$divs)")
    println("=" ^ 60)

    bp = make_helmduff_bp(divs)

    profile_scm = allocation_profile(label="SCM") do
        run_scm(bp)
    end

    profile_og = allocation_profile(label="Original") do
        populate_basins(bp)
    end

    return (scm=profile_scm, original=profile_og)
end

# =============================================================================
# TESTE 5: Escalabilidade por threads (requer julia -t N)
# =============================================================================

function test_thread_scaling(; divs=200, trials=3)
    max_t = Threads.nthreads()
    println("\n" * "=" ^ 60)
    println("TESTE 5: Escalabilidade por Threads (divs=$divs, max_threads=$max_t)")
    println("=" ^ 60)

    if max_t == 1
        println("  ⚠ Julia iniciada com 1 thread. Use `julia -t N` para testar scaling.")
        println("  Pulando teste.")
        return nothing
    end

    # Gera lista de threads: 1, 2, 4, ... até max_t
    counts = [1]
    t = 2
    while t <= max_t
        push!(counts, t)
        t *= 2
    end
    if last(counts) != max_t
        push!(counts, max_t)
    end

    results = thread_scaling(thread_counts=counts, trials=trials) do nt
        make_helmduff_bp(divs, threads=nt)
    end

    return results
end

# =============================================================================
# EXECUÇÃO
# =============================================================================

bench_scm  = test_benchmark_scm(divs=200, trials=3)
match      = test_compare_basins(divs=200)
ratios     = test_compare_benchmarks(divs=200, trials=3)
profiles   = test_allocation_profile(divs=200)
scaling    = test_thread_scaling(divs=200, trials=3)

println("\n" * "=" ^ 60)
println("Resumo Final")
println("=" ^ 60)
@printf("  match_fraction:  %.4f\n", match)
@printf("  SCM tempo médio: %.4f s\n", bench_scm.mean_time_s)
@printf("  Razão tempo:     %.2fx\n", ratios.time_ratio)
println("=" ^ 60)
println("Done.")