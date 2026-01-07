using Accessors
using DifferentialEquations
using Plots
using PlotThemes
using Plots.PlotMeasures
using DelimitedFiles
using LaTeXStrings
using Printf
using Colors
using Distributed
using ProgressMeter

include("structures.jl")
include("utils.jl")
include("milestone2.jl")
include("milestone4.jl")

# Teste
function report_attractors(result::BasinResult)
    # Cores usadas na função plotbasins
    colors_names = ["Black (Divergent)", "White", "Red", "Yellow", "Blue", 
                    "Green", "Cyan", "Purple", "Brown", "Pink"]
    
    for attr in result.attractors
        # O ID no grid é attr.number [cite: 33, 35]
        # O índice da cor no seu array 'colors' é id + 2 
        color_idx = attr.number + 2
        color_name = color_idx <= length(colors_names) ? colors_names[color_idx] : "Custom/Other"
        
        # Conta quantas células no grid pertencem a este atrator [cite: 41]
        num_cells = count(==(attr.number), result.cells)
        
        @printf("Atrator #%d:\n", attr.number)
        @printf("  - Tipo: %s\n", attr.kind) # [cite: 40]
        @printf("  - Cor no Mapa: %s\n", color_name)
        @printf("  - Células na Bacia: %d\n", num_cells)
        @printf("  - Pontos no Ciclo: %d\n", length(attr.points)) # [cite: 34]
        println("-"^20)
    end
end

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

# Visualização
function plotbasins(result::BasinResult)
    theme(:juno)

    colors = [colorant"black", colorant"white", colorant"red",
    colorant"yellow", colorant"blue", colorant"green", 
    colorant"cyan", colorant"purple", colorant"brown", 
    colorant"pink", colorant"aquamarine4", colorant"goldenrod",
    colorant"ivory4", colorant"salmon2", colorant"midnightblue",
    colorant"magenta", colorant"rosybrown4", colorant"teal"]

    x1 = collect(LinRange(result.region.range[1][1], result.region.range[1][2], result.region.elements[1]))
    x2 = collect(LinRange(result.region.range[2][1], result.region.range[2][2], result.region.elements[2]))

    basin = transpose(reshape(result.cells,(result.region.elements[1], result.region.elements[2])))
    colorgrid = map(c -> colors[c + 2], basin)

    xformatter(x) = string(round(x, digits=2))
    yformatter(y) = string(round(y, digits=2))

    plt = heatmap(
        x1, x2, colorgrid,
        size = (1100, 900),
        legend = :false,
        left_margin = [5mm 0mm],
        tickfontsize = 18, 
        framestyle = :box,
        xformatter = xformatter,
        yformatter = yformatter,
        yflip=false,
        aspect_ratio = :none,
    )

    for attr in result.attractors
        u1 = map(u -> u[1], attr.points)
        u2 = map(u -> u[2], attr.points)
        plot!(plt, u1, u2, seriestype = :scatter, 
                markershape = :cross, markercolor = :white, markersize = 16, legend = false)
    end

    xlabel!(L"x_1(0)", xguidefontsize = 18)
    ylabel!(L"x_2(0)", yguidefontsize = 18)

    return plt
end

# Benchmark functions for parallel testing
function bench_scm_helmduff_p(divs; e=0.077, verbose=false, threads=1)
    region = BasinRegion(
        [[-1.2, 1.5], [-1.5, 1.5]], 
        [divs, floor(Int64,(3.0/2.7)*divs)],
        [false, false],
        [[-2.0, 2.5],[-2.0, 2.0]],
    )

    bp = BasinProblem(
        helmholtz_duffing!, 
        [0.1, -1.2, -0.3, 2.0, e, 1.17],
        2*pi/1.17, 
        10, 
        1000, 
        20, 
        80, 
        region,
        threads,  # ← Parâmetro agora
    )

    println(">>> Running SCM with $(bp.region.elements[1])x$(bp.region.elements[2]) grid ($(threads) thread$(threads>1 ? "s" : ""))")
    scm = build_scmap(bp)
    @time result = find_attractors_from_scm(scm, bp)

    if verbose
        report_attractors(result)
    end

    return result
end

function bench_scm_spring_column_p(divs; omega=0.0, verbose=false, threads=1)
    region = BasinRegion(
        [[-3.14, 3.14], [-1.3, 1.3]], 
        [divs, divs],
        [true, false],
        [[-3.14, 3.14],[-2.0, 2.0]],
    )

    bp = BasinProblem(
        spring_column!, 
        [0.01, 0.05, 0.8, 0.01, 0.05, omega],
        7.853981634, 
        10, 
        1000, 
        20, 
        80, 
        region,
        threads,  # ← Parâmetro agora
    )

    println(">>> Running SCM with $(bp.region.elements[1])x$(bp.region.elements[2]) grid ($(threads) thread$(threads>1 ? "s" : ""))")
    scm = build_scmap(bp)
    @time result = find_attractors_from_scm(scm, bp)

    if verbose
        report_attractors(result)
    end

    return result
end

# Test functions for parallel testing
function test_scm_helmduff_p(divs; threads=1)
    GC.gc()
    bench_scm_helmduff_p(10, verbose=false, threads=threads) # pre-compile run
    GC.gc()
    result = bench_scm_helmduff_p(divs, threads=threads)
    plt = plotbasins(result)
    savefig(plt, "scm_basin_hd.png")
    display(plt)
    writedlm("scm_basin_hd.csv", result.cells, ',')
end

function test_scm_spring_column_p(divs; omega=0.0, threads=1)
    GC.gc()
    bench_scm_spring_column_p(10, omega=omega, verbose=false, threads=threads) # pre-compile run
    GC.gc()
    result = bench_scm_spring_column_p(divs, omega=omega, threads=threads)
    plt = plotbasins(result)
    savefig(plt, "scm_basin_sc.png")
    display(plt)
    writedlm("scm_basin_sc.csv", result.cells, ',')
end

#Test milestone4
function test_parallel_correctness(divs=100)
    println("\n" * "="^50)
    println("TESTE: Correção da Paralelização")
    println("="^50)
    
    # Serial
    println("\n[1 thread]")
    result_1t = bench_scm_helmduff_p(divs, verbose=false, threads=1)
    
    # Paralelo
    println("\n[4 threads]")
    result_4t = bench_scm_helmduff_p(divs, verbose=false, threads=4)
    
    # Comparar
    cells_match = result_1t.cells == result_4t.cells
    attractors_match = length(result_1t.attractors) == length(result_4t.attractors)
    
    println("\nResultado:")
    println("  Células idênticas: $(cells_match ? "✓" : "✗")")
    println("  Nº atratores igual: $(attractors_match ? "✓" : "✗")")
    
    if cells_match && attractors_match
        println("\n✅ PASSOU: Resultados idênticos!")
    else
        println("\n❌ FALHOU: Resultados diferentes!")
    end
    
    println("="^50)
    return cells_match && attractors_match
end

# Teste 2: Speedup
function test_parallel_speedup(divs=200)
    println("\n" * "="^50)
    println("TESTE: Speedup com Paralelização")
    println("="^50)
    
    max_threads = min(Threads.nthreads(), 8)
    println("\nThreads disponíveis: $(Threads.nthreads())")
    println("\nThreads | Tempo  | Speedup")
    println("-"^32)
    
    times = []
    for nt in 1:max_threads
        # Warm-up
        bench_scm_helmduff_p(10, verbose=false, threads=nt)
        GC.gc()
        
        # Medição
        t = @elapsed bench_scm_helmduff_p(divs, verbose=false, threads=nt)
        push!(times, t)
        
        speedup = times[1] / t
        @printf("%2d      | %.2fs | %.2fx\n", nt, t, speedup)
    end
    
    println("="^50)
    return times
end

# Teste completo Milestone 4
function test_milestone4(divs=100)
    println("\n🎯 TESTES MILESTONE 4 - PARALELIZAÇÃO 🎯\n")
    
    # Teste 1: Correção
    passed = test_parallel_correctness(divs)
    
    # Teste 2: Speedup
    test_parallel_speedup(divs)
    
    println("\n" * "="^50)
    if passed
        println("✅ Milestone 4 validada com sucesso!")
    else
        println("❌ Milestone 4 com problemas!")
    end
    println("="^50 * "\n")
end

test_scm_helmduff_p(200, threads=1)
test_scm_helmduff_p(200, threads=4)
test_scm_spring_column_p(200, threads=1)
test_scm_spring_column_p(200, threads=4)