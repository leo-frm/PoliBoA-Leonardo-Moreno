using Accessors
using DifferentialEquations
using Plots
using PlotThemes
using Plots.PlotMeasures
using LaTeXStrings
using Printf
using Colors
using Distributed
using ProgressMeter

include("utils.jl")
include("structures.jl")
include("milestone1.jl")
include("milestone2.jl")


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

# Benchmark functions
function bench_scm_helmduff(divs; e=0.077)
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
        1,
    )

    println(">>> Running SCM with $(bp.region.elements[1])x$(bp.region.elements[2]) grid")
    scm = build_scmap(bp)
    @time result = find_attractors_from_scm(scm, bp)
    return result
end

function bench_scm_spring_column(divs; omega=0.0)
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
        1,
    )

    println(">>> Running SCM with $(bp.region.elements[1])x$(bp.region.elements[2]) grid")
    scm = build_scmap(bp)
    @time result = find_attractors_from_scm(scm, bp)
    return result
end

# Test functions
function test_scm_helmduff(divs)
    GC.gc()
    bench_scm_helmduff(10) # pre-compile run
    GC.gc()
    result = bench_scm_helmduff(divs)
    plt = plotbasins(result)
    savefig(plt, "scm_basin_hd.png")
    display(plt)
end

function test_scm_spring_column(divs; omega=0.0)
    GC.gc()
    bench_scm_spring_column(10, omega=omega) # pre-compile run
    GC.gc()
    result = bench_scm_spring_column(divs, omega=omega)
    plt = plotbasins(result)
    savefig(plt, "scm_basin_sc.png")
    display(plt)
end

test_scm_helmduff(500)