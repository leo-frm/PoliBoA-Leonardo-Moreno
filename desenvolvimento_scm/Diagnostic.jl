using DifferentialEquations
using Plots
using Printf

include("structures.jl")
include("utils.jl")
include("milestone1.jl")
include("milestone2.jl")

# Helmholtz-Duffing
function helmholtz_duffing!(du, u, (a,b,c,d,e,ω), t)
    du[1] = u[2] 
    du[2] = -a*u[2] - b*u[1] - c*u[1]^2 - d*u[1]^3 + e*sin(ω*t)
end

# Matriz de distâncias entre centróides
function distance_matrix(result::BasinResult, bp::BasinProblem)
    ids = [attr.number for attr in result.attractors]
    n = length(ids)
    
    # Calcula centróide de cada atrator (integra e pega média dos últimos pontos)
    centroids = Dict{Int, Vector{Float64}}()
    
    for attr in result.attractors
        # Usa o primeiro ponto do atrator como condição inicial
        u0 = attr.points[1]
        
        prob = ODEProblem(bp.f, u0, (0.0, 100 * bp.period), bp.params)
        sol = solve(prob, Tsit5(), saveat=bp.period)
        
        # Média dos últimos 50 pontos
        last_points = sol.u[50:end]
        cx = sum(p[1] for p in last_points) / length(last_points)
        cy = sum(p[2] for p in last_points) / length(last_points)
        centroids[attr.number] = [cx, cy]
    end
    
    # Imprime matriz
    println("\nMatriz de distâncias entre centróides:")
    println("-" ^ 60)
    
    # Header
    print("       ")
    for id in ids
        @printf("  Attr%-3d", id)
    end
    println()
    
    # Linhas
    for id1 in ids
        @printf("Attr%-3d", id1)
        for id2 in ids
            dist = sqrt((centroids[id1][1] - centroids[id2][1])^2 + 
                       (centroids[id1][2] - centroids[id2][2])^2)
            @printf("  %6.4f", dist)
        end
        println()
    end
    
    return centroids
end

# Compara dois atratores específicos
function compare_attractors_timeseries(result::BasinResult, bp::BasinProblem, id1::Int, id2::Int)
    # Encontra os atratores
    attr1 = nothing
    attr2 = nothing
    for attr in result.attractors
        if attr.number == id1
            attr1 = attr
        end
        if attr.number == id2
            attr2 = attr
        end
    end
    
    if attr1 === nothing || attr2 === nothing
        println("Atrator alvo não encontrado!")
        return nothing
    end
    
    # Condições iniciais
    u0_1 = attr1.points[1]
    u0_2 = attr2.points[1]
    
    println("\nComparando atratores $id1 vs $id2")
    println("   CI atrator $id1: ($(round(u0_1[1], digits=4)), $(round(u0_1[2], digits=4)))")
    println("   CI atrator $id2: ($(round(u0_2[1], digits=4)), $(round(u0_2[2], digits=4)))")
    
    # Integra
    tempo = 100 * bp.period
    sol1 = solve(ODEProblem(bp.f, u0_1, (0.0, tempo), bp.params), Tsit5(), saveat=bp.period)
    sol2 = solve(ODEProblem(bp.f, u0_2, (0.0, tempo), bp.params), Tsit5(), saveat=bp.period)
    
    # Extrai dados após transiente
    x1_1 = [sol1.u[i][1] for i in 20:length(sol1.u)]
    x2_1 = [sol1.u[i][2] for i in 20:length(sol1.u)]
    x1_2 = [sol2.u[i][1] for i in 20:length(sol2.u)]
    x2_2 = [sol2.u[i][2] for i in 20:length(sol2.u)]
    
    # Calcula distância entre médias
    m1 = (sum(x1_1)/length(x1_1), sum(x2_1)/length(x2_1))
    m2 = (sum(x1_2)/length(x1_2), sum(x2_2)/length(x2_2))
    dist = sqrt((m1[1]-m2[1])^2 + (m1[2]-m2[2])^2)
    
    println("\nMédia atrator $id1: ($(round(m1[1], digits=4)), $(round(m1[2], digits=4)))")
    println("Média atrator $id2: ($(round(m2[1], digits=4)), $(round(m2[2], digits=4)))")
    println("Distância: $(round(dist, digits=6))")
    
    if dist < 0.01
        println("Atrator duplicado")
    elseif dist < 0.1
        println("Possivelmente o mesmo")
    else
        println("Atratores diferentes")
    end
    
    # Plot
    p1 = scatter(x1_1, x2_1, label="Atrator $id1", markersize=3, alpha=0.6,
                 title="Espaço de Fases", xlabel="x1", ylabel="x2")
    scatter!(p1, x1_2, x2_2, label="Atrator $id2", markersize=3, alpha=0.6)
    
    p2 = plot(x1_1, label="Atrator $id1", title="Série Temporal x1", xlabel="Período", ylabel="x1")
    plot!(p2, x1_2, label="Atrator $id2", linestyle=:dash)
    
    plt = plot(p1, p2, layout=(1,2), size=(1000, 400))
    savefig(plt, "comparacao_attr$(id1)_vs_attr$(id2).png")
    display(plt)
    
    return plt
end

function run_diagnostic(divs=200)
    region = BasinRegion(
        [[-1.2, 1.5], [-1.5, 1.5]], 
        [divs, floor(Int64,(3.0/2.7)*divs)],
        [false, false],
        [[-2.0, 2.5],[-2.0, 2.0]],
    )

    bp = BasinProblem(
        helmholtz_duffing!, 
        [0.1, -1.2, -0.3, 2.0, 0.077, 1.17],
        2*pi/1.17, 
        10, 
        1000,
        20, 
        80, 
        region,
        1,
    )

    println(">>> Construindo mapa...")
    scm = build_scmap(bp)
    
    println(">>> Encontrando atratores...")
    result = find_attractors_from_scm(scm, bp)
    
    report_attractors(result)
    centroids = distance_matrix(result, bp)
    
    return result, bp
end

# 1. Rodar o diagnóstico:
    result, bp = run_diagnostic(200)
#
# 2. Comparção de dois atratores:
    compare_attractors_timeseries(result, bp, 2, 4) 
    compare_attractors_timeseries(result, bp, 2, 5)  # Amarelo vs Ciano
# =============================================================================

println("Script carregado!")
println("Execute: result, bp = run_diagnostic()")