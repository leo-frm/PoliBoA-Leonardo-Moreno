# =============================================================================
# Diagnóstico do PoliBoA Original (populate_basins / search_for_attractor!)
# Uso: julia diagnostic_og.jl          (roda HD + SC)
# =============================================================================

include("poliboa.jl")

# ─── Relatório de atratores ──────────────────────────────────────────────────

function report_attractors(result::BasinResult)
    colors_names = ["Black (Divergent)", "White", "Red", "Yellow", "Blue", 
                    "Green", "Cyan", "Purple", "Brown", "Pink",
                    "Aquamarine", "Goldenrod", "Ivory", "Salmon",
                    "MidnightBlue", "Magenta", "RosyBrown", "Teal"]
    
    println("Detalhamento dos atratores")
    println("=" ^ 60)
    
    for attr in result.attractors
        color_idx = attr.number + 2
        color_name = color_idx <= length(colors_names) ? colors_names[color_idx] : "Custom"
        num_cells = count(==(attr.number), result.cells)
        
        println("Atrator #$(attr.number) ($(color_name))")
        println("-" ^ 40)
        @printf("  Tipo: %s\n", attr.kind)
        @printf("  Células na bacia: %d\n", num_cells)
        @printf("  Pontos no ciclo: %d\n", length(attr.points))
        println("  Pontos do atrator:")
        for (i, pt) in enumerate(attr.points)
            @printf("    [%d] x1 = %+.6f, x2 = %+.6f\n", i, pt[1], pt[2])
        end
        println("-" ^ 40)
    end

    # Resumo
    total = length(result.cells)
    n_divergent = count(==(-1), result.cells)
    n_classified = count(>(0), result.cells)
    println("\nResumo:")
    @printf("  Total de células:    %d\n", total)
    @printf("  Classificadas:       %d (%.1f%%)\n", n_classified, 100*n_classified/total)
    @printf("  Divergentes (preto): %d (%.1f%%)\n", n_divergent, 100*n_divergent/total)
    @printf("  Atratores únicos:    %d\n", length(result.attractors))
    println("=" ^ 60)
end

# ─── Matriz de distâncias entre centróides ───────────────────────────────────

function distance_matrix(result::BasinResult, f, params, period)
    ids = [attr.number for attr in result.attractors]
    n = length(ids)
    
    centroids = Dict{Int, Vector{Float64}}()
    
    for attr in result.attractors
        u0 = attr.points[1]
        
        prob = ODEProblem(f, u0, (0.0, 100 * period), params)
        sol = solve(prob, Tsit5(), saveat=period)
        
        last_points = sol.u[max(1, length(sol.u)-50):end]
        cx = sum(p[1] for p in last_points) / length(last_points)
        cy = sum(p[2] for p in last_points) / length(last_points)
        centroids[attr.number] = [cx, cy]
    end
    
    println("\nMatriz de distâncias entre centróides:")
    println("-" ^ (9 + 9 * n))
    
    print("       ")
    for id in ids
        @printf("  Attr%-3d", id)
    end
    println()
    
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

# ─── Comparação de séries temporais ──────────────────────────────────────────

function compare_attractors_timeseries(result::BasinResult, f, params, period, id1::Int, id2::Int)
    attr1 = findfirst(a -> a.number == id1, result.attractors)
    attr2 = findfirst(a -> a.number == id2, result.attractors)
    
    if attr1 === nothing || attr2 === nothing
        println("Atrator alvo não encontrado!")
        return nothing
    end
    
    a1 = result.attractors[attr1]
    a2 = result.attractors[attr2]
    u0_1 = a1.points[1]
    u0_2 = a2.points[1]
    
    println("\nComparando atratores $id1 vs $id2")
    println("   CI atrator $id1: ($(round(u0_1[1], digits=4)), $(round(u0_1[2], digits=4)))")
    println("   CI atrator $id2: ($(round(u0_2[1], digits=4)), $(round(u0_2[2], digits=4)))")
    
    tempo = 100 * period
    sol1 = solve(ODEProblem(f, u0_1, (0.0, tempo), params), Tsit5(), saveat=period)
    sol2 = solve(ODEProblem(f, u0_2, (0.0, tempo), params), Tsit5(), saveat=period)
    
    x1_1 = [sol1.u[i][1] for i in 20:length(sol1.u)]
    x2_1 = [sol1.u[i][2] for i in 20:length(sol1.u)]
    x1_2 = [sol2.u[i][1] for i in 20:length(sol2.u)]
    x2_2 = [sol2.u[i][2] for i in 20:length(sol2.u)]
    
    m1 = (sum(x1_1)/length(x1_1), sum(x2_1)/length(x2_1))
    m2 = (sum(x1_2)/length(x1_2), sum(x2_2)/length(x2_2))
    dist = sqrt((m1[1]-m2[1])^2 + (m1[2]-m2[2])^2)
    
    println("  Média atrator $id1: ($(round(m1[1], digits=4)), $(round(m1[2], digits=4)))")
    println("  Média atrator $id2: ($(round(m2[1], digits=4)), $(round(m2[2], digits=4)))")
    println("  Distância: $(round(dist, digits=6))")
    
    if dist < 0.01
        println("  → Atrator duplicado")
    elseif dist < 0.1
        println("  → Possivelmente o mesmo")
    else
        println("  → Atratores diferentes")
    end
    
    p1 = scatter(x1_1, x2_1, label="Atrator $id1", markersize=3, alpha=0.6,
                 title="Espaço de Fases", xlabel="x1", ylabel="x2")
    scatter!(p1, x1_2, x2_2, label="Atrator $id2", markersize=3, alpha=0.6)
    
    p2 = plot(x1_1, label="Atrator $id1", title="Série Temporal x1", xlabel="Período", ylabel="x1")
    plot!(p2, x1_2, label="Atrator $id2", linestyle=:dash)
    
    plt = plot(p1, p2, layout=(1,2), size=(1000, 400))
    savefig(plt, "og_comparacao_attr$(id1)_vs_attr$(id2).png")
    display(plt)
    
    return plt
end

# ─── Diagnóstico: Helmholtz-Duffing ──────────────────────────────────────────

function run_diagnostic_hd(divs=200; e=0.077)
    params_hd = [0.1, -1.2, -0.3, 2.0, e, 1.17]
    period_hd = 2*pi/1.17

    region = BasinRegion(
        [[-1.2, 1.5], [-1.5, 1.5]], 
        [divs, floor(Int64, (3.0/2.7)*divs)],
        [false, false],
        [[-2.0, 2.5], [-2.0, 2.0]],
    )

    bp = BasinProblem(
        helmholtz_duffing!, 
        params_hd,
        period_hd, 
        10, 1000, 20, 80, 
        region, Threads.nthreads(),
    )

    println("\n" * "=" ^ 60)
    println("  DIAGNÓSTICO ORIGINAL — HELMHOLTZ-DUFFING (e=$e, divs=$divs)")
    println("=" ^ 60)

    println(">>> [HD-OG] Rodando populate_basins...")
    @time result = populate_basins(bp)

    report_attractors(result)
    centroids = distance_matrix(result, helmholtz_duffing!, params_hd, period_hd)

    plt = plotbasins(result)
    savefig(plt, "og_basin_hd.png")
    display(plt)

    return result, params_hd, period_hd
end

# ─── Diagnóstico: Spring Column ──────────────────────────────────────────────

function run_diagnostic_sc(divs=1000; omega=0.0)
    params_sc = [0.01, 0.05, 0.8, 0.01, 0.05, omega]
    period_sc = 7.853981634

    region = BasinRegion(
        [[-3.14, 3.14], [-1.3, 1.3]], 
        [divs, divs],
        [true, false],
        [[-3.14, 3.14], [-2.0, 2.0]],
    )

    bp = BasinProblem(
        spring_column!, 
        params_sc,
        period_sc, 
        10, 1000, 20, 80, 
        region, Threads.nthreads(),
    )

    println("\n" * "=" ^ 60)
    println("  DIAGNÓSTICO ORIGINAL — SPRING COLUMN (ω=$omega, divs=$divs)")
    println("=" ^ 60)

    println(">>> [SC-OG] Rodando populate_basins...")
    @time result = populate_basins(bp)

    report_attractors(result)
    centroids = distance_matrix(result, spring_column!, params_sc, period_sc)

    plt = plotbasins(result)
    savefig(plt, "og_basin_sc.png")
    display(plt)

    return result, params_sc, period_sc
end

# =============================================================================
# EXECUÇÃO
# =============================================================================

# ── Helmholtz-Duffing ──
result_hd, params_hd, period_hd = run_diagnostic_hd(200)

ids_hd = [attr.number for attr in result_hd.attractors]
for i in 1:length(ids_hd)
    for j in i+1:length(ids_hd)
        compare_attractors_timeseries(result_hd, helmholtz_duffing!, params_hd, period_hd, ids_hd[i], ids_hd[j])
    end
end

# ── Spring Column ──
result_sc, params_sc, period_sc = run_diagnostic_sc(1000)

ids_sc = [attr.number for attr in result_sc.attractors]
for i in 1:length(ids_sc)
    for j in i+1:length(ids_sc)
        compare_attractors_timeseries(result_sc, spring_column!, params_sc, period_sc, ids_sc[i], ids_sc[j])
    end
end

# ── Resumo final ──
println("\n" * "▓" ^ 60)
println("▓  RESUMO FINAL — CÓDIGO ORIGINAL")
println("▓" ^ 60)
@printf("  [HD] Atratores: %d | Divergentes: %d (%.1f%%)\n", 
    length(result_hd.attractors), 
    count(==(-1), result_hd.cells),
    100*count(==(-1), result_hd.cells)/length(result_hd.cells))
@printf("  [SC] Atratores: %d | Divergentes: %d (%.1f%%)\n", 
    length(result_sc.attractors), 
    count(==(-1), result_sc.cells),
    100*count(==(-1), result_sc.cells)/length(result_sc.cells))
println("▓" ^ 60)