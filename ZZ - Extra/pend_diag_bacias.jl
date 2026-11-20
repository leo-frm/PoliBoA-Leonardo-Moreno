# Diagnóstico (fora do repositório, nada do repositório é alterado): bacia de referência
# do pêndulo por integração de alta precisão a partir do centro de cada célula, e
# concordância do original, do SCM e do GCM com essa referência. 200×200.
using Accessors, DifferentialEquations, ProgressMeter, Printf
const REPO = raw"C:\Users\Leonardo\Desktop\PoliBoA - Leonardo Moreno"
include(joinpath(REPO, "desenvolvimento_scm", "structures.jl"))
include(joinpath(REPO, "desenvolvimento_scm", "utils.jl"))
include(joinpath(REPO, "desenvolvimento_scm", "milestone2.jl"))
include(joinpath(REPO, "desenvolvimento_scm", "milestone4.jl"))
include(joinpath(REPO, "desenvolvimento_scm", "milestone5.jl"))
include(joinpath(REPO, "desenvolvimento_scm", "poliboa_original.jl"))
include(joinpath(REPO, "desenvolvimento_gcm", "structures_gcm.jl"))
include(joinpath(REPO, "desenvolvimento_gcm", "utils_gcm.jl"))
include(joinpath(REPO, "desenvolvimento_gcm", "milestone1_gcm.jl"))
include(joinpath(REPO, "desenvolvimento_gcm", "milestone2_gcm.jl"))
include(joinpath(REPO, "desenvolvimento_gcm", "milestone4_gcm.jl"))

function pendulum!(dx, x, (c, F, Ω), t)
    dx[1] = x[2]
    dx[2] = -c * x[2] - sin(x[1]) + F * cos(Ω * t)
end
# Mesma configuração de make_pendulum_bp (tests_milestone5.jl)
function make_pendulum_bp(divs; c=0.5, F=1.07, Ω=2/3, threads=1)
    region = BasinRegion([[-pi, pi], [-3.0, 3.0]], [divs, divs], [true, false], [[-pi, pi], [-5.0, 5.0]])
    return BasinProblem(pendulum!, [c, F, Ω], 2*pi/Ω, 10, 1000, 20, 80, region, threads)
end
const PAR = (0.5, 1.07, 2/3); const T = 2pi / PAR[3]
wrap(θ) = mod(θ + pi, 2pi) - pi
dist(a, b) = hypot(wrap(a[1] - b[1]), a[2] - b[2])
hp_integrator(u0) = init(ODEProblem(pendulum!, copy(u0), (0.0, T * 2000), PAR), Vern9();
    reltol=1e-10, abstol=1e-10, dense=false, save_everystep=false, save_start=false, maxiters=1e10)

# Atratores verdadeiros (período 2), a partir dos pontos do original, 400 períodos
true_attr = map([[-2.466, 1.155], [-0.141, 1.935]]) do u0
    integ = hp_integrator(u0); step!(integ, 400T, true)
    pts = Vector{Vector{Float64}}()
    for k in 1:2; step!(integ, T, true); push!(pts, [wrap(integ.u[1]), integ.u[2]]); end
    pts
end
println("Atratores de referência (Vern9 1e-10, período 2):")
for (j, a) in enumerate(true_attr)
    @printf("  R%d: %s\n", j, join([@sprintf("(%.4f, %.4f)", p[1], p[2]) for p in a], " "))
end

# Órbita vista pelo integrador do original (padrão, reltol 1e-3), últimos 12 pontos
println("\nÓrbita estroboscópica com o integrador do original, a partir de R1, períodos 389–400:")
let integ = init(ODEProblem(pendulum!, copy(true_attr[1][1]), (0.0, T * 2000), PAR);
                 dense=false, save_everystep=false, save_start=false, maxiters=1e10)
    step!(integ, 388T, true)
    for k in 1:12
        step!(integ, T, true); @printf("  (%.4f, %.4f)", wrap(integ.u[1]), integ.u[2]); k % 4 == 0 && println()
    end
end

# Bacia de referência: centro de cada célula, até ficar a < 1e-3 de um ponto de R1 ou R2
bp = make_pendulum_bp(200)
N = prod(bp.region.elements)
ref = zeros(Int64, N); nper = zeros(Int64, N)
t_ref = @elapsed Threads.@threads for c in 1:N
    u = zeros(2); store_cell_center!(u, c, bp.region)
    integ = hp_integrator(u)
    for k in 1:2000
        step!(integ, T, true)
        x = [wrap(integ.u[1]), integ.u[2]]
        j = findfirst(a -> any(dist(x, p) < 1e-3 for p in a), true_attr)
        if j !== nothing && k >= 10
            ref[c] = j; nper[c] = k; break
        end
    end
end
@printf("\nReferência: %.1f s; R1 %d células (%.4f), R2 %d (%.4f), sem convergir %d; períodos até convergir: mediana %d, máx %d\n",
    t_ref, count(==(1), ref), count(==(1), ref)/N, count(==(2), ref), count(==(2), ref)/N, count(==(0), ref),
    sort(nper[nper .> 0])[div(count(>(0), nper), 2)], maximum(nper))

# Células de fronteira na referência: algum vizinho 3×3 (x1 cíclico) com rótulo diferente
nx, ny = bp.region.elements
R = reshape(ref, nx, ny)
boundary = falses(nx, ny)
for i in 1:nx, j in 1:ny, di in -1:1, dj in -1:1
    jj = j + dj; (1 <= jj <= ny) || continue
    ii = mod1(i + di, nx)
    R[ii, jj] != R[i, j] && (boundary[i, j] = true)
end
boundary = vec(boundary)
@printf("Células de fronteira na referência (vizinhança 3×3): %d (%.2f%%)\n", count(boundary), 100count(boundary)/N)

# Mapeamento bijetivo guloso (mesma regra de compare_basins) de rótulos do método → referência
function best_map(cells, ref)
    pc = Dict{Tuple{Int,Int},Int}()
    for i in eachindex(cells); k = (cells[i], ref[i]); pc[k] = get(pc, k, 0) + 1; end
    cand = sort([(a, b, n) for ((a, b), n) in pc if a != -1 && b > 0], by = x -> -x[3])
    m = Dict{Int,Int}(); ua = Set{Int}(); ub = Set{Int}()
    for (a, b, _) in cand
        (a ∉ ua && b ∉ ub) && (m[a] = b; push!(ua, a); push!(ub, b))
    end
    return m
end
regioes = [("total", trues(N)), ("interior", .!boundary), ("fronteira", boundary)]
function report(nome, result)
    cells = vec(result.cells); m = best_map(cells, ref)
    hit = [get(m, cells[i], 0) == ref[i] for i in 1:N]
    fr = [count(==(a), cells)/N for a in sort(unique(cells))]
    @printf("  %-34s %d atr. | bacias %s | × ref: %s\n", nome, length(result.attractors),
        join([@sprintf("%.3f", f) for f in fr], "/"),
        join([@sprintf("%s %.2f%%", r, 100count(hit[s])/count(s)) for (r, s) in regioes], ", "))
    return cells
end

println("\nConcordância com a referência (compare_basins = mesma regra de mapeamento):")
res = Dict{String,Any}()
for nt in (1, 4)
    r = populate_basins(make_pendulum_bp(200; threads=nt)); res["og$nt"] = r
    report("Original, $nt thread(s)", r)
end
for nt in (1, 4)
    b = make_pendulum_bp(200; threads=nt)
    r = find_attractors_from_scm(build_scmap_parallel(b), b); res["scm$nt"] = r
    report("SCM, $nt thread(s)", r)
end
let b = make_pendulum_bp(200; threads=4)
    gp = GcmProblem(b, [10, 10], 10; sampling=:random)
    r = find_attractors_from_gcm(build_gcmap_parallel(gp), gp); res["gcm"] = r
    report("GCM, 4 threads", r)
end

# Variante só para diagnóstico (NÃO aplicada ao repositório): SCM com o integrador
# reiniciado por completo a cada célula, como map_point do GCM (reset_dt=true)
function compute_single_cell_mapping(cell_id::Int64, integrator, bp::BasinProblem) :: Int64
    u = zeros(Float64, length(bp.region.elements))
    store_cell_center!(u, cell_id, bp.region)
    return map_point(u, integrator, bp)
end
for nt in (1, 4)
    b = make_pendulum_bp(200; threads=nt)
    r = find_attractors_from_scm(build_scmap_parallel(b), b); res["scmR$nt"] = r
    report("SCM + reset_dt (variante), $nt thr", r)
end

println("\nReprodutibilidade (compare_basins entre execuções):")
for (a, b) in [("og1", "og4"), ("scm1", "scm4"), ("scmR1", "scmR4"), ("og1", "scm1"), ("og1", "gcm"), ("scm1", "gcm")]
    @printf("  %-6s × %-6s  %.4f\n", a, b, compare_basins(res[a], res[b]))
end
