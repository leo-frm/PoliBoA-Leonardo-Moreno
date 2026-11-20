# Diagnóstico: atratores verdadeiros do pêndulo (c=0.5, F=1.07, Ω=2/3) por integração
# de alta precisão (Vern9, tol 1e-12), independente de células. Só diagnóstico.
using DifferentialEquations, Printf, Random

function pendulum!(dx, x, (c, F, Ω), t)
    dx[1] = x[2]
    dx[2] = -c * x[2] - sin(x[1]) + F * cos(Ω * t)
end
const P = (0.5, 1.07, 2/3); const T = 2pi / P[3]
wrap(θ) = mod(θ + pi, 2pi) - pi

function strobe(u0, n_trans, n_keep; alg=Vern9(), tol=1e-12, default=false)
    prob = ODEProblem(pendulum!, copy(u0), (0.0, T * (n_trans + n_keep)), P)
    integ = default ? init(prob; dense=false, save_everystep=false, save_start=false, maxiters=1e10) :
                      init(prob, alg; reltol=tol, abstol=tol, dense=false, save_everystep=false, save_start=false, maxiters=1e10)
    step!(integ, T * n_trans, true)
    pts = Vector{Vector{Float64}}()
    for k in 1:n_keep
        step!(integ, T, true)
        push!(pts, [wrap(integ.u[1]), integ.u[2]])
    end
    return pts
end
dist(a, b) = hypot(wrap(a[1] - b[1]), a[2] - b[2])
function period_of(pts; tol=1e-6)
    n = length(pts)
    for k in 1:div(n, 2)
        all(dist(pts[n - j], pts[n - j - k]) < tol for j in 0:(n - 2k)) && return k
    end
    return 0
end

println("== A partir dos pontos dos atratores do original (pend_og_m1.log) ==")
for (nome, u0) in [("#1", [-2.466, 1.155]), ("#2", [-0.141, 1.935])]
    for (lab, kw) in [("Vern9 1e-12", (;)), ("integrador do original (padrão)", (; default=true))]
        pts = strobe(u0, 400, 16; kw...)
        k = period_of(pts; tol = haskey(kw, :default) ? 1e-3 : 1e-7)
        @printf("  %s, %-32s período %d; pontos: %s\n", nome, lab, k,
            join([@sprintf("(%.4f, %.4f)", p[1], p[2]) for p in pts[end-max(k,1)+1:end]], " "))
    end
end

println("\n== 1000 condições iniciais aleatórias em [-π,π]×[-3,3], Vern9 1e-12, 400 períodos ==")
rng = Xoshiro(2026)
reps = Vector{Vector{Vector{Float64}}}(); counts = Int[]; per = Int[]
for i in 1:1000
    u0 = [rand(rng) * 2pi - pi, rand(rng) * 6 - 3]
    pts = strobe(u0, 400, 16)
    k = period_of(pts; tol=1e-7)
    j = findfirst(r -> minimum(dist(pts[end], q) for q in r) < 1e-5, reps)
    if j === nothing
        push!(reps, pts[end-max(k,1)+1:end]); push!(counts, 1); push!(per, k)
    else
        counts[j] += 1
    end
end
for (j, r) in enumerate(reps)
    @printf("  atrator %d: período %d, %4d/1000 CIs; pontos: %s\n", j, per[j], counts[j],
        join([@sprintf("(%.4f, %.4f)", p[1], p[2]) for p in r], " "))
end
