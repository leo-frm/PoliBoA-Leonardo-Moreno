# Caso B do pêndulo com o PoliBoA original (incluído sem alteração), 200×200, 1 thread
include(raw"C:\Users\Leonardo\Desktop\PoliBoA - Leonardo Moreno\ZZ - Poli BoA og\poliboa.jl")

function pendulum!(dx, x, (c, F, Ω), t)
    dx[1] = x[2]
    dx[2] = -c * x[2] - sin(x[1]) + F * cos(Ω * t)
end

function pendulum_bp(c, F, Ω; divs=200)
    region = BasinRegion([[-pi, pi], [-3.0, 3.0]], [divs, divs], [true, false], [[-pi, pi], [-5.0, 5.0]])
    return BasinProblem(pendulum!, [c, F, Ω], 2*pi/Ω, 10, 1000, 20, 80, region, 1)
end

populate_basins(pendulum_bp(0.5, 1.07, 2/3; divs=10))   # pré-compilação
t = @elapsed r = populate_basins(pendulum_bp(0.5, 1.07, 2/3))
cells = vec(r.cells); N = length(cells)
@printf("Original — pêndulo c=0.5 F=1.07 Ω=2/3, 200x200: %d atratores, %.1f s\n", length(r.attractors), t)
for id in sort(unique(cells))
    @printf("    bacia %3d: %6d células (%.4f da área)\n", id, count(==(id), cells), count(==(id), cells) / N)
end
for a in r.attractors
    @printf("    #%d %-15s %2d células  pontos: %s\n", a.number, string(a.kind), length(a.points),
        join([@sprintf("(%.3f, %.3f)", p[1], p[2]) for p in a.points[1:min(end, 4)]], " "))
end
