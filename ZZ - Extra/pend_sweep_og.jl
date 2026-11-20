# Varredura com o PoliBoA original (ZZ - Poli BoA og/poliboa.jl, incluído sem alteração), grid 200×200
include(raw"C:\Users\Leonardo\Desktop\PoliBoA - Leonardo Moreno\ZZ - Poli BoA og\poliboa.jl")

const OUTDIR = @__DIR__

# Pêndulo amortecido forçado (servomecanismo do python-control, adimensionalizado):
# θ'' + c θ' + sin θ = F cos(Ω t)
function pendulum!(dx, x, (c, F, Ω), t)
    dx[1] = x[2]
    dx[2] = -c * x[2] - sin(x[1]) + F * cos(Ω * t)
end

function pendulum_bp(c, F, Ω; divs=200, x2max=3.0, x2ext=5.0, transient=10)
    region = BasinRegion(
        [[-pi, pi], [-x2max, x2max]],
        [divs, divs],
        [true, false],
        [[-pi, pi], [-x2ext, x2ext]],
    )
    return BasinProblem(pendulum!, [c, F, Ω], 2*pi/Ω, transient, 1000, 20, 80, region, 1)
end

function report(tag, result, t)
    cells = vec(result.cells); N = length(cells)
    @printf("\n%s: %d atratores, %.1f s\n", tag, length(result.attractors), t)
    for id in sort(unique(cells))
        @printf("    bacia %3d: %.4f da área\n", id, count(==(id), cells) / N)
    end
    for a in result.attractors
        @printf("    #%d %-15s %2d células  pontos: %s\n", a.number, string(a.kind), length(a.points),
            join([@sprintf("(%.3f, %.3f)", p[1], p[2]) for p in a.points[1:min(end, 4)]], " "))
    end
end

cases = [
    ("A_c0.10_F2.00_W1.000",  (0.1, 2.0, 1.0)),
    ("B_c0.50_F1.07_W0.667",  (0.5, 1.07, 2/3)),
    ("C_c0.10_F1.50_W1.000",  (0.1, 1.5, 1.0)),
]

populate_basins(pendulum_bp(0.1, 2.0, 1.0; divs=10))   # pré-compilação

tasks = map(cases) do (tag, (c, F, Ω))
    Threads.@spawn begin
        local r
        t = @elapsed r = populate_basins(pendulum_bp(c, F, Ω))
        (tag, r, t)
    end
end
for task in tasks
    tag, r, t = fetch(task)
    report(tag, r, t)
    savefig(plotbasins(r), joinpath(OUTDIR, "sweep_og_$(tag).png"))
end
