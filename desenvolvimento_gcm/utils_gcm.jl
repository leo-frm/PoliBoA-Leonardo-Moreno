# Funções auxiliares do GCM
# Depende de utils.jl (get_cell_number, store_cell_center!, is_inside_range,
# set_integrator!, adjust_cyclic)

using Random

# Número de amostras integradas por célula (m)
samples_per_cell(gp::GcmProblem) = gp.samples_per_cell

# Converte o número absoluto da célula no índice local dentro de cell_range
local_index(gcm::GeneralizedCellMap, cell::Int64) = cell - first(gcm.cell_range) + 1

function n_images(gcm::GeneralizedCellMap, cell::Int64) :: Int64 #I(z): número de imagens distintas
    k = local_index(gcm, cell)
    return Int64(gcm.offsets[k + 1] - gcm.offsets[k])
end

function images_of(gcm::GeneralizedCellMap, cell::Int64) #C(z,·): células-imagem, sem cópia
    k = local_index(gcm, cell)
    return view(gcm.images, gcm.offsets[k]:(gcm.offsets[k + 1] - 1))
end

function counts_of(gcm::GeneralizedCellMap, cell::Int64) #Contagens de amostras por imagem
    k = local_index(gcm, cell)
    return view(gcm.counts, gcm.offsets[k]:(gcm.offsets[k + 1] - 1))
end

function probabilities_of(gcm::GeneralizedCellMap, cell::Int64) :: Vector{Float64} #P(z,·); soma 1 por construção
    return counts_of(gcm, cell) ./ gcm.samples_per_cell
end

computed_cells(gcm::GeneralizedCellMap) = count(gcm.computed)

function memory_footprint(gcm::GeneralizedCellMap) :: Int64 #Bytes ocupados pelos arranjos do mapa
    return sizeof(gcm.offsets) + sizeof(gcm.images) + sizeof(gcm.counts) +
           cld(length(gcm.computed), 8)
end

# Escolhe quais subcélulas da célula cell_id serão integradas, sem reposição.
# As subcélulas são numeradas de 1 a K = prod(subgrid).
#   :random     → sorteio simples: as m primeiras de uma permutação aleatória
#   :stratified → divide 1..K em m faixas contíguas e sorteia uma em cada; como
#                 a numeração percorre a célula em ordem, as faixas são regiões
#                 vizinhas, e as amostras ficam espalhadas pela célula
# A semente depende da célula, então o sorteio é reprodutível e não depende da
# ordem de processamento nem do número de threads. Usar m maior com a mesma
# semente acrescenta subcélulas novas às já sorteadas (sem repetir), que é o
# processo de rodadas proposto: m = 10 e depois m = 20 equivale a duas rodadas de 10.
function select_samples(gp::GcmProblem, cell_id::Int64) :: Vector{Int64}
    K = prod(gp.subgrid)
    m = gp.samples_per_cell
    rng = Xoshiro(gp.seed * 1_000_000_000 + cell_id)

    if gp.sampling == :stratified
        return [rand(rng, (div((j - 1) * K, m) + 1):div(j * K, m)) for j in 1:m]
    else
        return randperm(rng, K)[1:m]
    end
end

# Escreve em u as coordenadas do centro da sample_id-ésima subcélula da célula
# cell_id. A célula é dividida em subgrid[i] partes na dimensão i (Hsu §11.12 e
# §12.1.3: amostras no interior da célula, nos centros das subcélulas). Quando o
# grid fino de referência é o grid grosso multiplicado pelo subgrid, cada amostra
# coincide com o centro de uma célula fina. Nenhuma amostra cai na fronteira
# entre células vizinhas.
function store_cell_sample!(
    u :: Vector{Float64},
    cell_id :: Int64,
    region :: BasinRegion,
    sample_id :: Int64,
    subgrid :: Vector{Int64},
)
    store_cell_center!(u, cell_id, region) # Parte do centro e desloca

    rest = sample_id - 1
    for i in eachindex(region.range)
        j = rest % subgrid[i] + 1 # Índice da subcélula na dimensão i
        rest = div(rest, subgrid[i])
        cell_length = (region.range[i][2] - region.range[i][1]) / region.elements[i]
        u[i] += ((j - 0.5) / subgrid[i] - 0.5) * cell_length
    end

    return u
end

# Integra uma condição inicial qualquer e devolve a célula de destino, ou -1 se
# divergir. Mesma lógica de compute_single_cell_mapping, generalizada para
# aceitar um ponto inicial arbitrário em vez do centro da célula; de fato
# compute_single_cell_mapping(c, ...) equivale a map_point(centro_de_c, ...).
# u0 não é modificado.
# O integrador é reiniciado por completo (reinit! com reset_dt) a cada ponto: o
# passo adaptativo não é herdado da integração anterior, então o resultado de um
# ponto não depende da ordem de processamento nem do número de threads.
function map_point(u0::Vector{Float64}, integrator, bp::BasinProblem) :: Int64
    reinit!(integrator, copy(u0); t0=0.0, reset_dt=true)
    step!(integrator, bp.period * bp.transient_cycles, true)

    set_integrator!(integrator; u=integrator.u, t=0.0)
    step!(integrator, bp.period, true)

    adjust_cyclic(integrator.u, bp.region.range, bp.region.is_cyclic)

    # Se o ponto final saiu da região estendida: divergência
    if !is_inside_range(integrator.u, bp.region.extended_range)
        return -1
    end

    # Se está dentro da estendida mas fora da região principal:
    # reintegra até maximum_extended_cycles tentando retornar ao grid
    if !is_inside_range(integrator.u, bp.region.range)
        for ext_cycle in 1:bp.maximum_extended_cycles
            set_integrator!(integrator; u=integrator.u, t=0.0)
            step!(integrator, bp.period, true)
            adjust_cyclic(integrator.u, bp.region.range, bp.region.is_cyclic)

            if !is_inside_range(integrator.u, bp.region.extended_range)
                return -1
            end

            if is_inside_range(integrator.u, bp.region.range)
                return get_cell_number(integrator.u, bp.region)
            end
        end
        return -1
    end

    return get_cell_number(integrator.u, bp.region)
end

# Instancia um integrador. Cada thread precisa do seu próprio, pois os
# integradores do DifferentialEquations.jl não são thread-safe.
function make_gcm_integrator(bp::BasinProblem)
    ode_problem = ODEProblem(
        bp.f,
        zeros(Float64, length(bp.region.elements)),
        (0.0, bp.period * bp.maximum_cycles),
        bp.params,
    )

    return init(
        ode_problem;
        dense=false,
        save_everystep=false,
        save_start=false,
        maxiters=1e10,
    )
end

# Imprime a distribuição de I(z), o número de imagens por célula.
# Células com I(z) = 1 comportam-se deterministicamente, como no SCM;
# células com I(z) > 1 são as de fronteira entre bacias.
function report_image_distribution(gcm::GeneralizedCellMap)
    done = computed_cells(gcm)
    if done == 0
        println("Mapa vazio.")
        return nothing
    end

    hist = Dict{Int64, Int64}()
    total_images = 0
    for cell in gcm.cell_range
        ni = n_images(gcm, cell)
        hist[ni] = get(hist, ni, 0) + 1
        total_images += ni
    end

    w = 52
    println("\n" * "=" ^ w)
    println("  DISTRIBUIÇÃO DE I(z) — imagens por célula")
    println("=" ^ w)
    @printf("  %-10s %12s %12s\n", "I(z)", "Células", "Fração")
    println("-" ^ w)
    for ni in sort(collect(keys(hist)))
        @printf("  %-10d %12d %11.2f%%\n", ni, hist[ni], 100.0 * hist[ni] / done)
    end
    println("-" ^ w)
    @printf("  Média de I(z)      : %.4f\n", total_images / done)
    @printf("  Máximo de I(z)     : %d  (limite teórico: %d)\n",
        maximum(keys(hist)), gcm.samples_per_cell)
    @printf("  Determinísticas    : %.2f%% das células (I(z) = 1)\n",
        100.0 * get(hist, 1, 0) / done)
    @printf("  Memória do mapa    : %.3f MB\n", memory_footprint(gcm) / 1024^2)
    println("=" ^ w * "\n")

    return nothing
end

# Converte uma imagem do mapa no índice de nó do grafo de transições: o
# sumidouro (-1) passa a ser o nó N+1, convenção de Hsu (§11.1) para a
# implementação computacional, em que o rótulo 0 não está disponível
gcm_node(img::Integer, N::Int64) :: Int64 = img < 0 ? N + 1 : Int64(img)

# Pré-imagens do nó z (R(z,·) de Hsu §11.2), sem cópia
function preimages_of(pre::PreImageArray, z::Int64)
    return view(pre.cells, pre.offsets[z]:(pre.offsets[z + 1] - 1))
end

# Imprime o resumo da análise do GCM: soluções periódicas do SCM dominante,
# grupos persistentes encontrados e fração de células de fronteira
function report_gcm_analysis(analysis::GcmAnalysis)
    N = length(analysis.group_of) - 1
    n_persistent = count(>(0), view(analysis.group_of, 1:N))
    n_single = count(==(1), view(analysis.domicile_count, 1:N))
    n_multiple = count(>(1), view(analysis.domicile_count, 1:N))

    w = 60
    println("\n" * "=" ^ w)
    println("  ANÁLISE DO GCM (Hsu, Cap. 11)")
    println("=" ^ w)
    @printf("  Soluções periódicas do SCM dominante : %d\n", analysis.n_periodic_solutions)
    @printf("  Grupos persistentes (com sumidouro)  : %d\n", length(analysis.groups))
    println("-" ^ w)
    for (g, members) in enumerate(analysis.groups)
        label = g == 1 ? "Sumidouro (divergência)" : "Grupo $g → Atrator #$(g - 1)"
        @printf("  %-28s %6d células, período %d\n", label, length(members), analysis.periods[g])
    end
    println("-" ^ w)
    @printf("  Células persistentes     : %6d (%.2f%%)\n", n_persistent, 100.0 * n_persistent / N)
    @printf("  Domicílio único          : %6d (%.2f%%)\n", n_single, 100.0 * n_single / N)
    @printf("  Domicílio múltiplo       : %6d (%.2f%%)  ← fronteira entre bacias\n",
        n_multiple, 100.0 * n_multiple / N)
    println("=" ^ w * "\n")

    return nothing
end