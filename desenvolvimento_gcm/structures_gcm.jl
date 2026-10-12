# Estruturas do Generalized Cell Mapping
# Reaproveita BasinRegion, BasinProblem, Attractor e BasinResult de structures.jl

mutable struct GcmProblem # Envolve um BasinProblem acrescentando os parâmetros de amostragem
    base :: BasinProblem #Problema base, idêntico ao usado no SCM
    subgrid :: Vector{Int64} #Subcélulas por dimensão dentro de cada célula (ex.: [10, 20] → 200 subcélulas)
    samples_per_cell :: Int64 #m: subcélulas sorteadas por célula, sem reposição
    seed :: Int64 #Semente do sorteio; a mesma semente reproduz exatamente o mesmo mapa
    sampling :: Symbol #:random (sorteio simples) ou :stratified (um sorteio em cada uma de m faixas)
end

# Construtor com valores padrão e verificações. Usar m igual ao total de
# subcélulas integra todas elas, como na grade fixa da versão anterior.
function GcmProblem(
    base :: BasinProblem,
    subgrid :: Vector{Int64},
    samples_per_cell :: Int64;
    seed :: Int64 = 1,
    sampling :: Symbol = :random,
)
    @assert length(subgrid) == length(base.region.elements) "subgrid precisa ter uma entrada por dimensão"
    @assert 1 <= samples_per_cell <= prod(subgrid) "m precisa estar entre 1 e o número de subcélulas ($(prod(subgrid)))"
    @assert sampling in (:random, :stratified) "sampling precisa ser :random ou :stratified"

    return GcmProblem(base, subgrid, samples_per_cell, seed, sampling)
end

# Mapa célula-a-célula generalizado, no formato dos três arranjos de Hsu (§11.1):
# I(z) número de imagens, C(z,i) i-ésima imagem, P(z,i) probabilidade de transição.
# Armazenado em CSR: as imagens da k-ésima célula da faixa ocupam
# images[offsets[k] : offsets[k+1]-1]. Guardam-se contagens inteiras em vez de
# probabilidades, já que o denominador (m amostras por célula) é constante.
mutable struct GeneralizedCellMap
    cell_range :: UnitRange{Int64} #Faixa de células coberta; 1:N no mapa completo, sub-faixa nos parciais
    offsets :: Vector{Int32} #Índices de início da lista de imagens de cada célula
    images :: Vector{Int32} #Células-destino concatenadas; -1 designa a célula sumidouro (divergência)
    counts :: Vector{Int32} #Número de amostras que chegou a cada célula-destino
    computed :: BitVector #Marcação de células já processadas
    samples_per_cell :: Int64 #m, denominador das probabilidades
    region :: BasinRegion
end

# Construtor: inicializa mapeamento vazio. As listas crescem conforme as células
# são processadas em ordem crescente dentro de cell_range.
function GeneralizedCellMap(
    region :: BasinRegion,
    n_samples :: Int64;
    cell_range :: UnitRange{Int64} = 1:prod(region.elements),
)
    offsets = Vector{Int32}(undef, length(cell_range) + 1)
    offsets[1] = 1

    return GeneralizedCellMap(
        cell_range,
        offsets,
        Int32[],
        Int32[],
        falses(length(cell_range)),
        n_samples,
        region,
    )
end

# Arranjo de pré-imagens (Hsu §11.2): R(z,j) é a j-ésima célula que tem z entre
# suas imagens, e J(z) o número de pré-imagens de z. Mesmo formato CSR do mapa:
# as pré-imagens do nó z ocupam cells[offsets[z] : offsets[z+1]-1].
# O nó N+1 representa a célula sumidouro (Hsu §11.1: N_c + 1 na implementação).
mutable struct PreImageArray
    offsets :: Vector{Int32} #Índices de início da lista de pré-imagens de cada nó
    cells :: Vector{Int32} #Pré-imagens concatenadas
end

# Resultado completo da análise do GCM (Hsu §11.3 a §11.10). Todos os vetores
# têm N+1 posições: as células 1..N do grid e o sumidouro na posição N+1.
mutable struct GcmAnalysis
    n_periodic_solutions :: Int64 #Nps: soluções periódicas do SCM dominante (§11.3)
    groups :: Vector{Vector{Int64}} #B_g: células de cada grupo persistente; groups[1] é o sumidouro
    periods :: Vector{Int64} #K(g): período de cada grupo persistente (§11.6)
    group_of :: Vector{Int64} #G(z): grupo persistente da célula, ou 0 se transiente
    domicile_count :: Vector{Int64} #Dm(z): número de domicílios da célula (§11.10)
    absorption :: Matrix{Float64} #Ag(z,g): probabilidade de absorção da célula z no grupo g (§11.9)
end