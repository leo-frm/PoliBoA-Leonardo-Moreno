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

#Funções base
mutable struct BasinRegion # Configura a região onde as bacias serão encontradas
    range :: Vector{Vector{Float64}} #Vetores que delimitam o espaço onde as bacias serão encontradas
    elements :: Vector{Int64} #Número de células (divisões) nesse espaço
    is_cyclic :: Vector{Bool}
    extended_range :: Vector{Vector{Float64}}
end 

mutable struct BasinProblem{F <: Function, P <: Union{Vector, Array}} #Configura a análise de bacias
    f :: F #Função que será integrada para encontrar as bacias
    params :: P #Parâmetros da função
    period :: Float64 #Período de integração
    transient_cycles :: Int64 #Número de ciclos transientes a serem ignorados
    maximum_cycles :: Int64 #Número máximo de ciclos a serem integrados
    maximum_extended_cycles :: Int64 #Número máximo de ciclos a serem integrados na região estendida
    attractor_tolerance :: Int64 #Número de vezes que uma célula deve ser visitada para ser considerada parte de um atrator
    region :: BasinRegion #Região onde as bacias serão encontradas
    threads :: Int64 #Número de threads a serem usadas na integração
end

@enum AttractorType begin #Classificação dos tipos de atratores encontrados
    none #Nenhum tipo de atrator identificado
    regular #Atrator periódico regular
    quasi_periodic #Atrator periódico quase regular
    divergent #Atrator divergente
end

mutable struct Attractor #Armazena informações sobre os atratores encontrados
    number :: Int64 #Identificador do atrator
    kind :: AttractorType #Tipo do atrator
    points :: Vector{Vector{Float64}} #Pontos que compõem o atrator
    finder_thread_id :: Int64 #ID da thread que encontrou o atrator
end

mutable struct BasinResult #Resultado da análise de bacias
    cells #Número do atrator de cada célula
    attractors :: Vector{Attractor} #Lista de atratores encontrados
    region :: BasinRegion #Região analisada
end

# Mapa célula‑a‑célula: para cada célula i, guarda-se o número da célula destino
mutable struct SimpleCellMap # Configura o mapa que armazena células e seus respectivos destinos
    target :: Vector{Int64} #target[i]: destino da célula i
    computed_cells :: Int64 #contador das células calculadas
    range :: BasinRegion
end

# Construtor: inicializa mapeamento vazio
function SimpleCellMap(region::BasinRegion)
    total_cells = prod(region.elements)
    return SimpleCellMap(
        zeros(Int64, total_cells),
        0,
        region,
    )
end