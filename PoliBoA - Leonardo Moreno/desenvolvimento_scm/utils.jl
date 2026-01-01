using DifferentialEquations

function get_cell_number(u, region::BasinRegion) :: Int64 #Retorna o número da célula correspondente às entradas: valores u e região da bacia
    # Initialize cell
    cell = 0

    # cell dimensions in each direction
    for i in eachindex(region.range) #Loop que itera de 1 até o número de dimensões da região
        length = (region.range[i][2] - region.range[i][1]) / region.elements[i] #Calcula o tamanho de cada célula na dimensão i
        index = clamp(ceil(Int64, (u[i] - region.range[i][1]) / length), 1, region.elements[i]) #u[i] - region.range[i][1] calcula a distância do ponto até a origem da região; dividido por length (tamanho de cada célula) dá a posição em "unidades de célula"; ceil arredonda para cima garantindo que todos os pontos dentro de uma célula tenham o mesmo número; clamp protege contra índices inválidos forçando o resultado a ficar entre 1 e o número máximo de células.
        cell = i == 1 ? index : cell + (index - 1) * prod(region.elements[1:i-1]) # Calcula o número da célula considerando a posição na dimensão i e as dimensões anteriores. O primeiro índice é simplesmente o índice calculado, enquanto os demais são ajustados para considerar a posição relativa dentro do espaço multidimensional.
    end

    # Return cell
    return cell
end

function store_cell_center!(u :: Vector{Float64}, cell::Int64, region::BasinRegion) #Retorna os valores exatos do centro da célula correspondente ao número da célula e região da bacia
    indexes = Vector{Int64}()

    pos = cell
    # First position
    pos_i = pos % region.elements[1] == 0 ? region.elements[1] : pos % region.elements[1] # Calcula a coordenada da célula na 1ª dimensão
    push!(indexes, pos_i) #Coloca a coordenada da célula na 1ª dimensão no vetor indexes

    for i = 2:lastindex(region.range)
        # Get position of generalized coordinate
        qi = pos - pos_i == 0 ? 1 : pos - pos_i # Variável intermediária qi que armazena a posição relativa da célula atual em relação à primeira dimensão
        # Update pos
        pos :: Int64 = floor(qi / region.elements[i - 1]) + 1 # Calcula a posição da célula em relação às linhas
        pos_i = pos % region.elements[i] == 0 ? region.elements[i] : pos % region.elements[i] #Calcula a coordenada
        push!(indexes, pos_i) # Adiciona a coordenada da célula na dimensão i ao vetor indexes
    end

    # cell dimensions in each direction
    for i in eachindex(region.range)
        length = (region.range[i][2] - region.range[i][1]) / region.elements[i] # Cálcula o tamanho de cada célula na dimensão i
        u[i] = (indexes[i] - 0.5) * length + region.range[i][1] # Calcula o centro da célula na dimensão i e armazena no vetor u
    end
end 

function is_inside_range(pos, range) #Verifica se os valores de pos estão dentro dos limites especificados em range
    # Declare in-line function for single dimension check
    is_inside(pos_i, range_i) = range_i[1] <= pos_i <= range_i[2]
    # Join the range checks. If all inside, returns true; returns false otherwise.
    return reduce(&, [is_inside(p, r) for (p, r) in zip(pos, range)])
end

function set_integrator!(integrator; u::Vector{Float64}=nothing, t::Float64=nothing)
    if u !== nothing
        set_u!(integrator, u)
        u_modified!(integrator, true)
    end

    t === nothing || set_t!(integrator, t)
end

function adjust_cyclic(u :: Vector{Float64}, range:: Vector{Vector{Float64}}, is_cyclic :: Vector{Bool})

    function adjust_cyclic_single(u :: Float64; maximum:: Float64 = pi, minimum :: Float64 = pi)
        # Number of turns
        turns = fld(u, 2 * pi)

        # Calculates the new value removing the number of turns
        adjusted_u = u - turns * 2 * pi

        # Checks if it is within the region and returns the corrected value
        if adjusted_u > maximum
            return adjusted_u - 2 * pi
        end

        # Checks if it is within the region and returns the corrected value
        if adjusted_u < minimum
            return adjusted_u + 2 * pi
        end

        # Returns the corrected value
        return adjusted_u
    end

    for i in eachindex(u)
        if is_cyclic[i]
            u[i] = adjust_cyclic_single(u[i], minimum = range[i][1], maximum = range[i][2])
        end
    end
end

function make_attractor_from_cells(cells; number::Int64, type::AttractorType, region::BasinRegion)
    points = map(cells) do c
        u =  zeros(Float64, length(region.elements))
        store_cell_center!(u, c, region)
        return u
    end

    return Attractor(number, type, points, Threads.threadid())
end
