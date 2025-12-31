using LinearAlgebra

function calcular_y(num_iteracoes, A, y0)
    y = Vector{Vector{Float64}}(undef, num_iteracoes + 1)
    y[1] = y0
    for n in 1:num_iteracoes
        y[n+1] = A * y[n]
    end
    return y
end

function y_analitico(n)
    return (0.5 + 0.5im)^(n + 1) * [0.0, 1.0, im] + (0.5 - 0.5im)^(n + 1) * [0.0, 1.0, -im] + 0.8^n * [1.0, 0.0, 0.0]
end

function comparar_metodos(num_iteracoes, A, y0)
    #método iterativo
    resultados_iterativo = calcular_y(num_iteracoes, A, y0)

    #método analítico
    resultados_analitico = Vector{Vector{ComplexF64}}(undef, num_iteracoes + 1)
    for n in 0:num_iteracoes
        resultados_analitico[n + 1] = y_analitico(n)
    end

    return resultados_iterativo, resultados_analitico
    iterativo, analitico

end

function main()
    num_iteracoes = 5
    A = [0.8 0.0 0.0;
         0.0   0.5 0.5;
         0.0  -0.5 0.5]
    y0 = [1.0, 1.0, -1.0]

    iterativo, analitico = comparar_metodos(num_iteracoes, A, y0)

    println("Resultados do método iterativo e analitico:")
    println("="^60)

    for i in eachindex(iterativo)
        n = i - 1
        println("n = $n:")
        println("  Iterativo: $(iterativo[i])")
        println("  Analítico: $(analitico[i])")

        diff = Complex.(iterativo[i]) - analitico[i]
        erro = norm(diff)
        println("  Erro (norma): $(erro)")
        println()
    end

    #return iterativo, analitico
end

main()