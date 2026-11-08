# Complementos para a Seção 8 (cerca de 5 a 10 minutos):
#   1. fração de passos no bloco 0 na Família 2 em função da sobreposição ω (o script 04
#      só imprime esses valores); grava resultados/estresse_bloco0.csv;
#   2. indicador L/L̄ das nove matrizes reais, com o método da potência corrigido
#      (a rodada do script 18a deu valores errados nas matrizes reais); grava
#      resultados/indicador_blocos_reais.csv.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)

function parte1()
    println("1. Família 2: passos no bloco 0 em função de ω")
    n, tb, r, N = 100, 5, 2, 10
    bl = blocos_consecutivos(n, tb)
    opc = Opcoes(ε=1e-2, maxiter=50_000)
    tab = Any[]
    for ω in (0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)
        P = linear_potencia_blocos(bl, r; p=1.5, sobreposicao=ω, semente=1)
        _, _, ests = multiplos_inicios(P, bl, N; opc=opc, semente=2, mostrar=false)
        conv = count(e -> e.status == :convergiu, ests)
        frac = sum(e.passos_bloco0 for e in ests) / max(1, sum(e.iter for e in ests))
        it = mediana([e.iter for e in ests]); b0 = mediana([e.passos_bloco0 for e in ests])
        @printf("   ω = %.2f  convergiram %2d/%d  fração no bloco 0 = %.4f  mediana de iterações = %.1f  mediana de passos no bloco 0 = %.1f\n",
                ω, conv, N, frac, it, b0)
        push!(tab, Any[ω, conv, frac, it, b0])
    end
    open(joinpath(pasta, "estresse_bloco0.csv"), "w") do io
        writedlm(io, ["omega" "convergiram" "fracao_bloco0" "mediana_iter" "mediana_passos_bloco0"], ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tab)
    end
end

function parte2()
    println("\n2. Indicador L/L̄ das matrizes reais")
    nomes = [("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
             ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
             ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10")]
    tab = Any[]
    for (g, nome) in nomes
        A = obter_matriz(g, nome); m, n = size(A)
        bl = blocos_consecutivos(n, cld(n, 10))
        v = indicador_blocos(A, bl)
        normas = [norm(A[:, l]) for l in 1:n]
        @printf("   %-10s L/L̄ = %6.2f   normas das colunas entre %.2f e %.2f\n", nome, v, minimum(normas), maximum(normas))
        push!(tab, Any[nome, m, n, length(bl), v, minimum(normas), maximum(normas)])
    end
    open(joinpath(pasta, "indicador_blocos_reais.csv"), "w") do io
        writedlm(io, ["matriz" "m" "n" "q" "L_sobre_media_Lj" "norma_coluna_min" "norma_coluna_max"], ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tab)
    end
end

parte1()
parte2()
println("\nResultados em ", abspath(pasta))
