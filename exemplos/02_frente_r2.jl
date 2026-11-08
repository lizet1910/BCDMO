# Família 1 com r = 2, construção (a): mesma matriz A, vetores b_1, b_2 distintos,
# p = (1.5, 1.5), isto é, β = 0.5. Frente de Pareto por múltiplos pontos iniciais.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO
include(joinpath(@__DIR__, "graficos.jl"))

m, n, tam_bloco = 200, 50, 5
N = 100                                   # número de pontos iniciais

A = matriz_sintetica(m, n; densidade=0.1, semente=1)
P = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
blocos = blocos_consecutivos(n, tam_bloco)
opc = Opcoes(ε=1e-2, maxiter=500_000)

X, F, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3)
resumo(ests)

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
salvar_estatisticas(joinpath(pasta, "frente_r2.csv"), ests, F)
savefig(plotar_frente(F; titulo="Família 1, r = 2"), joinpath(pasta, "frente_r2.pdf"))
savefig(plotar_frente(F; titulo="Família 1, r = 2"), joinpath(pasta, "frente_r2.png"))
println("Resultados gravados em ", abspath(pasta))
