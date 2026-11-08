# Família 1 com r = 3, construção (b): linhas de A repartidas em três grupos,
# expoentes p = (1.5, 1.7, 1.9), isto é, β_i = (0.5, 0.7, 0.9).
# Com p_i = 1.3 (β = 0.3) a convergência fica muito lenta; ver o LEIAME.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO
include(joinpath(@__DIR__, "graficos.jl"))

m, n, tam_bloco = 300, 50, 5
N = 200                                   # cerca de 15 minutos

A = matriz_sintetica(m, n; densidade=0.1, semente=1)
P = regressao_linhas(A, [1.5, 1.7, 1.9]; semente=2)
blocos = blocos_consecutivos(n, tam_bloco)
opc = Opcoes(ε=1e-2, maxiter=200_000)

X, F, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3)
resumo(ests)

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
salvar_estatisticas(joinpath(pasta, "frente_r3.csv"), ests, F)
fig = plotar_frente(F)
savefig(fig, joinpath(pasta, "frente_r3.pdf"))
savefig(fig, joinpath(pasta, "frente_r3.png"))
fig3d = plotar_frente_3d(F)
savefig(fig3d, joinpath(pasta, "frente_r3_3d.pdf"))
savefig(fig3d, joinpath(pasta, "frente_r3_3d.png"))
display(fig)
println("Resultados gravados em ", abspath(pasta))
