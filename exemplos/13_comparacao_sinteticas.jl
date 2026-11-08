# Mesma comparação do script 12 (blocos com modelo linear e de segunda ordem, só o
# bloco 0 e Fliege–Svaiter), nas matrizes sintéticas do script 09 com n = 2000,
# escalas iguais (c = 0) e diferentes entre os blocos (c = 1).
# Leva mais tempo que o script 12 (sobretudo só o bloco 0 e Fliege–Svaiter com c = 1).

include(joinpath(@__DIR__, "comum_comparacao.jl"))

N = 5
ε = 1e-2
m, n, dens, tb = 6000, 2000, 0.002, 100

A0 = matriz_sintetica(m, n; densidade=dens, semente=1)
blocos = blocos_consecutivos(n, tb)
cab = ["cenario" "c" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
tabela = Any[]
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
for c in [0.0, 1.0]
    local A, esc
    A, esc = escalar_blocos(A0, blocos, c)
    println("\nn = 2000, c = ", c, "  (escalas dos blocos de 1 a ", round(10.0^c; sigdigits=3), ")")
    cabecalho()
    for l in rodar(A, blocos, N, ε; escala=1 ./ esc)
        push!(tabela, Any["n = 2000", c, l...])
    end
    open(joinpath(pasta, "comparacao_sinteticas.csv"), "w") do io
        writedlm(io, cab, ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tabela)
    end
end
println("\nMedianas sobre os pontos iniciais. Tabela gravada em ", abspath(pasta))
