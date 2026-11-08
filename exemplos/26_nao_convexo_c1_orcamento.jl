# Problema não convexo (perda de Cauchy, construção (b)) na sintética com n = 2000 e c = 1,
# caso que não terminou num tempo razoável na rodada do script 19. Aqui, cada método tem o
# mesmo orçamento de tempo por execução; relatamos se atingiu G(x) ≤ ε, o G(x) final, as
# iterações e os gradientes equivalentes usados.
# Métodos: Algoritmo 1 com modelos de 1ª e 2ª ordem, só bloco 0 e Fliege–Svaiter BB.
# Duração máxima: 4 métodos × N pontos iniciais × orçamento (com N = 3 e 300 s, cerca de 1 hora).
# É preciso reiniciar o Julia antes, pois `src` mudou (campo tempo_max em Opcoes).
# Grava resultados/nao_convexo_c1_orcamento.csv (uma linha por execução).

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

orcamento = 300.0   # segundos por execução
N = 3
ε = 1e-2

A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1)
bl = blocos_consecutivos(2000, 100)
A, esc = escalar_blocos(A0, bl, 1.0)
escala = 1 ./ esc

metodos = [variantes[1], variantes[2], variantes[3], variante_fsbb]   # (rótulo, modelo, método, opções)
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
cab = ["metodo" "execucao" "status" "iter" "avalF" "grad_equiv" "G_final" "duplicacoes_p5" "tempo"]
tab = Any[]
println("Sintética n = 2000, c = 1, perda de Cauchy, construção (b); orçamento de ", orcamento, " s por execução")
for (rot, modelo, met, kw) in metodos
    P = problema(A, modelo; construcao=:b, perda=:cauchy)
    opc = Opcoes(; ε=ε, maxiter=100_000_000, maxaval=1_000_000_000, tempo_max=orcamento, kw...)
    _, _, ests = multiplos_inicios(P, bl, N; opc=opc, semente=3, escala=escala, mostrar=false, metodo=met)
    for (k, e) in enumerate(ests)
        @printf("  %-20s execução %d: %-10s iter %9d  grad. equiv %10.1f  G final %.3e  dupl. P5 %5d  tempo %7.1f s\n",
                rot, k, string(e.status), e.iter, e.grad_equiv, e.Gfinal, e.duplicacoes_p5, e.tempo)
        push!(tab, Any[rot, k, string(e.status), e.iter, e.avalF, e.grad_equiv, e.Gfinal, e.duplicacoes_p5, e.tempo])
    end
    open(joinpath(pasta, "nao_convexo_c1_orcamento.csv"), "w") do io
        writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
    end
end
println("\nResultados em ", abspath(joinpath(pasta, "nao_convexo_c1_orcamento.csv")))
