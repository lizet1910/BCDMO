# Compara o Algoritmo 1 (blocos próprios + bloco 0) com a variante que usa
# apenas o bloco 0 (análoga ao método de [B] com p = 1), e as regras de
# atualização de σ do Passo 6. Mesmos problemas e mesmos pontos iniciais para todos.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, Printf, DelimitedFiles

N = 20            # pontos iniciais por problema
ε = 1e-2

mediana(v) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

problemas = [
    ("Família 1, r = 2", regressao_mesma_matriz(matriz_sintetica(200, 50; densidade=0.1, semente=1), [1.5, 1.5]; semente=2), blocos_consecutivos(50, 5)),
    ("Família 1, r = 3", regressao_linhas(matriz_sintetica(300, 50; densidade=0.1, semente=1), [1.5, 1.7, 1.9]; semente=2), blocos_consecutivos(50, 5)),
    ("Família 2, sobrep. 0", linear_potencia_blocos(blocos_consecutivos(100, 5), 2; p=1.5, sobreposicao=0.0, semente=1), blocos_consecutivos(100, 5)),
    ("Família 2, sobrep. 1", linear_potencia_blocos(blocos_consecutivos(100, 5), 2; p=1.5, sobreposicao=1.0, semente=1), blocos_consecutivos(100, 5)),
]

variantes = [
    ("blocos, σ único ↑",  Opcoes(ε=ε, maxiter=200_000, sigma=:global)),
    ("blocos, σ único ↓",  Opcoes(ε=ε, maxiter=200_000, sigma=:decresce)),
    ("blocos, σ por bloco", Opcoes(ε=ε, maxiter=200_000, sigma=:porbloco)),
    ("só bloco 0",         Opcoes(ε=ε, maxiter=200_000, sigma=:decresce, blocos_proprios=false)),
]

cab = ["problema" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
tabela = Any[]
for (nome, P, blocos) in problemas
    println("\n", nome)
    @printf("  %-20s %6s %9s %9s %11s %9s %9s\n", "variante", "conv.", "iter", "aval. F", "grad. equiv", "bloco 0", "tempo(s)")
    for (rot, opc) in variantes
        local ests
        _, _, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3, mostrar=false)
        conv = count(e -> e.status == :convergiu, ests)
        lin = (mediana([e.iter for e in ests]), mediana([e.avalF for e in ests]),
               mediana([e.grad_equiv for e in ests]), mediana([e.passos_bloco0 for e in ests]),
               mediana([e.tempo for e in ests]))
        @printf("  %-20s %3d/%-2d %9.0f %9.0f %11.1f %9.0f %9.3f\n", rot, conv, N, lin...)
        push!(tabela, Any[nome, rot, conv, lin...])
    end
end

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
open(joinpath(pasta, "comparacao_bloco0.csv"), "w") do io
    writedlm(io, cab, ',')
    for l in tabela
        writedlm(io, permutedims(l), ',')
    end
end
println("\nValores: medianas sobre ", N, " pontos iniciais. Tabela gravada em ", abspath(pasta))
