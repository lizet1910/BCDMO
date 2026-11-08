# Perda de Cauchy, construção (b): valores de F nos pontos finais (rodada de cerca de 3 a 4 horas;
# pode ser interrompida e retomada: execuções já gravadas não são repetidas).
#   A. Matrizes reais (10 pontos iniciais) e sintética com n = 2000, c = 0 (5 pontos), os mesmos
#      da Tabela 7: método por blocos de 1ª ordem, variante só bloco 0 e Fliege–Svaiter BB.
#   B. Sintética com n = 2000, c = 1, orçamento de 300 s por execução, 3 pontos iniciais:
#      blocos de 1ª e de 2ª ordem, só bloco 0 e Fliege–Svaiter BB.
#   C. Método por blocos de 2ª ordem nas instâncias da parte A, exceto Franz8 (cerca de 16 min
#      por execução).
# Grava resultados/cauchy_valores_F.csv (uma linha por execução, com f_1 e f_2 no ponto final).

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

Sys.iswindows() && ccall((:SetThreadExecutionState, "kernel32"), UInt32, (UInt32,), 0x80000001)

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
arq = joinpath(pasta, "cauchy_valores_F.csv")
cab = ["instancia" "c" "variante" "execucao" "status" "iter" "avalF" "grad_equiv" "G_final" "f1" "f2" "tempo" "orcamento"]

tab = Any[]; feitas = Set{Tuple{String,Float64,String,Int}}()
if isfile(arq)
    d, _ = readdlm(arq, ','; header=true)
    for l in 1:size(d, 1)
        lin = Any[d[l, k] for k in 1:size(d, 2)]
        # execução com tempo muito acima do orçamento: provável suspensão do computador; refaz
        orc = Float64(lin[13]); Float64(lin[12]) > orc + 60 && continue
        push!(tab, lin); push!(feitas, (string(lin[1]), Float64(lin[2]), string(lin[3]), Int(lin[4])))
    end
end
gravar() = open(arq, "w") do io
    writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
end
pontos(n, N, escala) = (rng = MersenneTwister(3); [escala .* randn(rng, n) for _ in 1:N])

const V1  = ("blocos, 1ª ordem",  :linear, :bcdmo, (; sigma=:porbloco))
const V2  = ("blocos, 2ª ordem",  :hess,   :bcdmo, (; sigma=:porbloco))
const VB0 = ("só bloco 0",        :linear, :bcdmo, (; sigma=:decresce, blocos_proprios=false))
const VBB = ("Fliege–Svaiter BB", :linear, :fsbb,  (;))

function rodar!(nome, c, A, bl, escala, N, vars; orcamento=Inf)
    n = size(A, 2); xs = pontos(n, N, escala)
    probs = Dict{Symbol,Any}()
    for (rot, modelo, met, kw) in vars
        P = get!(() -> problema(A, modelo; construcao=:b, perda=:cauchy), probs, modelo)
        opc = Opcoes(; ε=1e-2, maxiter=(isinf(orcamento) ? 300_000 : 100_000_000),
                     maxaval=(isinf(orcamento) ? 3_000_000 : 1_000_000_000), tempo_max=orcamento, kw...)
        for (k, x0) in enumerate(xs)
            (nome, Float64(c), rot, k) in feitas && continue
            if met == :bcdmo
                _, Fx, e = bcdmo(P, x0, bl; opc=opc)
            else
                _, Fx, e = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval,
                                      passo_bb=true, memoria=10, tempo_max=opc.tempo_max)
            end
            push!(tab, Any[nome, c, rot, k, string(e.status), e.iter, e.avalF, e.grad_equiv, e.Gfinal,
                           Fx[1], Fx[2], e.tempo, isinf(orcamento) ? 1e9 : orcamento])
            push!(feitas, (nome, Float64(c), rot, k)); gravar()
            @printf("  %-9s c = %.0f %-18s %2d  %-9s grad. equiv %10.1f  F = (%.6e, %.6e)  G %.2e  tempo %7.1f s\n",
                    nome, c, rot, k, string(e.status), e.grad_equiv, Fx[1], Fx[2], e.Gfinal, e.tempo)
        end
    end
end

nomes_reais = [("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
               ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
               ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10")]
reais = [(nome, obter_matriz(g, nome)) for (g, nome) in nomes_reais]
A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1); blS = blocos_consecutivos(2000, 100)
AS0, esc0 = escalar_blocos(A0, blS, 0.0)
AS1, esc1 = escalar_blocos(A0, blS, 1.0)

println("A. Cauchy, 1ª ordem, valores de F nos pontos finais")
for (nome, A) in reais
    rodar!(nome, 0.0, A, blocos_consecutivos(size(A, 2), cld(size(A, 2), 10)), 1.0, 10, [V1, VB0, VBB])
end
rodar!("sintetica", 0.0, AS0, blS, 1 ./ esc0, 5, [V1, VB0, VBB])

println("\nB. Sintética com c = 1, orçamento de 300 s")
rodar!("sintetica", 1.0, AS1, blS, 1 ./ esc1, 3, [V1, V2, VB0, VBB]; orcamento=300.0)

println("\nC. Blocos de 2ª ordem")
for (nome, A) in reais
    nome == "Franz8" && continue
    rodar!(nome, 0.0, A, blocos_consecutivos(size(A, 2), cld(size(A, 2), 10)), 1.0, 10, [V2])
end
rodar!("sintetica", 0.0, AS0, blS, 1 ./ esc0, 5, [V2])

Sys.iswindows() && ccall((:SetThreadExecutionState, "kernel32"), UInt32, (UInt32,), 0x80000000)
println("\nResultados em ", abspath(arq))
