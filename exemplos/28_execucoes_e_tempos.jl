# Dados por execução para o Apêndice B e tempos com repetição (cerca de 2 horas).
#   A. Tabela 3, parte sintética (n = 2000, c ∈ {0, 1}), com os mesmos 10 pontos iniciais da
#      Tabela 2 (antes eram os 5 primeiros), para as cinco variantes da Tabela 3.
#   B. Tabela 1 e Tabela 2 com n = 50: cada execução é repetida 5 vezes e guarda-se o menor
#      tempo (tempos de milissegundos estão no limite da resolução do relógio);
#      Tabela 2 com n = 1000: uma execução de cada (tempos de segundos).
# Os iterados não mudam: as iterações e os gradientes equivalentes são os mesmos de antes.
# Grava, sem apagar os arquivos antigos:
#   resultados/tabela3_sinteticas_10.csv            (medianas)
#   resultados/tabela3_sinteticas_10_execucoes.csv  (uma linha por execução)
#   resultados/tabelas_1_2_execucoes.csv            (uma linha por execução)
# Pode ser rodado na mesma sessão do Julia.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
cab_exe = ["tabela" "problema" "c" "variante" "execucao" "status" "iter" "avalF" "grad_equiv" "passos_bloco0" "duplicacoes" "duplicacoes_p5" "tempo" "repeticoes"]
gravar(arq, cab, tab) = open(joinpath(pasta, arq), "w") do io
    writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
end

"Pontos iniciais como em multiplos_inicios: x0 = escala .* randn(n), semente 3."
pontos(n, N, escala) = (rng = MersenneTwister(3); [escala .* randn(rng, n) for _ in 1:N])

"Roda um método a partir de x0, `reps` vezes; devolve as estatísticas da primeira execução com o menor tempo."
function executar(P, x0, bl, met, opc; reps=1)
    melhor = nothing
    for _ in 1:reps
        if met == :bcdmo
            _, _, e = bcdmo(P, x0, bl; opc=opc)
        elseif met == :fs
            _, _, e = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval)
        else
            _, _, e = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval, passo_bb=true, memoria=10)
        end
        melhor === nothing ? (melhor = e) : (melhor.tempo = min(melhor.tempo, e.tempo))
    end
    return melhor
end

linha(tab_, prob, c, rot, k, e, reps) = Any[tab_, prob, c, rot, k, string(e.status), e.iter, e.avalF, e.grad_equiv,
                                            e.passos_bloco0, e.duplicacoes, e.duplicacoes_p5, e.tempo, reps]

# ---------------- B. Tabelas 1 e 2 ----------------
exe12 = Any[]
println("B. Tabela 1 (20 pontos iniciais, 5 repetições por execução)")
probs1 = [
    ("Família 1, r = 2", regressao_mesma_matriz(matriz_sintetica(200, 50; densidade=0.1, semente=1), [1.5, 1.5]; semente=2), blocos_consecutivos(50, 5)),
    ("Família 1, r = 3", regressao_linhas(matriz_sintetica(300, 50; densidade=0.1, semente=1), [1.5, 1.7, 1.9]; semente=2), blocos_consecutivos(50, 5)),
    ("Família 2, ω = 0", linear_potencia_blocos(blocos_consecutivos(100, 5), 2; p=1.5, sobreposicao=0.0, semente=1), blocos_consecutivos(100, 5)),
    ("Família 2, ω = 1", linear_potencia_blocos(blocos_consecutivos(100, 5), 2; p=1.5, sobreposicao=1.0, semente=1), blocos_consecutivos(100, 5)),
]
regras1 = [("(i)", (; sigma=:global)), ("(ii)", (; sigma=:decresce)), ("(iii)", (; sigma=:porbloco)),
           ("só bloco 0", (; sigma=:decresce, blocos_proprios=false))]
for (nome, P, bl) in probs1
    xs = pontos(dim(P), 20, 1.0)
    for (rot, kw) in regras1
        opc = Opcoes(; ε=1e-2, maxiter=200_000, kw...)
        for (k, x0) in enumerate(xs)
            push!(exe12, linha("1", nome, "", rot, k, executar(P, x0, bl, :bcdmo, opc; reps=5), 5))
        end
        @printf("   %-18s %-11s pronto\n", nome, rot)
        gravar("tabelas_1_2_execucoes.csv", cab_exe, exe12)
    end
end

println("B. Tabela 2 com n = 50 (20 pontos, 5 repetições) e n = 1000 (10 pontos, 1 execução)")
regras2 = [("regra por bloco", (; sigma=:porbloco)), ("regra (ii)", (; sigma=:decresce)),
           ("só bloco 0", (; sigma=:decresce, blocos_proprios=false))]
for (m, n, dens, tb, cs, N, reps) in ((200, 50, 0.1, 5, (0.0, 0.5, 1.0, 1.5), 20, 5),
                                      (3000, 1000, 0.005, 50, (0.0, 1.0), 10, 1))
    A0 = matriz_sintetica(m, n; densidade=dens, semente=1); bl = blocos_consecutivos(n, tb)
    for c in cs
        A, esc = escalar_blocos(A0, bl, c)
        P = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
        xs = pontos(n, N, 1 ./ esc)
        for (rot, kw) in regras2
            opc = Opcoes(; ε=1e-2, maxiter=200_000, kw...)
            for (k, x0) in enumerate(xs)
                push!(exe12, linha("2", "n = $n", c, rot, k, executar(P, x0, bl, :bcdmo, opc; reps=reps), reps))
            end
            @printf("   n = %-5d c = %.1f %-16s pronto\n", n, c, rot)
            gravar("tabelas_1_2_execucoes.csv", cab_exe, exe12)
        end
    end
end

# ---------------- A. Tabela 3, sintéticas com 10 pontos iniciais ----------------
println("\nA. Tabela 3, sintéticas, 10 pontos iniciais")
vars3 = [variantes; [variante_fsbb]]          # blocos 1ª, blocos 2ª, só bloco 0, Fliege–Svaiter, Fliege–Svaiter BB
A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1); bl = blocos_consecutivos(2000, 100)
exe3 = Any[]; med3 = Any[]
cab_med = ["instancia" "m" "n" "c" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
for c in (0.0, 1.0)
    A, esc = escalar_blocos(A0, bl, c)
    xs = pontos(2000, 10, 1 ./ esc)
    for (rot, modelo, met, kw) in vars3
        P = problema(A, modelo)
        opc = Opcoes(; ε=1e-2, maxiter=300_000, maxaval=3_000_000, kw...)
        es = [executar(P, x0, bl, met, opc) for x0 in xs]
        for (k, e) in enumerate(es)
            push!(exe3, linha("3", "sintetica", c, rot, k, e, 1))
        end
        med = Any["sintetica", 6000, 2000, c, rot, count(e -> e.status == :convergiu, es),
                  mediana([e.iter for e in es]), mediana([e.avalF for e in es]), mediana([e.grad_equiv for e in es]),
                  mediana([e.passos_bloco0 for e in es]), mediana([e.tempo for e in es])]
        push!(med3, med)
        @printf("   c = %.1f %-22s conv. %2d/10  grad. equiv %10.1f  tempo %8.3f\n", c, rot, med[6], med[9], med[11])
        gravar("tabela3_sinteticas_10_execucoes.csv", cab_exe, exe3)
        gravar("tabela3_sinteticas_10.csv", cab_med, med3)
    end
end
println("\nResultados em ", abspath(pasta))
