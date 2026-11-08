# Script curto (cerca de 30 a 40 minutos), antes da rodada grande (18b):
#   1. verificação rápida das peças novas: perda de Cauchy (não convexa), modelo de
#      segunda ordem com matrizes indefinidas, Fliege–Svaiter com passo BB e busca não
#      monótona, soma ponderada e métricas de frente;
#   2. sensibilidade da regra por bloco do Passo 6: (iii) σ_j = max{1, σ/2} a cada passo
#      aceito; (iii-a) redução à metade só após dois passos aceitos seguidos sem
#      duplicação; (iii-b) σ_j = max{1, σ/√2};
#   3. indicador L/L̄ (razão entre ‖A‖² e a média de ‖A U_j‖²) das matrizes das Tabelas 2 a 4.
# Resultados em resultados/sensibilidade_regra.csv e resultados/indicador_blocos.csv.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
gravar(arq, cab, tab) = open(joinpath(pasta, arq), "w") do io
    writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
end

# ---------------- 1. verificação rápida ----------------
println("1. Verificação rápida")
A = matriz_sintetica(120, 30; densidade=0.15, semente=1)
Pc = regressao_linhas(A, [1.5, 1.5]; semente=2, perda=:cauchy)
x = randn(MersenneTwister(5), 30)
g = jac_bloco(Pc, x, collect(1:30))
h = 1e-6; erro = 0.0
for l in 1:30
    e = zeros(30); e[l] = h
    df = (fvals(Pc, x .+ e) .- fvals(Pc, x .- e)) ./ (2h)
    global erro = max(erro, maximum(abs.(df .- g[l, :])))
end
@printf("   gradiente da perda de Cauchy, erro máximo por diferenças finitas: %.1e\n", erro)
blocos_t = blocos_consecutivos(30, 5)
for (rot, P) in (("Cauchy, 1ª ordem", Pc),
                 ("Cauchy, 2ª ordem", regressao_linhas(A, [1.5, 1.5]; semente=2, perda=:cauchy,
                                                       usar_hess=true, hess_comum=true)))
    _, _, ests = multiplos_inicios(P, blocos_t, 3; opc=Opcoes(ε=1e-2, maxiter=100_000), semente=3, mostrar=false)
    @printf("   %-18s convergiram %d/3, iterações %s\n", rot, count(e -> e.status == :convergiu, ests),
            string([e.iter for e in ests]))
end
Pq = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
for met in (:fs, :fsbb)
    _, _, ests = multiplos_inicios(Pq, blocos_t, 3; opc=Opcoes(ε=1e-2, maxiter=100_000), semente=3,
                                   mostrar=false, metodo=met)
    @printf("   %-18s convergiram %d/3, iterações %s\n", string(met), count(e -> e.status == :convergiu, ests),
            string([e.iter for e in ests]))
end
E = Escalarizado(Pq, [0.5, 0.5])
_, _, ests = multiplos_inicios(E, blocos_t, 3; opc=Opcoes(ε=1e-2, maxiter=100_000), semente=3, mostrar=false)
@printf("   %-18s convergiram %d/3\n", "soma ponderada", count(e -> e.status == :convergiu, ests))
F1 = [0.0 0.5 1.0; 1.0 0.5 0.0]; F2 = [0.0 1.0; 1.2 0.1]
for (k, m) in enumerate(metricas_frentes([F1, F2]))
    println("   métricas da frente de teste ", k, ": ", m)
end

# ---------------- 2. sensibilidade da regra por bloco ----------------
println("\n2. Sensibilidade da regra por bloco")
regras = [("(iii) σ/2 a cada sucesso", :porbloco), ("(iii-a) σ/2 após 2 sucessos", :porbloco2),
          ("(iii-b) σ/√2 a cada sucesso", :porbloco_raiz)]
casos = Any[
    ("Família 1, r = 2", regressao_mesma_matriz(matriz_sintetica(200, 50; densidade=0.1, semente=1), [1.5, 1.5]; semente=2), blocos_consecutivos(50, 5), 20, 1.0),
    ("Família 1, r = 3", regressao_linhas(matriz_sintetica(300, 50; densidade=0.1, semente=1), [1.5, 1.7, 1.9]; semente=2), blocos_consecutivos(50, 5), 20, 1.0),
    ("Família 2, ω = 0", linear_potencia_blocos(blocos_consecutivos(100, 5), 2; p=1.5, sobreposicao=0.0, semente=1), blocos_consecutivos(100, 5), 20, 1.0),
    ("Família 2, ω = 1", linear_potencia_blocos(blocos_consecutivos(100, 5), 2; p=1.5, sobreposicao=1.0, semente=1), blocos_consecutivos(100, 5), 20, 1.0),
]
A0 = matriz_sintetica(3000, 1000; densidade=0.005, semente=1); bl = blocos_consecutivos(1000, 50)
for c in (0.0, 1.0)
    Ac, esc = escalar_blocos(A0, bl, c)
    push!(casos, ("n = 1000, c = $c", regressao_mesma_matriz(Ac, [1.5, 1.5]; semente=2), bl, 10, 1 ./ esc))
end
for (grupo, nome) in (("JGD_Homology", "ch6-6-b3"), ("JGD_Franz", "Franz5"))
    Ar = obter_matriz(grupo, nome)
    push!(casos, (nome, regressao_mesma_matriz(Ar, [1.5, 1.5]; semente=2),
                  blocos_consecutivos(size(Ar, 2), cld(size(Ar, 2), 10)), 10, 1.0))
end
cabS = ["problema" "regra" "convergiram" "iter" "avalF" "grad_equiv" "tempo"]
tabS = Any[]
for (nome, P, bls, N, escala) in casos
    println("\n", nome)
    @printf("  %-30s %6s %9s %9s %11s %9s\n", "regra", "conv.", "iter", "aval. F", "grad. equiv", "tempo(s)")
    for (rot, regra) in regras
        _, _, ests = multiplos_inicios(P, bls, N; opc=Opcoes(ε=1e-2, maxiter=300_000, sigma=regra),
                                       semente=3, escala=escala, mostrar=false)
        conv = count(e -> e.status == :convergiu, ests)
        lin = (mediana([e.iter for e in ests]), mediana([e.avalF for e in ests]),
               mediana([e.grad_equiv for e in ests]), mediana([e.tempo for e in ests]))
        @printf("  %-30s %3d/%-2d %9.0f %9.0f %11.1f %9.3f\n", rot, conv, N, lin...)
        push!(tabS, Any[nome, rot, conv, lin...])
        gravar("sensibilidade_regra.csv", cabS, tabS)
    end
end

# ---------------- 3. indicador L/L̄ ----------------
println("\n3. Indicador L/L̄")
cabI = ["matriz" "m" "n" "q" "c" "L_sobre_media_Lj"]
tabI = Any[]
for (m, n, dens, tb, cs) in ((200, 50, 0.1, 5, (0.0, 0.5, 1.0, 1.5)), (3000, 1000, 0.005, 50, (0.0, 1.0)),
                             (6000, 2000, 0.002, 100, (0.0, 1.0)), (15000, 5000, 0.001, 100, (0.0, 1.0)))
    A0s = matriz_sintetica(m, n; densidade=dens, semente=1); bls = blocos_consecutivos(n, tb)
    for c in cs
        As, _ = escalar_blocos(A0s, bls, c)
        v = indicador_blocos(As, bls)
        @printf("  sintética n = %5d, c = %.1f: L/L̄ = %.2f\n", n, c, v)
        push!(tabI, Any["sintetica", m, n, length(bls), c, v])
        gravar("indicador_blocos.csv", cabI, tabI)
    end
end
for (grupo, nome) in (("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
                      ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
                      ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10"))
    Ar = obter_matriz(grupo, nome); m, n = size(Ar); bls = blocos_consecutivos(n, cld(n, 10))
    v = indicador_blocos(Ar, bls)
    @printf("  %-10s: L/L̄ = %.2f\n", nome, v)
    push!(tabI, Any[nome, m, n, length(bls), "", v])
    gravar("indicador_blocos.csv", cabI, tabI)
end
println("\nResultados gravados em ", abspath(pasta))
