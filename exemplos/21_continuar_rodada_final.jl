# Continuação do script 19, que parou na etapa 3 (problema não convexo nas sintéticas com
# c = 1, lenta demais). Roda as etapas 4 a 8 do script 19 e, no fim, o FF1 (script 20):
#   4. construção (b) com a perda convexa nas matrizes reais;
#   5. Tabela 4 nas matrizes reais com 20 pontos iniciais e perfis de desempenho;
#   6. Tabela 4 nas sintéticas;
#   7. Tabela 5 (ε = 1e-3) com 20 pontos iniciais e perfis;
#   8. experimentos das Tabelas 1 a 3, Figura 1 e Seções 8.3 e 8.4 com a regra nova;
#   e o FF1 (frente côncava).
# Tempo estimado: 7 a 9 horas. Deixe o computador na tomada e sem suspensão.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra, Plots

raiz = joinpath(@__DIR__, "..")
pasta_r = joinpath(raiz, "resultados"); mkpath(pasta_r)
pasta_f = joinpath(raiz, "figuras"); mkpath(pasta_f)
etapa(nome) = println("\n", "#"^70, "\n# ", nome, "  (", Libc.strftime("%H:%M", time()), ")\n", "#"^70)
function gravar(arq, cab, tab)
    open(joinpath(pasta_r, arq), "w") do io
        writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
    end
end

# ---------------- 0. cópia de segurança ----------------
for (o, d) in (("resultados", "resultados_sigma2"), ("figuras", "figuras_sigma2"))
    po = joinpath(raiz, o); pd = joinpath(raiz, d); mkpath(pd)
    for f in readdir(po)
        isfile(joinpath(pd, f)) || cp(joinpath(po, f), joinpath(pd, f))
    end
end

nomes_reais = [("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
               ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
               ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10")]
matrizes = Dict(nome => obter_matriz(g, nome) for (g, nome) in nomes_reais)
reais() = [(nome, matrizes[nome], blocos_consecutivos(size(matrizes[nome], 2), cld(size(matrizes[nome], 2), 10)), 1.0, "")
           for (_, nome) in nomes_reais]
function sinteticas(construcao_c)
    A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1); bl = blocos_consecutivos(2000, 100)
    lista = Any[]
    for c in construcao_c
        A, esc = escalar_blocos(A0, bl, c)
        push!(lista, ("sintetica", A, bl, 1 ./ esc, c))
    end
    return lista
end

cab_med = ["instancia" "m" "n" "c" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
cab_exe = ["instancia" "variante" "execucao" "status" "iter" "avalF" "grad_equiv" "duplicacoes" "duplicacoes_p5" "tempo"]

"Roda `vars` em cada instância de `lista`; grava medianas e execuções a cada instância."
function conjunto(arq, lista, N, ε; vars, construcao=:a, perda=:potencia)
    med = Any[]; exe = Any[]
    for (nome, A, bl, escala, c) in lista
        m, n = size(A)
        rot = c == "" ? nome : "$(nome), c = $(c)"
        println("\n", rot, "  (", m, " × ", n, ", ", length(bl), " blocos), ε = ", ε,
                ", construção ", construcao, ", perda ", perda)
        cabecalho()
        for l in rodar(A, bl, N, ε; escala=escala, vars=vars, construcao=construcao, perda=perda,
                       registro=exe, rotulo=rot)
            push!(med, Any[nome, m, n, c, l...])
        end
        gravar(arq * ".csv", cab_med, med)
        gravar(arq * "_execucoes.csv", cab_exe, exe)
    end
    return exe
end

"Perfis de desempenho (gradientes equivalentes e tempo) a partir das execuções."
function perfis(exe, nomes, arquivo)
    problemas = unique([(l[1], l[3]) for l in exe])
    for (col, medida) in ((7, "gradientes equivalentes"), (10, "tempo"))
        T = fill(Inf, length(problemas), length(nomes))
        for l in exe
            l[4] == "convergiu" || continue
            p = findfirst(==((l[1], l[3])), problemas); s = findfirst(==(l[2]), nomes)
            (p === nothing || s === nothing) && continue
            T[p, s] = Float64(l[col])
        end
        τ, ρ = perfil_desempenho(T)
        fig = plot(; xlabel="log₂ τ", ylabel="fração dos problemas", legend=:bottomright,
                   title=medida, ylims=(0, 1.02), size=(520, 400))
        for (s, nome) in enumerate(nomes)
            plot!(fig, vcat(0.0, log2.(τ)), vcat(ρ[1, s], ρ[:, s]); seriestype=:steppost, label=nome, lw=2)
        end
        savefig(fig, joinpath(pasta_f, "$(arquivo)_$(col == 7 ? "grad" : "tempo").pdf"))
    end
end

# ---------------- 1. verificação rápida ----------------
etapa("1. Verificação rápida")
Av = matriz_sintetica(120, 30; densidade=0.15, semente=1); blv = blocos_consecutivos(30, 5)
for (rot, P, met) in (("convexa, blocos", regressao_mesma_matriz(Av, [1.5, 1.5]; semente=2), :bcdmo),
                      ("Cauchy, 2ª ordem", regressao_linhas(Av, [1.5, 1.5]; semente=2, perda=:cauchy,
                                                            usar_hess=true, hess_comum=true), :bcdmo),
                      ("Fliege–Svaiter BB", regressao_mesma_matriz(Av, [1.5, 1.5]; semente=2), :fsbb),
                      ("soma ponderada", Escalarizado(regressao_mesma_matriz(Av, [1.5, 1.5]; semente=2), [0.3, 0.7]), :bcdmo))
    _, _, ests = multiplos_inicios(P, blv, 3; opc=Opcoes(ε=1e-2, maxiter=100_000), semente=3, mostrar=false, metodo=met)
    @printf("   %-20s convergiram %d/3, duplicações no Passo 5: %s\n", rot,
            count(e -> e.status == :convergiu, ests), string([e.duplicacoes_p5 for e in ests]))
end

vars_nc = [variantes[1], variantes[2], variantes[3], variante_fsbb]

# ---------------- 4. construção (b), perda convexa ----------------
etapa("4. Construção (b), perda convexa, matrizes reais")
conjunto("construcao_b_reais", reais(), 10, 1e-2; vars=vars_nc, construcao=:b)

# ---------------- 5 a 7. Tabelas 4 e 5 ----------------
vars_t4 = [variantes; [variante_fsbb, variante_b0_hess]]
nomes_t4 = [v[1] for v in vars_t4]
etapa("5. Tabela 4, matrizes reais, 20 pontos iniciais")
exe4 = conjunto("tabela4_reais", reais(), 20, 1e-2; vars=vars_t4)
perfis(exe4, nomes_t4, "perfil_tabela4")
etapa("6. Tabela 4, sintéticas")
conjunto("tabela4_sinteticas", sinteticas((0.0, 1.0)), 5, 1e-2; vars=[variantes; [variante_fsbb]])
etapa("7. Tabela 5 (ε = 1e-3), matrizes reais, 20 pontos iniciais")
vars_t5 = [variantes; [variante_fsbb]]
exe5 = conjunto("tabela5_reais", reais(), 20, 1e-3; vars=vars_t5)
perfis(exe5, [v[1] for v in vars_t5], "perfil_tabela5")

# ---------------- 8. experimentos anteriores com a regra nova ----------------
for (rot, arq) in (("01: exemplo da Observação 2.5", "01_exemplo_observacao.jl"),
                   ("07: Tabela 1", "07_comparacao_bloco0.jl"),
                   ("10: testes de aceitação", "10_testes_aceitacao.jl"),
                   ("04: uso do bloco 0 na Família 2", "04_estresse_bloco0.jl"),
                   ("02: frente com r = 2", "02_frente_r2.jl"),
                   ("03: frente com r = 3", "03_frente_r3.jl"),
                   ("11: Figura 1", "11_figuras_artigo.jl"),
                   ("05: Tabela 3", "05_matrizes_reais.jl"),
                   ("08: Tabela 2, n = 50 e n = 1000", "08_escalas_heterogeneas.jl"),
                   ("09: Tabela 2, n = 2000 e n = 5000", "09_matrizes_grandes.jl"))
    etapa("8. " * rot)
    include(joinpath(@__DIR__, arq))
end

etapa("FF1: frente côncava")
include(joinpath(@__DIR__, "20_frente_nao_convexa_ff1.jl"))

etapa("Fim")
println("Resultados em ", abspath(pasta_r), " e figuras em ", abspath(pasta_f))
