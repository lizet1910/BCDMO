# Refaz todos os experimentos do artigo com o teste de aceitação baseado em ε, que
# passou a ser o do Passo 6 (padrão de `Opcoes`). Antes, copia os resultados obtidos
# com o teste misto para `resultados_misto` e as figuras para `figuras_misto`.
#
# Ordem (do mais rápido ao mais lento; tudo é gravado ao fim de cada etapa):
#   01 exemplo da Observação 2.5; 07 Tabela 1; 04 bloco 0 na Família 2;
#   02, 03 e 11 frentes de Pareto e Figura 1; 05 Tabela 3; 08 e 09 Tabela 2;
#   12 e 13 Tabela 4; Tabela 5 (ε = 1e-3); variante só bloco 0 com 2ª ordem
#   (matrizes reais e, por último, as sintéticas, que levam cerca de 3 horas).
# Tempo total estimado: 7 a 8 horas.

raiz = joinpath(@__DIR__, "..")
for (origem, destino) in (("resultados", "resultados_misto"), ("figuras", "figuras_misto"))
    o = joinpath(raiz, origem); d = joinpath(raiz, destino); mkpath(d)
    isdir(o) || continue
    for f in readdir(o)
        isfile(joinpath(d, f)) || cp(joinpath(o, f), joinpath(d, f))
    end
end
println("Cópia dos resultados anteriores em resultados_misto e figuras_misto.")

etapa(nome) = println("\n", "#"^70, "\n# ", nome, "  (", Libc.strftime("%H:%M", time()), ")\n", "#"^70)

etapa("01: exemplo da Observação 2.5");          include(joinpath(@__DIR__, "01_exemplo_observacao.jl"))
etapa("07: Tabela 1");                           include(joinpath(@__DIR__, "07_comparacao_bloco0.jl"))
etapa("04: uso do bloco 0 na Família 2");        include(joinpath(@__DIR__, "04_estresse_bloco0.jl"))
etapa("02: frente com r = 2");                   include(joinpath(@__DIR__, "02_frente_r2.jl"))
etapa("03: frente com r = 3");                   include(joinpath(@__DIR__, "03_frente_r3.jl"))
etapa("11: Figura 1");                           include(joinpath(@__DIR__, "11_figuras_artigo.jl"))
etapa("05: Tabela 3");                           include(joinpath(@__DIR__, "05_matrizes_reais.jl"))
etapa("08: Tabela 2, n = 50 e n = 1000");        include(joinpath(@__DIR__, "08_escalas_heterogeneas.jl"))
etapa("09: Tabela 2, n = 2000 e n = 5000");      include(joinpath(@__DIR__, "09_matrizes_grandes.jl"))
etapa("12: Tabela 4, matrizes reais");           include(joinpath(@__DIR__, "12_comparacao_metodos.jl"))
etapa("13: Tabela 4, sintéticas");               include(joinpath(@__DIR__, "13_comparacao_sinteticas.jl"))

# ---------------- Tabela 5 (ε = 1e-3) e variante só bloco 0 com 2ª ordem ----------------
include(joinpath(@__DIR__, "comum_comparacao.jl"))
instancias_todas = [
    ("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
    ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
    ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10"),
]
pasta_r = joinpath(raiz, "resultados")
cab_r = ["instancia" "m" "n" "c" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
function gravar_r(arquivo, tabela)
    open(joinpath(pasta_r, arquivo), "w") do io
        writedlm(io, cab_r, ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tabela)
    end
end

etapa("Tabela 5: matrizes reais com ε = 1e-3")
tab5 = Any[]
for (grupo, nome) in instancias_todas
    local A, blocos
    A = obter_matriz(grupo, nome)
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    println("\n", nome, "  (", m, " × ", n, "), ε = 1e-3")
    cabecalho()
    for l in rodar(A, blocos, 5, 1e-3)
        push!(tab5, Any[nome, m, n, "", l...])
    end
    gravar_r("tabela5_eps1e-3.csv", tab5)
end

etapa("Tabela 4: só bloco 0 com 2ª ordem, matrizes reais")
tabQ3 = Any[]
for (grupo, nome) in instancias_todas
    local A, blocos
    A = obter_matriz(grupo, nome)
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    println("\n", nome, "  (", m, " × ", n, ")")
    cabecalho()
    for l in rodar(A, blocos, 5, 1e-2; vars=[variante_b0_hess])
        push!(tabQ3, Any[nome, m, n, "", l...])
    end
    gravar_r("revisao_Q3.csv", tabQ3)
end

etapa("Tabela 4: só bloco 0 com 2ª ordem, sintéticas (cerca de 3 horas)")
A0s = matriz_sintetica(6000, 2000; densidade=0.002, semente=1)
blocos_s = blocos_consecutivos(2000, 100)
for c in [0.0, 1.0]
    local A, esc
    A, esc = escalar_blocos(A0s, blocos_s, c)
    println("\nSintética n = 2000, c = ", c)
    cabecalho()
    for l in rodar(A, blocos_s, 5, 1e-2; escala=1 ./ esc, vars=[variante_b0_hess])
        push!(tabQ3, Any["sintetica", 6000, 2000, c, l...])
    end
    gravar_r("revisao_Q3.csv", tabQ3)
end

etapa("Fim")
println("Todos os resultados com o teste baseado em ε estão em ", abspath(pasta_r))
