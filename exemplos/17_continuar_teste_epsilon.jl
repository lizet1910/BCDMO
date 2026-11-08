# Continuação do script 16, a partir da Tabela 4, para quando a rodada anterior
# foi interrompida (por exemplo, porque o computador entrou em suspensão).
# Refaz, com o teste baseado em ε: Tabela 4 (scripts 12 e 13), Tabela 5 (ε = 1e-3)
# e a variante só bloco 0 com 2ª ordem nas matrizes reais. Não refaz as sintéticas
# dessa variante (cerca de 3 horas). Tempo estimado: 1h30 a 2h.
# Deixe o computador ligado na tomada e sem suspensão enquanto roda.

raiz = joinpath(@__DIR__, "..")
etapa(nome) = println("\n", "#"^70, "\n# ", nome, "  (", Libc.strftime("%H:%M", time()), ")\n", "#"^70)

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

etapa("Fim")
println("Resultados gravados em ", abspath(pasta_r))
