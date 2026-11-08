# Experimentos pedidos na revisão:
#   Q3: variante que usa apenas o bloco 0 com modelo de segunda ordem (matriz comum de
#       ordem n), a comparação justa para a coluna "2ª ordem" da Tabela 4; roda com
#       ε = 1e-2 nas matrizes reais com n ≤ 5000 e nas sintéticas com n = 2000.
#   Q4: as cinco variantes nas matrizes reais com ε = 1e-3.
# As tabelas são gravadas a cada instância em resultados/revisao_Q3.csv e
# resultados/revisao_Q4.csv. A parte Q4 é a mais demorada.

include(joinpath(@__DIR__, "comum_comparacao.jl"))

N = 5
instancias = [
    ("JGD_Homology", "ch6-6-b3"),
    ("JGD_Homology", "ch7-8-b2"),
    ("JGD_Homology", "mk10-b3"),
    ("JGD_Franz", "Franz5"),
    ("JGD_Franz", "Franz6"),
    ("JGD_Franz", "Franz7"),
    ("JGD_Franz", "Franz8"),
    ("JGD_Franz", "Franz9"),
    ("JGD_Franz", "Franz10"),
]
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
cab = ["instancia" "m" "n" "c" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]

function gravar(arquivo, tabela)
    open(joinpath(pasta, arquivo), "w") do io
        writedlm(io, cab, ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tabela)
    end
end

# verificação rápida (alguns segundos)
println("Verificação rápida (200 × 50, 5 blocos)")
cabecalho()
rodar(matriz_sintetica(200, 50; densidade=0.1, semente=1), blocos_consecutivos(50, 10), 2, 1e-2;
      vars=[variante_b0_hess])

# ---------------- Q3: só bloco 0 com segunda ordem, ε = 1e-2 ----------------
println("\n========== Q3: só bloco 0 com modelo de segunda ordem, ε = 1e-2 ==========")
tabQ3 = Any[]
for (grupo, nome) in instancias
    local A, blocos
    A = obter_matriz(grupo, nome)
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    println("\n", nome, "  (", m, " × ", n, ")")
    cabecalho()
    for l in rodar(A, blocos, N, 1e-2; vars=[variante_b0_hess])
        push!(tabQ3, Any[nome, m, n, "", l...])
    end
    gravar("revisao_Q3.csv", tabQ3)
end
A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1)
blocos_s = blocos_consecutivos(2000, 100)
for c in [0.0, 1.0]
    local A, esc
    A, esc = escalar_blocos(A0, blocos_s, c)
    println("\nSintética n = 2000, c = ", c)
    cabecalho()
    for l in rodar(A, blocos_s, N, 1e-2; escala=1 ./ esc, vars=[variante_b0_hess])
        push!(tabQ3, Any["sintetica", 6000, 2000, c, l...])
    end
    gravar("revisao_Q3.csv", tabQ3)
end

# ---------------- Q4: as cinco variantes, ε = 1e-3 ----------------
println("\n========== Q4: matrizes reais com ε = 1e-3 ==========")
tabQ4 = Any[]
for (grupo, nome) in instancias
    local A, blocos
    A = obter_matriz(grupo, nome)
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    println("\n", nome, "  (", m, " × ", n, ")")
    cabecalho()
    for l in rodar(A, blocos, N, 1e-3; vars=[variantes; [variante_b0_hess]])
        push!(tabQ4, Any[nome, m, n, "", l...])
    end
    gravar("revisao_Q4.csv", tabQ4)
end
println("\nMedianas sobre ", N, " pontos iniciais. Tabelas gravadas em ", abspath(pasta))
