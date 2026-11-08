# Comparação com outros métodos nas matrizes reais da coleção SuiteSparse:
#   (a) blocos, σ por bloco, modelo linear (B = 0), como nas tabelas anteriores;
#   (b) blocos, σ por bloco, modelo de segunda ordem com matriz comum
#       B_(j) = (1/r) Σ_i A_j' D̃_i A_j (Hessianas truncadas, como em [A]);
#       no bloco 0 usa-se B = 0, o que a Hipótese 2 permite;
#   (c) só o bloco 0, modelo linear;
#   (d) máxima descida de Fliege e Svaiter com busca de Armijo.
# Todos param com o mesmo critério G(x) ≤ ε e partem dos mesmos pontos iniciais.
# Antes das matrizes reais, uma verificação rápida num problema pequeno.
# As matrizes que faltarem em `dados` são baixadas e descompactadas.

include(joinpath(@__DIR__, "comum_comparacao.jl"))

N = 5             # pontos iniciais por instância
ε = 1e-2

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

# verificação rápida (alguns segundos): todas as variantes devem convergir
println("Verificação rápida (200 × 50, 5 blocos)")
cabecalho()
rodar(matriz_sintetica(200, 50; densidade=0.1, semente=1), blocos_consecutivos(50, 10), 2, ε)

cab = ["instancia" "m" "n" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
tabela = Any[]
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
for (grupo, nome) in instancias
    local A, blocos
    A = try
        obter_matriz(grupo, nome)
    catch e
        println("não foi possível obter ", nome, ": ", e)
        continue
    end
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    println("\n", nome, "  (", m, " × ", n, ", ", length(blocos), " blocos)")
    cabecalho()
    for l in rodar(A, blocos, N, ε)
        push!(tabela, Any[nome, m, n, l...])
    end
    # grava a tabela a cada instância, para não perder resultados se for interrompido
    open(joinpath(pasta, "comparacao_metodos.csv"), "w") do io
        writedlm(io, cab, ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tabela)
    end
end
println("\nMedianas sobre os pontos iniciais. Tabela gravada em ", abspath(pasta))
