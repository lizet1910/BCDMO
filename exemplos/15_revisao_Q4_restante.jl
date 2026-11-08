# Continuação da parte Q4 do script 14 (ε = 1e-3), sem a variante "só bloco 0, 2ª ordem",
# cuja matriz densa de ordem n torna cada iteração cara demais com essa tolerância.
# Roda as quatro variantes da Tabela 4 nas matrizes que faltaram (Franz5 a Franz10) e
# grava em resultados/revisao_Q4b.csv, sem apagar o revisao_Q4.csv do script 14.

include(joinpath(@__DIR__, "comum_comparacao.jl"))

N = 5
ε = 1e-3
instancias = [
    ("JGD_Franz", "Franz5"),
    ("JGD_Franz", "Franz6"),
    ("JGD_Franz", "Franz7"),
    ("JGD_Franz", "Franz8"),
    ("JGD_Franz", "Franz9"),
    ("JGD_Franz", "Franz10"),
]
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
cab = ["instancia" "m" "n" "c" "variante" "convergiram" "iter" "avalF" "grad_equiv" "passos_bloco0" "tempo"]
tabela = Any[]
for (grupo, nome) in instancias
    local A, blocos
    A = obter_matriz(grupo, nome)
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    println("\n", nome, "  (", m, " × ", n, "), ε = ", ε)
    cabecalho()
    for l in rodar(A, blocos, N, ε)
        push!(tabela, Any[nome, m, n, "", l...])
    end
    open(joinpath(pasta, "revisao_Q4b.csv"), "w") do io
        writedlm(io, cab, ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tabela)
    end
end
println("\nMedianas sobre ", N, " pontos iniciais. Tabela gravada em ", abspath(pasta))
