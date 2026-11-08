# Família 1 com matrizes reais da coleção SuiteSparse, as mesmas usadas em [A].
# O script baixa cada matriz (se ainda não estiver na pasta `dados`), descompacta
# com o comando `tar` do Windows e compara blocos com σ por bloco e só o bloco 0.
# Como em [A], os blocos têm cerca de 10% das variáveis (q = 10).

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, Printf, DelimitedFiles, Downloads

instancias = [
    ("JGD_Homology", "ch6-6-b3"),    # 5400 × 2400
    ("JGD_Homology", "ch7-8-b2"),    # 11760 × 1176
    ("JGD_Homology", "mk10-b3"),     # 4725 × 3150
    ("JGD_Franz",    "Franz5"),      # 7382 × 2882
]
N = 5             # pontos iniciais por instância
ε = 1e-2

mediana(v) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

dados = joinpath(@__DIR__, "..", "dados"); mkpath(dados)

function obter_matriz(grupo, nome)
    mtx = joinpath(dados, nome, nome * ".mtx")
    if !isfile(mtx)
        arq = joinpath(dados, nome * ".tar.gz")
        if !isfile(arq)
            # tenta o endereço antigo da Universidade da Flórida e, se falhar, o atual
            urls = ["https://www.cise.ufl.edu/research/sparse/MM/$grupo/$nome.tar.gz",
                    "https://suitesparse-collection-website.herokuapp.com/MM/$grupo/$nome.tar.gz"]
            for url in urls
                println("baixando ", url)
                try
                    Downloads.download(url, arq)
                    break
                catch e
                    println("  falhou: ", e)
                    isfile(arq) && rm(arq)
                end
            end
            isfile(arq) || error("nenhum endereço funcionou para $nome")
        end
        run(`tar -xzf $arq -C $dados`)
    end
    return ler_matriz(mtx)
end

variantes = [
    ("blocos, σ por bloco", (; sigma=:porbloco)),
    ("só bloco 0",          (; sigma=:decresce, blocos_proprios=false)),
]

cab = ["instancia" "m" "n" "variante" "convergiram" "iter" "grad_equiv" "passos_bloco0" "tempo"]
tabela = Any[]
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
for (grupo, nome) in instancias
    local A, P, blocos
    A = try
        obter_matriz(grupo, nome)
    catch e
        println("não foi possível obter ", nome, ": ", e)
        continue
    end
    m, n = size(A)
    blocos = blocos_consecutivos(n, cld(n, 10))
    P = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
    println("\n", nome, "  (", m, " × ", n, ", ", length(blocos), " blocos)")
    @printf("  %-20s %6s %9s %11s %9s %9s\n", "variante", "conv.", "iter", "grad. equiv", "bloco 0", "tempo(s)")
    for (rot, kw) in variantes
        local ests
        opc = Opcoes(; ε=ε, maxiter=300_000, kw...)
        _, _, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3, mostrar=false)
        conv = count(e -> e.status == :convergiu, ests)
        lin = (mediana([e.iter for e in ests]), mediana([e.grad_equiv for e in ests]),
               mediana([e.passos_bloco0 for e in ests]), mediana([e.tempo for e in ests]))
        @printf("  %-20s %3d/%-2d %9.0f %11.1f %9.0f %9.3f\n", rot, conv, N, lin...)
        push!(tabela, Any[nome, m, n, rot, conv, lin...])
        # grava a tabela a cada linha, para não perder resultados se for interrompido
        open(joinpath(pasta, "matrizes_reais.csv"), "w") do io
            writedlm(io, cab, ',')
            foreach(l -> writedlm(io, permutedims(l), ','), tabela)
        end
    end
end
println("\nMedianas sobre os pontos iniciais. Tabela gravada em ", abspath(pasta))
