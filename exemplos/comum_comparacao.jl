# Funções comuns aos scripts 12 e 13: obtenção das matrizes da SuiteSparse,
# variantes comparadas e execução de uma instância.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, Printf, DelimitedFiles, Downloads

mediana(v) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

dados = joinpath(@__DIR__, "..", "dados"); mkpath(dados)

function obter_matriz(grupo, nome)
    mtx = joinpath(dados, nome, nome * ".mtx")
    if !isfile(mtx)
        arq = joinpath(dados, nome * ".tar.gz")
        if !isfile(arq)
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

# Modelos usados pelas variantes:
#   :linear      B ≡ 0;
#   :hess        matriz comum (Hessianas truncadas) nos blocos próprios e B = 0 no bloco 0
#                quando n > 2000 (como na Tabela 4);
#   :hess_total  a mesma matriz também no bloco 0, de ordem n (usada só com n ≤ N_HESS_TOTAL).
const N_HESS_TOTAL = 5000

# variantes: (rótulo, modelo, método, opções)
variantes = [
    ("blocos, 1ª ordem", :linear, :bcdmo, (; sigma=:porbloco)),
    ("blocos, 2ª ordem", :hess,   :bcdmo, (; sigma=:porbloco)),
    ("só bloco 0",       :linear, :bcdmo, (; sigma=:decresce, blocos_proprios=false)),
    ("Fliege–Svaiter",   :linear, :fs,    (;)),
]
variante_b0_hess = ("só bloco 0, 2ª ordem", :hess_total, :bcdmo, (; sigma=:decresce, blocos_proprios=false))
variante_fsbb    = ("Fliege–Svaiter BB",    :linear,     :fsbb,  (;))

# construcao = :a (mesma matriz A para os dois objetivos) ou :b (linhas de A repartidas
# em dois grupos); perda = :potencia (convexa) ou :cauchy (não convexa)
function problema(A, modelo; semente_b=2, construcao=:a, perda=:potencia)
    if modelo == :linear
        kw = (; perda=perda)
    else
        kw = (; perda=perda, usar_hess=true, hess_comum=true, hess_max=(modelo == :hess ? 2000 : N_HESS_TOTAL))
    end
    construcao == :a && return regressao_mesma_matriz(A, [1.5, 1.5]; semente=semente_b, kw...)
    return regressao_linhas(A, [1.5, 1.5]; semente=semente_b, kw...)
end

# Roda as variantes `vars` a partir de N pontos iniciais e devolve uma linha de medianas
# por variante. Se `registro` for um vetor, acrescenta nele uma linha por execução:
# [rotulo, variante, k, status, iter, avalF, grad_equiv, duplicacoes, duplicacoes_p5, tempo].
function rodar(A, blocos, N, ε; semente_b=2, escala=1.0, vars=variantes, construcao=:a,
               perda=:potencia, registro=nothing, rotulo="")
    n = size(A, 2)
    probs = Dict{Symbol,Any}()
    linhas = []
    for (rot, modelo, met, kw) in vars
        if modelo == :hess_total && n > N_HESS_TOTAL
            @printf("  %-22s não rodada (n > %d)\n", rot, N_HESS_TOTAL)
            push!(linhas, Any[rot, "", "", "", "", "", ""])
            if registro !== nothing
                foreach(k -> push!(registro, Any[rotulo, rot, k, "nao_rodada", "", "", "", "", "", ""]), 1:N)
            end
            continue
        end
        P = get!(() -> problema(A, modelo; semente_b=semente_b, construcao=construcao, perda=perda), probs, modelo)
        opc = Opcoes(; ε=ε, maxiter=300_000, maxaval=3_000_000, kw...)
        _, F, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3,
                                       escala=escala, mostrar=false, metodo=met)
        conv = count(e -> e.status == :convergiu, ests)
        lin = (mediana([e.iter for e in ests]), mediana([e.avalF for e in ests]),
               mediana([e.grad_equiv for e in ests]), mediana([e.passos_bloco0 for e in ests]),
               mediana([e.tempo for e in ests]))
        @printf("  %-22s %3d/%-2d %9.0f %9.0f %11.1f %9.0f %9.3f\n", rot, conv, N, lin...)
        push!(linhas, Any[rot, conv, lin...])
        if registro !== nothing
            for (k, e) in enumerate(ests)
                push!(registro, Any[rotulo, rot, k, string(e.status), e.iter, e.avalF, e.grad_equiv,
                                    e.duplicacoes, e.duplicacoes_p5, e.tempo])
            end
        end
    end
    return linhas
end

cabecalho() = @printf("  %-22s %6s %9s %9s %11s %9s %9s\n",
                      "variante", "conv.", "iter", "aval. F", "grad. equiv", "bloco 0", "tempo(s)")

