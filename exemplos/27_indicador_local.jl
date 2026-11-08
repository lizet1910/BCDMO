# Indicador L/L̄ global e local (cerca de 5 a 10 minutos). Rodar depois do script 26.
#   1. L/L̄ da construção (b), com as matrizes A_1 e A_2 usadas nos experimentos;
#   2. indicador local, com a matriz W(x)^{1/2} A_i, W(x) = diag(min{1, |ϱ''((A_i x − b_i)_l)|}),
#      isto é, com os pesos da Hessiana de f_i truncados, no ponto inicial e no ponto final de
#      uma execução do Algoritmo 1 (1ª ordem), para a perda convexa nas construções (a) e (b)
#      e para a perda de Cauchy na construção (b); média sobre os objetivos;
#   3. gráfico de dispersão da razão "gradientes equivalentes, blocos / só bloco 0" (ambos com
#      modelo linear) contra L/L̄, para todas as instâncias das Tabelas 2, 3 e 6.
# Grava resultados/indicador_local.csv e figuras/dispersao_indicador_global.pdf e
# figuras/dispersao_indicador_local.pdf.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra, Plots
gr()

pasta_r = joinpath(@__DIR__, "..", "resultados"); pasta_f = joinpath(@__DIR__, "..", "figuras")

"Média, sobre os objetivos, do indicador com pesos da Hessiana truncada em x (ou sem pesos, se x = nothing)."
function indicador_medio(P::RegressaoLp, x, bl)
    v = 0.0
    for i in 1:nobj(P)
        Ai = P.A[i]
        if x === nothing
            M = Ai
        else
            r = Ai * x .- P.b[i]
            w = abs.(BCDMO._peso_hess.(r, P.p[i], 1.0, P.perda))
            M = Diagonal(sqrt.(w)) * Ai
        end
        v += indicador_blocos(M, bl)
    end
    return v / nobj(P)
end

instancias = Any[]
for (g, nome) in (("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
                  ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
                  ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10"))
    A = obter_matriz(g, nome)
    push!(instancias, (nome, A, blocos_consecutivos(size(A, 2), cld(size(A, 2), 10))))
end
let A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1)
    push!(instancias, ("sintetica", A0, blocos_consecutivos(2000, 100)))
end

cab = ["instancia" "global_a" "global_b" "local_a_x0" "local_a_final" "local_b_x0" "local_b_final" "local_cauchy_x0" "local_cauchy_final"]
tab = Any[]
opc = Opcoes(ε=1e-2, maxiter=300_000)
println("instância      global(a) global(b) | local (a): x0  final | local (b): x0  final | Cauchy: x0  final")
for (nome, A, bl) in instancias
    n = size(A, 2)
    x0 = randn(MersenneTwister(3), n)            # primeiro ponto inicial dos experimentos
    Pa = problema(A, :linear)
    Pb = problema(A, :linear; construcao=:b)
    Pc = problema(A, :linear; construcao=:b, perda=:cauchy)
    xa, _, _ = bcdmo(Pa, x0, bl; opc=opc)
    xb, _, _ = bcdmo(Pb, x0, bl; opc=opc)
    xc, _, _ = bcdmo(Pc, x0, bl; opc=opc)
    lin = Any[nome, indicador_blocos(A, bl), indicador_medio(Pb, nothing, bl),
              indicador_medio(Pa, x0, bl), indicador_medio(Pa, xa, bl),
              indicador_medio(Pb, x0, bl), indicador_medio(Pb, xb, bl),
              indicador_medio(Pc, x0, bl), indicador_medio(Pc, xc, bl)]
    @printf("%-13s %8.2f %8.2f | %10.2f %6.2f | %10.2f %6.2f | %8.2f %6.2f\n", lin...)
    push!(tab, lin)
    open(joinpath(pasta_r, "indicador_local.csv"), "w") do io
        writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
    end
end

# ---------------- gráfico de dispersão ----------------
function ler(arq)
    d, c = readdlm(joinpath(pasta_r, arq), ','; header=true)
    nomes = strip.(string.(vec(c)))
    return d, k -> findfirst(==(k), nomes)
end
"Razões gradientes(blocos)/gradientes(só bloco 0) por instância, a partir de um CSV de medianas."
function razoes(arq, rot_blocos, rot_b0; chave=(d, col, l) -> string(d[l, col("instancia")]))
    d, col = ler(arq)
    g = Dict{String,Dict{String,Float64}}()
    for l in 1:size(d, 1)
        v = string(d[l, col("variante")]); x = d[l, col("grad_equiv")]
        x isa Number || continue
        get!(g, chave(d, col, l), Dict{String,Float64}())[v] = Float64(x)
    end
    return Dict(k => v[rot_blocos] / v[rot_b0] for (k, v) in g if haskey(v, rot_blocos) && haskey(v, rot_b0))
end
ind = Dict(string(l[1]) => l for l in tab)

# sintéticas convexas (Tabela 2): indicador do script 18a
dI, cI = ler("indicador_blocos.csv")
ind_sint = Dict((Int(dI[l, cI("n")]), Float64(dI[l, cI("c")])) => Float64(dI[l, cI("L_sobre_media_Lj")])
                for l in 1:size(dI, 1) if string(dI[l, cI("matriz")]) == "sintetica")
pts = Dict{String,Vector{Tuple{Float64,Float64}}}()
for arq in ("escalas_heterogeneas.csv", "matrizes_grandes.csv")
    d, col = ler(arq)
    g = Dict{Tuple{Int,Float64},Dict{String,Float64}}()
    for l in 1:size(d, 1)
        nn = parse(Int, split(string(d[l, col("cenario")]), "= ")[2]); c = Float64(d[l, col("c")])
        get!(g, (nn, c), Dict{String,Float64}())[string(d[l, col("variante")])] = Float64(d[l, col("grad_equiv")])
    end
    for ((nn, c), v) in g
        fam = nn == 50 ? "sintéticas, n = 50" : "sintéticas, n ≥ 1000"
        push!(get!(pts, fam, Tuple{Float64,Float64}[]), (ind_sint[(nn, c)], v["blocos, σ por bloco"] / v["só bloco 0"]))
    end
end
ra = razoes("tabela4_reais.csv", "blocos, 1ª ordem", "só bloco 0")
rb = razoes("construcao_b_reais.csv", "blocos, 1ª ordem", "só bloco 0")
rc = razoes("nao_convexo_reais.csv", "blocos, 1ª ordem", "só bloco 0")
rs = razoes("nao_convexo_sinteticas.csv", "blocos, 1ª ordem", "só bloco 0")
pts["reais, construção (a)"] = [(ind[k][2], v) for (k, v) in ra]
pts["reais, construção (b)"] = [(ind[k][3], v) for (k, v) in rb]
pts["Cauchy, construção (b)"] = [(ind[k][3], v) for (k, v) in rc]
haskey(rs, "sintetica") && push!(pts["Cauchy, construção (b)"], (ind["sintetica"][3], rs["sintetica"]))
pts_local = deepcopy(pts)
pts_local["reais, construção (a)"] = [(ind[k][5], v) for (k, v) in ra]
pts_local["reais, construção (b)"] = [(ind[k][7], v) for (k, v) in rb]
pts_local["Cauchy, construção (b)"] = [(ind[k][9], v) for (k, v) in rc]
haskey(rs, "sintetica") && push!(pts_local["Cauchy, construção (b)"], (ind["sintetica"][9], rs["sintetica"]))

estilos = [("sintéticas, n = 50", "#999999", :circle), ("sintéticas, n ≥ 1000", "#0072B2", :circle),
           ("reais, construção (a)", "#009E73", :diamond), ("reais, construção (b)", "#E69F00", :utriangle),
           ("Cauchy, construção (b)", "#D55E00", :xcross)]
for (P_, arq, xl) in ((pts, "dispersao_indicador_global.pdf", "\$L/\\bar L\$ (matriz A)"),
                      (pts_local, "dispersao_indicador_local.pdf", "\$L/\\bar L\$ local (ponto final)"))
    fig = plot(; xlabel=xl, ylabel="gradientes: blocos / só bloco 0", yscale=:log10, legend=:topright,
               framestyle=:box, grid=false, size=(560, 420), guidefontsize=12, tickfontsize=10, legendfontsize=9)
    hline!(fig, [1.0]; color=:gray60, lw=1, label="")
    for (fam, cor, m) in estilos
        haskey(P_, fam) || continue
        v = P_[fam]
        scatter!(fig, first.(v), last.(v); label=fam, color=cor, marker=m, ms=6, markerstrokewidth=0.5)
    end
    savefig(fig, joinpath(pasta_f, arq))
end
println("\nResultados em ", abspath(joinpath(pasta_r, "indicador_local.csv")), "; figuras em ", abspath(pasta_f))
