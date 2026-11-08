# Redesenha as figuras do artigo a partir dos CSV de `resultados` (não roda nenhum
# experimento; leva poucos segundos). Grava em `figuras`, com os mesmos nomes usados no .tex:
#   Figura 1: frente_r2.pdf e frente_r3.pdf;
#   Figura 2: perfil_tabela4_grad.pdf e perfil_tabela4_tempo.pdf;
#   Figura 3: frentes_F6_FF1.pdf.
# Sem títulos dentro dos gráficos; cores da paleta de Okabe e Ito, distinguíveis também por
# daltônicos, e estilos de linha diferentes, de modo que as curvas se distinguem em preto e branco.
# Para a versão em inglês, troque `idioma = :pt` por `idioma = :en`.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, DelimitedFiles, Plots
gr()

idioma = :pt

const TEXTOS = Dict(
    :pt => Dict(:frac => "fração das execuções", :b1 => "blocos, 1ª ordem", :b2 => "blocos, 2ª ordem",
                :z1 => "só bloco 0", :z2 => "só bloco 0, 2ª ordem", :fs => "Fliege–Svaiter",
                :fsbb => "Fliege–Svaiter BB", :alg => "Algoritmo 1", :soma => "soma ponderada",
                :exata => "frente exata"),
    :en => Dict(:frac => "fraction of runs", :b1 => "blocks, first order", :b2 => "blocks, second order",
                :z1 => "block 0 only", :z2 => "block 0 only, second order", :fs => "Fliege–Svaiter",
                :fsbb => "Fliege–Svaiter BB", :alg => "Algorithm 1", :soma => "weighted sum",
                :exata => "exact front"),
)
tx(k) = TEXTOS[idioma][k]

# paleta de Okabe e Ito
const AZUL, VERMELHAO, VERDE, ROSA, LARANJA, PRETO =
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#000000"

rotulo(i) = "\$f_{$i}\$"
estilo = (guidefontsize=12, tickfontsize=10, legendfontsize=9, framestyle=:box, grid=false,
          size=(500, 400), left_margin=3Plots.mm, bottom_margin=3Plots.mm, right_margin=3Plots.mm)

pasta_res = joinpath(@__DIR__, "..", "resultados")
pasta_fig = joinpath(@__DIR__, "..", "figuras"); mkpath(pasta_fig)

"Lê um CSV com cabeçalho e devolve (dados, função que dá o índice de uma coluna pelo nome)."
function ler_csv(arquivo)
    dados, cab = readdlm(arquivo, ','; header=true)
    nomes = strip.(string.(vec(cab)))
    return dados, nome -> (c = findfirst(==(nome), nomes); c === nothing ? error("coluna $nome ausente em $arquivo") : c)
end

# ---------------- Figura 1: frentes da Família 1 ----------------
function ler_frente(arquivo, r)
    dados, col = ler_csv(arquivo)
    return permutedims(Float64.(dados[:, [col("f$i") for i in 1:r]]))     # r × N
end
function limites_cor(z)
    v = sort(collect(z)); m = length(v)
    lo = v[max(1, ceil(Int, 0.02 * m))]; hi = v[min(m, floor(Int, 0.98 * m) + 1)]
    return lo < hi ? (lo, hi) : (v[1], v[end])
end

F2 = ler_frente(joinpath(pasta_res, "frente_r2.csv"), 2); nd2 = nao_dominados(F2)
p = scatter(F2[1, nd2], F2[2, nd2]; ms=4, markerstrokewidth=0.3, color=AZUL, legend=false,
            xlabel=rotulo(1), ylabel=rotulo(2), estilo...)
savefig(p, joinpath(pasta_fig, "frente_r2.pdf"))

F3 = ler_frente(joinpath(pasta_res, "frente_r3.csv"), 3); nd3 = nao_dominados(F3)
p = scatter(F3[1, nd3], F3[2, nd3]; marker_z=F3[3, nd3], color=:viridis, clims=limites_cor(F3[3, nd3]),
            ms=4, markerstrokewidth=0.3, legend=false, colorbar=true,
            xlabel=rotulo(1), ylabel=rotulo(2), estilo...)
savefig(p, joinpath(pasta_fig, "frente_r3.pdf"))
println("Figura 1: ", length(nd2), " e ", length(nd3), " pontos não dominados")

# ---------------- Figura 2: perfis de desempenho da Tabela 3 ----------------
variantes = [("blocos, 1ª ordem", tx(:b1), AZUL, :solid),
             ("blocos, 2ª ordem", tx(:b2), VERMELHAO, :solid),
             ("só bloco 0", tx(:z1), VERDE, :dash),
             ("só bloco 0, 2ª ordem", tx(:z2), ROSA, :dash),
             ("Fliege–Svaiter", tx(:fs), LARANJA, :dot),
             ("Fliege–Svaiter BB", tx(:fsbb), PRETO, :dashdot)]

dados, col = ler_csv(joinpath(pasta_res, "tabela4_reais_execucoes.csv"))
ci, cv, ce, cs = col("instancia"), col("variante"), col("execucao"), col("status")
problemas = unique([(string(dados[l, ci]), string(dados[l, ce])) for l in 1:size(dados, 1)])
for (coluna, arquivo) in (("grad_equiv", "perfil_tabela4_grad.pdf"), ("tempo", "perfil_tabela4_tempo.pdf"))
    cm = col(coluna)
    Tm = fill(Inf, length(problemas), length(variantes))
    for l in 1:size(dados, 1)
        strip(string(dados[l, cs])) == "convergiu" || continue
        pr = findfirst(==((string(dados[l, ci]), string(dados[l, ce]))), problemas)
        s = findfirst(v -> v[1] == strip(string(dados[l, cv])), variantes)
        (pr === nothing || s === nothing) && continue
        Tm[pr, s] = Float64(dados[l, cm])
    end
    τ, ρ = perfil_desempenho(Tm)
    x = log2.(τ); xmax = maximum(x)
    fig = plot(; xlabel="\$\\log_2 \\tau\$", ylabel=tx(:frac), legend=:bottomright,
               xlims=(0, xmax), ylims=(0, 1.02), estilo...)
    for (s, (_, nome, cor, linha)) in enumerate(variantes)
        plot!(fig, vcat(x, xmax), vcat(ρ[:, s], ρ[end, s]); seriestype=:steppost,
              label=nome, color=cor, linestyle=linha, lw=2)
    end
    savefig(fig, joinpath(pasta_fig, arquivo))
end
println("Figura 2: ", length(problemas), " execuções")

# ---------------- Figura 3: FF1 ----------------
ts = range(-1, 1; length=400)
fig = plot(1 .- exp.(-(1 .- ts) .^ 2), 1 .- exp.(-(1 .+ ts) .^ 2); color=:gray60, lw=1.5,
           label=tx(:exata), xlabel=rotulo(1), ylabel=rotulo(2), legend=:bottomleft, estilo...)
# o Algoritmo 1 é desenhado primeiro, com círculos maiores e vazados, para continuar visível
# sob os pontos dos outros métodos, que praticamente coincidem com os dele
estilos_ff1 = [(1, tx(:alg),  (marker=:circle,    ms=7,   markercolor=:white, markerstrokecolor=AZUL, markerstrokewidth=1.3)),
               (2, tx(:z1),   (marker=:utriangle, ms=3.5, markercolor=VERDE,  markerstrokewidth=0)),
               (3, tx(:fsbb), (marker=:diamond,   ms=3,   markercolor=PRETO,  markerstrokewidth=0)),
               (4, tx(:soma), (marker=:xcross,    ms=7,   markercolor=VERMELHAO, markerstrokewidth=2))]
for (k, nome, est) in estilos_ff1
    F = ler_frente(joinpath(pasta_res, "frente_F6_$(k).csv"), 2); nd = nao_dominados(F)
    scatter!(fig, F[1, nd], F[2, nd]; label=nome, est...)
end
savefig(fig, joinpath(pasta_fig, "frentes_F6_FF1.pdf"))
println("Figura 3 gravada")

println("Figuras gravadas em ", abspath(pasta_fig))
