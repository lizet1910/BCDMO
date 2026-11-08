# Figuras finais para o artigo, redesenhadas a partir dos CSV de `resultados`
# (não roda o algoritmo de novo). Sem título, com rótulos f_1, f_2, f_3 e fontes
# maiores, no tamanho adequado para meia largura de texto. Gravadas em `figuras`.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, DelimitedFiles, Plots

# rótulo f_i: o GR interpreta textos entre cifrões como fórmulas do LaTeX
rotulo(i) = "\$f_{$i}\$"

function ler_frente(arquivo, r)
    dados, cab = readdlm(arquivo, ','; header=true)
    cols = [findfirst(==("f$i"), vec(cab)) for i in 1:r]
    return permutedims(Float64.(dados[:, cols]))       # r × N
end

function limites_cor(z)
    v = sort(collect(z)); m = length(v)
    lo = v[max(1, ceil(Int, 0.02 * m))]; hi = v[min(m, floor(Int, 0.98 * m) + 1)]
    return lo < hi ? (lo, hi) : (v[1], v[end])
end

estilo = (guidefontsize=12, tickfontsize=10, legend=false, framestyle=:box,
          size=(500, 400), left_margin=3Plots.mm, bottom_margin=3Plots.mm)

pasta_res = joinpath(@__DIR__, "..", "resultados")
pasta_fig = joinpath(@__DIR__, "..", "figuras"); mkpath(pasta_fig)

# r = 2
F2 = ler_frente(joinpath(pasta_res, "frente_r2.csv"), 2)
nd2 = nao_dominados(F2)
p2 = scatter(F2[1, nd2], F2[2, nd2]; ms=4, markerstrokewidth=0.3, color=:steelblue,
             xlabel=rotulo(1), ylabel=rotulo(2), estilo...)
savefig(p2, joinpath(pasta_fig, "frente_r2.pdf"))
println("r = 2: ", length(nd2), " não dominados de ", size(F2, 2))

# r = 3: plano (f_1, f_2), cor indicando f_3
F3 = ler_frente(joinpath(pasta_res, "frente_r3.csv"), 3)
nd3 = nao_dominados(F3)
p3 = scatter(F3[1, nd3], F3[2, nd3]; marker_z=F3[3, nd3], color=:turbo,
             clims=limites_cor(F3[3, nd3]), ms=4, markerstrokewidth=0.3,
             xlabel=rotulo(1), ylabel=rotulo(2), colorbar=true,
             right_margin=4Plots.mm, estilo...)   # a cor indica f_3 (dito na legenda do artigo)
savefig(p3, joinpath(pasta_fig, "frente_r3.pdf"))
println("r = 3: ", length(nd3), " não dominados de ", size(F3, 2))

# r = 3 em 3D (mesma escala de cores da projeção), em dois ângulos de visão
faixa(v) = (minimum(v) - 0.03 * (maximum(v) - minimum(v)), maximum(v) + 0.03 * (maximum(v) - minimum(v)))
for (nome, cam) in (("frente_r3_3d_a.pdf", (30, 30)), ("frente_r3_3d_b.pdf", (60, 20)))
    p3d = scatter3d(F3[1, nd3], F3[2, nd3], F3[3, nd3]; marker_z=F3[3, nd3], color=:turbo,
                    clims=limites_cor(F3[3, nd3]), ms=3, markerstrokewidth=0.2,
                    xlabel=rotulo(1), ylabel=rotulo(2), zlabel=rotulo(3), camera=cam,
                    xlims=faixa(F3[1, nd3]), ylims=faixa(F3[2, nd3]), zlims=faixa(F3[3, nd3]),
                    colorbar=false, guidefontsize=12, tickfontsize=9, legend=false, size=(500, 400))
    savefig(p3d, joinpath(pasta_fig, nome))
end

println("Figuras gravadas em ", abspath(pasta_fig))
