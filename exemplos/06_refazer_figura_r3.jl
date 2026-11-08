# Refaz as figuras da frente com r = 3 a partir de resultados/frente_r3.csv,
# sem rodar o algoritmo de novo. Troque `paleta` para testar outras escalas
# de cores (por exemplo, :turbo, :plasma, :viridis, :jet).

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, DelimitedFiles
include(joinpath(@__DIR__, "graficos.jl"))

paleta = :turbo

pasta = joinpath(@__DIR__, "..", "resultados")
dados, cab = readdlm(joinpath(pasta, "frente_r3.csv"), ','; header=true)
cols = [findfirst(==("f$i"), vec(cab)) for i in 1:3]
F = permutedims(Float64.(dados[:, cols]))      # 3 × N

fig = plotar_frente(F; paleta=paleta)
savefig(fig, joinpath(pasta, "frente_r3.pdf"))
savefig(fig, joinpath(pasta, "frente_r3.png"))
fig3d = plotar_frente_3d(F; paleta=paleta)
savefig(fig3d, joinpath(pasta, "frente_r3_3d.pdf"))
savefig(fig3d, joinpath(pasta, "frente_r3_3d.png"))
display(fig)
println("Figuras gravadas em ", abspath(pasta))
