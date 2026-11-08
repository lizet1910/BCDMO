# Família 2: f_i(x) = ⟨a_i, x⟩ + (1/p)Σ|x_l|^p, com cada bloco "pertencente" a um
# objetivo. O parâmetro `sobreposicao` controla quanto os a_i se misturam entre
# blocos. Mede-se a fração de passos no bloco 0 em função da sobreposição.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, Printf, Plots

n, tam_bloco, r = 100, 5, 2
blocos = blocos_consecutivos(n, tam_bloco)
opc = Opcoes(ε=1e-2, maxiter=50_000)
N = 10

println("Cada nível de sobreposição roda N execuções; a tabela é impressa no final.")
niveis = [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0]
linhas = String[]
fracoes = Float64[]
medianas = Int[]
for sob in niveis
    local P, ests
    P = linear_potencia_blocos(blocos, r; p=1.5, sobreposicao=sob, semente=1)
    _, _, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=2, mostrar=true)
    frac = sum(e.passos_bloco0 for e in ests) / max(1, sum(e.iter for e in ests))
    its = sort([e.iter for e in ests])[cld(N, 2)]
    push!(fracoes, frac); push!(medianas, its)
    push!(linhas, @sprintf("%12.2f | %27.4f | %20d", sob, frac, its))
end

println("\nsobreposição | fração de passos no bloco 0 | mediana de iterações")
foreach(println, linhas)

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
p1 = plot(niveis, 100 .* fracoes; marker=:circle, label="", xlabel="sobreposição",
          ylabel="passos no bloco 0 (%)", ylims=(0, 110 * maximum(fracoes)))
p2 = plot(niveis, medianas; marker=:circle, label="", xlabel="sobreposição",
          ylabel="mediana de iterações", ylims=(0, 1.1 * maximum(medianas)))
fig = plot(p1, p2; layout=(1, 2), size=(900, 350), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
savefig(fig, joinpath(pasta, "estresse_bloco0.pdf"))
savefig(fig, joinpath(pasta, "estresse_bloco0.png"))
display(fig)
println("Figura gravada em ", abspath(pasta))
