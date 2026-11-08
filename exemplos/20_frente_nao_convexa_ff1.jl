# Frente de Pareto não convexa (côncava): problema FF1 de Fonseca e Fleming,
#   f_1(x) = 1 − exp(−‖x − a‖²),  f_2(x) = 1 − exp(−‖x + a‖²),  a = (1, …, 1)/√n.
# O conjunto de Pareto é o segmento {t a : t ∈ [−1, 1]}, e a frente,
#   (1 − exp(−(1 − t)²), 1 − exp(−(1 + t)²)),  t ∈ [−1, 1],
# é côncava: a soma ponderada só alcança os extremos. Os gradientes são Lipschitz
# (Hipótese 1 com β = 1) e f_i ≥ 0 (Hipótese 3).
# Compara Algoritmo 1, só bloco 0, Fliege–Svaiter BB e soma ponderada, com n = 1000,
# blocos de 50 variáveis e 100 execuções por método. Cerca de 10 a 20 minutos.
# Pode ser rodado na mesma sessão do script 19 ou numa sessão nova.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra, Plots

struct FF1 <: BCDMO.ProblemaMO
    n::Int
    a::Vector{Float64}
end
FF1(n::Int) = FF1(n, fill(1 / sqrt(n), n))
BCDMO.nobj(::FF1) = 2
BCDMO.dim(P::FF1) = P.n
BCDMO.fvals(P::FF1, x) = [1 - exp(-sum(abs2, x .- P.a)), 1 - exp(-sum(abs2, x .+ P.a))]
function BCDMO.jac_bloco(P::FF1, x, idx)
    e1 = exp(-sum(abs2, x .- P.a)); e2 = exp(-sum(abs2, x .+ P.a))
    return hcat(2 .* (x[idx] .- P.a[idx]) .* e1, 2 .* (x[idx] .+ P.a[idx]) .* e2)
end

n = 1000; ε = 1e-3; N = 100
P = FF1(n); bl = blocos_consecutivos(n, 50)
escala = 1 / sqrt(n)                  # pontos iniciais com ‖x0‖ ≈ 1, perto do conjunto de Pareto
pesos = [[t, 1 - t] for t in range(0, 1; length=N)]
metodos = ["Algoritmo 1 (blocos)", "só bloco 0", "Fliege–Svaiter BB", "soma ponderada"]

res = Dict{String,Any}()
for (rot, kw, met) in (("Algoritmo 1 (blocos)", (; sigma=:porbloco), :bcdmo),
                       ("só bloco 0", (; sigma=:decresce, blocos_proprios=false), :bcdmo),
                       ("Fliege–Svaiter BB", (;), :fsbb))
    _, F, ests = multiplos_inicios(P, bl, N; opc=Opcoes(; ε=ε, maxiter=300_000, kw...), semente=3,
                                   escala=escala, mostrar=false, metodo=met)
    res[rot] = (F, sum(e.grad_equiv for e in ests), sum(e.tempo for e in ests),
                count(e -> e.status == :convergiu, ests))
end
let rng = MersenneTwister(3), F = zeros(2, N), ge = 0.0, tm = 0.0, cv = 0
    for k in 1:N
        x0 = escala .* randn(rng, n)
        x, _, est = bcdmo(Escalarizado(P, pesos[k]), x0, bl; opc=Opcoes(ε=ε, maxiter=300_000))
        F[:, k] = fvals(P, x); ge += est.grad_equiv; tm += est.tempo; cv += est.status == :convergiu
    end
    res["soma ponderada"] = (F, ge, tm, cv)
end

ms = metricas_frentes([res[m][1] for m in metodos])
pasta_r = joinpath(@__DIR__, "..", "resultados"); pasta_f = joinpath(@__DIR__, "..", "figuras")
println("\nFF1, n = $n, $(length(bl)) blocos, ε = $ε, $N execuções por método")
@printf("  %-22s %6s %7s %7s %7s %9s %11s %9s\n", "método", "conv.", "pureza", "Γ", "Δ", "hipervol.", "grad. equiv", "tempo(s)")
tab = Any[]
for (k, m) in enumerate(metodos)
    F, ge, tm, cv = res[m]
    @printf("  %-22s %6d %7.3f %7.3f %7.3f %9.4f %11.1f %9.2f\n", m, cv, ms[k].pureza, ms[k].Γ, ms[k].Δ, ms[k].hipervolume, ge, tm)
    push!(tab, Any["F6", "FF1, n = $n", m, cv, ms[k].pontos, ms[k].pureza, ms[k].Γ, ms[k].Δ, ms[k].hipervolume, ge, tm])
    open(joinpath(pasta_r, "frente_F6_$(k).csv"), "w") do io
        writedlm(io, permutedims(["f1", "f2"]), ','); writedlm(io, permutedims(F), ',')
    end
end
open(joinpath(pasta_r, "metricas_frente_FF1.csv"), "w") do io
    writedlm(io, ["caso" "descricao" "metodo" "convergiram" "pontos_nao_dominados" "pureza" "Gamma" "Delta" "hipervolume" "grad_equiv_total" "tempo_total"], ',')
    foreach(l -> writedlm(io, permutedims(l), ','), tab)
end

ts = range(-1, 1; length=400)
fig = plot(1 .- exp.(-(1 .- ts) .^ 2), 1 .- exp.(-(1 .+ ts) .^ 2); color=:black, lw=1, label="frente exata",
           xlabel="f₁", ylabel="f₂", legend=:bottomleft, size=(520, 400))
for (m, cor, marc) in zip(metodos, (:steelblue, :gray, :darkorange, :firebrick), (:circle, :utriangle, :diamond, :xcross))
    F = res[m][1]; nd = nao_dominados(F)
    scatter!(fig, F[1, nd], F[2, nd]; label=m, color=cor, marker=marc, ms=3, markerstrokewidth=0.2)
end
savefig(fig, joinpath(pasta_f, "frentes_F6_FF1.pdf"))
println("Resultados em ", abspath(pasta_r), "; figura em ", abspath(pasta_f))
