# Refaz o experimento FF1 do script 20 com pontos iniciais espalhados ao longo do
# conjunto de Pareto {t a : t ∈ [−1, 1]}:
#   x0_k = ξ_k a + ζ_k,  ξ_k igualmente espaçados em [−1,5; 1,5],  ζ_k = 0,3 z/√n, z ~ N(0, I),
# de modo que ‖ζ_k‖ ≈ 0,3. No script 20, x0 = randn/√n dá ⟨x0, a⟩ ≈ N(0, 1/n), e todos
# os métodos de descida terminam perto do meio da frente. Longe do conjunto de Pareto,
# F ≈ (1, 1) é quase constante (com n = 1000), e um ponto inicial distante seria
# aceito como ϵ-estacionário sem nenhuma iteração.
# Os mesmos pontos iniciais são usados por todos os métodos, inclusive pela soma
# ponderada, que usa o peso (t_k, 1 − t_k) com t_k igualmente espaçados em [0, 1].
# Grava resultados/metricas_frente_FF1.csv, frente_F6_k.csv e figuras/frentes_F6_FF1.pdf
# (sobrescrevendo os do script 20). Cerca de 2 a 5 minutos.

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

function rodar_ff1()
    n = 1000; ε = 1e-3; N = 100
    P = FF1(n); bl = blocos_consecutivos(n, 50)
    rng = MersenneTwister(3)
    ξ = range(-1.5, 1.5; length=N)
    X0 = hcat([ξ[k] .* P.a .+ (0.3 / sqrt(n)) .* randn(rng, n) for k in 1:N]...)
    pesos = [[t, 1 - t] for t in range(0, 1; length=N)]
    metodos = ["Algoritmo 1 (blocos)", "só bloco 0", "Fliege–Svaiter BB", "soma ponderada"]

    res = Dict{String,Any}()
    for (rot, kw, met) in (("Algoritmo 1 (blocos)", (; sigma=:porbloco), :bcdmo),
                           ("só bloco 0", (; sigma=:decresce, blocos_proprios=false), :bcdmo),
                           ("Fliege–Svaiter BB", (;), :fsbb))
        _, Fm, ests = multiplos_inicios(P, bl, N; opc=Opcoes(; ε=ε, maxiter=300_000, kw...),
                                        pontos=X0, mostrar=false, metodo=met)
        res[rot] = (Fm, sum(e.grad_equiv for e in ests), sum(e.tempo for e in ests),
                    count(e -> e.status == :convergiu, ests), sum(e.iter for e in ests))
    end
    Fs = zeros(2, N); ge = 0.0; tm = 0.0; cv = 0; it = 0
    for k in 1:N
        x, _, est = bcdmo(Escalarizado(P, pesos[k]), X0[:, k], bl; opc=Opcoes(ε=ε, maxiter=300_000))
        Fs[:, k] = fvals(P, x); ge += est.grad_equiv; tm += est.tempo
        cv += est.status == :convergiu; it += est.iter
    end
    res["soma ponderada"] = (Fs, ge, tm, cv, it)

    ms = metricas_frentes([res[m][1] for m in metodos])
    pasta_r = joinpath(@__DIR__, "..", "resultados"); pasta_f = joinpath(@__DIR__, "..", "figuras")
    println("\nFF1, n = $n, $(length(bl)) blocos, ε = $ε, $N execuções por método (pontos iniciais espalhados)")
    @printf("  %-22s %6s %7s %7s %7s %9s %11s %9s\n", "método", "conv.", "pureza", "Γ", "Δ", "hipervol.", "grad. equiv", "tempo(s)")
    tab = Any[]
    for (k, m) in enumerate(metodos)
        Fm, ge, tm, cv, it = res[m]
        @printf("  %-22s %6d %7.3f %7.3f %7.3f %9.4f %11.1f %9.2f\n", m, cv, ms[k].pureza, ms[k].Γ, ms[k].Δ, ms[k].hipervolume, ge, tm)
        push!(tab, Any["F6", "FF1, n = $n", m, cv, ms[k].pontos, ms[k].pureza, ms[k].Γ, ms[k].Δ, ms[k].hipervolume, ge, tm, it])
        open(joinpath(pasta_r, "frente_F6_$(k).csv"), "w") do io
            writedlm(io, permutedims(["f1", "f2"]), ','); writedlm(io, permutedims(Fm), ',')
        end
    end
    open(joinpath(pasta_r, "metricas_frente_FF1.csv"), "w") do io
        writedlm(io, ["caso" "descricao" "metodo" "convergiram" "pontos_nao_dominados" "pureza" "Gamma" "Delta" "hipervolume" "grad_equiv_total" "tempo_total" "iter_total"], ',')
        foreach(l -> writedlm(io, permutedims(l), ','), tab)
    end

    ts = range(-1, 1; length=400)
    fig = plot(1 .- exp.(-(1 .- ts) .^ 2), 1 .- exp.(-(1 .+ ts) .^ 2); color=:black, lw=1, label="frente exata",
               xlabel="f₁", ylabel="f₂", legend=:bottomleft, size=(520, 400))
    for (m, cor, marc) in zip(metodos, (:steelblue, :gray, :darkorange, :firebrick), (:circle, :utriangle, :diamond, :xcross))
        Fm = res[m][1]; nd = nao_dominados(Fm)
        scatter!(fig, Fm[1, nd], Fm[2, nd]; label=m, color=cor, marker=marc, ms=3, markerstrokewidth=0.2)
    end
    savefig(fig, joinpath(pasta_f, "frentes_F6_FF1.pdf"))
    println("Resultados em ", abspath(pasta_r), "; figura em ", abspath(pasta_f))
end

rodar_ff1()
