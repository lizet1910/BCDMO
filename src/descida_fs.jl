# ---------------------------------------------------------------------------
# Métodos de comparação que atualizam todas as variáveis: máxima descida
# multiobjetivo de Fliege e Svaiter [FS00] com busca de Armijo por retrocesso.
#
#   d(x) = −Σ_i γ_i ∇f_i(x),  γ minimiza ‖Σ_i γ_i ∇f_i(x)‖² no simplexo.
#
# Duas versões da busca:
#   passo_bb = false, memoria = 1: t parte de 1 e o teste é monótono,
#       f_i(x + t d) ≤ f_i(x) + c t ⟨∇f_i(x), d⟩  para todo i  (versão original);
#   passo_bb = true, memoria = M: t parte do passo de Barzilai–Borwein
#       t_BB = ⟨s, s⟩/⟨s, y⟩, com s = x_k − x_{k−1} e y = (V_k − V_{k−1})γ_k,
#       truncado em [1e−10, 1e10] (t = 1 se ⟨s, y⟩ ≤ 0), e o teste é o não monótono
#       do tipo máximo de Mita, Fukuda e Yamashita (2019),
#       f_i(x + t d) ≤ max_{0≤l<M} f_i(x_{k−l}) + c t ⟨∇f_i(x), d⟩  para todo i.
#
# Pela otimalidade de γ, ⟨∇f_i(x), d⟩ ≤ −‖d‖² < 0 para todo i, de modo que a busca
# termina. O critério de parada é o mesmo do Algoritmo 1, G(x) ≤ ε.
# ---------------------------------------------------------------------------

"""
    descida_fs(P, x0; ε, c_armijo, maxiter, maxaval, passo_bb, memoria, tempo_max)

Devolve `(x, F(x), estatisticas)`. Em `estatisticas`, `iter` conta as iterações
(cada uma com um gradiente completo), `avalF` as avaliações de F e `duplicacoes`
as reduções do passo na busca de Armijo.
"""
function descida_fs(P::ProblemaMO, x0::AbstractVector; ε::Real=1e-4, c_armijo::Real=1e-4,
                    maxiter::Int=100_000, maxaval::Int=1_000_000,
                    passo_bb::Bool=false, memoria::Int=1, tempo_max::Real=Inf)
    t0 = time()
    n = dim(P); r = nobj(P)
    @assert length(x0) == n && memoria >= 1
    todos = collect(1:n)
    est = Estatisticas()
    x = Vector{Float64}(x0)
    Fx = fvals(P, x); est.avalF = 1
    historico = [copy(Fx)]                    # valores de F nos últimos `memoria` iterados
    xant = Float64[]; Vant = zeros(0, 0)
    while true
        V = jac_bloco(P, x, todos); est.grad_equiv += 1
        G0, _ = medida_G(V); est.avalG0 += 1
        est.Gfinal = G0
        if G0 <= ε
            est.status = :convergiu
            break
        elseif est.iter >= maxiter
            est.status = :maxiter
            break
        elseif est.avalF >= maxaval
            est.status = :maxaval
            break
        elseif time() - t0 >= tempo_max
            est.status = :tempo
            break
        end
        γ = r == 1 ? [1.0] : qp_simplex(Symmetric(Matrix(V' * V)))
        d = -(V * γ)
        inc = V' * d                          # ⟨∇f_i(x), d⟩, todos negativos
        t = 1.0
        if passo_bb && !isempty(xant)
            s = x .- xant
            y = (V .- Vant) * γ
            sy = dot(s, y)
            t = sy > 0 ? clamp(dot(s, s) / sy, 1e-10, 1e10) : 1.0
        end
        Fref = reduce((a, b) -> max.(a, b), historico)   # máximo componente a componente
        aceito = false
        while t >= 1e-20
            xt = x .+ t .* d
            Ft = fvals(P, xt); est.avalF += 1
            if all(Ft .<= Fref .+ c_armijo * t .* inc)
                xant = x; Vant = V
                x = xt; Fx = Ft
                aceito = true
                break
            end
            t /= 2; est.duplicacoes += 1
        end
        if !aceito                            # só por erro de arredondamento
            est.status = :passo_minimo
            break
        end
        push!(historico, copy(Fx))
        length(historico) > memoria && popfirst!(historico)
        est.iter += 1
        est.passos_bloco0 += 1
    end
    est.σfinal = NaN
    est.tempo = time() - t0
    return x, Fx, est
end
