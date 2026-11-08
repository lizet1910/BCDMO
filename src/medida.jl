# ---------------------------------------------------------------------------
# Medida G_j(x) = min_{γ ∈ Δ^{r−1}} ‖ Σ_i γ_i ∇_(j) f_i(x) ‖_∞   (Definição 2.2)
#
# Entrada: V (n_j × r), colunas ∇_(j) f_i(x). Saída: (valor, γ).
#   r = 1: norma do máximo.
#   r = 2: busca da razão áurea em t ∈ [0,1] (a função é convexa).
#   r = 3: busca da razão áurea encaixada no triângulo (a minimização parcial
#          de uma função convexa é convexa).
#   r ≥ 4: problema de programação linear (Observação 2.5), resolvido pelo HiGHS.
# ---------------------------------------------------------------------------

function medida_G(V::AbstractMatrix{<:Real})
    r = size(V, 2)
    if r == 1
        return norm(view(V, :, 1), Inf), [1.0]
    elseif r == 2
        return _medida_G_r2(view(V, :, 1), view(V, :, 2))
    elseif r == 3
        return _medida_G_r3(view(V, :, 1), view(V, :, 2), view(V, :, 3))
    else
        return _medida_G_lp(V)
    end
end

function _medida_G_r2(a, b)
    h(t) = maximum(abs(t * a[l] + (1 - t) * b[l]) for l in eachindex(a))
    ρ = (sqrt(5) - 1) / 2
    lo, hi = 0.0, 1.0
    t1 = hi - ρ * (hi - lo); h1 = h(t1)
    t2 = lo + ρ * (hi - lo); h2 = h(t2)
    for _ in 1:80
        if h1 <= h2
            hi = t2; t2 = t1; h2 = h1
            t1 = hi - ρ * (hi - lo); h1 = h(t1)
        else
            lo = t1; t1 = t2; h1 = h2
            t2 = lo + ρ * (hi - lo); h2 = h(t2)
        end
    end
    cands = [0.0, 1.0, (lo + hi) / 2]
    vals = h.(cands)
    k = argmin(vals)
    return vals[k], [cands[k], 1 - cands[k]]
end

"Minimiza a função convexa h em [lo, hi] pela razão áurea; devolve (t, h(t))."
function _aurea(h, lo, hi; iter::Int=60)
    ρ = (sqrt(5) - 1) / 2
    t1 = hi - ρ * (hi - lo); h1 = h(t1)
    t2 = lo + ρ * (hi - lo); h2 = h(t2)
    for _ in 1:iter
        if h1 <= h2
            hi = t2; t2 = t1; h2 = h1
            t1 = hi - ρ * (hi - lo); h1 = h(t1)
        else
            lo = t1; t1 = t2; h1 = h2
            t2 = lo + ρ * (hi - lo); h2 = h(t2)
        end
    end
    cands = [lo, hi, (lo + hi) / 2]
    vals = h.(cands)
    k = argmin(vals)
    return cands[k], vals[k]
end

function _medida_G_r3(a, b, c)
    h(s, t) = maximum(abs(s * a[l] + t * b[l] + (1 - s - t) * c[l]) for l in eachindex(a))
    interno(s) = _aurea(t -> h(s, t), 0.0, 1.0 - s)[2]
    s, _ = _aurea(interno, 0.0, 1.0)
    t, v = _aurea(t -> h(s, t), 0.0, 1.0 - s)
    γ = max.([s, t, 1 - s - t], 0.0)
    return v, γ ./ sum(γ)
end

function _medida_G_lp(V)
    r = size(V, 2)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, γ[1:r] >= 0)
    @variable(model, t)
    @constraint(model, sum(γ) == 1)
    e = V * γ
    @constraint(model, e .<= t)
    @constraint(model, e .>= -t)
    @objective(model, Min, t)
    optimize!(model)
    γv = max.(value.(γ), 0.0)
    γv ./= sum(γv)
    return norm(V * γv, Inf), γv
end
