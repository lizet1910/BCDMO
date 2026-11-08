# ---------------------------------------------------------------------------
# Subproblema por bloco (Seção 4):
#   min_z max_i φ_i(z),  φ_i(z) = ⟨g_i, z⟩ + ½⟨B_i z, z⟩ + σ‖z‖²
#
# V (n_j × r) tem colunas g_i = ∇_(j) f_i(x^k). Bs pode ser
#   nothing         → B_i = 0 (modelo linear)
#   uma matriz      → B_i = B comum (Observação 4.3: dual quadrático em γ)
#   vetor de matrizes → B_i por objetivo (problema com restrições quadráticas)
# ---------------------------------------------------------------------------

struct SolucaoSub
    s::Vector{Float64}
    γ::Vector{Float64}
    residuo::Float64   # ‖Σ_i γ_i ∇φ_i(s)‖  (condição (6) da Definição 4.1)
    φmax::Float64      # max_i φ_i(s)       (condição (5): deve ser ≤ 0)
end

_Bmul(::Nothing, i, s) = zero(s)
_Bmul(B::AbstractMatrix, i, s) = B * s
_Bmul(Bs::AbstractVector, i, s) = Bs[i] * s

function resolver_subproblema(V::AbstractMatrix{<:Real}, Bs, σ::Real)
    r = size(V, 2)
    if Bs === nothing
        γ = r == 1 ? [1.0] : qp_simplex(Symmetric(Matrix(V' * V)))
        s = -(V * γ) ./ (2σ)
    elseif Bs isa AbstractMatrix
        C = cholesky(Symmetric(Matrix(Bs) + 2σ * I))   # falha se não for positiva definida
        W = C \ Matrix(V)
        γ = r == 1 ? [1.0] : qp_simplex(Symmetric(Matrix(V' * W)))
        s = -(W * γ)
    else
        s, γ = _qcqp(V, Bs, σ)
    end
    # verificação das condições de par aceitável
    gradΣ = zeros(length(s))
    φmax = -Inf
    for i in 1:r
        Bs_i = _Bmul(Bs, i, s)
        gradΣ .+= γ[i] .* (view(V, :, i) .+ Bs_i .+ 2σ .* s)
        φmax = max(φmax, dot(view(V, :, i), s) + 0.5 * dot(s, Bs_i) + σ * dot(s, s))
    end
    return SolucaoSub(s, γ, norm(gradΣ), φmax)
end

"""Resolve min ½ γ'Mγ sobre o simplexo: fórmula fechada para r = 2, enumeração
exata das faces do simplexo para 3 ≤ r ≤ 6 e Clarabel para r ≥ 7."""
function qp_simplex(M::AbstractMatrix)
    r = size(M, 1)
    if r == 2
        a = M[1, 1] - 2 * M[1, 2] + M[2, 2]
        q(t) = 0.5 * (M[1, 1] * t^2 + 2 * M[1, 2] * t * (1 - t) + M[2, 2] * (1 - t)^2)
        t = a > 0 ? clamp((M[2, 2] - M[1, 2]) / a, 0.0, 1.0) : (q(0.0) <= q(1.0) ? 0.0 : 1.0)
        return [t, 1 - t]
    elseif r <= 6
        return _qp_simplex_faces(M)
    end
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, γ[1:r] >= 0)
    @constraint(model, sum(γ) == 1)
    @objective(model, Min, 0.5 * sum(M[i, l] * γ[i] * γ[l] for i in 1:r, l in 1:r))
    optimize!(model)
    γv = max.(value.(γ), 0.0)
    return γv ./ sum(γv)
end

function _qcqp(V, Bs, σ)
    nj, r = size(V)
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, z[1:nj])
    @variable(model, t)
    cons = [begin
                H = Matrix(Symmetric(Bs[i] + 2σ * I))
                @constraint(model, sum(V[l, i] * z[l] for l in 1:nj) + 0.5 * (z' * H * z) <= t)
            end for i in 1:r]
    @objective(model, Min, t)
    optimize!(model)
    s = value.(z)
    γ = max.(-dual.(cons), 0.0)
    γ = sum(γ) > 0 ? γ ./ sum(γ) : fill(1 / r, r)
    return s, γ
end

# Para cada face S do simplexo, resolve o sistema KKT de min ½γ_S'M_SSγ_S com
# Σγ_S = 1; entre as soluções com γ_S ≥ 0, escolhe a de menor valor. Os vértices
# são sempre viáveis, de modo que há ao menos um candidato.
function _qp_simplex_faces(M::AbstractMatrix)
    r = size(M, 1)
    melhor_v = Inf
    melhor_γ = zeros(r)
    for mascara in 1:(2^r - 1)
        S = [i for i in 1:r if (mascara >> (i - 1)) & 1 == 1]
        k = length(S)
        K = zeros(k + 1, k + 1)
        K[1:k, 1:k] = M[S, S]
        K[1:k, k+1] .= 1.0
        K[k+1, 1:k] .= 1.0
        rhs = vcat(zeros(k), 1.0)
        sol = try
            K \ rhs
        catch
            continue
        end
        (all(isfinite, sol) && norm(K * sol - rhs) <= 1e-9 * (1 + norm(K))) || continue
        γS = sol[1:k]
        minimum(γS) >= -1e-12 || continue
        γ = zeros(r)
        γ[S] = max.(γS, 0.0)
        γ ./= sum(γ)
        v = 0.5 * dot(γ, M * γ)
        if v < melhor_v
            melhor_v = v
            melhor_γ = γ
        end
    end
    return melhor_γ
end
