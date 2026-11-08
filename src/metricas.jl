# ---------------------------------------------------------------------------
# Métricas para comparar aproximações da frente de Pareto (r = 2 ou 3), perfis
# de desempenho de Dolan e Moré e o indicador L/L̄ de heterogeneidade dos blocos.
#
# Frentes: cada frente é uma matriz r × N de valores de F. As métricas usam os
# valores normalizados por objetivo, (f − mín)/(máx − mín), com mín e máx tomados
# na união de todas as frentes comparadas. Seguem Custódio, Madeira, Vaz e Vicente
# (SIAM J. Optim., 2011):
#   pureza   fração dos pontos não dominados do método que continuam não dominados
#            na união de todas as frentes;
#   Γ        maior distância entre pontos consecutivos ao longo de cada objetivo,
#            incluindo os extremos da frente de referência (menor é melhor);
#   Δ        uniformidade da distribuição dos pontos (menor é melhor);
#   hipervolume  volume dominado em relação ao ponto (1.1, …, 1.1) (maior é melhor).
# ---------------------------------------------------------------------------

function _normalizar(frentes::Vector{<:AbstractMatrix})
    U = reduce(hcat, frentes)
    lo = vec(minimum(U; dims=2)); hi = vec(maximum(U; dims=2))
    esc = max.(hi .- lo, eps())
    return [(F .- lo) ./ esc for F in frentes]
end

"Hipervolume (minimização) dos pontos em F (2 × N) em relação a R."
function _hv2(F::AbstractMatrix, R)
    idx = sortperm(F[1, :])
    hv = 0.0; prev = R[2]
    for k in idx
        f1, f2 = F[1, k], F[2, k]
        (f1 >= R[1] || f2 >= R[2]) && continue
        if f2 < prev
            hv += (R[1] - f1) * (prev - f2)
            prev = f2
        end
    end
    return hv
end

"Hipervolume em 2 ou 3 dimensões (em 3D, por fatias em f₃)."
function hipervolume(F::AbstractMatrix, R::AbstractVector)
    r = size(F, 1)
    r == 2 && return _hv2(F, R)
    @assert r == 3 "hipervolume implementado para r = 2 e r = 3"
    dentro = [k for k in 1:size(F, 2) if all(F[:, k] .< R)]
    G = F[:, dentro]
    ordem = sortperm(G[3, :])
    z = vcat(G[3, ordem], R[3])
    hv = 0.0
    for k in 1:length(ordem)
        altura = z[k+1] - z[k]
        altura > 0 || continue
        hv += altura * _hv2(G[1:2, ordem[1:k]], R[1:2])
    end
    return hv
end

"Γ e Δ de Custódio et al. para a frente F (já normalizada), com extremos de referência lo e hi."
function espalhamento(F::AbstractMatrix, lo::AbstractVector, hi::AbstractVector)
    r, N = size(F)
    Γ = 0.0; Δ = 0.0
    for j in 1:r
        v = vcat(lo[j], sort(F[j, :]), hi[j])
        δ = diff(v)                               # δ_0, …, δ_N
        Γ = max(Γ, maximum(δ))
        if N >= 2
            δint = δ[2:end-1]
            δm = sum(δint) / length(δint)
            num = δ[1] + δ[end] + sum(abs.(δint .- δm))
            den = δ[1] + δ[end] + (N - 1) * δm
            Δ = max(Δ, den > 0 ? num / den : 0.0)
        else
            Δ = max(Δ, 1.0)
        end
    end
    return Γ, Δ
end

"""
    metricas_frentes(frentes)

Recebe um vetor de frentes (r × N_s, uma por método) e devolve, para cada método,
uma tupla (pontos não dominados, pureza, Γ, Δ, hipervolume).
"""
function metricas_frentes(frentes::Vector{<:AbstractMatrix})
    Fn = _normalizar(frentes)
    proprias = [F[:, nao_dominados(F)] for F in Fn]
    U = reduce(hcat, proprias)
    ref = U[:, nao_dominados(U)]
    lo = vec(minimum(ref; dims=2)); hi = vec(maximum(ref; dims=2))
    R = fill(1.1, size(U, 1))
    res = NamedTuple[]
    for F in proprias
        na_ref = count(k -> any(l -> F[:, k] == ref[:, l], 1:size(ref, 2)), 1:size(F, 2))
        Γ, Δ = espalhamento(F, lo, hi)
        push!(res, (pontos=size(F, 2), pureza=na_ref / size(F, 2), Γ=Γ, Δ=Δ,
                    hipervolume=hipervolume(F, R)))
    end
    return res
end

"""
    perfil_desempenho(T)

T: matriz problemas × métodos com o custo de cada execução (Inf para falha).
Devolve (τ, ρ), com ρ[k, s] = fração dos problemas em que o método s tem razão de
desempenho ≤ τ[k] (Dolan e Moré, 2002).
"""
function perfil_desempenho(T::AbstractMatrix)
    melhor = vec(minimum(T; dims=2))
    validos = findall(isfinite, melhor)
    Rz = T[validos, :] ./ melhor[validos]
    τ = sort(unique(filter(isfinite, vec(Rz))))
    ρ = [count(<=(t), Rz[:, s]) / length(validos) for t in τ, s in 1:size(T, 2)]
    return τ, ρ
end

# λ_max(MᵀM) pelo método da potência (matrizes esparsas grandes)
function _lmax(M; iter::Int=300)
    v = randn(MersenneTwister(1), size(M, 2))     # vetor inicial aleatório (A·1 = 0 nas matrizes de homologia)
    v ./= norm(v)
    λ = 0.0
    for _ in 1:iter
        w = M' * (M * v)
        λ = norm(w)
        λ == 0 && return 0.0
        v = w ./ λ
    end
    return λ
end

"""
    indicador_blocos(A, blocos)

Razão L/L̄ entre ‖A‖² e a média de ‖A U_j‖² sobre os blocos, a quantidade
clássica da análise da descida por blocos (Nesterov, 2012; Wright, 2015).
"""
function indicador_blocos(A, blocos)
    L = _lmax(A)
    Lj = [length(b) <= 1500 ? eigmax(Symmetric(Matrix(A[:, b]' * A[:, b]))) : _lmax(A[:, b])
          for b in blocos]
    return L / (sum(Lj) / length(Lj))
end
