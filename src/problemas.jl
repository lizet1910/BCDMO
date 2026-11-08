# ---------------------------------------------------------------------------
# Interface dos problemas
#
# Todo problema P <: ProblemaMO deve implementar
#   nobj(P)               número de objetivos r
#   dim(P)                dimensão n
#   fvals(P, x)           vetor (f_1(x), ..., f_r(x))
#   jac_bloco(P, x, idx)  matriz n_j × r cujas colunas são ∇_(j) f_i(x),
#                         onde idx são as coordenadas do bloco j
# e, opcionalmente,
#   matrizes_B(P, x, idx) nothing (modelo linear), uma matriz comum a todos
#                         os objetivos, ou um vetor com uma matriz por objetivo
#   cota_B(P)             constante B da Hipótese 2
#   B_semidefinida(P)     true se as matrizes B_{i,(j)} são sempre semidefinidas
#                         positivas (o teste do Passo 5 reduz-se a 2σ ≥ θ)
# ---------------------------------------------------------------------------

abstract type ProblemaMO end

cota_B(::ProblemaMO) = 0.0
B_semidefinida(::ProblemaMO) = false
matrizes_B(::ProblemaMO, x, idx) = nothing

"""
    fvals_passo(P, x, xt, idx, s)

Valores F(xt) no ponto tentativo xt = x + U_j s. Por padrão chama `fvals(P, xt)`;
problemas podem especializar para aproveitar que só o bloco `idx` mudou.
"""
fvals_passo(P::ProblemaMO, x, xt, idx, s) = fvals(P, xt)

# ---------------------------------------------------------------------------
# Família 1: f_i(x) = (1/p_i) ‖A_i x − b_i‖_{p_i}^{p_i},  1 < p_i ≤ 2
# O gradiente é Hölder com expoente β_i = p_i − 1.
# ---------------------------------------------------------------------------

# Estado da Família 1 no iterado corrente x: resíduos r_i = A_i x − b_i, pesos
# w_i = sign(r_i)|r_i|^(p_i−1) e valores f_i. Um passo no bloco j altera apenas
# as linhas em que A_i[:, idx] tem elementos não nulos; o ponto tentativo guarda
# só essas linhas, e a aceitação atualiza o estado nelas. Assim, o custo de um
# passo de bloco é proporcional a nnz(A_i[:, idx]), como em [A].
const _BlocoCompacto = Tuple{SparseMatrixCSC{Float64,Int},Vector{Int}}

# Contador de "passadas pela matriz": cada produto com as colunas idx (gradiente parcial
# ou avaliação de F num ponto tentativo) soma Σ_i nnz(A_i[:, idx]) / Σ_i nnz(A_i); um
# recálculo completo do resíduo soma 1. A cópia das linhas alteradas na aceitação de um
# passo não é contada. Zerar com `BCDMO.PASSADAS[] = 0.0` antes de cada execução.
const PASSADAS = Ref(0.0)
_nnz_total(P) = sum(nnz, P.A)

mutable struct _CacheResiduo
    x::Vector{Float64}
    r::Vector{Vector{Float64}}
    w::Vector{Vector{Float64}}
    f::Vector{Float64}
    xt::Vector{Float64}
    linhas_t::Vector{Vector{Int}}
    rt::Vector{Vector{Float64}}
    ft::Vector{Float64}
    blocos::Dict{UInt64,Vector{_BlocoCompacto}}
    aceitos::Int        # atualizações incrementais desde o último recálculo completo
end
_CacheResiduo() = _CacheResiduo(Float64[], Vector{Float64}[], Vector{Float64}[], Float64[],
                                Float64[], Vector{Int}[], Vector{Float64}[], Float64[],
                                Dict{UInt64,Vector{_BlocoCompacto}}(), 0)

struct RegressaoLp <: ProblemaMO
    A::Vector{SparseMatrixCSC{Float64,Int}}
    b::Vector{Vector{Float64}}
    p::Vector{Float64}
    usar_hess::Bool   # true: B_{i,(j)} = A_j' D̃_i A_j truncada, como em [A]
    hess_comum::Bool  # true: matriz comum B_{(j)} = (1/r) Σ_i A_j' D̃_i A_j (subproblema com Cholesky)
    hess_max::Int     # blocos com mais coordenadas usam B = 0 (em geral, só o bloco 0)
    perda::Symbol     # :potencia, ρ(t) = |t|^p/p (convexa), ou :cauchy, ρ(t) = log(1 + |t|^p) (não convexa)
    D::Float64        # truncamento dos pesos da Hessiana
    B::Float64        # cota da Hipótese 2
    cache::_CacheResiduo
end

function RegressaoLp(A::AbstractVector, b::AbstractVector, p::AbstractVector;
                     usar_hess::Bool=false, hess_comum::Bool=false,
                     hess_max::Int=2000, D::Real=1.0, perda::Symbol=:potencia)
    @assert perda in (:potencia, :cauchy) "perda deve ser :potencia ou :cauchy"
    r = length(A)
    @assert length(b) == r && length(p) == r
    n = size(A[1], 2)
    @assert all(size(Ai, 2) == n for Ai in A) "todas as A_i devem ter n colunas"
    @assert all(1 < pe <= 2 for pe in p) "exige 1 < p_i ≤ 2"
    Asp = [SparseMatrixCSC{Float64,Int}(Ai) for Ai in A]
    B = usar_hess ? Float64(D) * maximum(norma2_estimada(Ai)^2 for Ai in Asp) : 0.0
    return RegressaoLp(Asp, [Vector{Float64}(bi) for bi in b], Vector{Float64}(p),
                       usar_hess, hess_comum, hess_max, perda, Float64(D), B, _CacheResiduo())
end

nobj(P::RegressaoLp) = length(P.A)
dim(P::RegressaoLp) = size(P.A[1], 2)
cota_B(P::RegressaoLp) = P.B
B_semidefinida(P::RegressaoLp) = P.perda === :potencia   # com :cauchy, os pesos podem ser negativos

# Perda por resíduo e suas derivadas. Com :potencia, f_i = Σ_l |r_l|^p / p; com :cauchy,
# f_i = Σ_l log(1 + |r_l|^p), cuja derivada p|t|^(p−1)sign(t)/(1 + |t|^p) também é
# (p−1)-Hölder, mas cuja segunda derivada p|t|^(p−2)(p − 1 − |t|^p)/(1 + |t|^p)^2 muda de sinal.
_rho(t, p, perda) = perda === :potencia ? abs(t)^p / p : log1p(abs(t)^p)
function _peso(t, p, perda=:potencia)
    perda === :potencia && return sign(t) * abs(t)^(p - 1)
    return p * sign(t) * abs(t)^(p - 1) / (1 + abs(t)^p)
end

# Para cada objetivo, A_i[:, idx] restrita às linhas não nulas, e essas linhas.
function _bloco_compacto(P::RegressaoLp, idx)
    c = P.cache
    return get!(c.blocos, hash(idx)) do
        [begin
             Aj = P.A[i][:, idx]
             linhas = sort!(unique(rowvals(Aj)))
             (Aj[linhas, :], linhas)
         end for i in 1:nobj(P)]
    end
end

function _recalcular!(P::RegressaoLp, x)
    c = P.cache
    c.x = copy(x)
    c.r = [P.A[i] * x .- P.b[i] for i in 1:nobj(P)]
    c.w = [_peso.(c.r[i], P.p[i], P.perda) for i in 1:nobj(P)]
    c.f = [sum(_rho.(c.r[i], P.p[i], P.perda)) for i in 1:nobj(P)]
    c.xt = Float64[]
    c.aceitos = 0
    PASSADAS[] += 1.0
    return c
end

# Deixa o estado em x: nada a fazer se x é o iterado guardado; atualização
# incremental se x é o último ponto tentativo; recálculo completo caso contrário
# (e a cada 1000 atualizações incrementais, para evitar acúmulo de arredondamento).
function _sincronizar!(P::RegressaoLp, x)
    c = P.cache
    if length(c.x) == length(x) && c.x == x
        return c
    elseif length(c.xt) == length(x) && c.xt == x && c.aceitos < 1000
        for i in 1:nobj(P)
            l = c.linhas_t[i]
            c.r[i][l] .= c.rt[i]
            c.w[i][l] .= _peso.(c.rt[i], P.p[i], P.perda)
        end
        c.f = copy(c.ft)
        c.x = c.xt
        c.xt = Float64[]
        c.aceitos += 1
        return c
    end
    return _recalcular!(P, x)
end

fvals(P::RegressaoLp, x) = copy(_sincronizar!(P, x).f)

function fvals_passo(P::RegressaoLp, x, xt, idx, s)
    c = _sincronizar!(P, x)
    blk = _bloco_compacto(P, idx)
    r = nobj(P)
    c.linhas_t = Vector{Vector{Int}}(undef, r)
    c.rt = Vector{Vector{Float64}}(undef, r)
    c.ft = zeros(r)
    PASSADAS[] += sum(nnz(blk[i][1]) for i in 1:r) / _nnz_total(P)
    for i in 1:r
        Ajr, linhas = blk[i]
        pe = P.p[i]
        rvelho = c.r[i][linhas]
        rnovo = rvelho .+ Ajr * s
        c.ft[i] = c.f[i] + (sum(_rho.(rnovo, pe, P.perda)) - sum(_rho.(rvelho, pe, P.perda)))
        c.linhas_t[i] = linhas
        c.rt[i] = rnovo
    end
    c.xt = copy(xt)
    return copy(c.ft)
end

function jac_bloco(P::RegressaoLp, x, idx)
    c = _sincronizar!(P, x)
    blk = _bloco_compacto(P, idx)
    V = zeros(length(idx), nobj(P))
    PASSADAS[] += sum(nnz(blk[i][1]) for i in 1:nobj(P)) / _nnz_total(P)
    for i in 1:nobj(P)
        Ajr, linhas = blk[i]
        V[:, i] = Ajr' * c.w[i][linhas]
    end
    return V
end

# peso da Hessiana truncada: (p − 1)|t|^(p−2), limitado por D
function _peso_hess(t, pe, D, perda=:potencia)
    if perda === :potencia
        return pe == 2 ? min(D, 1.0) : (t == 0 ? D : min(D, (pe - 1) * abs(t)^(pe - 2)))
    end
    t == 0 && return pe == 2 ? clamp(2.0, -D, D) : D
    a = abs(t); u = a^pe
    return clamp(pe * a^(pe - 2) * (pe - 1 - u) / (1 + u)^2, -D, D)   # pode ser negativo
end

function matrizes_B(P::RegressaoLp, x, idx)
    (P.usar_hess && length(idx) <= P.hess_max) || return nothing
    c = _sincronizar!(P, x)
    blk = _bloco_compacto(P, idx)
    r = nobj(P)
    termo(i) = (Ajr = blk[i][1]; d = _peso_hess.(c.r[i][blk[i][2]], P.p[i], P.D, P.perda); Ajr' * Diagonal(d) * Ajr)
    if P.hess_comum
        # soma esparsa e uma única conversão para matriz densa (economiza memória no bloco 0)
        S = termo(1)
        for i in 2:r
            S += termo(i)
        end
        return Matrix(S) ./ r
    end
    return [Matrix(termo(i)) for i in 1:r]
end

"""Cota superior garantida para ‖A‖₂: √(‖A‖₁‖A‖_∞). Usada na constante B da
Hipótese 2 quando `usar_hess = true`."""
function norma2_estimada(A)
    Aabs = abs.(A)
    n1 = maximum(sum(Aabs; dims=1))      # ‖A‖₁: maior soma de coluna
    ninf = maximum(sum(Aabs; dims=2))    # ‖A‖_∞: maior soma de linha
    return sqrt(n1 * ninf)
end

"Matriz esparsa aleatória m × n, sem colunas nulas."
function matriz_sintetica(m::Int, n::Int; densidade::Real=0.05, semente::Int=1)
    rng = MersenneTwister(semente)
    A = sprandn(rng, m, n, densidade)
    A += sparse(rand(rng, 1:m, n), collect(1:n), randn(rng, n), m, n)
    return A
end

"Lê uma matriz no formato Matrix Market (.mtx), como as da coleção SuiteSparse."
ler_matriz(caminho::AbstractString) = SparseMatrixCSC{Float64,Int}(MatrixMarket.mmread(caminho))

"""
    escalar_blocos(A, blocos, c)

Multiplica as colunas do bloco j por 10^(c·(j−1)/(q−1)), de modo que as escalas
dos blocos variem de 1 a 10^c. Devolve (A escalada, vetor de fatores por coluna).
"""
function escalar_blocos(A, blocos::Vector{Vector{Int}}, c::Real)
    q = length(blocos)
    esc = ones(size(A, 2))
    for (j, bl) in enumerate(blocos)
        esc[bl] .= 10.0^(c * (j - 1) / max(q - 1, 1))
    end
    return A * Diagonal(esc), esc
end

"Construção (a): mesma matriz A para todos os objetivos, vetores b_i aleatórios."
function regressao_mesma_matriz(A, p::AbstractVector; semente::Int=1, kwargs...)
    rng = MersenneTwister(semente)
    r = length(p)
    m = size(A, 1)
    b = [randn(rng, m) for _ in 1:r]
    return RegressaoLp([A for _ in 1:r], b, p; kwargs...)
end

"Construção (b): as linhas de A são repartidas aleatoriamente em r grupos."
function regressao_linhas(A, p::AbstractVector; semente::Int=1, kwargs...)
    rng = MersenneTwister(semente)
    r = length(p)
    m = size(A, 1)
    perm = randperm(rng, m)
    cortes = round.(Int, range(0, m, length=r + 1))
    As = [A[perm[cortes[i]+1:cortes[i+1]], :] for i in 1:r]
    b = [randn(rng, size(Ai, 1)) for Ai in As]
    return RegressaoLp(As, b, p; kwargs...)
end

# ---------------------------------------------------------------------------
# Família 2: f_i(x) = ⟨a_i, x⟩ + (1/p) Σ_l |x_l|^p,  1 < p ≤ 2
# Com a = I e p = 2, é o exemplo da Observação 2.5.
# ---------------------------------------------------------------------------

struct LinearPotencia <: ProblemaMO
    a::Matrix{Float64}   # n × r, coluna i = a_i
    p::Float64
end

nobj(P::LinearPotencia) = size(P.a, 2)
dim(P::LinearPotencia) = size(P.a, 1)

function fvals(P::LinearPotencia, x)
    s = sum(abs.(x) .^ P.p) / P.p
    return [dot(view(P.a, :, i), x) + s for i in 1:nobj(P)]
end

function jac_bloco(P::LinearPotencia, x, idx)
    xb = x[idx]
    g = sign.(xb) .* abs.(xb) .^ (P.p - 1)
    return P.a[idx, :] .+ g
end

"""
Família 2 com conflito controlado entre blocos: cada bloco é atribuído a um
objetivo (em rodízio), e a_i é aleatório nos seus blocos. Nos demais blocos,
a_i recebe ruído de intensidade `sobreposicao`. Com sobreposicao = 0 e x = 0,
tem-se G_j(0) = 0 para todo bloco próprio e G(0) > 0.
"""
function linear_potencia_blocos(blocos::Vector{Vector{Int}}, r::Int;
                                p::Real=1.5, sobreposicao::Real=0.0, semente::Int=1)
    rng = MersenneTwister(semente)
    n = sum(length, blocos)
    a = zeros(n, r)
    for (k, bl) in enumerate(blocos)
        dono = mod1(k, r)
        for i in 1:r
            a[bl, i] .= i == dono ? randn(rng, length(bl)) : sobreposicao .* randn(rng, length(bl))
        end
    end
    return LinearPotencia(a, Float64(p))
end

# ---------------------------------------------------------------------------
# Família 3: JOS1, f_1 = ‖x‖²/n, f_2 = ‖x − 2·1‖²/n  (caso Lipschitz, β = 1)
# ---------------------------------------------------------------------------

struct JOS1 <: ProblemaMO
    n::Int
end

nobj(::JOS1) = 2
dim(P::JOS1) = P.n
fvals(P::JOS1, x) = [sum(abs2, x) / P.n, sum(abs2, x .- 2) / P.n]
jac_bloco(P::JOS1, x, idx) = hcat(2 .* x[idx] ./ P.n, 2 .* (x[idx] .- 2) ./ P.n)

"Partição de 1:n em blocos consecutivos de tamanho `tam` (o último pode ser menor)."
blocos_consecutivos(n::Int, tam::Int) = [collect(i:min(i + tam - 1, n)) for i in 1:tam:n]

# ---------------------------------------------------------------------------
# Escalarização por soma ponderada: f_w(x) = Σ_i w_i f_i(x), como problema com r = 1.
# Com r = 1, o Algoritmo 1 é o método escalar por blocos de [A] (Observação 7.7).
# Usa o modelo linear (B = 0) e herda o cache do problema original.
# ---------------------------------------------------------------------------

struct Escalarizado{T<:ProblemaMO} <: ProblemaMO
    P::T
    w::Vector{Float64}
end

nobj(::Escalarizado) = 1
dim(E::Escalarizado) = dim(E.P)
fvals(E::Escalarizado, x) = [dot(E.w, fvals(E.P, x))]
fvals_passo(E::Escalarizado, x, xt, idx, s) = [dot(E.w, fvals_passo(E.P, x, xt, idx, s))]
jac_bloco(E::Escalarizado, x, idx) = reshape(jac_bloco(E.P, x, idx) * E.w, :, 1)
