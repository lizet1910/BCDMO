# ---------------------------------------------------------------------------
# Algoritmo 1 (BCDQR-MO)
# ---------------------------------------------------------------------------

Base.@kwdef struct Opcoes
    α::Float64 = 1e-4
    θ::Float64 = 1.0
    ε::Float64 = 1e-4
    σ0::Float64 = 1.0            # Dados do Algoritmo 1: σ0 ≥ 1
    teste::Symbol = :epsilon     # :epsilon (Passo 6), :misto (teste de [A]) ou :sigma (testes alternativos do Teorema 7.6)
    selecao::Symbol = :ciclica   # :ciclica ou :aleatoria
    blocos_proprios::Bool = true # false: variante que usa apenas o bloco 0
    sigma::Symbol = :porbloco    # regra do Passo 6 após um passo aceito com parâmetro σ no bloco j:
                                 #   :porbloco → σ_j = max(1, σ/reducao), demais blocos inalterados (Teorema 7.4(b))
                                 #   :global   → parâmetro único, σ_ℓ = σ para todo ℓ, como em [A] (Teorema 7.4(c))
                                 #   :decresce → parâmetro único, σ_ℓ = max(1, σ/reducao) para todo ℓ (Teorema 7.4(c))
                                 #   :porbloco2 → como :porbloco, mas reduz à metade só após dois passos
                                 #                aceitos seguidos no bloco j sem duplicação; senão mantém σ
                                 #   :porbloco_raiz → σ_j = max(1, σ/√2), demais blocos inalterados
                                 #   :reinicia → parâmetro único reiniciado no piso, σ_ℓ = 1 para todo ℓ
                                 #               (Teorema 7.4(a); não satisfaz σ_j ≥ σ/2)
                                 # Todas satisfazem σ_j ≥ σ/2 e estão cobertas pelo Teorema 7.4(b).
                                 # Com teste = :sigma, o piso 1 passa a max{1, (B+θ)/2} (Teorema 7.6).
    reducao::Float64 = sqrt(2)   # fator de redução das regras :porbloco e :decresce; com 1 < reducao ≤ 2,
                                 # σ_j ≥ σ/2 e as cotas do Teorema 7.4 valem (reducao = 2: versão anterior)
    maxiter::Int = 100_000       # iterações externas
    maxaval::Int = 1_000_000     # avaliações de F
    tempo_max::Float64 = Inf     # limite de tempo (s); ao atingi-lo, status = :tempo
    semente::Int = 1
    verbose::Bool = false
    guardar_historico::Bool = false
end

mutable struct Estatisticas
    iter::Int            # iterações externas (passos aceitos)
    avalF::Int           # avaliações de F
    avalGj::Int          # cálculos de G_j, j ≥ 1
    avalG0::Int          # cálculos de G_0
    passos_bloco0::Int   # passos aceitos no bloco 0
    duplicacoes::Int     # duplicações de σ (Passo 5 e rejeições do Passo 6); em descida_fs, reduções do passo
    grad_equiv::Float64  # gradientes completos equivalentes: Σ n_j / n
    σfinal::Float64
    Gfinal::Float64      # último G_0 calculado
    residuo_max::Float64 # maior resíduo da condição (6) entre os passos aceitos
    status::Symbol       # :convergiu, :maxiter, :maxaval ou :tempo
    tempo::Float64       # segundos
    historico::Vector{Vector{Float64}}
    duplicacoes_p5::Int  # das duplicações, as feitas no Passo 5 (λ_min(B + 2σI) < θ)
end

Estatisticas() = Estatisticas(0, 0, 0, 0, 0, 0, 0.0, NaN, NaN, 0.0, :executando, 0.0, Vector{Float64}[], 0)

function _escolhe_bloco(ativo, ultimo, q, selecao, rng)
    if selecao == :aleatoria
        return rand(rng, findall(ativo))
    end
    for d in 1:q
        j = mod1(ultimo + d, q)
        ativo[j] && return j
    end
    error("nenhum bloco ativo")
end

function _decrescimo(opc::Opcoes, σ, s)
    ns2 = sum(abs2, s)
    if opc.teste == :misto
        return min(opc.α * opc.ε^2 / (16σ), opc.α * ns2)
    elseif opc.teste == :epsilon
        return opc.α * opc.ε^2 / (16σ)
    elseif opc.teste == :sigma
        return opc.α * σ * ns2
    else
        error("teste desconhecido: $(opc.teste)")
    end
end

# Passo 5: verifica λ_min(B_{i,(j)} + 2σI) ≥ θ para todo i. Sem matrizes, ou com
# matrizes semidefinidas positivas, basta 2σ ≥ θ; caso contrário, testa se
# B_{i,(j)} + (2σ − θ)I admite fatoração de Cholesky.
_pd(B, σ, θ) = issuccess(cholesky(Symmetric(Matrix(B) + (2σ - θ) * I); check=false))
function _autovalor_ok(P, Bs, σ, θ)
    (Bs === nothing || B_semidefinida(P)) && return 2σ >= θ
    Bs isa AbstractMatrix && return _pd(Bs, σ, θ)
    return all(_pd(B, σ, θ) for B in Bs)
end

"""
    bcdmo(P, x0, blocos; opc = Opcoes())

Executa o Algoritmo 1 a partir de `x0`, com a partição `blocos` de 1:n.
Devolve `(x, F(x), estatisticas)`.
"""
function bcdmo(P::ProblemaMO, x0::AbstractVector, blocos::Vector{Vector{Int}};
               opc::Opcoes=Opcoes())
    t0 = time()
    @assert opc.ε > 0
    rng = MersenneTwister(opc.semente)
    n = dim(P); q = length(blocos)
    @assert length(x0) == n
    @assert sort(reduce(vcat, blocos)) == collect(1:n) "os blocos devem particionar 1:n"

    @assert opc.σ0 >= 1 "o Algoritmo 1 exige σ0 ≥ 1"
    @assert 1 < opc.reducao <= 2 "o fator de redução deve estar em (1, 2]"
    # piso das regras do Passo 6: 1; o teste com σ‖s‖² (Teorema 7.6) exige (B+θ)/2
    σmin = opc.teste == :sigma ? max(1.0, (max(1.0, cota_B(P)) + opc.θ) / 2) : 1.0
    σk = max(opc.σ0, σmin)
    σbloco = fill(σk, q + 1)                         # usado só com as regras por bloco
    sucessos = zeros(Int, q + 1)                     # passos aceitos seguidos sem duplicação (:porbloco2)
    ε = opc.ε
    est = Estatisticas()

    x = Vector{Float64}(x0)
    Fx = fvals(P, x); est.avalF = 1
    opc.guardar_historico && push!(est.historico, copy(Fx))
    todos = collect(1:n)
    ativo = trues(q)
    ultimo = 0

    while true
        if est.iter >= opc.maxiter
            est.status = :maxiter
            break
        elseif time() - t0 >= opc.tempo_max
            est.status = :tempo
            break
        end
        fill!(ativo, opc.blocos_proprios)            # Passo 1: C ← {1,…,q}
        nativos = opc.blocos_proprios ? q : 0
        aceito = false; parar = false
        while !aceito && !parar
            if nativos == 0                          # Passo 7: teste de parada com G_0
                V = jac_bloco(P, x, todos); est.grad_equiv += 1
                G0, _ = medida_G(V); est.avalG0 += 1
                est.Gfinal = G0
                if G0 <= ε
                    est.status = :convergiu
                    parar = true
                    break
                end
                j = 0; idx = todos
            else                                     # Passos 2 e 3: escolhe j ∈ C e testa G_j
                j = _escolhe_bloco(ativo, ultimo, q, opc.selecao, rng)
                idx = blocos[j]
                V = jac_bloco(P, x, idx); est.grad_equiv += length(idx) / n
                Gj, _ = medida_G(V); est.avalGj += 1
                if Gj <= ε
                    ativo[j] = false; nativos -= 1
                    continue
                end
            end
            Bs = matrizes_B(P, x, idx)
            porbloco = opc.sigma in (:porbloco, :porbloco2, :porbloco_raiz)
            σ = porbloco ? σbloco[j+1] : σk                  # Passo 4
            σinicial = σ
            while true                               # Passos 5 e 6: passo e teste de aceitação
                if !_autovalor_ok(P, Bs, σ, opc.θ)   # Passo 5: λ_min(B_{i,(j)} + 2σI) < θ para algum i
                    σ *= 2; est.duplicacoes += 1; est.duplicacoes_p5 += 1
                    continue
                end
                sol = resolver_subproblema(V, Bs, σ)
                xt = copy(x); xt[idx] .+= sol.s
                Ft = fvals_passo(P, x, xt, idx, sol.s); est.avalF += 1
                if all(Ft .<= Fx .- _decrescimo(opc, σ, sol.s))
                    x = xt; Fx = Ft
                    if opc.sigma == :global
                        σk = σ
                    elseif opc.sigma == :decresce
                        σk = max(σmin, σ / opc.reducao)
                    elseif opc.sigma == :reinicia
                        σk = σmin
                    elseif opc.sigma == :porbloco
                        σbloco[j+1] = max(σmin, σ / opc.reducao)
                    elseif opc.sigma == :porbloco_raiz
                        σbloco[j+1] = max(σmin, σ / sqrt(2))
                    elseif opc.sigma == :porbloco2
                        if σ == σinicial
                            sucessos[j+1] += 1
                        else
                            sucessos[j+1] = 0
                        end
                        if sucessos[j+1] >= 2
                            σbloco[j+1] = max(σmin, σ / 2); sucessos[j+1] = 0
                        else
                            σbloco[j+1] = σ
                        end
                    else
                        error("opção sigma desconhecida: $(opc.sigma)")
                    end
                    est.iter += 1
                    est.residuo_max = max(est.residuo_max, sol.residuo)
                    if j == 0
                        est.passos_bloco0 += 1
                    else
                        ultimo = j
                    end
                    opc.guardar_historico && push!(est.historico, copy(Fx))
                    aceito = true
                    break
                end
                σ *= 2; est.duplicacoes += 1
                if est.avalF >= opc.maxaval
                    est.status = :maxaval
                    parar = true
                    break
                elseif time() - t0 >= opc.tempo_max
                    est.status = :tempo
                    parar = true
                    break
                end
            end
        end
        parar && break
        if opc.verbose && est.iter % 100 == 0
            @printf("iter %6d  σ = %.2e  F = %s\n", est.iter, σk, string(round.(Fx, sigdigits=6)))
        end
    end
    if est.status != :convergiu
        # G no ponto final (não contado em grad_equiv), para relatar execuções interrompidas
        est.Gfinal = medida_G(jac_bloco(P, x, todos))[1]
    end
    est.σfinal = opc.sigma in (:porbloco, :porbloco2, :porbloco_raiz) ? maximum(σbloco) : σk
    est.tempo = time() - t0
    return x, Fx, est
end
