# ---------------------------------------------------------------------------
# Frente de Pareto: múltiplos pontos iniciais, filtro de não dominados,
# resumo e gravação das estatísticas
# ---------------------------------------------------------------------------

_domina(u, v) = all(u .<= v) && any(u .< v)

"Índices das colunas não dominadas de F (r × N)."
function nao_dominados(F::AbstractMatrix)
    N = size(F, 2)
    manter = trues(N)
    for a in 1:N
        for b in 1:N
            if b != a && _domina(view(F, :, b), view(F, :, a))
                manter[a] = false
                break
            end
        end
    end
    return findall(manter)
end

"""
    multiplos_inicios(P, blocos, N; opc, escala, semente, metodo)

Roda o Algoritmo 1 (metodo = :bcdmo) ou a máxima descida de Fliege e Svaiter
(metodo = :fs, ou :fsbb com passo de Barzilai–Borwein e busca não monótona;
ε, maxiter e maxaval vêm de `opc`) a partir de N pontos iniciais
x0 = escala·randn(n), ou a partir das colunas de `pontos` (n × N), se fornecida;
os pontos iniciais não dependem do método.
Devolve (X, F, estatisticas), com X (n × N) e F (r × N).
"""
function multiplos_inicios(P::ProblemaMO, blocos, N::Int;
                           opc::Opcoes=Opcoes(), escala=1.0, semente::Int=1,
                           mostrar::Bool=true, metodo::Symbol=:bcdmo,
                           pontos::Union{Nothing,AbstractMatrix}=nothing)
    rng = MersenneTwister(semente)
    n = dim(P); r = nobj(P)
    pontos === nothing || @assert size(pontos) == (n, N)
    X = zeros(n, N); F = zeros(r, N)
    ests = Estatisticas[]
    for k in 1:N
        x0 = pontos === nothing ? escala .* randn(rng, n) : Vector{Float64}(pontos[:, k])
        if metodo == :fs
            x, Fx, est = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval, tempo_max=opc.tempo_max)
        elseif metodo == :fsbb
            x, Fx, est = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval,
                                    passo_bb=true, memoria=10, tempo_max=opc.tempo_max)
        else
            x, Fx, est = bcdmo(P, x0, blocos; opc=opc)
        end
        X[:, k] = x; F[:, k] = Fx
        push!(ests, est)
        mostrar && @printf("execução %4d/%d  %-10s iter = %6d  bloco0 = %5d  tempo = %.2fs\n",
                           k, N, string(est.status), est.iter, est.passos_bloco0, est.tempo)
    end
    return X, F, ests
end

_mediana(v) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

"Imprime medianas das estatísticas de várias execuções."
function resumo(ests::Vector{Estatisticas})
    conv = count(e -> e.status == :convergiu, ests)
    @printf("convergiram: %d de %d\n", conv, length(ests))
    for (nome, f) in (("iterações externas", e -> e.iter), ("avaliações de F", e -> e.avalF),
                      ("cálculos de G_j", e -> e.avalGj), ("cálculos de G_0", e -> e.avalG0), ("gradientes equivalentes", e -> e.grad_equiv),
                      ("passos no bloco 0", e -> e.passos_bloco0), ("tempo (s)", e -> e.tempo))
        @printf("  mediana de %-20s %12.4g\n", nome * ":", _mediana(Float64.(map(f, ests))))
    end
    frac = sum(e.passos_bloco0 for e in ests) / max(1, sum(e.iter for e in ests))
    @printf("  fração de passos no bloco 0: %.4f\n", frac)
end

"Grava as estatísticas (uma linha por execução) em CSV."
function salvar_estatisticas(arquivo::AbstractString, ests::Vector{Estatisticas}, F::AbstractMatrix)
    r = size(F, 1)
    cab = ["execucao" "status" "iter" "avalF" "avalGj" "avalG0" "passos_bloco0" "duplicacoes" "grad_equiv" "sigma_final" "G_final" "residuo_max" "tempo" ["f$i" for i in 1:r]...]
    linhas = [Any[k, string(e.status), e.iter, e.avalF, e.avalGj, e.avalG0, e.passos_bloco0,
                  e.duplicacoes, e.grad_equiv, e.σfinal, e.Gfinal, e.residuo_max, e.tempo, F[:, k]...]
              for (k, e) in enumerate(ests)]
    open(arquivo, "w") do io
        writedlm(io, cab, ',')
        for l in linhas
            writedlm(io, permutedims(l), ',')
        end
    end
end
