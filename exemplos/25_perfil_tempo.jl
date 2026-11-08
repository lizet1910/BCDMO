# Perfil do tempo por componente (cerca de 10 a 15 minutos). Não altera nada em `src`.
# Reproduz, com cronômetros, o laço do Algoritmo 1 (regra por bloco, teste do Passo 6) e o
# de Fliege–Svaiter BB, chamando as mesmas funções de `src`, e mede o tempo gasto em:
#   gradiente     jac_bloco (gradientes parciais ou completos);
#   G             medida_G (G_j ou G_0, razão áurea);
#   subproblema   matrizes B, teste do Passo 5 e solução do subproblema (no FS-BB: γ e passo BB);
#   F             avaliação de F no ponto tentativo;
#   cópia         construção do ponto tentativo, xt = copy(x); xt[idx] .+= s;
#   sincronização atualização do resíduo guardado após um passo aceito e comparações de vetores;
#   restante      laço, escolha do bloco, testes etc.
# Antes, confere que o laço cronometrado gera as mesmas iterações que `bcdmo`.
# Grava resultados/perfil_tempo.csv.

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

const COMP = (:grad, :G, :sub, :F, :copia, :sinc)

macro cron(t, k, ex)
    quote
        local a = time_ns()
        local v = $(esc(ex))
        $(esc(t))[$(esc(k))] += (time_ns() - a) / 1e9
        v
    end
end

"Laço do Algoritmo 1 com cronômetros (regra por bloco com ρ = √2 e teste do Passo 6)."
function bcdmo_cron(P, x0, blocos; ε=1e-2, α=1e-4, θ=1.0, maxiter=300_000, blocos_proprios=true)
    t = Dict(k => 0.0 for k in COMP)
    n = dim(P); q = length(blocos); todos = collect(1:n)
    σb = fill(1.0, q + 1); ativo = trues(q); ultimo = 0; iter = 0; ge = 0.0
    rng = MersenneTwister(1)
    x = Vector{Float64}(x0); Fx = fvals(P, x)
    t0 = time_ns()
    convergiu = false
    while iter < maxiter && !convergiu
        fill!(ativo, blocos_proprios); nat = blocos_proprios ? q : 0
        aceito = false
        while !aceito
            @cron t :sinc BCDMO._sincronizar!(P, x)
            if nat == 0
                V = @cron t :grad jac_bloco(P, x, todos); ge += 1
                G0 = @cron t :G medida_G(V)[1]
                if G0 <= ε
                    convergiu = true; break
                end
                j = 0; idx = todos
            else
                j = BCDMO._escolhe_bloco(ativo, ultimo, q, :ciclica, rng); idx = blocos[j]
                V = @cron t :grad jac_bloco(P, x, idx); ge += length(idx) / n
                Gj = @cron t :G medida_G(V)[1]
                if Gj <= ε
                    ativo[j] = false; nat -= 1; continue
                end
            end
            σ = σb[j+1]
            while true
                a = time_ns()
                Bs = matrizes_B(P, x, idx)
                while !BCDMO._autovalor_ok(P, Bs, σ, θ)
                    σ *= 2
                end
                sol = resolver_subproblema(V, Bs, σ)
                t[:sub] += (time_ns() - a) / 1e9
                xt = @cron t :copia (y = copy(x); y[idx] .+= sol.s; y)
                Ft = @cron t :F fvals_passo(P, x, xt, idx, sol.s)
                if all(Ft .<= Fx .- α * ε^2 / (16σ))
                    x = xt; Fx = Ft; σb[j+1] = max(1.0, σ / sqrt(2)); iter += 1
                    j > 0 && (ultimo = j)
                    aceito = true; break
                end
                σ *= 2
            end
        end
    end
    total = (time_ns() - t0) / 1e9
    return (iter=iter, ge=ge, total=total, t=t, convergiu=convergiu)
end

"Laço de Fliege–Svaiter BB (memória 10) com cronômetros, como em descida_fs."
function fsbb_cron(P, x0; ε=1e-2, c=1e-4, memoria=10, maxiter=300_000)
    t = Dict(k => 0.0 for k in COMP)
    n = dim(P); todos = collect(1:n)
    x = Vector{Float64}(x0); Fx = fvals(P, x); hist = [copy(Fx)]
    xant = Float64[]; Vant = zeros(0, 0); iter = 0
    t0 = time_ns(); convergiu = false
    while iter < maxiter
        V = @cron t :grad jac_bloco(P, x, todos)
        G0 = @cron t :G medida_G(V)[1]
        if G0 <= ε
            convergiu = true; break
        end
        a = time_ns()
        γ = BCDMO.qp_simplex(Symmetric(Matrix(V' * V)))
        d = -(V * γ); inc = V' * d; τ = 1.0
        if !isempty(xant)
            s = x .- xant; y = (V .- Vant) * γ; sy = dot(s, y)
            τ = sy > 0 ? clamp(dot(s, s) / sy, 1e-10, 1e10) : 1.0
        end
        Fref = reduce((u, v) -> max.(u, v), hist)
        t[:sub] += (time_ns() - a) / 1e9
        while true
            xt = @cron t :copia (x .+ τ .* d)
            Ft = @cron t :F fvals(P, xt)
            if all(Ft .<= Fref .+ c * τ .* inc)
                xant = x; Vant = V; x = xt; Fx = Ft; break
            end
            τ /= 2
        end
        push!(hist, copy(Fx)); length(hist) > memoria && popfirst!(hist)
        iter += 1
    end
    total = (time_ns() - t0) / 1e9
    return (iter=iter, ge=float(iter + 1), total=total, t=t, convergiu=convergiu)
end

instancias = Any[]
for (g, nome) in (("JGD_Homology", "ch7-8-b2"), ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz9"))
    A = obter_matriz(g, nome)
    push!(instancias, (nome, A, blocos_consecutivos(size(A, 2), cld(size(A, 2), 10)), 1.0, 3))
end
let A0 = matriz_sintetica(6000, 2000; densidade=0.002, semente=1), bl = blocos_consecutivos(2000, 100)
    A, esc = escalar_blocos(A0, bl, 1.0)
    push!(instancias, ("sintética, c = 1", A, bl, 1 ./ esc, 1))
end

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
cab = ["instancia" "metodo" "execucao" "iter" "grad_equiv" "tempo_total" "grad" "G" "subproblema" "F" "copia" "sincronizacao" "restante"]
tab = Any[]
for (nome, A, bl, escala, N) in instancias
    P = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
    n = size(A, 2)
    rng = MersenneTwister(3)
    println("\n", nome, " (n = ", n, ", ", length(bl), " blocos)")
    for k in 1:N
        x0 = escala .* randn(rng, n)
        # conferência: o laço cronometrado deve reproduzir o bcdmo
        if k == 1
            _, _, e = bcdmo(P, x0, bl; opc=Opcoes(ε=1e-2, maxiter=300_000))
            r = bcdmo_cron(P, x0, bl)
            @printf("  conferência: bcdmo %d iterações, laço cronometrado %d iterações%s\n", e.iter, r.iter,
                    e.iter == r.iter ? " (iguais)" : "  *** DIFERENTES ***")
        end
        for (met, f) in (("Algoritmo 1", () -> bcdmo_cron(P, x0, bl)),
                         ("só bloco 0", () -> bcdmo_cron(P, x0, bl; blocos_proprios=false)),
                         ("Fliege–Svaiter BB", () -> fsbb_cron(P, x0)))
            r = f()
            soma = sum(values(r.t)); resto = r.total - soma
            fr = [r.t[c] / r.total for c in COMP]
            @printf("  %-18s it %6d  GE %8.1f  tempo %7.3f s | grad %4.1f%%  G %4.1f%%  sub %4.1f%%  F %4.1f%%  cópia %4.1f%%  sinc %4.1f%%  resto %4.1f%%\n",
                    met, r.iter, r.ge, r.total, (100 .* fr)..., 100 * resto / r.total)
            push!(tab, Any[nome, met, k, r.iter, r.ge, r.total, [r.t[c] for c in COMP]..., resto])
            open(joinpath(pasta, "perfil_tempo.csv"), "w") do io
                writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
            end
        end
    end
    # custos de referência: um gradiente completo, uma avaliação completa de F e operações de tamanho n
    x = randn(MersenneTwister(5), n); todos = collect(1:n)
    BCDMO._recalcular!(P, x)
    tg = minimum(@elapsed(jac_bloco(P, x, todos)) for _ in 1:20)
    tf = minimum(@elapsed(BCDMO._recalcular!(P, x .+ 1e-9)) for _ in 1:20)
    y = copy(x)
    tc = minimum(@elapsed(copy(x)) for _ in 1:200)
    tq = minimum(@elapsed(x == y) for _ in 1:200)
    tG = minimum(@elapsed(medida_G(jac_bloco(P, x, todos))) for _ in 1:20) - tg
    @printf("  referência: gradiente completo %.2e s, F completa %.2e s, G_0 %.2e s, copy(x) %.2e s, x == y %.2e s\n",
            tg, tf, tG, tc, tq)
end
println("\nResultados em ", abspath(joinpath(pasta, "perfil_tempo.csv")))
