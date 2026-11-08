# Tabela 2 com todos os concorrentes e passadas pela matriz (rodada longa; pode ser
# interrompida e retomada: execuções já gravadas não são repetidas).
# É preciso reiniciar o Julia antes, pois `src` mudou (contador BCDMO.PASSADAS e regra
# σ = :reinicia).
#
#   A. Sintéticas da Tabela 2, construção (a), n ∈ {1000, 2000, 5000}, c ∈ {0, 1}, os mesmos
#      10 pontos iniciais de antes, todas na mesma sessão:
#        blocos, σ por bloco (ρ = √2);  só bloco 0 com ρ = √2, ρ = 2 e σ reiniciado em 1;
#        Fliege–Svaiter (FS);  Fliege–Svaiter BB (FS-BB).
#      Com n = 5000, das três versões da variante só bloco 0 roda-se só ρ = √2 e a que usou
#      menos gradientes equivalentes com n = 2000 e c = 1 (para economizar tempo).
#   B. Matrizes reais da Tabela 3: método por blocos de 1ª e de 2ª ordem, 20 pontos
#      iniciais, só para medir as passadas pela matriz (nos métodos que atualizam todas as
#      variáveis, passadas = gradientes equivalentes + avaliações de F, que já estão gravados).
# Cada execução tem limite de tempo `orcamento` (status "tempo" se for atingido).
# Ordem: n = 1000; n = 2000; parte B; n = 5000, c = 0; n = 5000, c = 1 (a mais demorada).
#
# Grava, sem apagar os arquivos antigos:
#   resultados/tabela2_concorrentes_execucoes.csv   (uma linha por execução)
#   resultados/tabela2_concorrentes.csv             (medianas e quartis)
#   resultados/passadas_reais_execucoes.csv         (parte B, uma linha por execução)

include(joinpath(@__DIR__, "comum_comparacao.jl"))
using Random, LinearAlgebra

# No Windows, impede que o computador entre em suspensão enquanto o script roda
# (ES_CONTINUOUS | ES_SYSTEM_REQUIRED); a tela pode desligar normalmente. Uma suspensão
# no meio de uma execução invalida o tempo medido.
Sys.iswindows() && ccall((:SetThreadExecutionState, "kernel32"), UInt32, (UInt32,), 0x80000001)

orcamento = 1200.0     # segundos por execução
N = 10                 # pontos iniciais nas sintéticas
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
arq_exe = joinpath(pasta, "tabela2_concorrentes_execucoes.csv")
arq_med = joinpath(pasta, "tabela2_concorrentes.csv")
arq_reais = joinpath(pasta, "passadas_reais_execucoes.csv")
cab = ["instancia" "n" "c" "variante" "execucao" "status" "iter" "avalF" "grad_equiv" "passadas" "passos_bloco0" "G_final" "tempo"]

"Lê um CSV de execuções (se existir) e devolve as linhas e o conjunto de chaves (instancia, c, variante, execucao)."
function ler_execucoes(arq)
    linhas = Any[]; feitas = Set{Tuple{String,Float64,String,Int}}()
    isfile(arq) || return linhas, feitas
    d, _ = readdlm(arq, ','; header=true)
    for l in 1:size(d, 1)
        lin = Any[d[l, k] for k in 1:size(d, 2)]
        push!(linhas, lin)
        cc = lin[3] isa Number ? Float64(lin[3]) : 0.0
        push!(feitas, (string(lin[1]), cc, string(lin[4]), Int(lin[5])))
    end
    return linhas, feitas
end
gravar(arq, tab) = open(arq, "w") do io
    writedlm(io, cab, ','); foreach(l -> writedlm(io, permutedims(l), ','), tab)
end

"Pontos iniciais como em multiplos_inicios: x0 = escala .* randn(n), semente 3."
pontos(n, N, escala) = (rng = MersenneTwister(3); [escala .* randn(rng, n) for _ in 1:N])

"Executa um método a partir de x0 e devolve (estatísticas, passadas)."
function executar(P, x0, bl, met, opc)
    BCDMO.PASSADAS[] = 0.0
    if met == :bcdmo
        _, _, e = bcdmo(P, x0, bl; opc=opc)
    elseif met == :fs
        _, _, e = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval, tempo_max=opc.tempo_max)
    else
        _, _, e = descida_fs(P, x0; ε=opc.ε, maxiter=opc.maxiter, maxaval=opc.maxaval, tempo_max=opc.tempo_max,
                             passo_bb=true, memoria=10)
    end
    return e, BCDMO.PASSADAS[]
end

const V_BLOCOS = ("blocos, σ por bloco", :bcdmo, (; sigma=:porbloco))
const V_B0     = [("só bloco 0, ρ = √2",       :bcdmo, (; sigma=:decresce, blocos_proprios=false)),
                  ("só bloco 0, ρ = 2",        :bcdmo, (; sigma=:decresce, reducao=2.0, blocos_proprios=false)),
                  ("só bloco 0, σ reiniciado", :bcdmo, (; sigma=:reinicia, blocos_proprios=false))]
const V_FS     = [("Fliege–Svaiter", :fs, (;)), ("Fliege–Svaiter BB", :fsbb, (;))]

exe, feitas = ler_execucoes(arq_exe)
# Execuções com tempo maior que o orçamento só ocorrem se o computador suspendeu no meio;
# são descartadas e refeitas.
let ruins = [l for l in exe if Float64(l[13]) > orcamento + 60]
    for l in ruins
        println("refazendo (tempo inválido, provável suspensão): ", l[1], ", c = ", l[3], ", ", l[4], ", execução ", l[5])
        delete!(feitas, (string(l[1]), Float64(l[3]), string(l[4]), Int(l[5])))
    end
    filter!(l -> Float64(l[13]) <= orcamento + 60, exe)
end

function rodar_sintetica!(exe, feitas, n, c, vars)
    m, dens, tb, maxit = n == 1000 ? (3000, 0.005, 50, 200_000) :
                         n == 2000 ? (6000, 0.002, 100, 300_000) : (15000, 0.001, 100, 300_000)
    A0 = matriz_sintetica(m, n; densidade=dens, semente=1); bl = blocos_consecutivos(n, tb)
    A, esc = escalar_blocos(A0, bl, c)
    P = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
    xs = pontos(n, N, 1 ./ esc)
    inst = "n = $n"
    println("\n", inst, ", c = ", c)
    for (rot, met, kw) in vars
        opc = Opcoes(; ε=1e-2, maxiter=maxit, maxaval=50_000_000, tempo_max=orcamento, kw...)
        for (k, x0) in enumerate(xs)
            (inst, Float64(c), rot, k) in feitas && continue
            e, ps = executar(P, x0, bl, met, opc)
            push!(exe, Any[inst, n, c, rot, k, string(e.status), e.iter, e.avalF, e.grad_equiv, ps,
                           e.passos_bloco0, e.Gfinal, e.tempo])
            push!(feitas, (inst, Float64(c), rot, k))
            gravar(arq_exe, exe)
            @printf("  %-26s %2d  %-10s grad. equiv %10.1f  passadas %10.1f  tempo %8.2f s\n",
                    rot, k, string(e.status), e.grad_equiv, ps, e.tempo)
        end
    end
end

"Variante só bloco 0 com menos gradientes equivalentes (mediana) com n = 2000 e c = 1."
function melhor_b0(exe)
    melhor = V_B0[1]; gm = Inf
    for v in V_B0
        g = [Float64(l[9]) for l in exe if string(l[1]) == "n = 2000" && Float64(l[3]) == 1.0 && string(l[4]) == v[1]]
        isempty(g) && continue
        mg = mediana(g)
        mg < gm && (gm = mg; melhor = v)
    end
    return melhor
end

# ---------------- A. n = 1000 e n = 2000 ----------------
for n in (1000, 2000), c in (0.0, 1.0)
    rodar_sintetica!(exe, feitas, n, c, [V_BLOCOS; V_B0; V_FS])
end

# ---------------- B. matrizes reais: passadas do método por blocos ----------------
println("\nB. Matrizes reais, método por blocos (passadas pela matriz)")
exe_r, feitas_r = ler_execucoes(arq_reais)
for (g, nome) in (("JGD_Homology", "ch6-6-b3"), ("JGD_Homology", "ch7-8-b2"), ("JGD_Homology", "mk10-b3"),
                  ("JGD_Franz", "Franz5"), ("JGD_Franz", "Franz6"), ("JGD_Franz", "Franz7"),
                  ("JGD_Franz", "Franz8"), ("JGD_Franz", "Franz9"), ("JGD_Franz", "Franz10"))
    A = obter_matriz(g, nome); n = size(A, 2)
    bl = blocos_consecutivos(n, cld(n, 10))
    xs = pontos(n, 20, 1.0)
    for (rot, modelo) in (("blocos, 1ª ordem", :linear), ("blocos, 2ª ordem", :hess))
        P = problema(A, modelo)
        opc = Opcoes(; ε=1e-2, maxiter=300_000, maxaval=3_000_000, sigma=:porbloco)
        for (k, x0) in enumerate(xs)
            (nome, 0.0, rot, k) in feitas_r && continue
            e, ps = executar(P, x0, bl, :bcdmo, opc)
            push!(exe_r, Any[nome, n, "", rot, k, string(e.status), e.iter, e.avalF, e.grad_equiv, ps,
                             e.passos_bloco0, e.Gfinal, e.tempo])
            push!(feitas_r, (nome, 0.0, rot, k))
            gravar(arq_reais, exe_r)
        end
        g_ = [Float64(l[9]) for l in exe_r if string(l[1]) == nome && string(l[4]) == rot]
        p_ = [Float64(l[10]) for l in exe_r if string(l[1]) == nome && string(l[4]) == rot]
        @printf("  %-9s %-18s grad. equiv %8.1f  passadas %8.1f\n", nome, rot, mediana(g_), mediana(p_))
    end
end

# ---------------- A. n = 5000 ----------------
vb0 = melhor_b0(exe)
vars5000 = vb0 === V_B0[1] ? [V_BLOCOS; V_B0[1]; V_FS] : [V_BLOCOS; V_B0[1]; vb0; V_FS]
println("\nCom n = 5000, variante só bloco 0 adicional: ", vb0 === V_B0[1] ? "nenhuma" : vb0[1])
for c in (0.0, 1.0)
    rodar_sintetica!(exe, feitas, 5000, c, vars5000)
end

# ---------------- resumo: medianas e quartis ----------------
quartil(v, t) = (s = sort(v); h = (length(s) - 1) * t + 1; i = floor(Int, h); i == length(s) ? s[i] : s[i] + (h - i) * (s[i+1] - s[i]))
cab_med = ["instancia" "c" "variante" "convergiram" "execucoes" "grad_equiv" "grad_q1" "grad_q3" "passadas" "passadas_q1" "passadas_q3" "avalF" "iter" "tempo" "tempo_q1" "tempo_q3" "G_final_max"]
med = Any[]
println("\nResumo (medianas [quartis])")
for inst in ("n = 1000", "n = 2000", "n = 5000"), c in (0.0, 1.0)
    for rot in unique(string(l[4]) for l in exe)
        ls = [l for l in exe if string(l[1]) == inst && Float64(l[3]) == c && string(l[4]) == rot]
        isempty(ls) && continue
        g = Float64.(getindex.(ls, 9)); p = Float64.(getindex.(ls, 10)); t = Float64.(getindex.(ls, 13))
        conv = count(l -> string(l[6]) == "convergiu", ls)
        push!(med, Any[inst, c, rot, conv, length(ls), mediana(g), quartil(g, 0.25), quartil(g, 0.75),
                       mediana(p), quartil(p, 0.25), quartil(p, 0.75), mediana(Float64.(getindex.(ls, 8))),
                       mediana(Float64.(getindex.(ls, 7))), mediana(t), quartil(t, 0.25), quartil(t, 0.75),
                       maximum(Float64.(getindex.(ls, 12)))])
        @printf("  %-9s c = %.0f  %-26s conv. %2d/%-2d  grad. %9.0f [%9.0f; %9.0f]  passadas %9.0f  tempo %8.2f\n",
                inst, c, rot, conv, length(ls), med[end][6], med[end][7], med[end][8], med[end][9], med[end][14])
    end
end
open(arq_med, "w") do io
    writedlm(io, cab_med, ','); foreach(l -> writedlm(io, permutedims(l), ','), med)
end
println("\nResultados em ", abspath(pasta))
Sys.iswindows() && ccall((:SetThreadExecutionState, "kernel32"), UInt32, (UInt32,), 0x80000000)
