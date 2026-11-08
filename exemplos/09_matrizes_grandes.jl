# Matrizes sintéticas grandes e esparsas (regime em que a descida por blocos deve
# compensar): compara blocos com σ por bloco e só o bloco 0, com escalas iguais
# (c = 0) e diferentes entre os blocos (c = 1).
# Com N = 10 pontos iniciais, a rodada completa leva cerca de 2 horas; a tabela
# é gravada a cada linha, de modo que uma interrupção não perde o que já rodou.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, Printf, DelimitedFiles

mediana(v) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

variantes = [
    ("blocos, σ por bloco", (; sigma=:porbloco)),
    ("só bloco 0",          (; sigma=:decresce, blocos_proprios=false)),
]

# (nome, m, n, densidade, tamanho do bloco, valores de c, N)
cenarios = [
    ("n = 2000", 6000,  2000, 0.002, 100, [0.0, 1.0], 10),
    ("n = 5000", 15000, 5000, 0.001, 100, [0.0, 1.0], 10),
]

cab = ["cenario" "c" "variante" "convergiram" "iter" "grad_equiv" "tempo"]
tabela = Any[]
pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
for (nome, m, n, dens, tb, cs, N) in cenarios
    A0 = matriz_sintetica(m, n; densidade=dens, semente=1)
    blocos = blocos_consecutivos(n, tb)
    for c in cs
        local A, esc, P
        A, esc = escalar_blocos(A0, blocos, c)
        P = regressao_mesma_matriz(A, [1.5, 1.5]; semente=2)
        println("\n", nome, ", c = ", c, "  (escalas dos blocos de 1 a ", round(10.0^c; sigdigits=3), ")")
        @printf("  %-20s %6s %9s %11s %9s\n", "variante", "conv.", "iter", "grad. equiv", "tempo(s)")
        for (rot, kw) in variantes
            local ests
            opc = Opcoes(; ε=1e-2, maxiter=300_000, kw...)
            _, _, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3, escala=1 ./ esc, mostrar=false)
            conv = count(e -> e.status == :convergiu, ests)
            it = mediana([e.iter for e in ests]); ge = mediana([e.grad_equiv for e in ests])
            tm = mediana([e.tempo for e in ests])
            @printf("  %-20s %3d/%-2d %9.0f %11.1f %9.3f\n", rot, conv, N, it, ge, tm)
            push!(tabela, Any[nome, c, rot, conv, it, ge, tm])
            # grava a tabela a cada linha, para não perder resultados se for interrompido
            open(joinpath(pasta, "matrizes_grandes.csv"), "w") do io
                writedlm(io, cab, ',')
                foreach(l -> writedlm(io, permutedims(l), ','), tabela)
            end
        end
    end
end
println("\nMedianas sobre os pontos iniciais. Tabela gravada em ", abspath(pasta))
