# Onde os blocos podem ganhar do passo conjunto: escalas diferentes entre os
# blocos (as colunas do bloco j são multiplicadas por 10^(c(j−1)/(q−1))) e n maior.
# Compara: blocos com σ por bloco, blocos com σ único (reduzido à metade) e só o bloco 0.

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO, Printf, DelimitedFiles

mediana(v) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

variantes = [
    ("blocos, σ por bloco", (; sigma=:porbloco)),
    ("blocos, σ único ↓",   (; sigma=:decresce)),
    ("só bloco 0",          (; sigma=:decresce, blocos_proprios=false)),
]

# (nome, m, n, densidade, tamanho do bloco, valores de c, N)
cenarios = [
    ("n = 50",   200,  50,  0.1,   5, [0.0, 0.5, 1.0, 1.5], 20),
    ("n = 1000", 3000, 1000, 0.005, 50, [0.0, 1.0], 10),
]

cab = ["cenario" "c" "variante" "convergiram" "iter" "grad_equiv" "tempo"]
tabela = Any[]
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
            opc = Opcoes(; ε=1e-2, maxiter=200_000, kw...)
            _, _, ests = multiplos_inicios(P, blocos, N; opc=opc, semente=3, escala=1 ./ esc, mostrar=false)
            conv = count(e -> e.status == :convergiu, ests)
            it = mediana([e.iter for e in ests]); ge = mediana([e.grad_equiv for e in ests])
            tm = mediana([e.tempo for e in ests])
            @printf("  %-20s %3d/%-2d %9.0f %11.1f %9.3f\n", rot, conv, N, it, ge, tm)
            push!(tabela, Any[nome, c, rot, conv, it, ge, tm])
        end
    end
end

pasta = joinpath(@__DIR__, "..", "resultados"); mkpath(pasta)
open(joinpath(pasta, "escalas_heterogeneas.csv"), "w") do io
    writedlm(io, cab, ',')
    for l in tabela
        writedlm(io, permutedims(l), ',')
    end
end
println("\nMedianas sobre os pontos iniciais. Tabela gravada em ", abspath(pasta))
