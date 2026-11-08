# Exemplo da Observação 2.5: f_i(x) = x_i + ½‖x‖², blocos {x1} e {x2}, x0 = 0.
# Em x = 0, G_1 = G_2 = 0 e G = 1/2: o algoritmo deve recorrer ao bloco 0.
# Os pontos Pareto-estacionários são x = −γ, γ no simplexo; com x0 = 0, espera-se
# convergência a (−1/2, −1/2).

isdefined(Main, :BCDMO) || include(joinpath(@__DIR__, "..", "src", "BCDMO.jl"))
using .BCDMO

P = LinearPotencia([1.0 0.0; 0.0 1.0], 2.0)
blocos = [[1], [2]]

println("G_1(0) = ", medida_G(jac_bloco(P, zeros(2), [1]))[1])
println("G_2(0) = ", medida_G(jac_bloco(P, zeros(2), [2]))[1])
println("G(0)   = ", medida_G(jac_bloco(P, zeros(2), [1, 2]))[1])

x, F, est = bcdmo(P, zeros(2), blocos; opc=Opcoes(ε=1e-6))
println("\nstatus: ", est.status)
println("x final: ", x)
println("iterações: ", est.iter, "  (no bloco 0: ", est.passos_bloco0, ")")
println("cálculos de G_j: ", est.avalGj, "   de G_0: ", est.avalG0)
