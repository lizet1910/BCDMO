"""
BCDMO: descida coordenada por blocos com regularização quadrática para
otimização multiobjetivo irrestrita (Algoritmo 1 do artigo).
"""
module BCDMO

using LinearAlgebra, SparseArrays, Random, Printf, DelimitedFiles
using JuMP, HiGHS, Clarabel, MatrixMarket

include("problemas.jl")
include("medida.jl")
include("subproblema.jl")
include("algoritmo.jl")
include("descida_fs.jl")
include("pareto.jl")
include("metricas.jl")

export ProblemaMO, RegressaoLp, LinearPotencia, JOS1, Escalarizado,
       nobj, dim, fvals, fvals_passo, jac_bloco, matrizes_B, cota_B, B_semidefinida,
       matriz_sintetica, ler_matriz, regressao_mesma_matriz, regressao_linhas,
       linear_potencia_blocos, blocos_consecutivos, escalar_blocos,
       medida_G, resolver_subproblema,
       Opcoes, Estatisticas, bcdmo, descida_fs,
       nao_dominados, multiplos_inicios, resumo, salvar_estatisticas,
       hipervolume, espalhamento, metricas_frentes, perfil_desempenho, indicador_blocos

end # module
