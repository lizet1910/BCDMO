# BCDMO

Julia code for the paper

> *A block coordinate descent method with a certified joint step for multiobjective optimization under Hölder continuity*

It implements Algorithm 1 of the paper (block coordinate descent with quadratic regularization, a per-block regularization parameter and a joint step in all variables), the methods it is compared with, the test problems, and the scripts that produce every table and figure of the paper and of its supplementary material (Online Resource 1).

## Requirements

- Julia 1.13.1 (the version used for the experiments).
- Packages: JuMP, HiGHS, Clarabel, MatrixMarket, Plots (the exact versions are fixed in `Project.toml` and `Manifest.toml`).
- Internet access the first time the real matrices are used: the scripts download them from the SuiteSparse Matrix Collection into `dados/` and unpack them with `tar`.

To install the packages, start Julia in the repository folder and run

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Layout

| Folder / file | Contents |
|---|---|
| `src/BCDMO.jl` | Main module. |
| `src/problemas.jl` | Test problems: Family 1 (ℓp regression, constructions (a) and (b), convex loss and Cauchy loss), Family 2 (linear part plus power term, overlap parameter ω), FF1, synthetic matrices with block scales between 1 and 10^c. |
| `src/medida.jl` | Block measures G_j (golden-section search for r = 2 and r = 3, linear programming for r ≥ 4). |
| `src/subproblema.jl` | Block subproblem and the acceptable-pair conditions. |
| `src/algoritmo.jl` | Algorithm 1 (`bcdmo`) and its options (`Opcoes`). |
| `src/descida_fs.jl` | Fliege–Svaiter steepest descent with Armijo search (FS) and with Barzilai–Borwein step and nonmonotone search (FS-BB). |
| `src/pareto.jl`, `src/metricas.jl` | Multiple starting points, nondominated filter, front metrics, performance profiles, the L/L̄ indicator. |
| `exemplos/` | Experiment scripts (see below). |
| `resultados/` | CSV files with the results: medians per instance, and one row per run in the files ending in `_execucoes.csv`. |
| `figuras/` | Figures of the paper and of the supplementary material. |

## Main options (`Opcoes`)

- `ε`: tolerance of the stopping test G(x) ≤ ε.
- `sigma`: update rule in Step 6: `:porbloco` (per-block parameter, default), `:global` (single nondecreasing parameter, as in Amaral et al.) or `:decresce` (single parameter, reduced after each accepted step).
- `reducao`: reduction factor ρ ∈ (1, 2] of the per-block rule (default √2).
- `teste`: acceptance test: `:epsilon` (Step 6, default), `:misto` and `:sigma` (the alternative tests of Theorem 7.6).
- `selecao`: block selection, `:ciclica` (default) or `:aleatoria`.
- Family 1: `perda=:cauchy` replaces |t|^p/p by log(1 + |t|^p); `usar_hess=true` and `hess_comum=true` give the second-order model of Sections 8.5 and 8.6.

## Reproducing the results

Each script is run from the Julia REPL, for example

```julia
include("exemplos/29_tabela2_concorrentes.jl")
```

All random data and starting points are generated with fixed seeds. The scripts are numbered in the order in which they were written, and several later scripts rerun earlier ones with the final settings (per-block rule with ρ = √2 and the acceptance test of Step 6). The tables below give, for each result, the script that produced the final version and the file where it is stored. Times depend on the machine; the other quantities (iterations, gradient equivalents, matrix passes) do not.

### Paper

| Result | Script | Output in `resultados/` |
|---|---|---|
| Example of Remark 2.5 (Section 8.3) | `01_exemplo_observacao.jl` | printed on screen |
| Table 1 | `07_comparacao_bloco0.jl` (run by `21`, step 8); times: `28_execucoes_e_tempos.jl` | `comparacao_bloco0.csv`, `tabelas_1_2_execucoes.csv` |
| Family 2 with several values of ω (Section 8.3) | `23_complementos_secao8.jl` | `estresse_bloco0.csv` |
| Table 2, n = 50 and 1000 | `08_escalas_heterogeneas.jl` (run by `21`, step 8); times for n = 50: `28` | `escalas_heterogeneas.csv`, `tabelas_1_2_execucoes.csv` |
| Table 2, n = 2000 and 5000 | `09_matrizes_grandes.jl` (run by `21`, step 8); times: `29` | `matrizes_grandes.csv`, `tabela2_concorrentes.csv` |
| Tables 3 and 4 | `29_tabela2_concorrentes.jl` | `tabela2_concorrentes.csv`, `tabela2_concorrentes_execucoes.csv` |
| Tables 5 and 6, real matrices | `21_continuar_rodada_final.jl`, step 5 | `tabela4_reais.csv`, `tabela4_reais_execucoes.csv` |
| Tables 5 and 6, synthetic matrices | `28_execucoes_e_tempos.jl`, part A | `tabela3_sinteticas_10.csv`, `tabela3_sinteticas_10_execucoes.csv` |
| Tables 7 and 8, convex loss | `21_continuar_rodada_final.jl`, step 4 | `construcao_b_reais.csv`, `construcao_b_reais_execucoes.csv` |
| Tables 7 and 8, Cauchy loss | `19_rodada_final.jl`, step 3 | `nao_convexo_reais.csv`, `nao_convexo_sinteticas.csv` and the `_execucoes` files |
| Values of F at the final points, Cauchy loss (Section 8.6) | `30_cauchy_valores_F.jl`, parts A and C | `cauchy_valores_F.csv` |
| Synthetic problem with c = 1 and a 300 s budget (Section 8.6) | `26_nao_convexo_c1_orcamento.jl`, `30_cauchy_valores_F.jl` part B | `nao_convexo_c1_orcamento.csv`, `cauchy_valores_F.csv` |
| Figure 1 | `27_indicador_local.jl` | `figuras/dispersao_indicador_global.pdf` |

### Supplementary material (Online Resource 1)

| Result | Script | Output in `resultados/` |
|---|---|---|
| S1.1, computation of G_j | `src/medida.jl` | — |
| S1.2, choice of ρ | `18a_sensibilidade_e_indicador.jl` | `sensibilidade_regra.csv` |
| S1.3, alternative acceptance tests | `10_testes_aceitacao.jl` (run by `21`, step 8) | `testes_aceitacao.csv` |
| S1.4, composition of the time; Table S9 | `25_perfil_tempo.jl`; `21`, step 5 | `perfil_tempo.csv`, `tabela4_reais.csv` |
| S1.5, Family 2 | `23_complementos_secao8.jl` | `estresse_bloco0.csv` |
| S2.1, other rules for the variant with block 0 only | `29_tabela2_concorrentes.jl` | `tabela2_concorrentes.csv` |
| S2.2, Barzilai–Borwein step | `21`, step 5 | `tabela4_reais.csv` |
| S2.3, performance profiles (Figure S1) | `21`, step 5; drawn by `24_figuras_finais.jl` | `tabela4_reais_execucoes.csv`, `figuras/perfil_tabela4_*.pdf` |
| S2.4, ε = 10⁻³ (Table S1) | `21`, step 7 | `tabela5_reais.csv`, `tabela5_reais_execucoes.csv` |
| S2.5, indicator L/L̄ | `23_complementos_secao8.jl`, `27_indicador_local.jl` | `indicador_blocos.csv`, `indicador_blocos_reais.csv`, `indicador_local.csv` |
| S3.1, fronts of Algorithm 1 (Figure S2) | `02_frente_r2.jl`, `03_frente_r3.jl` (run by `21`, step 8); drawn by `24` | `frente_r2.csv`, `frente_r3.csv` |
| S3.2, comparison with the weighted sum (Table S2) | `19_rodada_final.jl`, step 2 | `metricas_frentes.csv`, `frente_F*_*.csv` |
| S3.2, problem FF1 (Figure S3) | `22_ff1_pontos_espalhados.jl`; drawn by `24` | `metricas_frente_FF1.csv`, `frente_F6_*.csv` |
| S4, quartiles and counts (Tables S3–S10) | the scripts above | the `_execucoes.csv` files; matrix passes on the real matrices: `passadas_reais_execucoes.csv` (`29`) |

### Running times

The experiments were run on a laptop with an AMD Ryzen 5 5500U processor, 16 GB of memory and Windows, with Julia 1.13.1 on a single thread. On this machine, the longest scripts took about: `19` and `21` together, 9–11 h; `29`, several hours; `30`, 3–4 h; `09` and `28`, about 2 h each; `26`, about 1 h. The others take from a few seconds to about 40 minutes.

## License

MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the paper (reference to be completed after publication).
