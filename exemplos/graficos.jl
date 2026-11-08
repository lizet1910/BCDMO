# Funções de gráfico (dependem do Plots). Usadas pelos scripts de exemplo.
using Plots

# Limites da escala de cores entre os quantis de 2% e 98%, para que poucos valores
# extremos não comprimam as cores dos demais pontos.
function _limites_cor(z)
    v = sort(collect(z)); m = length(v)
    lo = v[max(1, ceil(Int, 0.02 * m))]; hi = v[min(m, floor(Int, 0.98 * m) + 1)]
    return lo < hi ? (lo, hi) : (v[1], v[end])
end

"""
    plotar_frente(F; titulo = "")

F é r × N (uma coluna por execução). Para r = 2, desenha todos os pontos e
destaca os não dominados; para r = 3, desenha os não dominados no plano (f1, f2),
com a cor indicando f3 (para a versão em 3D, use `plotar_frente_3d`).
"""
function plotar_frente(F::AbstractMatrix; titulo::AbstractString="", paleta=:turbo)
    r = size(F, 1)
    nd = nao_dominados(F)
    if r == 2
        p = scatter(F[1, :], F[2, :]; label="pontos obtidos", ms=3, alpha=0.35,
                    xlabel="f1", ylabel="f2", title=titulo, legend=:topright)
        scatter!(p, F[1, nd], F[2, nd]; label="não dominados", ms=4)
        return p
    elseif r == 3
        return scatter(F[1, nd], F[2, nd]; marker_z=F[3, nd], color=paleta, clims=_limites_cor(F[3, nd]),
                       ms=5, markerstrokewidth=0.3, label="", xlabel="f1", ylabel="f2",
                       colorbar_title="f3", title=titulo, size=(650, 500))
    else
        error("plotar_frente: use r = 2 ou r = 3")
    end
end

"Frente com r = 3 em 3D (apenas pontos não dominados), com a cor indicando f3."
function plotar_frente_3d(F::AbstractMatrix; titulo::AbstractString="", paleta=:turbo)
    nd = nao_dominados(F)
    return scatter3d(F[1, nd], F[2, nd], F[3, nd]; marker_z=F[3, nd], color=paleta,
                     clims=_limites_cor(F[3, nd]),
                     ms=3, label="", xlabel="f1", ylabel="f2", zlabel="f3",
                     title=titulo, size=(650, 550), colorbar=false)
end
