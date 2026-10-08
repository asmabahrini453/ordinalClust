# coclustOrdinal

Co-clustering de données ordinales par modèle à blocs latents (LBM) avec une loi CUB dans chaque bloc.

Projet 2026–2027, M2 MALIA, Université Lumière Lyon 2.
Auteures : Martine Ouedraogo et Asma Bahrini.

## Présentation

Le package regroupe **simultanément les lignes et les colonnes** d'une matrice de notes ordinales (1 à m). L'inférence se fait par un algorithme **SEM-Gibbs** lancé depuis plusieurs partitions aléatoires. Le meilleur essai est choisi par le critère **ICL**, qui sert aussi à comparer plusieurs valeurs de (K, L).

## Installation

```r
install.packages("coclustOrdinal_0.1.0.tar.gz", repos = NULL, type = "source")
library(coclustOrdinal)
```

## Exemple

```r
set.seed(1)
sim <- simulate_ordinal_lbm(
  n = 40, d = 30, m = 5,
  alpha = c(0.5, 0.5), beta = c(0.5, 0.5),
  xi = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2),
  pi = matrix(0.9, 2, 2)
)

res <- ordinal_coclust(sim$X, K = 2, L = 2, n_init = 5)
res$ICL
table(res$row_class, sim$row_class)
```

## Fonctions exportées

| Fonction | Rôle |
| --- | --- |
| `ordinal_coclust(X, K, L, n_init = 20, ...)` | co-clustering : lance `n_init` essais et renvoie le meilleur (ICL maximal) |
| `simulate_ordinal_lbm(n, d, m, alpha, beta, xi, pi)` | simulation de données selon le modèle |

`ordinal_coclust()` renvoie : `row_prob`, `col_prob` (probabilités d'appartenance), `row_class`, `col_class` (partitions), `parameters` (`alpha`, `beta`, `xi`, `pi`), `ICL`.

Aide : `?ordinal_coclust`. Vignette : `vignette("coclustOrdinal")`.

## Notations du projet

(garde ici ton tableau actuel, sans changement)

Attention : `pi` désigne le paramètre CUB π et masque la constante `pi` de R.