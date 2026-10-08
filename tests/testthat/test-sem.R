# FICHIER : test-sem.R
# RÔLE : vérifier sem_lbm() (boucle SEM-Gibbs, cours partie 3 p35-36).
# CE QU'ON TESTE :
#   1) sur des données simulées, le SEM retrouve la partition et les paramètres
#   2) il fonctionne avec une initialisation qui ne contient que les partitions
#   3) il refuse un burn_in trop grand (burn_in >= max_iter)

# Petit simulateur local (même rôle que dans test-mstep.R).
simuler_test_sem <- function(row_class, col_class, xi, pi, m) {
  X <- matrix(0, length(row_class), length(col_class))
  for (i in seq_along(row_class)) {
    for (j in seq_along(col_class)) {
      k <- row_class[i]; l <- col_class[j]
      X[i, j] <- sample(1:m, 1, prob = cub_probability(1:m, m, xi[k, l], pi[k, l]))
    }
  }
  X
}

# Mesure d'accord entre deux partitions en 2 classes, sans tenir compte des numéros :
# 1 = identiques (même si les étiquettes sont inversées), 0.5 = aucun lien.
concordance <- function(a, b) max(mean(a == b), mean(a != b))

# 1) Cas complet : on simule un tableau avec une structure connue, on part d'une partition
#    ALÉATOIRE, et on vérifie que le SEM retrouve les groupes et les paramètres.
test_that("sem_lbm retrouve la partition et les parametres", {
  set.seed(1)
  m <- 5
  vrai_row <- rep(1:2, each = 20)
  vrai_col <- rep(1:2, each = 15)
  xi_vrai <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
  pi_vrai <- matrix(0.9, nrow = 2, ncol = 2)
  X <- simuler_test_sem(vrai_row, vrai_col, xi_vrai, pi_vrai, m)
  
  set.seed(2)
  init <- list(
    row_class = sample(1:2, 40, replace = TRUE),
    col_class = sample(1:2, 30, replace = TRUE)
  )
  res <- sem_lbm(X, K = 2, L = 2, m = m, init = init, max_iter = 100, burn_in = 20)
  
  # la sortie contient tout ce qu'attendent estimate_partition et ordinal_coclust
  expect_named(res, c("alpha", "beta", "xi", "pi", "row_class", "col_class",
                      "iterations", "burn_in"))
  expect_equal(sum(res$alpha), 1)
  expect_equal(sum(res$beta), 1)
  # la partition trouvée coïncide avec la vraie (à plus de 95 %)
  expect_gt(concordance(res$row_class, vrai_row), 0.95)
  expect_gt(concordance(res$col_class, vrai_col), 0.95)
  # xi : proche du vrai, ou du vrai avec les classes inversées (changement d'étiquettes)
  xi_ok <- all(abs(res$xi - xi_vrai) < 0.15) ||
    all(abs(res$xi - xi_vrai[2:1, 2:1]) < 0.15)
  expect_true(xi_ok)
})

# 2) Branchement avec initialize_lbm : l'initialisation ne contient QUE row_class et
#    col_class ; sem_lbm doit calculer lui-même les paramètres de départ (étape M).
#    Les données sont sans structure : on vérifie seulement que tout tourne et que les
#    sorties sont cohérentes.
test_that("sem_lbm accepte une initialisation sans parametres", {
  set.seed(3)
  X <- matrix(sample(1:5, 40 * 30, replace = TRUE), 40, 30)
  init <- initialize_lbm(X, K = 3, L = 2, m = 5)
  res <- sem_lbm(X, K = 3, L = 2, m = 5, init = init, max_iter = 20, burn_in = 5)
  
  expect_equal(sum(res$alpha), 1)
  expect_equal(sum(res$beta), 1)
  expect_length(res$row_class, 40)
  expect_length(res$col_class, 30)
})

# 3) Garde-fou : avec burn_in = max_iter, il ne resterait aucune itération à moyenner.
test_that("sem_lbm refuse un burn_in trop grand", {
  X <- matrix(sample(1:5, 20 * 10, replace = TRUE), 20, 10)
  init <- initialize_lbm(X, K = 2, L = 2, m = 5)
  expect_error(sem_lbm(X, 2, 2, 5, init, max_iter = 10, burn_in = 10), "burn_in")
})