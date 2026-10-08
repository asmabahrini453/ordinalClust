# FICHIER : test-ordinal_coclust.R
# CE QU'ON TESTE : la fonction principale ordinal_coclust
#   1. la sortie a la bonne forme (noms, tailles, probabilités qui somment à 1)
#   2. sur des données simulées très structurées, on retrouve les vrais groupes
#   3. une matrice ET un data.frame sont acceptés
#   4. les entrées invalides donnent une erreur claire
#   5. un run qui plante est sauté, avec un avertissement
#   6. si tous les runs plantent, erreur claire
#   7. si tous les runs ont un groupe vide, on garde le meilleur avec un avertissement

# Données de test : 40 lignes, 30 colonnes, 2 x 2 groupes, m = 5 (blocs bien différents)
xi_vrai <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2)
pi_vrai <- matrix(0.9, nrow = 2, ncol = 2)
set.seed(42)
sim <- simulate_ordinal_lbm(n = 40, d = 30, m = 5, alpha = c(0.5, 0.5), beta = c(0.5, 0.5),
                            xi = xi_vrai, pi = pi_vrai)

# Comparaison "au renommage près" (les numéros de groupes sont arbitraires)
concordance <- function(a, b) max(mean(a == b), mean(a != b))

# Un seul calcul, réutilisé par les tests 1 à 2 (réglages courts pour que ce soit rapide)
res <- ordinal_coclust(sim$X, K = 2, L = 2, n_init = 3, max_iter = 50, burn_in = 10)

test_that("1. la sortie a la bonne forme", {
  expect_named(res, c("row_prob", "col_prob", "row_class", "col_class", "parameters",
                      "ICL", "K", "L", "iterations", "best_init"))
  expect_named(res$parameters, c("alpha", "beta", "xi", "pi"))
  expect_equal(dim(res$row_prob), c(40, 2))
  expect_equal(dim(res$col_prob), c(30, 2))
  expect_equal(rowSums(res$row_prob), rep(1, 40))
  expect_equal(rowSums(res$col_prob), rep(1, 30))
  expect_length(res$row_class, 40)
  expect_length(res$col_class, 30)
  expect_true(all(res$row_class %in% 1:2) && all(res$col_class %in% 1:2))
  expect_equal(sum(res$parameters$alpha), 1)
  expect_equal(sum(res$parameters$beta), 1)
  expect_equal(dim(res$parameters$xi), c(2, 2))
  expect_true(is.finite(res$ICL))
  expect_true(res$best_init %in% 1:3)
})

test_that("2. on retrouve les vrais groupes", {
  expect_gt(concordance(res$row_class, sim$row_class), 0.9)
  expect_gt(concordance(res$col_class, sim$col_class), 0.9)
})

test_that("3. un data.frame est accepté", {
  res_df <- ordinal_coclust(as.data.frame(sim$X), K = 2, L = 2, n_init = 1,
                            max_iter = 30, burn_in = 5)
  expect_equal(dim(res_df$row_prob), c(40, 2))
})

test_that("4. les entrées invalides donnent une erreur", {
  X <- sim$X
  expect_error(ordinal_coclust(X, K = 2, L = 2, n_init = 0), "n_init")
  expect_error(ordinal_coclust(X, K = 0, L = 2), "K")
  expect_error(ordinal_coclust(X, K = 41, L = 2), "K")          # plus de groupes que de lignes
  expect_error(ordinal_coclust(X, K = 2, L = 31), "L")          # plus de groupes que de colonnes
  expect_error(ordinal_coclust(X, K = 2, L = 2, max_iter = 10, burn_in = 10), "burn_in")
  expect_error(ordinal_coclust(X, K = 2, L = 2, m = 3), "supérieures à m")
  X_na <- X; X_na[1, 1] <- NA
  expect_error(ordinal_coclust(X_na, K = 2, L = 2), "manquantes")
  X_dec <- X + 0.5
  expect_error(ordinal_coclust(X_dec, K = 2, L = 2), "entiers")
  expect_error(ordinal_coclust("abc", K = 2, L = 2), "nombres")
})

# Pour les tests 5, 6 et 7, on remplace temporairement une fonction du package par une version
# qui plante (ou qui produit un groupe vide), afin de vérifier que ordinal_coclust réagit bien.
vrai_sem <- sem_lbm
vrai_estimate <- estimate_partition

test_that("5. un run qui plante est sauté, avec un avertissement", {
  etat <- new.env(); etat$appels <- 0
  testthat::local_mocked_bindings(
    sem_lbm = function(...) {
      etat$appels <- etat$appels + 1
      if (etat$appels == 1) stop("panne simulée")     # seul le 1er run plante
      vrai_sem(...)
    }
  )
  expect_warning(
    res5 <- ordinal_coclust(sim$X, K = 2, L = 2, n_init = 3, max_iter = 30, burn_in = 5),
    "1 initialisation"
  )
  expect_true(res5$best_init %in% 2:3)     # le meilleur run n'est pas le 1er (planté)
})

test_that("6. si tous les runs plantent, erreur claire", {
  testthat::local_mocked_bindings(sem_lbm = function(...) stop("panne simulée"))
  expect_error(ordinal_coclust(sim$X, K = 2, L = 2, n_init = 2, max_iter = 30, burn_in = 5),
               "Aucune initialisation")
})

test_that("7. si tous les runs ont un groupe vide, on garde le meilleur avec un avertissement", {
  # on force estimate_partition à mettre toutes les lignes dans le groupe 1 : le groupe 2 est vide
  testthat::local_mocked_bindings(
    estimate_partition = function(...) {
      p <- vrai_estimate(...)
      p$row_class[] <- 1L
      p
    }
  )
  expect_warning(
    res7 <- ordinal_coclust(sim$X, K = 2, L = 2, n_init = 2, max_iter = 30, burn_in = 5),
    "groupe vide"
  )
  expect_true(is.finite(res7$ICL))
})