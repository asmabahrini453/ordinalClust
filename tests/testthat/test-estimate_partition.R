# FICHIER : test-estimate_partition.R
# CE QU'ON TESTE : la fonction estimate_partition (partition finale + probabilités)
#   1. la forme de la sortie (tailles des vecteurs et des matrices)
#   2. les probabilités sont bien des probabilités (entre 0 et 1, somme 1 par ligne)
#   3. les groupes sont des entiers entre 1 et K (ou 1 et L), et ce sont les plus probables
#   4. avec des données très structurées et les vrais paramètres, on retrouve les vrais groupes
#   5. la fonction marche sans partition de départ, et avec n_gibbs = 1

# Données de test communes : 40 lignes, 30 colonnes, 2 x 2 groupes, m = 5 modalités.
# Les blocs sont très différents (xi = 0.8 ou 0.2, pi = 0.9), donc faciles à retrouver.
alpha_vrai <- c(0.5, 0.5)
beta_vrai  <- c(0.5, 0.5)
xi_vrai    <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2)
pi_vrai    <- matrix(0.9, nrow = 2, ncol = 2)

set.seed(123)
sim <- simulate_ordinal_lbm(n = 40, d = 30, m = 5, alpha = alpha_vrai, beta = beta_vrai,
                            xi = xi_vrai, pi = pi_vrai)

# Fonction de comparaison "au renommage près" : les numéros de groupes sont arbitraires
# (le groupe 1 de l'un peut être le groupe 2 de l'autre), donc on prend le meilleur des deux
concordance <- function(a, b) max(mean(a == b), mean(a != b))

res <- estimate_partition(sim$X, alpha_vrai, beta_vrai, xi_vrai, pi_vrai, K = 2, L = 2, m = 5)

test_that("1. la sortie a les bonnes tailles", {
  expect_named(res, c("row_class", "col_class", "row_prob", "col_prob"))
  expect_length(res$row_class, 40)
  expect_length(res$col_class, 30)
  expect_equal(dim(res$row_prob), c(40, 2))
  expect_equal(dim(res$col_prob), c(30, 2))
})

test_that("2. les probabilités sont entre 0 et 1 et somment à 1", {
  expect_true(all(res$row_prob >= 0 & res$row_prob <= 1))
  expect_true(all(res$col_prob >= 0 & res$col_prob <= 1))
  expect_equal(rowSums(res$row_prob), rep(1, 40))
  expect_equal(rowSums(res$col_prob), rep(1, 30))
})

test_that("3. les groupes sont dans 1..K (1..L) et correspondent à la plus grande probabilité", {
  expect_true(all(res$row_class %in% 1:2))
  expect_true(all(res$col_class %in% 1:2))
  # le groupe choisi a une probabilité maximale
  expect_equal(res$row_prob[cbind(1:40, res$row_class)], apply(res$row_prob, 1, max))
  expect_equal(res$col_prob[cbind(1:30, res$col_class)], apply(res$col_prob, 1, max))
})

test_that("4. avec les vrais paramètres, on retrouve les vrais groupes", {
  expect_gt(concordance(res$row_class, sim$row_class), 0.9)
  expect_gt(concordance(res$col_class, sim$col_class), 0.9)
})

test_that("5. marche sans partition de départ, avec une partition donnée, et avec n_gibbs = 1", {
  # sans partition de départ (c'est le cas du res ci-dessus) : déjà testé.
  # avec une partition de départ donnée (la vraie) :
  res2 <- estimate_partition(sim$X, alpha_vrai, beta_vrai, xi_vrai, pi_vrai, K = 2, L = 2, m = 5,
                             row_class = sim$row_class, col_class = sim$col_class)
  expect_gt(concordance(res2$row_class, sim$row_class), 0.9)
  # un seul tirage gardé : les probabilités ne valent que 0 ou 1
  res3 <- estimate_partition(sim$X, alpha_vrai, beta_vrai, xi_vrai, pi_vrai, K = 2, L = 2, m = 5,
                             n_gibbs = 1, n_burn = 0)
  expect_true(all(res3$row_prob %in% c(0, 1)))
})