# Tests de simulate_ordinal_lbm.R
# - dimensions du tableau, observations entières dans 1..m
# - aucun groupe vide, classes de la bonne longueur
# - les proportions de groupes s'approchent de alpha et beta
# - dans chaque bloc, cub_em retrouve les vrais (xi, pi)
# - erreurs sur les entrées invalides, avertissement si m < 4, reproductibilité

alpha_test <- c(0.4, 0.6)
beta_test <- c(0.5, 0.5)
xi_test <- matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2)   # xi[k, l]
pi_test <- matrix(c(0.9, 0.8, 0.85, 0.9), nrow = 2)  # pi[k, l]

test_that("simulate_ordinal_lbm : dimensions et valeurs du tableau", {
  set.seed(1)
  sim <- simulate_ordinal_lbm(60, 20, 5, alpha_test, beta_test, xi_test, pi_test)
  expect_equal(dim(sim$X), c(60, 20))
  expect_true(all(sim$X %in% 1:5))
  expect_false(anyNA(sim$X))
  expect_length(sim$row_class, 60)
  expect_length(sim$col_class, 20)
  expect_setequal(names(sim), c("X", "row_class", "col_class", "alpha", "beta", "xi", "pi"))
})

test_that("simulate_ordinal_lbm : aucun groupe vide, même avec peu de lignes", {
  set.seed(2)
  for (essai in 1:50){
    sim <- simulate_ordinal_lbm(6, 4, 5, c(0.1, 0.9), c(0.1, 0.9), xi_test, pi_test)
    expect_setequal(sim$row_class, 1:2)
    expect_setequal(sim$col_class, 1:2)
  }
})

test_that("simulate_ordinal_lbm : les proportions de groupes s'approchent de alpha et beta", {
  set.seed(3)
  sim <- simulate_ordinal_lbm(5000, 5000, 5, alpha_test, beta_test, xi_test, pi_test)
  expect_equal(as.numeric(table(sim$row_class)) / 5000, alpha_test, tolerance = 0.03)
  expect_equal(as.numeric(table(sim$col_class)) / 5000, beta_test, tolerance = 0.03)
})

test_that("simulate_ordinal_lbm : dans chaque bloc, cub_em retrouve les vrais (xi, pi)", {
  set.seed(4)
  sim <- simulate_ordinal_lbm(400, 100, 5, alpha_test, beta_test, xi_test, pi_test)
  for (k in 1:2) for (l in 1:2){
    bloc <- sim$X[sim$row_class == k, sim$col_class == l]
    est <- cub_em(as.vector(bloc), m = 5)
    expect_lt(abs(est$xi - xi_test[k, l]), 0.05)
    expect_lt(abs(est$pi - pi_test[k, l]), 0.05)
  }
})

test_that("simulate_ordinal_lbm : même graine, même résultat", {
  set.seed(5); a <- simulate_ordinal_lbm(30, 10, 5, alpha_test, beta_test, xi_test, pi_test)
  set.seed(5); b <- simulate_ordinal_lbm(30, 10, 5, alpha_test, beta_test, xi_test, pi_test)
  expect_identical(a, b)
})

test_that("simulate_ordinal_lbm : erreurs et avertissement", {
  expect_error(simulate_ordinal_lbm(60, 20, 5, c(0.5, 0.6), beta_test, xi_test, pi_test))    # alpha de somme 1,1
  expect_error(simulate_ordinal_lbm(60, 20, 5, alpha_test, beta_test, matrix(0.5, 3, 3), pi_test))  # mauvaise taille
  expect_error(simulate_ordinal_lbm(60, 20, 5, alpha_test, beta_test, xi_test, matrix(0, 2, 2)))    # pi = 0 interdit
  expect_error(simulate_ordinal_lbm(60, 20, 5, alpha_test, beta_test, matrix(1.5, 2, 2), pi_test))  # xi > 1
  expect_error(simulate_ordinal_lbm(1, 20, 5, alpha_test, beta_test, xi_test, pi_test))      # n < K
  expect_warning(simulate_ordinal_lbm(60, 20, 3, alpha_test, beta_test, xi_test, pi_test))   # m < 4
})