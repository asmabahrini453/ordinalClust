# Tests de cub_probability.R
# - la somme des m probabilités vaut 1
# - toutes les probabilités sont positives
# - pi = 1 et xi quelconque donne la binomiale décalée ; pi proche de 0 donne presque l'uniforme
# - erreur si pi est hors de (0, 1] (protège contre la constante 3.14159 de R)

test_that("cub_probability : la somme des m probabilités vaut 1", {
  for (m in c(3, 5, 7, 10)) for (pi_cub in c(0.05, 0.3, 0.7, 1)) for (xi_cub in c(0, 0.2, 0.5, 0.8, 1)){
    expect_equal(sum(cub_probability(1:m, m, xi_cub, pi_cub)), 1)
  }
})

test_that("cub_probability : toutes les probabilités sont positives", {
  # pi < 1 : le terme "hasard" (1 - pi)/m est > 0, donc toutes les probabilités sont > 0
  for (pi_cub in c(0.05, 0.3, 0.7, 0.99)) for (xi_cub in c(0, 0.2, 0.5, 0.8, 1)){
    expect_true(all(cub_probability(1:5, 5, xi_cub, pi_cub) > 0))
  }
  # pi = 1 : plus de hasard, certaines probabilités peuvent valoir exactement 0 (jamais < 0)
  for (xi_cub in c(0, 0.2, 0.5, 0.8, 1)){
    expect_true(all(cub_probability(1:5, 5, xi_cub, 1) >= 0))
  }
})

test_that("cub_probability : pi = 1 donne la binomiale décalée, pi proche de 0 presque l'uniforme", {
  m <- 6
  # pi = 1 : P(x) = PB(x) = dbinom(x - 1, m - 1, 1 - xi), quel que soit xi
  for (xi_cub in c(0.1, 0.35, 0.5, 0.9)){
    expect_equal(cub_probability(1:m, m, xi_cub, 1), dbinom(0:(m - 1), m - 1, 1 - xi_cub))
  }
  # pi proche de 0 : P(x) est presque 1/m pour toutes les valeurs
  expect_equal(cub_probability(1:m, m, 0.3, 1e-6), rep(1 / m, m), tolerance = 1e-4)
})

test_that("cub_probability : erreur si pi est hors de (0, 1]", {
  expect_error(cub_probability(3, 5, 0.2, 0))      # pi = 0 est interdit
  expect_error(cub_probability(3, 5, 0.2, -0.3))   # pi négatif
  expect_error(cub_probability(3, 5, 0.2, 1.2))    # pi > 1
  # protège contre la constante 3.14159 de R : si on oublie de définir pi, R utilise pi = 3.14159...
  expect_error(cub_probability(3, 5, 0.2, pi))
})