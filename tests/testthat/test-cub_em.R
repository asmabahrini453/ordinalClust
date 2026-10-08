# Tests de cub_em.R
# - sur des observations simulées avec (xi, pi) connus, l'estimation s'en approche
# - la log-vraisemblance ne diminue pas d'une itération à l'autre
# - cas limites : bloc vide, une seule valeur

# Simuler n observations CUB (même idée que pejSim de ordinalClust : sample() avec les
# probabilités du modèle)
simuler_cub <- function(n, m, xi, pi){
  return(sample(1:m, size = n, replace = TRUE, prob = cub_probability(1:m, m, xi, pi)))
}

test_that("cub_em : sur des observations simulées, l'estimation s'approche de (xi, pi)", {
  set.seed(1)
  cas <- list(c(0.3, 0.8, 5), c(0.7, 0.5, 7), c(0.1, 0.9, 10), c(0.5, 0.6, 4))   # (xi, pi, m)
  for (cc in cas){
    x <- simuler_cub(20000, cc[3], cc[1], cc[2])
    est <- cub_em(x, m = cc[3])
    expect_lt(abs(est$xi - cc[1]), 0.03)
    expect_lt(abs(est$pi - cc[2]), 0.03)
  }
})

test_that("cub_em : la log-vraisemblance ne diminue pas d'une itération à l'autre", {
  set.seed(1)
  x <- simuler_cub(500, 5, 0.4, 0.6)
  # tol = 0 : on force exactement k itérations, et on relève la log-vraisemblance obtenue
  ll <- sapply(1:40, function(k) cub_em(x, m = 5, max_iter = k, tol = 0)$loglik)
  expect_true(all(diff(ll) > -1e-9))
})

test_that("cub_em : la sortie contient xi, pi, loglik, iterations et converged", {
  set.seed(1)
  x <- simuler_cub(200, 5, 0.4, 0.6)
  expect_setequal(names(cub_em(x, m = 5)), c("xi", "pi", "loglik", "iterations", "converged"))
})

test_that("cub_em : cas limite du bloc vide (valeurs neutres, sans erreur)", {
  est <- cub_em(numeric(0), m = 5)
  expect_equal(est$iterations, 0)
  expect_false(est$converged)
  expect_equal(c(est$xi, est$pi), c(0.5, 0.5))
})

test_that("cub_em : cas limite d'une seule valeur (paramètres bornés, pas d'erreur)", {
  # toutes les observations égales
  est <- cub_em(rep(5, 10), m = 5)
  expect_true(est$xi >= 0.001 && est$xi <= 0.999)
  expect_true(est$pi >= 0.001 && est$pi <= 0.999)
  expect_true(est$converged)
  expect_true(is.finite(est$loglik))
  # une seule observation
  est <- cub_em(4, m = 5)
  expect_true(est$xi >= 0.001 && est$xi <= 0.999)
  expect_true(est$pi >= 0.001 && est$pi <= 0.999)
  expect_true(is.finite(est$loglik))
})

test_that("cub_em : erreur si une observation sort de 1..m", {
  expect_error(cub_em(c(1, 2, 7), m = 5))
})