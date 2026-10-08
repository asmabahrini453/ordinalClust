# Tests de compute_icl.R
# - renvoie un seul nombre fini
# - le vrai (K, L) a un ICL plus grand que des (K, L) trop petits sur un faux tableau facile

# Ajuste (alpha, beta, xi, pi) pour des partitions données : proportions observées, puis cub_em
# sur chaque bloc (c'est ce que fait l'étape M du SEM, mais ici sans dépendre de mstep_lbm)
ajuster_blocs <- function(X, row_class, col_class, K, L, m){
  alpha <- tabulate(row_class, K) / nrow(X)
  beta <- tabulate(col_class, L) / ncol(X)
  xi <- matrix(0.5, K, L)
  pi <- matrix(0.5, K, L)
  for (k in 1:K) for (l in 1:L){
    obs <- as.vector(X[row_class == k, col_class == l])
    if (length(obs) > 0){
      est <- cub_em(obs, m)
      xi[k, l] <- est$xi
      pi[k, l] <- est$pi
    }
  }
  return(list(alpha = alpha, beta = beta, xi = xi, pi = pi))
}

# ICL d'un couple de partitions (ajustement des paramètres puis compute_icl)
icl_de <- function(X, row_class, col_class, K, L, m){
  a <- ajuster_blocs(X, row_class, col_class, K, L, m)
  return(compute_icl(X, row_class, col_class, a$alpha, a$beta, a$xi, a$pi, K, L, m))
}

# Faux tableau facile : 2 x 2 blocs aux paramètres très différents
xi_facile <- matrix(c(0.15, 0.85, 0.85, 0.15), nrow = 2)
pi_facile <- matrix(0.9, nrow = 2, ncol = 2)

test_that("compute_icl : renvoie un seul nombre fini", {
  set.seed(1)
  sim <- simulate_ordinal_lbm(100, 40, 5, c(0.5, 0.5), c(0.5, 0.5), xi_facile, pi_facile)
  icl <- compute_icl(sim$X, sim$row_class, sim$col_class, sim$alpha, sim$beta, sim$xi, sim$pi, 2, 2, 5)
  expect_type(icl, "double")
  expect_length(icl, 1)
  expect_true(is.finite(icl))
})

test_that("compute_icl : le vrai (K, L) bat des (K, L) trop petits sur un tableau facile", {
  n <- 100; d <- 40; m <- 5
  for (graine in 1:5){
    set.seed(graine)
    sim <- simulate_ordinal_lbm(n, d, m, c(0.5, 0.5), c(0.5, 0.5), xi_facile, pi_facile)
    icl_vrai <- icl_de(sim$X, sim$row_class, sim$col_class, 2, 2, m)
    icl_11 <- icl_de(sim$X, rep(1, n), rep(1, d), 1, 1, m)                  # tout dans un seul bloc
    icl_12 <- icl_de(sim$X, rep(1, n), sim$col_class, 1, 2, m)              # lignes fusionnées
    icl_21 <- icl_de(sim$X, sim$row_class, rep(1, d), 2, 1, m)              # colonnes fusionnées
    expect_gt(icl_vrai, icl_11)
    expect_gt(icl_vrai, icl_12)
    expect_gt(icl_vrai, icl_21)
  }
})

test_that("compute_icl : le vrai (K, L) bat aussi un (K, L) trop grand (la pénalité joue son rôle)", {
  n <- 100; d <- 40; m <- 5
  for (graine in 1:5){
    set.seed(graine)
    sim <- simulate_ordinal_lbm(n, d, m, c(0.5, 0.5), c(0.5, 0.5), xi_facile, pi_facile)
    # on coupe au hasard le groupe 1 de lignes et celui de colonnes en deux sous-groupes : (3, 3)
    lignes3 <- sim$row_class; lignes3[sim$row_class == 1 & runif(n) < 0.5] <- 3
    colonnes3 <- sim$col_class; colonnes3[sim$col_class == 1 & runif(d) < 0.5] <- 3
    expect_gt(icl_de(sim$X, sim$row_class, sim$col_class, 2, 2, m),
              icl_de(sim$X, lignes3, colonnes3, 3, 3, m))
  }
})

test_that("compute_icl : valeur calculée à la main pour K = L = 1", {
  # Un seul bloc : ICL = somme des ln P(x) + 0 (ln alpha = ln beta = 0) - (nu/2) ln(n d), avec nu = 2
  X <- matrix(c(1, 2, 3, 4, 5, 3), nrow = 2)   # n = 2, d = 3, m = 5
  attendu <- sum(log(cub_probability(as.vector(X), 5, 0.4, 0.8))) - (2 / 2) * log(2 * 3)
  obtenu <- compute_icl(X, c(1, 1), c(1, 1, 1), 1, 1, matrix(0.4), matrix(0.8), 1, 1, 5)
  expect_equal(obtenu, attendu)
})

test_that("compute_icl : un groupe vide ne donne ni NaN ni erreur", {
  X <- matrix(c(1, 2, 3, 4, 5, 3), nrow = 2)
  icl <- compute_icl(X, c(1, 1), c(1, 1, 1), c(1, 0), 1, matrix(c(0.4, 0.5), 2), matrix(c(0.8, 0.5), 2), 2, 1, 5)
  expect_true(is.finite(icl))
})

test_that("compute_icl : erreur si les tailles ne collent pas", {
  X <- matrix(c(1, 2, 3, 4, 5, 3), nrow = 2)
  xi <- matrix(0.4); pi <- matrix(0.8)
  expect_error(compute_icl(X, c(1, 1, 1), c(1, 1, 1), 1, 1, xi, pi, 1, 1, 5))   # row_class trop long
  expect_error(compute_icl(X, c(1, 1), c(1, 1), 1, 1, xi, pi, 1, 1, 5))         # col_class trop court
  expect_error(compute_icl(X, c(1, 2), c(1, 1, 1), 1, 1, xi, pi, 1, 1, 5))      # groupe 2 > K
  expect_error(compute_icl(X, c(1, 1), c(1, 1, 1), c(0.5, 0.5), 1, xi, pi, 1, 1, 5))   # alpha de mauvaise taille
  expect_error(compute_icl(X, c(1, 1), c(1, 1, 1), 1, 1, matrix(0.4, 2, 2), pi, 1, 1, 5))  # xi de mauvaise taille
})