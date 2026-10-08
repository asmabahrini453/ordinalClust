# FICHIER : test-mstep.R
# RÔLE : vérifier mstep_lbm() (étape M, cours partie 3 p36).
# CE QU'ON TESTE :
#   1) les paramètres sortent bien formés (somme des proportions = 1, xi et pi dans ]0,1[)
#   2) si on donne la VRAIE partition de données simulées, on retrouve les vrais paramètres
#   3) une classe vide ne fait pas planter la fonction (proportion plancher, bloc neutre)

# Petit simulateur local : tire chaque note d'un tableau selon la loi CUB de son bloc.
# (Le vrai simulateur du package viendra plus tard ; celui-ci sert seulement ici.)
simuler_test <- function(row_class, col_class, xi, pi, m) {
  X <- matrix(0, length(row_class), length(col_class))
  for (i in seq_along(row_class)) {
    for (j in seq_along(col_class)) {
      k <- row_class[i]; l <- col_class[j]
      X[i, j] <- sample(1:m, 1, prob = cub_probability(1:m, m, xi[k, l], pi[k, l]))
    }
  }
  X
}

# 1) Forme des sorties : 2 classes de lignes, 2 classes de colonnes.
test_that("mstep_lbm renvoie des parametres bien formes", {
  set.seed(1)
  m <- 5
  row_class <- rep(1:2, each = 20)    # 40 lignes : 20 dans la classe 1, 20 dans la classe 2
  col_class <- rep(1:2, each = 15)    # 30 colonnes : 15 et 15
  xi_vrai <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
  pi_vrai <- matrix(0.9, nrow = 2, ncol = 2)
  X <- simuler_test(row_class, col_class, xi_vrai, pi_vrai, m)
  
  res <- mstep_lbm(X, row_class, col_class, K = 2, L = 2, m = m)
  
  expect_equal(sum(res$alpha), 1)                       # proportions de lignes : somme 1
  expect_equal(sum(res$beta), 1)                        # proportions de colonnes : somme 1
  expect_equal(dim(res$xi), c(2, 2))                    # un xi par bloc
  expect_equal(dim(res$pi), c(2, 2))                    # un pi par bloc
  expect_true(all(res$xi > 0 & res$xi < 1))
  expect_true(all(res$pi > 0 & res$pi < 1))
})

# 2) Estimation : avec la vraie partition, l'étape M doit retrouver les vrais paramètres
#    (à 0.1 près, car chaque bloc ne contient que 300 notes).
test_that("mstep_lbm retrouve les vrais parametres avec la vraie partition", {
  set.seed(1)
  m <- 5
  row_class <- rep(1:2, each = 20)
  col_class <- rep(1:2, each = 15)
  xi_vrai <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
  pi_vrai <- matrix(0.9, nrow = 2, ncol = 2)
  X <- simuler_test(row_class, col_class, xi_vrai, pi_vrai, m)
  
  res <- mstep_lbm(X, row_class, col_class, K = 2, L = 2, m = m)
  
  expect_equal(res$alpha, c(0.5, 0.5))                  # 20 lignes sur 40 dans chaque classe
  expect_equal(res$beta, c(0.5, 0.5))                   # 15 colonnes sur 30 dans chaque classe
  expect_true(all(abs(res$xi - xi_vrai) < 0.1))
  expect_true(all(abs(res$pi - pi_vrai) < 0.1))
})

# 3) Classe vide : on met toutes les lignes dans la classe 1, la classe 2 est donc vide.
#    La fonction ne doit pas planter : la proportion de la classe 2 reste > 0 (plancher
#    de 0.001) et son bloc reçoit les valeurs neutres xi = pi = 0.5.
test_that("mstep_lbm gere une classe vide sans erreur", {
  set.seed(1)
  m <- 5
  row_class <- rep(1:2, each = 20)
  col_class <- rep(1:2, each = 15)
  xi_vrai <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
  pi_vrai <- matrix(0.9, nrow = 2, ncol = 2)
  X <- simuler_test(row_class, col_class, xi_vrai, pi_vrai, m)
  
  res <- mstep_lbm(X, rep(1, 40), col_class, K = 2, L = 2, m = m)
  
  expect_true(all(res$alpha > 0))                        # aucune proportion nulle
  expect_equal(sum(res$alpha), 1)                        # la somme reste 1
  expect_equal(as.vector(res$xi[2, ]), c(0.5, 0.5))      # bloc de la classe vide : neutre
})