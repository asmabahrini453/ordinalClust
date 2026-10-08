# FICHIER : test-gibbs.R
# RÔLE : vérifier gibbs_row(), gibbs_col() et gibbs_lbm() (étape SE-Gibbs, cours partie 3 p36).
# CE QU'ON TESTE :
#   1) les sorties ont la bonne forme et les probabilités somment à 1
#   2) quand les données ont une structure nette, Gibbs retrouve les bonnes classes
#   3) gibbs_lbm alterne lignes puis colonnes et corrige une partition de départ fausse
#
# LES DONNÉES DE TEST : un petit tableau 4 x 4 construit pour être sans ambiguïté.
#   lignes 1-2 avec colonnes 1-2 -> notes basses (1, 2)
#   lignes 3-4 avec colonnes 3-4 -> notes basses (1, 2)
#   les deux autres blocs        -> notes hautes (4, 5)
#   xi grand (0.8) = notes basses ; xi petit (0.2) = notes hautes

X_test <- matrix(
  c(1, 1, 5, 5,
    2, 1, 4, 5,
    5, 4, 1, 2,
    5, 5, 2, 1),
  nrow = 4, byrow = TRUE
)
xi_test    <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
pi_test    <- matrix(0.9, nrow = 2, ncol = 2)
alpha_test <- c(0.5, 0.5)
beta_test  <- c(0.5, 0.5)

# 1) gibbs_row : forme des sorties (4 lignes, 2 classes) et probabilités qui somment à 1.
test_that("gibbs_row renvoie des sorties bien formees", {
  set.seed(1)
  res <- gibbs_row(X_test, c(1, 1, 2, 2), alpha_test, xi_test, pi_test, K = 2, m = 5)
  
  expect_length(res$row_class, 4)                 # une classe par ligne
  expect_equal(dim(res$row_prob), c(4, 2))        # matrice n x K
  expect_true(all(res$row_class %in% 1:2))        # classes entre 1 et K
  expect_equal(as.vector(rowSums(res$row_prob)), rep(1, 4), tolerance = 1e-8)  # somme = 1
})

# 2) gibbs_row : sachant les bonnes classes de colonnes, il doit placer chaque ligne
#    dans sa vraie classe (probabilité > 0.9 pour la bonne classe).
test_that("gibbs_row retrouve les classes de lignes", {
  set.seed(1)
  res <- gibbs_row(X_test, c(1, 1, 2, 2), alpha_test, xi_test, pi_test, K = 2, m = 5)
  
  expect_gt(res$row_prob[1, 1], 0.9)   # ligne 1 -> classe 1
  expect_gt(res$row_prob[2, 1], 0.9)   # ligne 2 -> classe 1
  expect_gt(res$row_prob[3, 2], 0.9)   # ligne 3 -> classe 2
  expect_gt(res$row_prob[4, 2], 0.9)   # ligne 4 -> classe 2
  expect_equal(res$row_class, c(1, 1, 2, 2))   # et le tirage tombe sur la bonne classe
})

# 3) gibbs_col : même vérification de forme que pour les lignes (4 colonnes, 2 classes).
test_that("gibbs_col renvoie des sorties bien formees", {
  set.seed(1)
  res <- gibbs_col(X_test, c(1, 1, 2, 2), beta_test, xi_test, pi_test, L = 2, m = 5)
  
  expect_length(res$col_class, 4)
  expect_equal(dim(res$col_prob), c(4, 2))        # matrice d x L
  expect_true(all(res$col_class %in% 1:2))
  expect_equal(as.vector(rowSums(res$col_prob)), rep(1, 4), tolerance = 1e-8)
})

# 4) gibbs_col : sachant les bonnes classes de lignes, il doit retrouver les classes de colonnes.
test_that("gibbs_col retrouve les classes de colonnes", {
  set.seed(1)
  res <- gibbs_col(X_test, c(1, 1, 2, 2), beta_test, xi_test, pi_test, L = 2, m = 5)
  
  expect_gt(res$col_prob[1, 1], 0.9)   # colonne 1 -> classe 1
  expect_gt(res$col_prob[2, 1], 0.9)   # colonne 2 -> classe 1
  expect_gt(res$col_prob[3, 2], 0.9)   # colonne 3 -> classe 2
  expect_gt(res$col_prob[4, 2], 0.9)   # colonne 4 -> classe 2
  expect_equal(res$col_class, c(1, 1, 2, 2))
})

# 5) gibbs_lbm : renvoie bien les quatre éléments annoncés, et si on part de la bonne
#    partition, il la conserve.
test_that("gibbs_lbm renvoie les quatre sorties et garde une bonne partition", {
  set.seed(1)
  res <- gibbs_lbm(X_test, c(1, 1, 2, 2), c(1, 1, 2, 2),
                   alpha_test, beta_test, xi_test, pi_test, K = 2, L = 2, m = 5)
  
  expect_named(res, c("row_class", "col_class", "row_prob", "col_prob"))
  expect_equal(res$row_class, c(1, 1, 2, 2))
  expect_equal(res$col_class, c(1, 1, 2, 2))
})

# 6) L'alternance lignes/colonnes : on part d'une partition FAUSSE (1 2 1 2) et on fait
#    5 passages. Gibbs doit retrouver la bonne structure. On compare "à l'étiquetage près"
#    car les numéros de classes pourraient être échangés (1 1 2 2 ou 2 2 1 1).
test_that("gibbs_lbm corrige une partition de depart fausse", {
  set.seed(1)
  row_class <- c(1, 2, 1, 2)
  col_class <- c(1, 2, 1, 2)
  for (it in 1:5) {
    res <- gibbs_lbm(X_test, row_class, col_class,
                     alpha_test, beta_test, xi_test, pi_test, K = 2, L = 2, m = 5)
    row_class <- res$row_class
    col_class <- res$col_class
  }
  concordance <- function(a, b) max(mean(a == b), mean(a != b))
  expect_equal(concordance(row_class, c(1, 1, 2, 2)), 1)   # 1 = partition identique
  expect_equal(concordance(col_class, c(1, 1, 2, 2)), 1)
})