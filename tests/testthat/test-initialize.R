# FICHIER : test-initialize.R
# RÔLE : vérifier initialize_lbm(), la partition de départ aléatoire.
# CE QU'ON TESTE :
#   1) les partitions ont la bonne longueur et des classes entre 1 et K (ou L)
#   2) aucune classe n'est vide au départ (sinon mstep_lbm partirait d'un bloc vide)
#   3) la fonction est bien aléatoire (deux appels donnent deux partitions différentes)
#   4) un cas impossible (plus de classes que de lignes) donne une erreur claire

# 1) Forme des sorties : 40 lignes, 30 colonnes, K = 3 classes de lignes, L = 2 classes de colonnes.
test_that("initialize_lbm renvoie des partitions de la bonne longueur", {
  X <- matrix(sample(1:5, 40 * 30, replace = TRUE), 40, 30)
  init <- initialize_lbm(X, K = 3, L = 2, m = 5)
  
  expect_length(init$row_class, 40)           # une classe par ligne
  expect_length(init$col_class, 30)           # une classe par colonne
  expect_true(all(init$row_class %in% 1:3))   # classes de lignes dans 1..K
  expect_true(all(init$col_class %in% 1:2))   # classes de colonnes dans 1..L
})

# 2) Aucune classe vide : on répète 20 fois pour être sûr (le tirage est aléatoire).
test_that("initialize_lbm n'a aucune classe vide", {
  X <- matrix(sample(1:5, 40 * 30, replace = TRUE), 40, 30)
  for (i in 1:20) {
    init <- initialize_lbm(X, K = 5, L = 4, m = 5)
    expect_equal(length(unique(init$row_class)), 5)   # les 5 classes sont présentes
    expect_equal(length(unique(init$col_class)), 4)   # les 4 classes sont présentes
  }
})

# 3) Le hasard : deux appels successifs ne doivent pas donner la même partition.
test_that("initialize_lbm est aleatoire", {
  X <- matrix(sample(1:5, 40 * 30, replace = TRUE), 40, 30)
  a <- initialize_lbm(X, K = 3, L = 2, m = 5)
  b <- initialize_lbm(X, K = 3, L = 2, m = 5)
  expect_false(identical(a$row_class, b$row_class))
})

# 4) Cas limite : 50 classes pour 40 lignes est impossible, on attend un message d'erreur.
test_that("initialize_lbm refuse plus de classes que de lignes", {
  X <- matrix(sample(1:5, 40 * 30, replace = TRUE), 40, 30)
  expect_error(initialize_lbm(X, K = 50, L = 2, m = 5), "au moins autant")
})