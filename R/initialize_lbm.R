# FICHIER : initialize_lbm.R
# RÔLE : fabriquer une partition de départ aléatoire (lignes et colonnes) pour le SEM-Gibbs.
# COURS : partie 1 p50 (plusieurs initialisations aléatoires pour éviter les mauvais optimums) ;
#         partie 3 p35-36 (le SEM démarre d'une partition (v, w)) ;
#         article ordinalClust p5 (Algorithm 1, initialisation).
# FONCTION : initialize_lbm(X, K, L, m)
#   X : tableau n x d de notes       K, L : nombre de classes de lignes / colonnes
#   m : nombre de modalités (non utilisé ici, gardé pour avoir la même signature que le schéma)
#   sortie : list(row_class, col_class)
#     row_class : vecteur de longueur n, valeurs dans 1..K
#     col_class : vecteur de longueur d, valeurs dans 1..L
# À RETENIR :
#   - chaque classe reçoit au moins un élément : sinon mstep_lbm partirait d'une classe vide
#   - ensuite on mélange au hasard, pour que la partition ne dépende pas de l'ordre des lignes
#   - la fonction est aléatoire : chaque appel donne une partition différente
#     (ordinal_coclust l'appelle n_init fois, 20 par défaut)
# UTILISÉE PAR : ordinal_coclust.R, puis sem_lbm.R (argument init)

##########################################################################################""
# CE QU'ON FAIT ICI 
# Le SEM-Gibbs a besoin d'une partition de départ, mais on n'en connaît aucune.
# On en tire une au hasard. Comme le résultat dépend du point de départ (le SEM peut
# s'arrêter sur un mauvais optimum), on recommence plusieurs fois et on garde la
# meilleure solution (partie 1 p50).
#
# Pour ne pas tomber sur une classe vide dès le départ, on procède en deux temps :
#   1) on donne un élément à chaque classe : 1, 2, ..., K
#   2) on répartit les autres éléments au hasard dans les K classes
#   3) on mélange le tout


# Petite fonction interne : tire "taille" classes parmi 1..nb_classes,
# avec au moins un élément dans chaque classe.
tirer_classes <- function(taille, nb_classes) {
  
  if (taille < nb_classes) {
    stop("Il faut au moins autant de lignes (ou de colonnes) que de classes.")
  }
  
  # chaque classe apparaît au moins une fois, le reste est tiré au hasard
  classes <- c(
    seq_len(nb_classes),
    sample(seq_len(nb_classes), size = taille - nb_classes, replace = TRUE)
  )
  
  # on mélange pour que l'ordre ne donne aucune information
  classes[sample.int(taille)]
}


initialize_lbm <- function(X, K, L, m) {
  
  n <- nrow(X)    # nombre de lignes
  d <- ncol(X)    # nombre de colonnes
  
  row_class <- tirer_classes(n, K)
  col_class <- tirer_classes(d, L)
  
  return(
    list(
      row_class = row_class,
      col_class = col_class
    )
  )
}