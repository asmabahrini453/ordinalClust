# FICHIER : cub_probability.R        RESPONSABLE (proposition) : Martine, à écrire en premier
# RÔLE : probabilité P(X = x) sous la loi CUB, pour un couple (xi, pi) et m modalités.
# COURS : partie 2 p58 (loi CUB = pi * binomiale décalée + (1 - pi) * uniforme).
# FONCTION : cub_probability(x, m, xi, pi)
#   x  : valeur(s) dans 1..m          m  : nombre de modalités
#   xi : paramètre feeling dans [0,1] pi : paramètre de mélange dans (0,1]
#   sortie : prob, une probabilité par valeur de x
# À RETENIR :
#   - la binomiale décalée est dbinom(x - 1, m - 1, 1 - xi) (voir la démonstration 3)
#   - pour la simulation et le calcul des blocs, on l'appelle avec x = 1:m (somme = 1)
#   - PRÉCAUTION : pi est un ARGUMENT, jamais la constante 3.14159 de R
#     -> vérifier pi dans (0, 1] et xi dans [0, 1], sinon erreur
# UTILISÉE PAR : cub_em.R, gibbs.R, compute_icl.R, simulate_ordinal_lbm.R
# ----------------------------------------------------------------------


# CE QU'ON FAIT ICI (explication pour l'oral)
# Le modèle CUB dit qu'une personne qui donne une note x entre 1 et m
#   - soit exprime son vrai ressenti ("feeling")      avec probabilité pi
#   - soit répond au hasard ("incertitude")           avec probabilité 1 - pi
#
# Donc  P(X = x) = pi * PB(x) + (1 - pi) * PU(x)      (cours, partie 2 p58)
#
#   PB(x) = binomiale décalée = probabilité "feeling" de la note x
#         = choose(m-1, x-1) * xi^(m-x) * (1-xi)^(x-1)
#         = dbinom(x - 1, m - 1, 1 - xi)     <- c'est la même chose, écrite avec R
#           (x - 1 va de 0 à m-1 : c'est ce que fait une binomiale ; "décalée" = le -1)
#   PU(x) = loi uniforme = 1/m                (toutes les notes ont la même chance)
#
# Effet de xi : petit xi -> notes hautes ; grand xi -> notes basses
#   (pour x = m, PB vaut (1 - xi)^(m-1), grand quand xi est petit).
# Effet de pi : proche de 1 -> les gens sont sûrs d'eux ; proche de 0 -> hasard pur.

cub_probability <- function(x, m, xi, pi){
  
  if (m < 2 || m != round(m)) stop("m doit etre un entier >= 2")
  if (pi <= 0 || pi > 1) stop("pi doit etre dans (0, 1]")
  if (xi < 0 || xi > 1) stop("xi doit etre dans [0, 1]")
  if (any(x < 1 | x > m | x != round(x))) stop("x doit contenir des entiers entre 1 et m")
  
  # ce sont les deux ingrédients du mélange
  binomiale_decalee <- dbinom(x - 1, m - 1, 1 - xi)   # PB(x) est le "feeling"
  uniforme <- 1 / m                                   # PU(x) est le "hasard"
  
  # le mélange : pi * feeling + (1 - pi) * hasard
  prob <- pi * binomiale_decalee + (1 - pi) * uniforme
  
  return(prob)
}