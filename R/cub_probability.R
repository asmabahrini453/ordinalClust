# FICHIER : cub_probability.R      
# RÔLE :  Le modèle CUB dit qu'une personne qui donne une observation x entre 1 et m (où m est le nombre de modalités)
#   - soit exprime son vrai ressenti ("feeling")      avec probabilité pi
#   - soit répond au hasard ("incertitude")           avec probabilité 1 - pi
#
# Donc la formule de CUB est :  P(X = x) = pi * PB(x) + (1 - pi) * PU(x)      (cours, partie 2 p58)
#
#   PB(x) = binomiale décalée = probabilité "feeling" de l'observation x
#         = (m-1, x-1) * xi^(m-x) * (1-xi)^(x-1)
#         = dbinom(x - 1, m - 1, 1 - xi)  
#           (x - 1 va de 0 à m-1 : c'est ce que fait une binomiale ; "décalée" = le -1)
#   PU(x) = loi uniforme = 1/m                (toutes les modalités ont la même chance)
#Précautions:
# Effet de xi : petit xi -> observations hautes ; grand xi -> observations basses
#   (pour x = m, PB vaut (1 - xi)^(m-1), grand quand xi est petit).
# Effet de pi : proche de 1 -> on est sur ; proche de 0 -> hasard pur.

# UTILISÉE PAR : cub_em.R, gibbs.R, compute_icl.R, simulate_ordinal_lbm.R
#########################################################################################

#' @importFrom stats dbinom

cub_probability <- function(x, m, xi, pi){
  
  if (m < 2 || m != round(m)) stop("m doit etre un entier >= 2")
  if (pi <= 0 || pi > 1) stop("pi doit etre dans (0, 1]")
  if (xi < 0 || xi > 1) stop("xi doit etre dans [0, 1]")
  if (any(x < 1 | x > m | x != round(x))) stop("x doit contenir des entiers entre 1 et m")
  
  # ce sont les deux ingrédients du mélange
  binomiale_decalee <- dbinom(x - 1, m - 1, 1 - xi)   # PB(x) est le "feeling"
  uniforme <- 1 / m                                   # PU(x) est le "hasard"
  
  # le mélange :c'est la formule de cub 
  prob <- pi * binomiale_decalee + (1 - pi) * uniforme
  
  return(prob)
}
