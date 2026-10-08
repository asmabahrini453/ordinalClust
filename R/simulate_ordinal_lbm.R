# FICHIER : simulate_ordinal_lbm.R
# RÔLE : fabriquer un faux tableau dont on connaît les vrais groupes et paramètres
#        (pour les tests, la page d'aide et la vignette).
# COURS : partie 3 p29-31 (hypothèses 1 et 2 du LBM) ; partie 2 p58 ;
#         manuel ordinalClust : jeu Msimulated ; article ordinalClust p4 (figure 2).
# FONCTION : simulate_ordinal_lbm(n, d, m, alpha, beta, xi, pi)
#   n : nombre de lignes     d : nombre de colonnes     m : nombre de modalités
#   alpha : proportions des K groupes de lignes (somme = 1)
#   beta  : proportions des L groupes de colonnes (somme = 1)
#   xi, pi : matrices K x L, les paramètres CUB de chaque bloc (k, l)
#   sortie : list(X, row_class, col_class, alpha, beta, xi, pi)   (vrais groupes et vrais paramètres)
# À RETENIR :
#   - tirer row_class avec alpha, col_class avec beta (hypothèse 1)
#   - tirer chaque observation dans cub_probability(1:m, m, xi[k, l], pi[k, l]) (hypothèse 2)
#   - prendre m >= 4 (avec m = 3 la CUB s'identifie mal, à vérifier avec le prof)
#   - aucun groupe vide
###########################################################################################
# CE QU'ON FAIT ICI 
# C'est le LBM lu "à l'envers" : au lieu d'estimer les groupes à partir des données,
# on choisit les groupes et les paramètres, puis on fabrique les données. On sait donc
# ce que l'algorithme devra retrouver.
#
#   Hypothèse 1 (les groupes) : chaque ligne i reçoit un groupe v_i = k avec probabilité alpha[k],
#     chaque colonne j reçoit un groupe w_j = l avec probabilité beta[l], indépendamment.
#     En R : row_class[i] = k et col_class[j] = l.
#
#   Hypothèse 2 (les observations) : sachant les groupes, l'observation x_ij ne dépend que
#     de son bloc (k, l) = (row_class[i], col_class[j]), et suit la loi CUB de ce bloc :
#         P(x_ij = x) = cub_probability(x, m, xi[k, l], pi[k, l])
#     Les observations sont tirées indépendamment les unes des autres.
#
# Pour un bloc (k, l) : on calcule les m probabilités avec cub_probability(1:m, ...)
# (leur somme vaut 1), puis on tire avec sample(1:m, prob = ...).
#
# Aucun groupe vide : si le tirage laisse un groupe sans ligne (ou sans colonne), on
# recommence. Sinon le bloc correspondant n'aurait aucune observation et on ne pourrait
# pas le retrouver. On tire donc "conditionnellement à n'avoir aucun groupe vide".


#' Simuler un tableau de notes ordinales selon le modèle à blocs latents
#'
#' Tire les groupes de lignes et de colonnes, puis une note par case selon la loi CUB
#' du bloc correspondant. Sert aux exemples, aux tests et à la vignette : on connaît la
#' vérité, donc on peut vérifier que \code{\link{ordinal_coclust}} la retrouve.
#'
#' @param n nombre de lignes.
#' @param d nombre de colonnes.
#' @param m nombre de modalités (prendre m >= 4).
#' @param alpha proportions des groupes de lignes (somme = 1).
#' @param beta proportions des groupes de colonnes (somme = 1).
#' @param xi,pi matrices K x L des paramètres CUB de chaque bloc.
#'
#' @return Une liste : \code{X} (le tableau), \code{row_class} et \code{col_class} (vrais
#'   groupes), \code{alpha}, \code{beta}, \code{xi}, \code{pi} (vrais paramètres).
#'
#' @examples
#' xi <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2)
#' pi <- matrix(0.9, nrow = 2, ncol = 2)
#' sim <- simulate_ordinal_lbm(20, 15, 5, c(.5, .5), c(.5, .5), xi, pi)
#' dim(sim$X)
#'
#' @export

simulate_ordinal_lbm <- function(n, d, m, alpha, beta, xi, pi){
  
  K <- length(alpha)    # nombre de groupes de lignes : il y a un alpha[k] par groupe
  L <- length(beta)     # nombre de groupes de colonnes : il y a un beta[l] par groupe
  
  # Vérifications des entrées (une erreur claire vaut mieux qu'un résultat faux)
  #    "n != round(n)" vérifie que n est entier : round(60) = 60 mais round(60.5) = 61
  if (n < K || n != round(n)) stop("n doit etre un entier >= K (nombre de groupes de lignes)")      # il faut au moins une ligne par groupe
  if (d < L || d != round(d)) stop("d doit etre un entier >= L (nombre de groupes de colonnes)")    # il faut au moins une colonne par groupe
  if (m < 2 || m != round(m)) stop("m doit etre un entier >= 2")
  if (m < 4) warning("m < 4 : avec peu de modalites la CUB s'identifie mal")                        # consigne de l'en-tête : simple avertissement
  #    alpha et beta sont des probabilités : positives et de somme 1. On tolère 1e-8 d'écart,
  #    car l'ordinateur arrondit (0.1 + 0.2 ne fait pas exactement 0.3)
  if (any(alpha < 0) || abs(sum(alpha) - 1) > 1e-8) stop("alpha doit etre positif et de somme 1")
  if (any(beta < 0) || abs(sum(beta) - 1) > 1e-8) stop("beta doit etre positif et de somme 1")
  #    xi et pi : une valeur par bloc (k, l), donc une matrice K x L (K lignes, L colonnes)
  if (!is.matrix(xi) || any(dim(xi) != c(K, L))) stop("xi doit etre une matrice K x L")
  if (!is.matrix(pi) || any(dim(pi) != c(K, L))) stop("pi doit etre une matrice K x L")
  if (any(xi < 0 | xi > 1)) stop("xi doit etre dans [0, 1]")      # mêmes domaines que dans cub_probability
  if (any(pi <= 0 | pi > 1)) stop("pi doit etre dans (0, 1]")
  
  # Hypothèse 1 (partie 3 p29-31) : tirer les groupes des lignes et des colonnes
  # Petite fonction interne : tire "taille" groupes avec les probabilités "proportions"
  # et recommence tant qu'un groupe est vide.
  tirer_groupes <- function(taille, proportions){
    nb_groupes <- length(proportions)
    for (essai in 1:1000){                                  # 1000 essais maximum, pour ne jamais tourner indéfiniment
      classe <- sample(1:nb_groupes, size = taille,         # un numéro de groupe (1..nb_groupes) par ligne (ou colonne)
                       replace = TRUE,                      # avec remise : plusieurs lignes peuvent avoir le même groupe
                       prob = proportions)                  # le groupe k est tiré avec la probabilité alpha[k] (ou beta[l])
      if (length(unique(classe)) == nb_groupes) return(classe)   # tous les groupes sont présents : on garde ce tirage
    }
    stop("impossible d'obtenir des groupes tous non vides : augmenter n (ou d) ou changer les proportions")
  }
  row_class <- tirer_groupes(n, alpha)    # row_class[i] = groupe de la ligne i     (v_i dans le cours)
  col_class <- tirer_groupes(d, beta)     # col_class[j] = groupe de la colonne j   (w_j dans le cours)
  
  # Hypothèse 2 (partie 3 p29-31) : tirer les observations, bloc par bloc
  X <- matrix(NA_integer_, nrow = n, ncol = d)   # tableau vide n x d ; NA permet de repérer une case oubliée
  for (k in 1:K){                                # on parcourt les K x L blocs
    for (l in 1:L){
      lignes <- which(row_class == k)            # numéros des lignes du groupe k
      colonnes <- which(col_class == l)          # numéros des colonnes du groupe l
      # P(x_ij = x) pour x = 1..m : la loi CUB du bloc (k, l)  (partie 2 p58). Somme = 1.
      probas <- cub_probability(1:m, m, xi[k, l], pi[k, l])
      # une observation par case du bloc : length(lignes) * length(colonnes) cases,
      # tirées indépendamment (replace = TRUE) dans 1..m avec les probabilités "probas"
      X[lignes, colonnes] <- sample(1:m, size = length(lignes) * length(colonnes),
                                    replace = TRUE, prob = probas)
    }
  }
  
  # le Resultat : le faux tableau, et la vérité (groupes et paramètres) pour pouvoir comparer
  return(list(X = X, row_class = row_class, col_class = col_class,
              alpha = alpha, beta = beta, xi = xi, pi = pi))
}
