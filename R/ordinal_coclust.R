# FICHIER : ordinal_coclust.R
# RÔLE : fonction principale du package. Co-clustering de données ordinales :
#        modèle à blocs latents (LBM), loi CUB dans chaque bloc, inférence SEM-Gibbs.
# COURS : énoncé MBL-Projet-2027 ; partie 1 p50 (plusieurs départs) ; partie 3 p35-37.
# FONCTION : ordinal_coclust(X, K, L, n_init = 20, max_iter = 100, burn_in = 20, m = NULL)
#   X      : tableau n x d de notes entières entre 1 et m (matrice ou data.frame)
#   K, L   : nombre de classes de lignes / de colonnes
#   n_init : nombre d'initialisations aléatoires (20 par défaut)
#   max_iter, burn_in : réglages du SEM (voir sem_lbm)
#   m      : nombre de modalités ; par défaut le maximum de X
#   sortie : list(row_prob, col_prob, row_class, col_class,
#                 parameters = list(alpha, beta, xi, pi), ICL, K, L, iterations, best_init)
# À RETENIR :
#   - pour chaque départ : initialize_lbm -> sem_lbm -> estimate_partition -> compute_icl
#   - on garde le départ dont l'ICL est le plus GRAND (c'est un critère à maximiser)
#   - un départ qui plante est écarté (avec un avertissement) ; un départ sans groupe vide
#     passe toujours avant un départ avec groupe vide
#   - les paramètres utilisés pour l'ICL sont les theta estimés par le SEM,
#     la partition est celle de estimate_partition
# UTILISE : initialize_lbm, sem_lbm, estimate_partition, compute_icl

##########################################################################################
# CE QU'ON FAIT ICI (explication pour l'oral)
# Le SEM-Gibbs dépend de sa partition de départ : il peut s'arrêter sur une mauvaise
# solution. On le lance donc n_init fois à partir de départs aléatoires différents
# (partie 1 p50), et on garde la meilleure solution. Pour comparer, on utilise l'ICL
# (partie 3 p37), qui mesure l'ajustement du modèle en pénalisant sa complexité :
# plus l'ICL est grand, meilleur est le résultat.
#
# Un run complet :
#   1) initialize_lbm     : partition de départ au hasard
#   2) sem_lbm            : estime theta (moyenne après burn-in)
#   3) estimate_partition : partition finale et probabilités, theta fixé
#   4) compute_icl        : note du run
#
# Deux sécurités :
#   - si un run plante, on le saute et on passe au suivant (tryCatch), mais on prévient
#     avec un avertissement, pour ne pas cacher un vrai problème ;
#   - si un run finit avec un groupe vide, on le garde en réserve : un run "complet" (sans
#     groupe vide) passe toujours avant. Si tous les runs ont un groupe vide, on garde le
#     meilleur ICL avec un avertissement (K ou L est sans doute trop grand). L'ICL d'un tel
#     run reste calculable et il est pénalisé pour ses K (ou L) groupes, donc plus bas :
#     c'est ce qui permet de comparer plusieurs valeurs de K et L.
#   On ne s'arrête avec une erreur que si TOUS les runs ont planté.


#' Co-clustering de données ordinales (modèle à blocs latents et loi CUB)
#'
#' Regroupe simultanément les lignes et les colonnes d'un tableau de notes ordinales.
#' Dans chaque bloc (groupe de lignes, groupe de colonnes), les notes suivent une loi CUB.
#' Les paramètres sont estimés par un algorithme SEM-Gibbs lancé à partir de plusieurs
#' partitions aléatoires ; on garde la solution qui a le plus grand critère ICL.
#'
#' @param X matrice (ou data.frame) n x d de notes entières entre 1 et m, sans valeur manquante.
#' @param K nombre de groupes de lignes.
#' @param L nombre de groupes de colonnes.
#' @param n_init nombre d'initialisations aléatoires (20 par défaut).
#' @param max_iter nombre d'itérations du SEM-Gibbs (100 par défaut).
#' @param burn_in nombre d'itérations de chauffe, ignorées dans la moyenne des paramètres
#'   (20 par défaut ; doit être strictement inférieur à \code{max_iter}).
#' @param m nombre de modalités. Par défaut, la plus grande note observée dans \code{X}.
#'
#' @return Une liste avec :
#' \describe{
#'   \item{row_prob}{matrice n x K des probabilités d'appartenance des lignes aux groupes.}
#'   \item{col_prob}{matrice d x L des probabilités d'appartenance des colonnes aux groupes.}
#'   \item{row_class}{groupe de chaque ligne (le plus probable).}
#'   \item{col_class}{groupe de chaque colonne (le plus probable).}
#'   \item{parameters}{liste des paramètres estimés : \code{alpha} et \code{beta}
#'     (proportions des groupes), \code{xi} et \code{pi} (matrices K x L des paramètres CUB).}
#'   \item{ICL}{valeur du critère ICL de la meilleure initialisation (plus grand = meilleur).}
#'   \item{K, L}{nombres de groupes demandés.}
#'   \item{iterations}{nombre d'itérations du SEM.}
#'   \item{best_init}{numéro de l'initialisation retenue.}
#' }
#'
#' @details Les numéros de groupes sont arbitraires (\emph{label switching}) : le groupe 1
#'   d'un résultat peut correspondre au groupe 2 d'un autre. Un avertissement est émis si des
#'   initialisations ont échoué, ou si tous les essais ont un groupe vide (K ou L trop grand).
#'
#' @examples
#' xi <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2)
#' pi <- matrix(0.9, nrow = 2, ncol = 2)
#' set.seed(1)
#' sim <- simulate_ordinal_lbm(n = 40, d = 30, m = 5, alpha = c(.5, .5),
#'                             beta = c(.5, .5), xi = xi, pi = pi)
#' res <- ordinal_coclust(sim$X, K = 2, L = 2, n_init = 3, max_iter = 40, burn_in = 10)
#' res$ICL
#' table(res$row_class, sim$row_class)
#'
#' @export

ordinal_coclust <- function(X, K, L, n_init = 20, max_iter = 100, burn_in = 20, m = NULL) {
  
  # ---------- vérification des entrées ----------
  
  X <- as.matrix(X)    # accepte aussi un data.frame
  
  if (!is.numeric(X)) {
    stop("X doit contenir des nombres.")
  }
  if (anyNA(X)) {
    stop("X ne doit pas contenir de valeurs manquantes (NA).")
  }
  if (any(X != round(X)) || any(X < 1)) {
    stop("X doit contenir des entiers sup\u00e9rieurs ou \u00e9gaux \u00e0 1.")
  }
  
  if (is.null(m)) {
    m <- max(X)                          # par défaut : plus grande note observée
  }
  if (m < 2 || m != round(m)) {
    stop("m doit \u00eatre un entier sup\u00e9rieur ou \u00e9gal \u00e0 2.")
  }
  if (max(X) > m) {
    stop("X contient des notes sup\u00e9rieures \u00e0 m.")
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Petite fonction : TRUE si x est UN seul entier compris entre lower et upper
  entier_valide <- function(x, lower, upper) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) &&
      x == floor(x) && x >= lower && x <= upper
  }
  
  if (!entier_valide(K, 1, n)) {
    stop("K doit \u00eatre un entier entre 1 et le nombre de lignes de X.")
  }
  if (!entier_valide(L, 1, d)) {
    stop("L doit \u00eatre un entier entre 1 et le nombre de colonnes de X.")
  }
  if (!entier_valide(n_init, 1, Inf)) {
    stop("n_init doit \u00eatre un entier sup\u00e9rieur ou \u00e9gal \u00e0 1.")
  }
  if (!entier_valide(max_iter, 1, Inf)) {
    stop("max_iter doit \u00eatre un entier sup\u00e9rieur ou \u00e9gal \u00e0 1.")
  }
  if (!entier_valide(burn_in, 0, max_iter - 1)) {
    stop("burn_in doit \u00eatre un entier entre 0 et max_iter - 1.")
  }
  
  # ---------- n_init essais, on garde le meilleur ----------
  
  best <- NULL
  best_init <- NA_integer_
  n_echecs <- 0               # nombre de runs qui ont planté
  derniere_raison <- ""       # raison de la dernière panne (pour le message)
  
  for (r in seq_len(n_init)) {
    
    # Un run complet. En cas d'erreur, tryCatch renvoie le TEXTE de l'erreur
    # (sinon il renvoie la liste des résultats du run).
    run <- tryCatch({
      
      # 1) partition de départ aléatoire
      init <- initialize_lbm(X, K, L, m)
      
      # 2) SEM-Gibbs : estimation de theta
      sem <- sem_lbm(X, K, L, m, init, max_iter, burn_in)
      
      # 3) partition finale et probabilités, avec theta fixé
      part <- estimate_partition(
        X = X,
        row_class = sem$row_class,
        col_class = sem$col_class,
        alpha = sem$alpha,
        beta = sem$beta,
        xi = sem$xi,
        pi = sem$pi,
        K = K,
        L = L,
        m = m
      )
      
      # un groupe final vide ? On le note, mais on ne plante pas
      complet <- all(tabulate(part$row_class, nbins = K) > 0) &&
        all(tabulate(part$col_class, nbins = L) > 0)
      
      # 4) note du run (ICL : plus grand = meilleur)
      icl <- compute_icl(
        X = X,
        row_class = part$row_class,
        col_class = part$col_class,
        alpha = sem$alpha,
        beta = sem$beta,
        xi = sem$xi,
        pi = sem$pi,
        K = K,
        L = L,
        m = m
      )
      
      if (!is.finite(icl)) {
        stop("ICL non calculable")
      }
      
      list(sem = sem, part = part, icl = icl, complet = complet)
      
    }, error = function(e) conditionMessage(e))
    
    # run planté : on le compte et on passe au suivant
    if (is.character(run)) {
      n_echecs <- n_echecs + 1
      derniere_raison <- run
      next
    }
    
    # On garde ce run s'il est meilleur. Un run complet (sans groupe vide) passe toujours
    # avant un run avec groupe vide ; à complétude égale, le plus grand ICL gagne.
    if (is.null(best) ||
        (run$complet && !best$complet) ||
        (run$complet == best$complet && run$icl > best$icl)) {
      best <- run
      best_init <- r
    }
  }
  
  # ---------- bilan ----------
  
  if (is.null(best)) {
    stop("Aucune initialisation n'a abouti (derni\u00e8re raison : ", derniere_raison, ").")
  }
  
  if (n_echecs > 0) {
    warning(n_echecs, " initialisation(s) sur ", n_init,
            " \u00e9cart\u00e9e(s) (derni\u00e8re raison : ", derniere_raison, ").")
  }
  
  if (!best$complet) {
    warning("Tous les essais ont au moins un groupe vide dans la partition finale : ",
            "K ou L est sans doute trop grand pour ces donn\u00e9es.")
  }
  
  return(
    list(
      row_prob   = best$part$row_prob,
      col_prob   = best$part$col_prob,
      row_class  = best$part$row_class,
      col_class  = best$part$col_class,
      parameters = list(
        alpha = best$sem$alpha,
        beta  = best$sem$beta,
        xi    = best$sem$xi,
        pi    = best$sem$pi
      ),
      ICL        = best$icl,
      K          = K,
      L          = L,
      iterations = best$sem$iterations,
      best_init  = best_init
    )
  )
}
