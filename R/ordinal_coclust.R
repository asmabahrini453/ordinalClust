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
    stop("X doit contenir des entiers supérieurs ou égaux à 1.")
  }
  
  if (is.null(m)) {
    m <- max(X)                          # par défaut : plus grande note observée
  }
  if (m < 2 || m != round(m)) {
    stop("m doit être un entier supérieur ou égal à 2.")
  }
  if (max(X) > m) {
    stop("X contient des notes supérieures à m.")
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Petite fonction : TRUE si x est UN seul entier compris entre lower et upper
  entier_valide <- function(x, lower, upper) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) &&
      x == floor(x) && x >= lower && x <= upper
  }
  
  if (!entier_valide(K, 1, n)) {
    stop("K doit être un entier entre 1 et le nombre de lignes de X.")
  }
  if (!entier_valide(L, 1, d)) {
    stop("L doit être un entier entre 1 et le nombre de colonnes de X.")
  }
  if (!entier_valide(n_init, 1, Inf)) {
    stop("n_init doit être un entier supérieur ou égal à 1.")
  }
  if (!entier_valide(max_iter, 1, Inf)) {
    stop("max_iter doit être un entier supérieur ou égal à 1.")
  }
  if (!entier_valide(burn_in, 0, max_iter - 1)) {
    stop("burn_in doit être un entier entre 0 et max_iter - 1.")
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
    stop("Aucune initialisation n'a abouti (dernière raison : ", derniere_raison, ").")
  }
  
  if (n_echecs > 0) {
    warning(n_echecs, " initialisation(s) sur ", n_init,
            " écartée(s) (dernière raison : ", derniere_raison, ").")
  }
  
  if (!best$complet) {
    warning("Tous les essais ont au moins un groupe vide dans la partition finale : ",
            "K ou L est sans doute trop grand pour ces données.")
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