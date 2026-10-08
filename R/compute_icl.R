# FICHIER : compute_icl.R
# RÔLE : calculer le critère ICL d'un co-clustering (pour choisir K et L : plus l'ICL est grand,
#        meilleur est le modèle).
# COURS : partie 3 p37 (ICL du LBM : vraisemblance complète moins une pénalité) ;
#         partie 3 p29-31 (hypothèses 1 et 2 : de quoi est faite la vraisemblance complète) ;
#         partie 2 p58 (loi CUB).
# FONCTION : compute_icl(X, row_class, col_class, alpha, beta, xi, pi, K, L, m)
#   X : matrice n x d des données ordinales       row_class, col_class : partitions finales
#   alpha, beta : proportions estimées des K groupes de lignes et des L groupes de colonnes
#   xi, pi : matrices K x L des paramètres CUB estimés      K, L, m : nombres de groupes et de modalités
#   sortie : ICL, un seul nombre
# À RETENIR :
#   - n <- nrow(X) et d <- ncol(X) se récupèrent dans la fonction
#   - ICL = ln p(x, v, w ; theta) - (K-1)/2 log n - (L-1)/2 log d - K L nu/2 log(n d), avec nu = 2 pour la CUB
#   - ln p(x, v, w ; theta) est la log-vraisemblance COMPLÈTE : elle utilise les groupes (row_class, col_class)
#   - on veut le plus grand ICL (et non le plus petit comme pour d'autres critères)
# UTILISÉE PAR : ordinal_coclust.R (pour comparer les n_init essais et les couples (K, L))

#################################
# CE QU'ON FAIT ICI 
# Le problème : la vraisemblance observée du LBM est impossible à calculer (il faudrait sommer
# sur les K^n x L^d partitions possibles), donc BIC aussi. Mais la vraisemblance COMPLÈTE,
# celle où l'on connaît les groupes, se calcule facilement. On l'évalue avec les groupes
# estimés (v, w) : c'est le principe de l'ICL.
#
# Vraisemblance complète (hypothèses 1 et 2 du LBM, partie 3 p29-31) :
#   p(x, v, w ; theta) = [produit sur les lignes de alpha_{v_i}] x [produit sur les colonnes de beta_{w_j}]
#                        x [produit sur les cases de la loi CUB du bloc (v_i, w_j)]
# donc, en passant au logarithme :
#   ln p = somme_k n_k ln(alpha_k) + somme_l d_l ln(beta_l)
#          + somme_{blocs (k,l)} somme_{observations du bloc} ln P(x | xi_kl, pi_kl)
#   avec n_k = nombre de lignes du groupe k et d_l = nombre de colonnes du groupe l.
#
# La pénalité (partie 3 p37) punit les modèles trop compliqués :
#   (K-1)/2 log n      : K - 1 proportions alpha libres, estimées avec n lignes
#   (L-1)/2 log d      : L - 1 proportions beta libres, estimées avec d colonnes
#   K L nu/2 log(n d)  : K x L blocs, nu paramètres par bloc (nu = 2 : xi et pi), estimés avec n x d cases
# Sans pénalité, plus de groupes donnerait toujours une meilleure vraisemblance.
###########################################################

compute_icl <- function(X, row_class, col_class, alpha, beta, xi, pi, K, L, m){
  
  n <- nrow(X)    # nombre de lignes
  d <- ncol(X)    # nombre de colonnes
  
  #vérif des entrées (des tailles qui ne collent pas donneraient un résultat faux)
  if (length(row_class) != n) stop("row_class doit avoir une valeur par ligne de X")
  if (length(col_class) != d) stop("col_class doit avoir une valeur par colonne de X")
  if (any(row_class < 1 | row_class > K)) stop("row_class doit contenir des groupes entre 1 et K")
  if (any(col_class < 1 | col_class > L)) stop("col_class doit contenir des groupes entre 1 et L")
  if (length(alpha) != K) stop("alpha doit avoir K valeurs")
  if (length(beta) != L) stop("beta doit avoir L valeurs")
  if (!is.matrix(xi) || any(dim(xi) != c(K, L))) stop("xi doit etre une matrice K x L")
  if (!is.matrix(pi) || any(dim(pi) != c(K, L))) stop("pi doit etre une matrice K x L")
  
  # Partie "groupes" de la log-vraisemblance complète : somme_k n_k ln(alpha_k) + somme_l d_l ln(beta_l)
  n_k <- tabulate(row_class, nbins = K)     # n_k[k] = nombre de lignes du groupe k
  d_l <- tabulate(col_class, nbins = L)     # d_l[l] = nombre de colonnes du groupe l
  #    On ne garde que les groupes non vides : un groupe vide a un effectif 0, et 0 * ln(0) donnerait NaN
  k_pleins <- n_k > 0
  l_pleins <- d_l > 0
  loglik_groupes <- sum(n_k[k_pleins] * log(alpha[k_pleins])) + sum(d_l[l_pleins] * log(beta[l_pleins]))
  
  #Partie "observations" : somme, sur chaque bloc, des ln P(x | xi_kl, pi_kl)
  loglik_blocs <- 0
  for (k in 1:K){
    for (l in 1:L){
      obs <- X[row_class == k, col_class == l]     # les observations du bloc (k, l)
      if (length(obs) > 0){                        # un bloc vide n'apporte rien
        effectif <- tabulate(obs, nbins = m)       # effectif[x] = combien d'observations valent x dans ce bloc
        proba <- cub_probability(1:m, m, xi[k, l], pi[k, l])   # P(x) pour x = 1..m, loi CUB du bloc
        utiles <- effectif > 0                     # on ignore les valeurs absentes (évite 0 * ln(0))
        loglik_blocs <- loglik_blocs + sum(effectif[utiles] * log(proba[utiles]))
      }
    }
  }
  loglik_complete <- loglik_groupes + loglik_blocs   # ln p(x, v, w ; theta)
  
  #pénalité (partie 3 p37), avec nu = 2 paramètres par bloc (xi et pi)
  nu <- 2
  penalite <- (K - 1) / 2 * log(n) + (L - 1) / 2 * log(d) + K * L * nu / 2 * log(n * d)
  
  #ICL = vraisemblance complète - pénalité (plus grand = meilleur)
  ICL <- loglik_complete - penalite
  return(ICL)
}