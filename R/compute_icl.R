# FICHIER : compute_icl.R
#
# RÔLE : calculer le critère ICL d'un modèle de co-clustering ordinal.
#        L'ICL combine la log-vraisemblance complète, calculée à partir
#        des groupes estimés des lignes et des colonnes, et une pénalité
#        qui tient compte de la complexité du modèle.
#
# COURS : partie 3 p29-31 (vraisemblance complète du LBM : hypothèses 1 et 2) ;
#         partie 3 p37 (ICL : vraisemblance complète moins pénalité) ;
#         partie 2 p58 (loi CUB utilisée dans chaque bloc).
#
# FONCTION : compute_icl(X, row_class, col_class, alpha, beta, xi, pi, K, L, m)
#   X : matrice n x d des données ordinales
#   row_class : partition des n lignes, donnant le groupe de chaque ligne
#   col_class : partition des d colonnes, donnant le groupe de chaque colonne
#   alpha : proportions des K groupes de lignes
#   beta : proportions des L groupes de colonnes
#   xi, pi : paramètres de la loi CUB pour chacun des K x L blocs
#   K, L : nombres de groupes de lignes et de colonnes
#   m : nombre de modalités ordinales
#   sortie : valeur du critère ICL
#
# À RETENIR :
#   - n <- nrow(X) et d <- ncol(X)
#
#   - La vraisemblance complète du LBM est :
#       p(x, v, w ; theta) =
#         [produit_i alpha_{v_i}]
#         [produit_j beta_{w_j}]
#         [produit_{i,j} P(x_ij | xi_{v_i,w_j}, pi_{v_i,w_j})]
#
#   - Son logarithme est donc :
#       ln p(x, v, w ; theta) =
#         somme_k n_k ln(alpha_k)
#         + somme_l d_l ln(beta_l)
#         + somme_{k,l} somme_{i,j dans le bloc (k,l)}
#           ln P(x_ij | xi_kl, pi_kl)
#
#   - L'ICL est obtenu en retirant à cette log-vraisemblance une pénalité :
#       ICL = ln p(x, v, w ; theta)
#             - (K-1)/2 log(n)
#             - (L-1)/2 log(d)
#             - K L nu/2 log(n d)
#
#   - Les trois termes de pénalité correspondent respectivement :
#       (K-1) paramètres pour les proportions alpha des lignes ;
#       (L-1) paramètres pour les proportions beta des colonnes ;
#       K x L blocs contenant nu paramètres chacun.
#
#   - Pour la loi CUB, nu = 2 car chaque bloc possède deux paramètres :
#     xi et pi.
#
#   - Le meilleur modèle est celui qui maximise l'ICL :
#     une grande vraisemblance complète est recherchée, tout en pénalisant
#     les modèles trop complexes avec trop de groupes.
#
# UTILISÉE PAR : ordinal_coclust.R
#                pour comparer les différentes initialisations et les
#                différents couples (K, L), puis retenir le modèle
#                ayant le plus grand ICL.
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
