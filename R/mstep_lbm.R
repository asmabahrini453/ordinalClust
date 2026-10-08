# FICHIER : mstep_lbm.R
# RÔLE : étape M de l'algorithme : mettre à jour les paramètres du modèle
#        (alpha, beta, xi, pi) à partir des partitions tirées par Gibbs.
# COURS : partie 3 p35-36 (M step : on maximise la vraisemblance complète) ;
#         partie 3 p29-33 (hypothèses du LBM) ; partie 1 p49 (p_k = n_k / n).
# FONCTION : mstep_lbm(X, row_class, col_class, K, L, m)
#   sortie : list(alpha, beta, xi, pi)
#     alpha : vecteur de longueur K     beta : vecteur de longueur L
#     xi, pi : matrices K x L (un couple CUB par bloc)
# À RETENIR :
#   - la vraisemblance complète se découpe en trois morceaux indépendants,
#     donc on maximise chacun séparément
#   - alpha et beta : simples proportions
#   - xi et pi : on appelle cub_em() sur les notes de chaque bloc
#   - cas limite : une classe vide. On garde une proportion minimale (0.001)
#     pour que Gibbs puisse encore la retirer ; le bloc vide reçoit xi = pi = 0.5
#     (valeurs neutres renvoyées par cub_em). Choix pratique, absent du cours.
# UTILISE : cub_em() (cub_em.R)
# UTILISÉE PAR : sem_lbm.R

##########################################################################################""
# CE QU'ON FAIT ICI
# Quand on connaît les partitions v (lignes) et w (colonnes), la log-vraisemblance
# complète se découpe en trois parties (partie 3, hypothèses) :
#
#   ln p(x, v, w) = somme sur i de ln alpha_{v_i}                      (proportions des lignes)
#                 + somme sur j de ln beta_{w_j}                       (proportions des colonnes)
#                 + somme sur i,j de ln CUB( x_ij ; xi[v_i, w_j], pi[v_i, w_j] )   (notes de chaque bloc)
#
# Chaque partie ne dépend que de ses propres paramètres, on les maximise donc une par une :
#
#   alpha_k = (nombre de lignes dans la classe k) / n          (comme p_k = n_k / n, partie)
#   beta_l  = (nombre de colonnes dans la classe l) / d
#   (xi[k, l], pi[k, l]) = estimation CUB sur les notes du bloc (k, l)
#
# Pour les blocs, il n'y a pas de formule directe (la CUB est un mélange) : on utilise
# donc l'EM de cub_em() sur les notes du bloc, une fois pour chaque couple (k, l).
#
# Classe vide : si Gibbs vide une classe, sa proportion vaudrait 0 et log(0) = -Inf dans
# gibbs_row / gibbs_col : elle ne pourrait plus jamais être retirée. On impose donc une
# proportion minimale de 0.001, puis on renormalise pour que la somme reste égale à 1.


mstep_lbm <- function(X, row_class, col_class, K, L, m) {
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Proportions des classes de lignes et de colonnes
  # (tabulate compte aussi les classes vides : leur effectif vaut 0)
  alpha <- tabulate(row_class, nbins = K) / n
  beta  <- tabulate(col_class, nbins = L) / d
  
  # Cas limite : classe vide. On remonte les proportions nulles à 0.001,
  # puis on renormalise pour que alpha et beta somment toujours à 1.
  alpha <- pmax(alpha, 0.001)
  beta  <- pmax(beta,  0.001)
  alpha <- alpha / sum(alpha)
  beta  <- beta / sum(beta)
  
  # Matrices des paramètres CUB (une valeur par bloc)
  xi <- matrix(NA_real_, nrow = K, ncol = L)
  pi <- matrix(NA_real_, nrow = K, ncol = L)
  
  # Estimation des paramètres CUB pour chacun des K x L blocs
  for (k in seq_len(K)) {
    
    row_index <- which(row_class == k)        # lignes de la classe k
    
    for (l in seq_len(L)) {
      
      col_index <- which(col_class == l)      # colonnes de la classe l
      
      # bloc (k, l) : drop = FALSE garde une matrice même avec une seule ligne ou colonne
      # (si la classe est vide, le bloc est vide et cub_em renvoie 0.5 / 0.5)
      block <- X[row_index, col_index, drop = FALSE]
      block_values <- as.vector(block)
      
      # EM sur les notes de ce bloc
      cub_result <- cub_em(x = block_values, m = m)
      
      xi[k, l] <- cub_result$xi
      pi[k, l] <- cub_result$pi
    }
  }
  
  return(
    list(
      alpha = alpha,
      beta  = beta,
      xi    = xi,
      pi    = pi
    )
  )
}