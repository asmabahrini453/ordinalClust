# FICHIER : estimate_partition.R
# RÔLE : une fois les paramètres estimés par sem_lbm, trouver la partition finale des lignes et
#        des colonnes, et les probabilités d'appartenance à chaque groupe.
# COURS : partie 3 p35-36 (SEM-Gibbs : la partition finale se tire avec theta fixé) ;
#         partie 3 p29-33 (hypothèses du LBM) ; partie 2 p58 (loi CUB).
# FONCTION : estimate_partition(X, alpha, beta, xi, pi, K, L, m,
#                               n_gibbs = 50, n_burn = 10, row_class = NULL, col_class = NULL)
#   X : matrice n x d des données ordinales
#   alpha, beta : proportions estimées (longueurs K et L)
#   xi, pi : matrices K x L des paramètres CUB estimés
#   K, L, m : nombres de groupes de lignes, de colonnes, et de modalités
#   n_gibbs : nombre de tirages gardés       n_burn : nombre de tirages de chauffe jetés
#   row_class, col_class : partition de départ (optionnelle, par exemple celle de sem_lbm)
#   sortie : list(row_class, col_class, row_prob, col_prob)
# À RETENIR :
#   - theta est FIXÉ pendant toute la fonction : on ne refait plus de M-step
#   - on enchaîne gibbs_lbm (tirer v | w, puis w | v) et on COMPTE les groupes tirés
#   - fréquence d'un groupe = probabilité d'appartenance (row_prob, col_prob)
#   - groupe final = celui qui a la plus grande probabilité
# UTILISÉE PAR : ordinal_coclust.R (une fois par initialisation, après sem_lbm)

#################################
# CE QU'ON FAIT ICI
# Après sem_lbm on a une bonne estimation theta = (alpha, beta, xi, pi), mais la partition
# n'est pas "estimée" par l'algorithme : à chaque itération elle est tirée au hasard.
# Le cours propose donc de la fixer ainsi (partie 3 p36) :
#   1. on garde theta fixe ;
#   2. on lance Gibbs : on tire les groupes des lignes sachant les colonnes, puis les groupes
#      des colonnes sachant les lignes, et on recommence ;
#   3. on jette les premiers tirages (chauffe, "burn-in"), car la chaîne part d'une partition
#      arbitraire et n'est pas encore stable ;
#   4. sur les tirages gardés, on compte combien de fois chaque ligne est tombée dans chaque
#      groupe : cette fréquence est la probabilité d'appartenance ;
#   5. le groupe final d'une ligne est celui où elle est tombée le plus souvent.
# Même principe pour les colonnes.
###########################################################

estimate_partition <- function(X, alpha, beta, xi, pi, K, L, m,
                               n_gibbs = 50, n_burn = 10,
                               row_class = NULL, col_class = NULL){
  
  n <- nrow(X)    # nombre de lignes
  d <- ncol(X)    # nombre de colonnes
  
  # Partition de départ : si on n'en donne pas, on en tire une au hasard
  # (rep(..., length.out = n) met chaque groupe au moins une fois, sample mélange)
  if (is.null(row_class)) row_class <- sample(rep(seq_len(K), length.out = n))
  if (is.null(col_class)) col_class <- sample(rep(seq_len(L), length.out = d))
  
  # Compteurs : row_counts[i, k] = nombre de fois où la ligne i est tombée dans le groupe k
  row_counts <- matrix(0, nrow = n, ncol = K)
  col_counts <- matrix(0, nrow = d, ncol = L)
  
  # Boucle de Gibbs : n_burn tirages jetés, puis n_gibbs tirages comptés
  for (iter in seq_len(n_burn + n_gibbs)){
    
    # un tour de Gibbs : v | w puis w | v, avec theta fixé
    gibbs_result <- gibbs_lbm(X = X, row_class = row_class, col_class = col_class,
                              alpha = alpha, beta = beta, xi = xi, pi = pi,
                              K = K, L = L, m = m)
    row_class <- gibbs_result$row_class
    col_class <- gibbs_result$col_class
    
    # après la chauffe, on compte les groupes tirés
    if (iter > n_burn){
      for (i in seq_len(n)) row_counts[i, row_class[i]] <- row_counts[i, row_class[i]] + 1
      for (j in seq_len(d)) col_counts[j, col_class[j]] <- col_counts[j, col_class[j]] + 1
    }
  }
  
  # Fréquences = probabilités d'appartenance (chaque ligne de row_prob somme à 1)
  row_prob <- row_counts / n_gibbs
  col_prob <- col_counts / n_gibbs
  
  # Groupe final = groupe le plus fréquent (en cas d'égalité, le premier)
  list(row_class = max.col(row_prob, ties.method = "first"),
       col_class = max.col(col_prob, ties.method = "first"),
       row_prob  = row_prob,
       col_prob  = col_prob)
}