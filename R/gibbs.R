# FICHIER : gibbs.R
# RÔLE : étape SE-Gibbs de l'algorithme : tirer au hasard les classes des lignes,
#        puis celles des colonnes, sachant les paramètres du moment.
# COURS : partie 3 p35-36 (Stochastic EM within Gibbs) ; partie 3 p29-33 (hypothèses du LBM) ;
#         partie 1 p48 (formule de Bayes de l'étape E) ; article ordinalClust p5 (Algorithm 1).
# FONCTIONS :
#   gibbs_row(X, col_class, alpha, xi, pi, K, m)
#       sortie : list(row_class, row_prob)    row_prob : matrice n x K
#   gibbs_col(X, row_class, beta, xi, pi, L, m)
#       sortie : list(col_class, col_prob)    col_prob : matrice d x L
#   gibbs_lbm(X, row_class, col_class, alpha, beta, xi, pi, K, L, m)
#       sortie : list(row_class, col_class, row_prob, col_prob)
# NOTATIONS :
#   X : tableau n x d de notes entre 1 et m        K, L : nombre de classes de lignes / colonnes
#   row_class, col_class : classe de chaque ligne / colonne (v et w dans le cours)
#   alpha, beta : proportions des classes (longueurs K et L)
#   xi, pi : matrices K x L, un couple CUB par bloc (k, l)
# À RETENIR :
#   - on travaille en logarithmes : un produit de probabilités devient une somme
#   - on soustrait le maximum avant l'exponentielle pour ne pas obtenir des 0
#   - on TIRE la classe au hasard (sample), on ne prend pas la plus probable : c'est le "S" de SEM
#   - les paramètres ne changent pas ici : c'est mstep_lbm qui les met à jour
# UTILISE : cub_probability() (cub_probability.R)
# UTILISÉE PAR : sem_lbm.R, estimate_partition.R

##########################################################################################""
# CE QU'ON FAIT ICI 
# Dans le LBM, on ne connaît ni la partition des lignes v ni celle des colonnes w.
# La loi conjointe de (v, w) est impossible à calculer (K^n x L^d possibilités, partie 3),
# alors on les génère avec un échantillonneur de Gibbs : on tire v en gardant w fixé,
# puis w en gardant v fixé.
#
# Grâce aux hypothèses du modèle (lignes indépendantes, notes indépendantes dans un bloc),
# la probabilité qu'une ligne i soit dans la classe k, sachant w et les paramètres, vaut :
#
#   P(v_i = k | x, w) proportionnelle à  alpha_k  *  produit sur j de CUB( x_ij ; xi[k, w_j], pi[k, w_j] )
#
# et de même pour une colonne j dans la classe l, sachant v :
#
#   P(w_j = l | x, v) proportionnelle à  beta_l   *  produit sur i de CUB( x_ij ; xi[v_i, l], pi[v_i, l] )
#
# Une passe de gibbs_lbm :
#   1) score de chaque classe : log(alpha_k) + somme des log CUB    (logs : on évite les nombres trop petits)
#   2) on retire le maximum, on prend l'exponentielle, on divise par la somme -> probabilités (somme = 1)
#   3) on tire la classe au hasard avec ces probabilités
#   4) on fait pareil pour les colonnes, avec les NOUVELLES classes de lignes


# ------------------------------------------------------------------------------------------
# gibbs_row : tire la classe de chaque ligne, les colonnes étant fixées
# ------------------------------------------------------------------------------------------
gibbs_row <- function(X, col_class, alpha, xi, pi, K, m) {
  
  n <- nrow(X)
  d <- ncol(X)
  
  row_class <- integer(n)                      # classe tirée pour chaque ligne
  row_prob  <- matrix(0, nrow = n, ncol = K)   # probabilités d'appartenance (n x K)
  
  for (i in seq_len(n)) {
    
    log_score <- numeric(K)                    # un score (en log) par classe possible
    
    for (k in seq_len(K)) {
      
      # on part de la proportion de la classe k
      log_score[k] <- log(alpha[k])
      
      # on ajoute le log de la probabilité de chaque note de la ligne
      for (j in seq_len(d)) {
        l <- col_class[j]                      # classe de la colonne j (connue)
        
        p_x <- cub_probability(
          x  = X[i, j],
          m  = m,
          xi = xi[k, l],
          pi = pi[k, l]
        )
        
        log_score[k] <- log_score[k] + log(p_x)
      }
    }
    
    # on retire le plus grand score avant exp() pour éviter d'obtenir des 0
    log_score <- log_score - max(log_score)
    
    # retour aux probabilités : somme égale à 1
    prob <- exp(log_score)
    prob <- prob / sum(prob)
    
    row_prob[i, ] <- prob
    
    # tirage au hasard de la classe, selon ces probabilités
    row_class[i] <- sample(seq_len(K), size = 1, prob = prob)
  }
  
  return(
    list(
      row_class = row_class,
      row_prob  = row_prob
    )
  )
}


# ------------------------------------------------------------------------------------------
# gibbs_col : tire la classe de chaque colonne, les lignes étant fixées
# ------------------------------------------------------------------------------------------
# Même calcul que gibbs_row, mais cette fois c'est la classe k de la ligne i qui est
# connue, et on teste chaque classe l possible pour la colonne j.
gibbs_col <- function(X, row_class, beta, xi, pi, L, m) {
  
  n <- nrow(X)
  d <- ncol(X)
  
  col_class <- integer(d)                      # classe tirée pour chaque colonne
  col_prob  <- matrix(0, nrow = d, ncol = L)   # probabilités d'appartenance (d x L)
  
  for (j in seq_len(d)) {
    
    log_score <- numeric(L)                    # un score (en log) par classe possible
    
    for (l in seq_len(L)) {
      
      # on part de la proportion de la classe l
      log_score[l] <- log(beta[l])
      
      # on ajoute le log de la probabilité de chaque note de la colonne
      for (i in seq_len(n)) {
        k <- row_class[i]                      # classe de la ligne i (connue)
        
        p_x <- cub_probability(
          x  = X[i, j],
          m  = m,
          xi = xi[k, l],
          pi = pi[k, l]
        )
        
        log_score[l] <- log_score[l] + log(p_x)
      }
    }
    
    # même normalisation que pour les lignes
    log_score <- log_score - max(log_score)
    prob <- exp(log_score)
    prob <- prob / sum(prob)
    
    col_prob[j, ] <- prob
    col_class[j]  <- sample(seq_len(L), size = 1, prob = prob)
  }
  
  return(
    list(
      col_class = col_class,
      col_prob  = col_prob
    )
  )
}


# ------------------------------------------------------------------------------------------
# gibbs_lbm : un passage complet = lignes, puis colonnes
# ------------------------------------------------------------------------------------------
# On alterne : d'abord les lignes (avec les colonnes actuelles), puis les colonnes en
# utilisant les NOUVELLES classes de lignes (c'est ce qui fait un échantillonneur de Gibbs).
# L'argument row_class n'est pas utilisé dans le calcul : on le garde pour que la fonction
# ait les mêmes entrées que dans le schéma du projet.
gibbs_lbm <- function(X, row_class, col_class, alpha, beta, xi, pi, K, L, m) {
  
  # 1) classes des lignes
  row_result <- gibbs_row(
    X = X,
    col_class = col_class,
    alpha = alpha,
    xi = xi,
    pi = pi,
    K = K,
    m = m
  )
  
  # 2) classes des colonnes, avec les nouvelles classes de lignes
  col_result <- gibbs_col(
    X = X,
    row_class = row_result$row_class,
    beta = beta,
    xi = xi,
    pi = pi,
    L = L,
    m = m
  )
  
  return(
    list(
      row_class = row_result$row_class,
      col_class = col_result$col_class,
      row_prob  = row_result$row_prob,
      col_prob  = col_result$col_prob
    )
  )
}