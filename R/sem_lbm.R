# FICHIER : sem_lbm.R
# RÔLE : boucle SEM-Gibbs : alterne tirage des partitions et mise à jour des paramètres,
#        puis renvoie les paramètres estimés (moyenne après le burn-in).
# COURS : partie 3 p35-36 (Stochastic EM within Gibbs) ; article ordinalClust p5 (Algorithm 1).
# FONCTION : sem_lbm(X, K, L, m, init, max_iter = 100, burn_in = 20)
#   init : list(row_class, col_class), partition de départ (renvoyée par initialize_lbm).
#          Si init contient aussi alpha, beta, xi, pi, ils servent de paramètres de départ.
#   sortie : list(alpha, beta, xi, pi, row_class, col_class, iterations, burn_in)
#     alpha, beta, xi, pi : paramètres estimés θ̂ (moyenne des itérations après le burn-in)
#     row_class, col_class : dernière partition tirée (point de départ de estimate_partition)
# À RETENIR :
#   - pas de critère d'arrêt tol : le SEM ne converge pas vers un point (la chaîne fluctue)
#   - on jette les burn_in premières itérations, puis on moyenne les paramètres
#   - une classe vide est gérée par mstep_lbm (proportion plancher de 0.001)
# UTILISE : gibbs_lbm() (gibbs.R), mstep_lbm() (mstep_lbm.R)
# UTILISÉE PAR : ordinal_coclust.R (une fois par initialisation aléatoire)

##########################################################################################""
# CE QU'ON FAIT ICI 
# On part d'une partition (v, w) donnée par initialize_lbm. Une première étape M
# donne les paramètres de départ θ(0). Puis, à chaque itération q :
#
#   SE-Gibbs : on tire v puis w avec les paramètres θ(q)           -> gibbs_lbm
#   M        : on met à jour θ(q) en θ(q+1) à partir de (v, w)    -> mstep_lbm
#
# Les burn_in premières itérations servent à "laisser la chaîne se mettre en place" :
# on ne les utilise pas. Pour les suivantes, on stocke les paramètres, et à la fin
#
#   θ̂ = moyenne des θ(q) pour q > burn_in
#
# C'est l'estimateur du cours ("θ̂ est obtenu à partir de la distribution
# d'échantillonnage après un burn-in", partie 3 p36).


sem_lbm <- function(
    X,
    K,
    L,
    m,
    init,
    max_iter = 100,
    burn_in = 20
) {
  
  if (burn_in >= max_iter) {
    stop("burn_in doit \u00eatre strictement inf\u00e9rieur \u00e0 max_iter.")
  }
  
  # Initialisation : partition de départ
  row_class <- init$row_class
  col_class <- init$col_class
  
  # Paramètres de départ θ(0) : si initialize_lbm ne les fournit pas,
  # on les calcule par une étape M sur la partition initiale
  if (is.null(init$alpha)) {
    theta0 <- mstep_lbm(X, row_class, col_class, K, L, m)
  } else {
    theta0 <- init
  }
  
  alpha <- theta0$alpha
  beta  <- theta0$beta
  xi    <- theta0$xi
  pi    <- theta0$pi
  
  # Nombre d'itérations conservées
  n_keep <- max_iter - burn_in
  
  # Stockage des paramètres après burn-in
  alpha_samples <- matrix(NA_real_, nrow = n_keep, ncol = K)
  beta_samples  <- matrix(NA_real_, nrow = n_keep, ncol = L)
  xi_samples    <- array(NA_real_, dim = c(n_keep, K, L))
  pi_samples    <- array(NA_real_, dim = c(n_keep, K, L))
  
  keep_index <- 0
  
  # Boucle SEM
  for (iter in seq_len(max_iter)) {
    
    # ----- SE-Gibbs : on tire les partitions avec les paramètres actuels -----
    
    gibbs_result <- gibbs_lbm(
      X = X,
      row_class = row_class,
      col_class = col_class,
      alpha = alpha,
      beta = beta,
      xi = xi,
      pi = pi,
      K = K,
      L = L,
      m = m
    )
    
    row_class <- gibbs_result$row_class
    col_class <- gibbs_result$col_class
    
    # ----- M-step : on met à jour les paramètres à partir de ces partitions -----
    
    mstep_result <- mstep_lbm(
      X = X,
      row_class = row_class,
      col_class = col_class,
      K = K,
      L = L,
      m = m
    )
    
    alpha <- mstep_result$alpha
    beta  <- mstep_result$beta
    xi    <- mstep_result$xi
    pi    <- mstep_result$pi
    
    # ----- Stockage après le burn-in -----
    
    if (iter > burn_in) {
      
      keep_index <- keep_index + 1
      
      alpha_samples[keep_index, ] <- alpha
      beta_samples[keep_index, ]  <- beta
      
      xi_samples[keep_index, , ] <- xi
      pi_samples[keep_index, , ] <- pi
    }
  }
  
  # Estimation finale des paramètres : moyenne des itérations conservées
  
  alpha_hat <- colMeans(alpha_samples)
  beta_hat  <- colMeans(beta_samples)
  
  xi_hat <- apply(xi_samples, c(2, 3), mean)
  pi_hat <- apply(pi_samples, c(2, 3), mean)
  
  return(
    list(
      alpha = alpha_hat,
      beta = beta_hat,
      xi = xi_hat,
      pi = pi_hat,
      row_class = row_class,       # dernière partition tirée (pour estimate_partition)
      col_class = col_class,
      iterations = max_iter,
      burn_in = burn_in
    )
  )
}
