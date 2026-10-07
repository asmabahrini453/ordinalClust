# FICHIER : sem_lbm.R                RESPONSABLE (proposition) : Martine
# RÔLE : boucle SEM-Gibbs d'un run (étapes 2 à 5 du schéma).
# COURS : partie 3 p35-36 ; article ordinalClust p5 (burn-in, moyenne des paramètres).
# FONCTION : sem_lbm(X, K, L, m, init, max_iter = 100, burn_in = 20)
#   init : liste produite par initialize_lbm()
#   sortie : list(row_class, col_class, alpha, beta, xi, pi, iterations)
# À RETENIR (corrections du schéma) :
#   - SEM ne converge pas vers un point : on fait max_iter itérations, sans critère d'arrêt sur tol
#   - les paramètres finaux (theta chapeau) sont la MOYENNE des itérations après le burn-in
#   - à chaque itération : gibbs_lbm(), puis mstep_lbm()
#   - garder l'état final (row_class, col_class) pour estimate_partition()
