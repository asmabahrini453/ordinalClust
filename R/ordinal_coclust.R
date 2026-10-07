# FICHIER : ordinal_coclust.R        RESPONSABLE : les deux (assemblage vendredi)
# RÔLE : fonction principale du package (interface utilisateur), étapes 1 à 8.
# COURS : énoncé MBL-Projet-2027 ; partie 1 p50 (plusieurs départs) ; partie 3 p35-37.
# FONCTION : ordinal_coclust(X, K, L, n_init = 20, max_iter = 100, burn_in = 20, m = NULL)
#   m = NULL : par défaut max(X) ; vérifier que X contient des entiers de 1 à m
#   pour chaque départ : initialize_lbm -> sem_lbm -> estimate_partition -> compute_icl
#   garder le départ dont l'ICL est le plus grand
#   sortie : list(row_prob, col_prob, row_class, col_class,
#                 parameters = list(alpha, beta, xi, pi), ICL, K, L, iterations, best_init)
# À RETENIR :
#   - les 20 départs sont lancés à l'intérieur de la fonction (énoncé)
#   - vérifier les entrées dans utils.R avant de commencer
