# FICHIER : cub_em.R             
# RÔLE : estimer xi et pi à partir des notes d'un seul bloc (petit EM à deux composantes).
# COURS : partie 2 p58 (la CUB est un mélange estimé par EM) ; partie 1 p47-49 (étapes E puis M) ;
#         démonstration 3 de la feuille de démonstrations.
# FONCTION : cub_em(x, m, max_iter = 100, tol = 1e-6)
#   x : vecteur des notes du bloc      m : nombre de modalités
#   sortie : list(xi, pi, loglik, iterations, converged)
# À RETENIR :
#   - ici tol et converged ont un sens : c'est un EM déterministe (pas SEM)
#   - le mu du cours (partie 1) est la moyenne gaussienne : sans rapport avec xi
#   - cas limites : bloc vide ou une seule valeur, somme des tau = 0
#   - borner xi et pi à [0.001, 0.999] pour éviter les log(0)
# UTILISÉE PAR : mstep_lbm.R (appelée K x L fois par itération)
