# FICHIER : initialize_lbm.R         RESPONSABLE (proposition) : Martine
# RÔLE : point de départ aléatoire d'un run (étape 1 du schéma).
# COURS : partie 1 p50 (plusieurs départs aléatoires) ; article ordinalClust p5 (init "random").
# FONCTION : initialize_lbm(X, K, L, m)
#   sortie : list(row_class, col_class, alpha, beta, xi, pi)
#   dimensions : row_class n ; col_class d ; alpha K ; beta L ; xi K x L ; pi K x L
# À RETENIR :
#   - tirer row_class et col_class au hasard, puis calculer les paramètres de départ (mstep_lbm)
#   - GROUPES VIDES : refuser ou retirer un départ où un groupe est vide
#     (l'article utilise "randomBurnin" et un minimum de cases par bloc)
