# FICHIER : simulate_ordinal_lbm.R   RESPONSABLE : Martine
# RÔLE : fabriquer un faux tableau dont on connaît les vrais groupes et paramètres
#        (pour les tests, la page d'aide et la vignette).
# COURS : partie 3 p29-31 (hypothèses 1 et 2 du LBM) ; partie 2 p58 ;
#         manuel ordinalClust : jeu Msimulated ; article ordinalClust p4 (figure 2).
# FONCTION : simulate_ordinal_lbm(n, d, m, alpha, beta, xi, pi)
#   sortie : list(X, row_class, col_class, alpha, beta, xi, pi)   (vrais groupes et vrais paramètres)
# À RETENIR :
#   - tirer row_class avec alpha, col_class avec beta (hypothèse 1)
#   - tirer chaque note dans cub_probability(1:m, m, xi[k, l], pi[k, l]) (hypothèse 2)
#   - prendre m >= 4 (avec m = 3 la CUB s'identifie mal, à vérifier avec le prof)
#   - aucun groupe vide
