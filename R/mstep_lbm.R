# FICHIER : mstep_lbm.R              RESPONSABLE (proposition) : la binôme
# RÔLE : étape 4 : recalculer les paramètres à partir des classes tirées.
# COURS : démonstration 1 (alpha et beta, Lagrange) ; démonstration 3 (xi et pi par bloc) ;
#         partie 3 p36 (étape M) ; article ordinalClust p5 (2.1 et 2.2).
# FONCTION : mstep_lbm(X, row_class, col_class, K, L, m)
#   sortie : list(alpha, beta, xi, pi)   alpha K ; beta L ; xi K x L ; pi K x L
# À RETENIR :
#   - alpha = effectifs des groupes de lignes / n ; beta = effectifs des groupes de colonnes / d
#   - pour chaque bloc (k, l) : notes de X[row_class == k, col_class == l], puis cub_em()
#   - bloc vide ou groupe vide : prévoir un comportement (valeurs par défaut ou erreur claire)
