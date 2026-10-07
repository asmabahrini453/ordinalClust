# FICHIER : estimate_partition.R     RESPONSABLE (proposition) : la binôme
# RÔLE : étape 6 : partition finale et probabilités d'appartenance, avec theta chapeau FIXÉ.
# COURS : partie 3 p36 (partition obtenue par Gibbs avec theta chapeau) ; article ordinalClust p5.
# FONCTION : estimate_partition(X, row_class, col_class, alpha, beta, xi, pi, K, L, m, n_gibbs = 50)
#   (signature complétée par rapport au dictionnaire : départ + nombre de tirages ; à valider à deux)
#   sortie : list(row_class, col_class, row_prob, col_prob)
# À RETENIR :
#   - enchaîner n_gibbs fois gibbs_lbm() avec les paramètres fixés
#   - row_prob et col_prob = fréquence d'affectation à chaque groupe sur ces tirages
#   - row_class et col_class = groupe le plus fréquent
