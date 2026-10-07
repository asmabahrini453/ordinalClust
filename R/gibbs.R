# FICHIER : gibbs.R                  RESPONSABLE (proposition) : Martine
# RÔLE : étapes 2 et 3 du schéma : tirer les classes des lignes, puis celles des colonnes.
# COURS : partie 3 p36 (SE-Gibbs) ; démonstration 2 ; article ordinalClust p5 (Algorithm 1, 1.1 et 1.2).
# FONCTIONS :
#   gibbs_row(X, col_class, alpha, xi, pi, K, m)  -> list(row_class, row_prob)   row_prob : n x K
#   gibbs_col(X, row_class, beta, xi, pi, L, m)   -> list(col_class, col_prob)   col_prob : d x L
#   gibbs_lbm(X, row_class, col_class, alpha, beta, xi, pi, K, L, m)
#        gibbs_row, puis gibbs_col avec le NOUVEAU row_class
#        -> list(row_class, col_class, row_prob, col_prob)
# À RETENIR :
#   - P(v_i = k) proportionnelle à alpha_k * produit sur les colonnes h de CUB(x_ih ; xi[k, w_h], pi[k, w_h])
#   - travailler en logarithmes et soustraire le maximum avant l'exponentielle (log-sum-exp)
#   - tirer un groupe par ligne : sample(1:K, 1, prob = ...)
