# FICHIER : compute_icl.R            RESPONSABLE (proposition) : la binôme
# RÔLE : étape 7 : critère ICL d'un run.
# COURS : partie 3 p37 (formule à trois pénalités) ; partie 2 p28-31 ; démonstration 4.
# FONCTION : compute_icl(X, row_class, col_class, alpha, beta, xi, pi, K, L, m)
#   sortie : ICL, un seul nombre (le plus grand est le meilleur)
# À RETENIR :
#   - vraisemblance complète = somme des log alpha[k] + somme des log beta[l] + log CUB de chaque note
#   - pénalité = (K-1)/2 * log(n) + (L-1)/2 * log(d) + K*L*2/2 * log(n*d), car nu = 2 par bloc
#   - logarithme naturel (log)
