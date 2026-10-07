# FICHIER : cub_probability.R        RESPONSABLE (proposition) : Martine, à écrire en premier
# RÔLE : probabilité P(X = x) sous la loi CUB, pour un couple (xi, pi) et m modalités.
# COURS : partie 2 p58 (loi CUB = pi * binomiale décalée + (1 - pi) * uniforme).
# FONCTION : cub_probability(x, m, xi, pi)
#   x  : valeur(s) dans 1..m          m  : nombre de modalités
#   xi : paramètre feeling dans [0,1] pi : paramètre de mélange dans (0,1]
#   sortie : prob, une probabilité par valeur de x
# À RETENIR :
#   - la binomiale décalée est dbinom(x - 1, m - 1, 1 - xi) (voir la démonstration 3)
#   - pour la simulation et le calcul des blocs, on l'appelle avec x = 1:m (somme = 1)
#   - PRÉCAUTION : pi est un ARGUMENT, jamais la constante 3.14159 de R
#     -> vérifier pi dans (0, 1] et xi dans [0, 1], sinon erreur
# UTILISÉE PAR : cub_em.R, gibbs.R, compute_icl.R, simulate_ordinal_lbm.R
