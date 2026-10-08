# FICHIER : essai_bout_en_bout.R
# RÔLE : essai à la main de tout l'algorithme (script de travail, pas une fonction du package).
#        On fabrique un faux tableau dont on connaît la vérité, on lance ordinal_coclust(),
#        puis on compare avec la vérité et on regarde si l'ICL choisit le bon (K, L).
# COURS : partie 3 p37 (ICL) ; partie 3 p29-31 (LBM).
# À LANCER : depuis la racine du projet, après devtools::load_all()
# DURÉE : environ une minute

devtools::load_all()

# Indice de Rand ajusté (ARI) : compare deux partitions, au renommage des groupes près.
# 1 = partitions identiques, proche de 0 = aucun lien.
ari <- function(a, b){
  tab <- table(a, b)                            # tableau croisé des deux partitions
  n <- sum(tab)
  paires <- function(x) sum(choose(x, 2))       # nombre de paires d'éléments rangés ensemble
  somme_lignes <- paires(rowSums(tab))
  somme_colonnes <- paires(colSums(tab))
  somme_cases <- paires(tab)
  attendu <- somme_lignes * somme_colonnes / choose(n, 2)   # valeur attendue si le hasard seul
  return((somme_cases - attendu) / ((somme_lignes + somme_colonnes) / 2 - attendu))
}

# 1/Construire un faux tableau : 60 lignes, 20 colonnes, m = 5, 2 x 2 blocs aux paramètres très différents
xi_vrai <- matrix(c(0.15, 0.85, 0.85, 0.15), nrow = 2)
pi_vrai <- matrix(0.9, nrow = 2, ncol = 2)
set.seed(1)
sim <- simulate_ordinal_lbm(n = 60, d = 20, m = 5, alpha = c(0.5, 0.5), beta = c(0.5, 0.5),
                            xi = xi_vrai, pi = pi_vrai)

# 2/Lancer l'algorithme avec les vrais (K, L)
res <- ordinal_coclust(sim$X, K = 2, L = 2, n_init = 5)

# 3/ Comparer avec la vérité
cat("ARI des lignes    :", round(ari(sim$row_class, res$row_class), 3), "\n")   # 1 = groupes retrouvés
cat("ARI des colonnes  :", round(ari(sim$col_class, res$col_class), 3), "\n")
cat("ICL               :", round(res$ICL, 1), "\n")
print(table(vrai = sim$row_class, estime = res$row_class))   # une diagonale (ou anti-diagonale) pleine
cat("xi estimé (les numéros de groupes peuvent être échangés) :\n"); print(round(res$parameters$xi, 2))
cat("xi vrai :\n"); print(xi_vrai)

# 4/ L'ICL choisit-il le bon (K, L) ? On essaie plusieurs couples, on garde le plus grand ICL.
couples <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2), c(3, 3))
icl <- sapply(couples, function(kl){
  suppressWarnings(ordinal_coclust(sim$X, K = kl[1], L = kl[2], n_init = 5)$ICL)
})
names(icl) <- sapply(couples, function(kl) paste0(kl[1], "x", kl[2]))
print(round(icl, 1))
cat("Meilleur (K x L) selon l'ICL :", names(icl)[which.max(icl)], "(le vrai modèle est 2x2)\n")