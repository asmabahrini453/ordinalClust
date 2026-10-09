# FICHIER : essai_dataqol.R        (script de travail, pas une fonction du package)
# RÔLE : tester ordinal_coclust sur un VRAI jeu de données ordinales (dataqol).
#
# ----------------------------------------------------------------------
# LE JEU DE DONNÉES
# ----------------------------------------------------------------------
# dataqol vient du package R "ordinalClust" (package de référence pour le
# co-clustering de données ordinales). Il contient les réponses de patients à
# un questionnaire de qualité de vie ("qol" = quality of life).
#   - lignes   : les patients (les mesures)
#   - colonnes : les questions q1 à q28
#   - valeurs  : réponse sur une échelle ordinale de 1 à m (m = 4)
#
# Ce qu'on retire du tableau d'origine (31 colonnes, 117 lignes) :
#   - Id        : c'est le numéro du patient, pas une note
#   - q29, q30  : ils ont une autre échelle (7 modalités) ; notre modèle
#                 suppose le même m pour toutes les questions
#   - les lignes avec des NA (le package ne les gère pas encore)
# Remarque : chaque patient apparaît 3 fois (3 mesures). Ces lignes ne sont pas
# strictement indépendantes ; on peut ne garder qu'une ligne par patient
# (option une_ligne_par_patient plus bas).
#
# ----------------------------------------------------------------------
# CE QU'ON FAIT ICI
# ----------------------------------------------------------------------
# Ici on ne connaît PAS la vérité : on regarde si les groupes trouvés ont un sens.
#   1/ Charger dataqol, garder q1 à q28, vérifier les valeurs, retirer les NA.
#   2/ Lancer ordinal_coclust pour plusieurs (K, L) et garder le plus grand ICL.
#   3/ Lire le résultat (effectifs, groupes de questions, xi et pi par bloc).
#   4/ graphique: la matrice réorganisée par bloc
#
# Lecture d'un bloc (k, l) : xi petit -> notes hautes, xi grand -> notes basses ;
# pi proche de 1 -> réponses nettes, pi proche de 0 -> hasard.
# Attention : vérifier avec ?dataqol si 1 veut dire "pas du tout" ou "beaucoup".
#
# À RETENIR :
#   - chronométrer un premier essai avant de comparer tous les (K, L)
#   - numéros de groupes arbitraires (label switching)

# ----------------------------------------------------------------------

devtools::load_all()          # NOTRE package : ordinal_coclust vient de lui

# ---- Réglages -----------------------------------------------------------
rapide <- FALSE                  # TRUE : essai rapide ; FALSE : comparaison complète
une_ligne_par_patient <- TRUE  # TRUE : garde seulement la 1re mesure de chaque patient

# ---- 1/ Charger les données ---------------------------------------------
# Le package ordinalClust ne sert qu'à récupérer les données (pas ses fonctions).
data(dataqol, package = "ordinalClust")
if (une_ligne_par_patient) dataqol <- dataqol[!duplicated(dataqol$Id), ]
X <- as.matrix(dataqol[, paste0("q", 1:28)])       # on garde q1 à q28 seulement
cat("Dimensions :", dim(X), "\n")

# ---- 2/ Vérifier les valeurs --------------------------------------------
print(table(X, useNA = "ifany"))                   # on attend 1, 2, 3, 4 et des NA

# ---- 3/ Enlever les lignes avec des NA ----------------------------------
n_avant <- nrow(X)
X <- X[complete.cases(X), , drop = FALSE]
cat("Lignes gardées :", nrow(X), "sur", n_avant, "\n")

# ---- 4/ Vérifier que les notes vont de 1 à m ----------------------------
storage.mode(X) <- "integer"
m <- max(X)
cat("m =", m, "modalités\n")
stopifnot(min(X) == 1, all(X >= 1 & X <= m))

# ---- 5/ Un premier essai chronométré ------------------------------------
set.seed(1)
temps <- system.time(res <- ordinal_coclust(X, K = 2, L = 2, n_init = 2))
print(temps)
cat("ICL (2x2) :", round(res$ICL, 1), "\n")

# ---- 6/ Choisir (K, L) avec l'ICL ---------------------------------------
if (rapide) {
  couples <- list(c(1,1), c(2,2), c(3,2), c(2,3))
  n_init_icl <- 2; n_init_final <- 3
} else {
  couples <- list(c(1,1), c(2,1), c(1,2), c(2,2), c(3,2), c(2,3), c(3,3))
  n_init_icl <- 5; n_init_final <- 10
}
set.seed(1)
icl <- sapply(couples, function(kl){
  suppressWarnings(ordinal_coclust(X, K = kl[1], L = kl[2], n_init = n_init_icl)$ICL)
})
names(icl) <- sapply(couples, function(kl) paste0(kl[1], "x", kl[2]))
print(round(icl, 1))
best <- couples[[which.max(icl)]]
cat("Meilleur (K x L) :", names(icl)[which.max(icl)], "\n")

# ---- 7/ Résultat final avec le meilleur couple --------------------------
set.seed(1)
fin <- ordinal_coclust(X, K = best[1], L = best[2], n_init = n_init_final)
cat("Effectifs des groupes de lignes   :", table(fin$row_class), "\n")
cat("Effectifs des groupes de colonnes :", table(fin$col_class), "\n")
cat("xi estimé :\n");  print(round(fin$parameters$xi, 2))
cat("pi estimé :\n");  print(round(fin$parameters$pi, 2))

# ---- 8/ Quelles questions dans quel groupe ? ----------------------------
print(split(colnames(X), fin$col_class))


# ---- 9/ Graphique : la matrice avant / après réorganisation -------------
palette_notes <- hcl.colors(m, "YlOrRd", rev = TRUE)

dessiner_matrice <- function(M, titre, bornes_l = NULL, bornes_c = NULL){
  nl <- nrow(M); nc <- ncol(M)
  image(x = 1:nc, y = 1:nl, z = t(M[nl:1, , drop = FALSE]),
        col = palette_notes, zlim = c(1, m),
        xlab = "Questions", ylab = "Patients", main = titre, axes = FALSE)
  box()
  if (!is.null(bornes_c)) abline(v = bornes_c + 0.5, lwd = 2)
  if (!is.null(bornes_l)) abline(h = nl - bornes_l + 0.5, lwd = 2)
}

graphique_matrice <- function(){
  X_ord <- X[order(fin$row_class), order(fin$col_class)]
  bornes_l <- cumsum(table(fin$row_class)); bornes_l <- bornes_l[-length(bornes_l)]
  bornes_c <- cumsum(table(fin$col_class)); bornes_c <- bornes_c[-length(bornes_c)]
  par(mfrow = c(1, 2))
  dessiner_matrice(X, "Tableau d'origine")
  dessiner_matrice(X_ord, paste0("Tableau réorganisé (", best[1], "x", best[2], ")"),
                   bornes_l, bornes_c)
  par(mfrow = c(1, 1))
}

graphique_matrice()

# ---- 10/ Sauvegarder le graphique dans dev/ -----------------------------
dir.create("dev", showWarnings = FALSE)
png("dev/matrice_dataqol.png", width = 1100, height = 550); graphique_matrice(); dev.off()