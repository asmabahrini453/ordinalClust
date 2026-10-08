# FICHIER : testthat.R
# RÔLE : lanceur des tests. Il ne teste rien lui-même : il dit à R d'exécuter tous les
#        fichiers de tests/testthat/ (c'est ce que fait R CMD check). On n'y touche plus.

library(testthat)
library(cubcoclust)

test_check("cubcoclust")