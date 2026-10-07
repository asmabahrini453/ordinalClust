# ordinalClust
# Package R pour le Co-clustering de données ordinales

## Présentation

L'objectif est de développer un **package R permettant de réaliser un co-clustering de données ordinales**, en s'appuyant sur un modèle de blocs latents. Le modèle utilisé repose sur le **modèle CUB pour données ordinales**, avec une inférence réalisée à l'aide d'un **algorithme SEM** et une sélection du nombre de blocs basée sur le **critère ICL**.

Le package permettra notamment de :
- regrouper simultanément les **lignes et les colonnes** d'une matrice de données ordinales ;
- estimer les probabilités d'appartenance aux différents clusters ;
- déterminer l'affectation des lignes et des colonnes aux clusters ;
- calculer le **critère ICL** ;
- retourner les paramètres estimés du modèle.

##  Objectif

L'objectif final est de fournir un **package R installable et documenté**, accompagné d'une aide pour la fonction principale ainsi que d'une vignette présentant son utilisation.

---
**Projet 2026–2027 — M2 MALIA**  
Université Lumière Lyon 2

## Notations du projet

Les mêmes noms sont utilisés partout dans le code et la documentation.

| Nom R | Math | Dimension | Signification |
| --- | --- | --- | --- |
| X | X | n x d | matrice de données ordinales |
| n | n | 1 | nombre de lignes |
| d | d | 1 | nombre de colonnes |
| m | m | 1 | nombre de modalités ordinales |
| K | K | 1 | nombre de clusters de lignes |
| L | L | 1 | nombre de clusters de colonnes |
| x | x | variable | une ou plusieurs observations ordinales |
| row_class | v | n | cluster de chaque ligne |
| col_class | w | d | cluster de chaque colonne |
| row_prob | | n x K | probabilités d'appartenance des lignes |
| col_prob | | d x L | probabilités d'appartenance des colonnes |
| alpha | α | K | proportions des clusters de lignes |
| beta | β | L | proportions des clusters de colonnes |
| xi | ξ | K x L | paramètre feeling de la CUB, par bloc |
| pi | π | K x L | paramètre de mélange de la CUB, par bloc |
| theta | θ | liste | ensemble des paramètres du modèle |
| loglik | ℓ | 1 | log-vraisemblance |
| max_iter | | 1 | nombre d'itérations (EM de la CUB, SEM) |
| burn_in | | 1 | itérations de burn-in du SEM |
| tol | ε | 1 | seuil de convergence (EM de la CUB) |
| n_init | | 1 | nombre de départs aléatoires (20 par défaut) |
| ICL | ICL | 1 | critère de sélection (le plus grand est le meilleur) |

Attention : `pi` désigne le paramètre CUB π et masque la constante `pi` de R. Il doit toujours être un argument des fonctions qui l'utilisent.

K et L sont les nombres de clusters ; dans les boucles, k = 1..K et l = 1..L désignent un cluster particulier.
