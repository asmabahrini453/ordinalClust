# FICHIER : cub_em.R             
# RÔLE : estimer xi et pi à partir des notes d'un seul bloc.On a les notes d'un bloc
#et on cherche les (xi, pi) qui les expliquent le mieux (maximum de vraisemblance). 
#On ne peut pas le faire directement car, pour chaque  note, on ne sait pas si elle vient 
# du "feeling" ou du "hasard" (on ne peut pas maximiser la vraissemblance) : c'est une
# VARIABLE LATENTE, comme le groupe z dans un mélange. D'où l'algorithme EM.

# COURS : partie 2 p58 (la CUB est un mélange estimé par EM) ; partie 1 p47-49 (étapes E puis M) ;
#         démonstration 3 de la feuille de démonstrations.

# FONCTION : cub_em(x, m, max_iter = 100, tol = 1e-6)
#  Entrées:  x : vecteur des observations du bloc      m : nombre de modalités ,
#           max_iter : nombre maximal d'itérations de l'EM CUB,  tol : seuil utilisé pour décider la convergence
#  Sorties : list(xi, pi, loglik, iterations, converged)

# principe : 
#   groupe k (composante)         ->  "feeling" ou "hasard"   (2 composantes)
#   proportion p_k                ->  pi (feeling) et 1 - pi (hasard)
#   densité f_k(x)                ->  PB(x) (feeling) et 1/m (hasard)
#   responsabilité t_k(x_i)       ->  tau(x) = P(feeling | note x): la proba. d'appartenance de xi au grp K 
# 
# Une itération :
#   Étape E :  t_k(x_i) = pi * PB(x) / ( pi * PB(x) + (1 - pi)/m )      (Bayes)
#   Étape M : pi = moyenne des tau                                    (comme p_k = n_k / n, p49)
#             xi = (m - moyenne des x pondérée par tau) / (m - 1)
#             (parmi les notes "feeling", X - 1 suit une binomiale(m-1, 1-xi) ;
#              le maximum de vraisemblance d'une binomiale est la moyenne divisée par
#              m-1, donc 1 - xi = (moyenne de x - 1)/(m-1))
#   Arrêt   : |loglik(q+1) - loglik(q)| < tol                     (le critère de convergence).
## La preuve de convergence du cours garantit que loglik ne baisse jamais.

# À RETENIR :
#   - ici tol et converged ont un sens : c'est un EM déterministe (pas SEM)
#   - cas limites : bloc vide ou une seule valeur, somme des tau = 0
#   - borner xi et pi à [0.001, 0.999] pour éviter les log(0)

# UTILISÉE PAR : mstep_lbm.R (appelée K x L fois par itération)

##########################################################################################""

cub_em <- function(x, m, max_iter = 100, tol = 1e-6){
  
  # le Cas limite : bloc vide. On renvoie des valeurs neutres, sans converger.
  #    (mstep_lbm peut ainsi appeler cub_em sans planter si un bloc se vide)
  if (length(x) == 0){
    return(list(xi = 0.5, pi = 0.5, loglik = 0, iterations = 0, converged = FALSE))
  }
  if (any(x < 1 | x > m | x != round(x))) stop("x doit contenir des entiers entre 1 et m")
  
  #  les notes ne prennent que m valeurs, donc on compte
  #    combien de fois chaque note apparaît (beaucoup plus rapide que de boucler sur x)
  notes <- 1:m
  effectif <- tabulate(x, nbins = m)   # effectif[j] = nombre de notes égales à j
  n <- sum(effectif)                   # nombre total de notes
  
  # Initialisation (theta^(0) dans le cours) : pi = 0.5 ; xi par la méthode des
  #    moments, car la moyenne d'un CUB vaut à peu près m - xi * (m - 1)
  borner <- function(v) min(max(v, 0.001), 0.999)   # garde les paramètres dans [0.001, 0.999]
  pi <- 0.5
  xi <- borner((m - mean(x)) / (m - 1))
  
  # ALGO EM
  loglik_old <- -Inf
  converged <- FALSE
  for (iter in 1:max_iter){
    
    # Étape E : probabilité que chaque note vienne du feeling
    pb <- dbinom(notes - 1, m - 1, 1 - xi)             # PB(x) pour x = 1..m
    tau <- pi * pb / (pi * pb + (1 - pi) / m)
    
    # Étape M : on remet à jour pi puis xi
    poids_feeling <- sum(effectif * tau)               # "nombre" de notes issues du feeling
    pi <- borner(poids_feeling / n)
    if (poids_feeling > 1e-12){                        # cas limite : somme des tau = 0
      xi <- borner((m - sum(effectif * tau * notes) / poids_feeling) / (m - 1))
    }
    
    # Log-vraisemblance avec les nouveaux paramètres : somme des ln P(x_i)
    loglik <- sum(effectif * log(cub_probability(notes, m, xi, pi)))
    
    # Critère d'arrêt : la log-vraisemblance ne bouge presque plus
    if (abs(loglik - loglik_old) < tol){
      converged <- TRUE
      break
    }
    loglik_old <- loglik
  }
  
  return(list(xi = xi, pi = pi, loglik = loglik, iterations = iter, converged = converged))
}
