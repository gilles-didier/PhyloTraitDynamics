# PhyloTraitDynamics

Petit package R compagnon d'article pour les calculs et figures autour de :

1. arbres phylogénétiques fixés ;
2. distributions d'âge de MRCA sous birth-death généralisé ;
3. espérance de la variance empirique sous birth-death + Brownien ;
4. variance de la moyenne empirique sous birth-death + Brownien.

Le package n'est pas conçu pour le CRAN. Il sert à reproduire et organiser des scripts de calcul déjà validés.

## Principe de construction

Les scripts initiaux sont conservés dans :

```r
inst/scripts/original/
```

Les fonctions publiques sont des wrappers conservateurs. Elles ne modifient pas les formules, les intégrandes ou les boucles internes des scripts source. Les noms publics cherchent seulement à rendre l'API plus cohérente.

## Installation locale

Depuis le répertoire parent :

```r
install.packages("PhyloTraitDynamics_0.1.0.tar.gz", repos = NULL, type = "source")
```

ou, depuis le dossier source :

```r
remotes::install_local("PhyloTraitDynamics")
```

## Fonctions publiques principales

### Arbres fixes

```r
fixed_tree_simulate_brownian_realization()
fixed_tree_plot_brownian_realization()

fixed_tree_compute_theoretical_summary()
fixed_tree_plot_theoretical_summary()
```

### MRCA

```r
mrca_age_compute_distribution()
mrca_age_compute_dynamics()
mrca_age_plot_dynamics()
```

### Espérance de la variance empirique

```r
empirical_variance_compute_expectation()
empirical_variance_simulate()
empirical_variance_summarise_simulation()
empirical_variance_plot_expectation()
```

Convention préservée : la variance empirique vaut 0 quand il y a moins de deux lignées vivantes.

Les conditionnements disponibles sont :

```r
conditioning = "none"
conditioning = "survival"   # conditionnement à N(t) > 0
```

Il n'y a pas de conditionnement à `N(t) >= 2` ni à la survie à un temps final dans les scripts sources.

### Variance de la moyenne empirique

```r
empirical_mean_compute_variance()
empirical_mean_simulate()
empirical_mean_summarise_simulation()
empirical_mean_plot_variance()
empirical_mean_plot_paths()
```

La théorie actuellement exposée correspond à :

```r
Var(empirical mean at time t | N(t) > 0)
```

## Taux birth-death

Les arguments publics sont toujours :

```r
birth
death
```

Ils peuvent être :

- des constantes numériques positives ou nulles ;
- des fonctions R du temps.

Si une fonction de taux n'est pas vectorisée, le wrapper tente une évaluation point par point. Si cela échoue, une erreur explicite est renvoyée.

## Arguments temporels

Les fonctions publiques utilisent en général :

```r
time_start = 0
time_end
time_step
```

Une grille interne est construite avec `seq(time_start, time_end, by = time_step)`, en ajoutant `time_end` si nécessaire.

Pour les simulations birth-death, le processus reste simulé depuis le temps 0 jusqu'à `time_end`, puis les temps demandés sont extraits.

## Exemple minimal : espérance de la variance empirique

```r
library(PhyloTraitDynamics)

th <- empirical_variance_compute_expectation(
  birth = 0.8,
  death = 0.2,
  sigma2 = 1,
  time_step = 0.02,
  time_end = 3,
  conditioning = "none",
  method = "auto"
)

plot(th$time, th$empirical_variance_expectation, type = "l")
```

## Exemple minimal : variance de la moyenne empirique

```r
th_mean <- empirical_mean_compute_variance(
  birth = function(t) 0.6 + 0.1 * t,
  death = 0.2,
  sigma2 = 1,
  time_end = 3,
  time_step = 0.02
)

plot(th_mean$time, th_mean$empirical_mean_variance, type = "l")
```

## Note conservatrice

Le package expose les calculs des scripts initiaux. Il ne tente pas encore de factoriser les duplications entre scripts, ni de transformer le code en bibliothèque générale.

En particulier, `fixed_tree_compute_theoretical_summary()` reprend la sortie de `compute_emp_mean_var_timeseries_general()` pour la variance de la moyenne empirique, sans modifier le calcul interne.
