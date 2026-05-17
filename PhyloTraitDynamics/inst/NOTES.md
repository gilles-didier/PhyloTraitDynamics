# Notes de construction

- Les scripts originaux sont conservés tels quels dans `inst/scripts/original/`.
- Les fonctions publiques sont des wrappers.
- Les fonctions auxiliaires originales ne sont pas exportées.
- Les noms publics emploient `birth` et `death` plutôt que `lambda_fun` et `mu_fun`.
- La convention `S^2(t) = 0` quand moins de deux valeurs sont présentes est préservée.
- Les conditionnements exposés correspondent uniquement à ceux présents dans les scripts initiaux.
- Aucun conditionnement à `N(t) >= 2` ni à la survie à un temps final n'a été ajouté.
