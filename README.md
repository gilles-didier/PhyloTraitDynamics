# PhyloTraitDynamics

R code accompanying the manuscript:

[**Phylogenetic dynamics of MRCA ages and empirical moments of a Brownian trait**](https://arxiv.org/abs/2605.29736)  
Gilles Didier

This repository contains R functions and scripts used to compute and illustrate
the theoretical results presented in the manuscript. The focus is on the temporal
dynamics of empirical moments of Brownian traits evolving on phylogenetic trees,
and on MRCA ages in generalized birth-death trees.

## Overview

The package provides tools for computing:

- the distributional dynamics of the empirical moments of Brownian traits evolving on a given fixed tree;
- the distributional dynamics of the MRCA age of two uniformly sampled lineages in generalized birth-death trees;
- the dynamics of the variance of the empirical mean and of the expected empirical variance of a Brownian trait on generalized birth-death trees;
- simulation-based estimates of the variance of the empirical mean and of the expected empirical variance of a Brownian trait on generalized birth-death trees;
- various functions to plot the results;
- numerical reproduction of the figures in the manuscript.

The code is intended primarily as research code accompanying the paper. It is
not currently available on CRAN.

## Installation

The package can be installed directly from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "gilles-didier/PhyloTraitDynamics",
  subdir = "PhyloTraitDynamics"
)
```

Then load it with:

```r
library(PhyloTraitDynamics)
```

## Main contents

The repository contains:

- the R package implementing the functionalities above;
- its manual, `PhyloTraitDynamics-manual.pdf`;
- the R script used to generate the figures of the manuscript, `scriptFigures.R`.

## Citation

If you use this code, please cite the associated manuscript:

```bibtex
@misc{didier2026phylogeneticdynamicsmoments,
      title={Phylogenetic dynamics of MRCA ages and empirical moments of a Brownian trait}, 
      author={Gilles Didier},
      year={2026},
      eprint={2605.29736},
      archivePrefix={arXiv},
      primaryClass={q-bio.PE},
      url={https://arxiv.org/abs/2605.29736}, 
}
```

## License

Please see the `LICENSE` file for licensing information.

## Author

Gilles Didier  
Institut Montpelliérain Alexander Grothendieck  
Université de Montpellier  
Montpellier, France
