# ARZIMM

Statistical modeling and inference of microbial interaction and stability for microbial dynamical systems

This package is developed to model microbial dynamical systems from longitudinal microbiome data and infer microbial interaction and stability. ARZIMM models the excess zero abundance and the non-zero abundances separately; and use a random effect model to borrow strength across subjects.

NeedsCompilation: No

Depends: R(>= 3.6.2)

Imports: phyloseq, glmnet, stringr, expm, lme4, ggplot2

License: GPL-2

The manual file is "ARZIMM-manual.pdf".

Installation of ARZIMM in R:
```r
library("devtools")
install_github("Hlch1992/ARZIMM",force=T)
```

A remark should be make that there requires some patience to install the dependent packages. 
