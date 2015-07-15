## Pipelines for GxE Analysis using Markers & Env. Covariates

This repository present a set of modules for analysis of data from multi-environment trials using genotype (SNP) and environmental covariates. The methods presented here are largley based on [Jarquin et al. (2014)](http://link.springer.com/article/10.1007%2Fs00122-013-2243-1#page-1). The project has ben funded by [Arvalis](http://www.arvalisinstitutduvegetal.fr/index.html) (PI: G. de los Campos).

**Contact**:  Gustavo de los Campos [ gdeloscampos@gmail.com ] &  Marco Lopez-Cruz [ malctony@hotmail.com ]

#####Before you start

The modules are based on R (>=3.1.0 [download](http://cran.r-project.org)). Before you start, be sure to install the packages [BGLR](http://cran.r-project.org/web/packages/BGLR/index.html) and [lme4](http://cran.r-project.org/web/packages/lme4/index.html).

```R
  install.packages("BGLR",repos="http://cran.r-project.org")
  install.packages("lme4",repos="http://cran.r-project.org")
```


#####Modules:

* [Module 1: Variance Components (lmer)](https://github.com/gdlc/ARVALIS/blob/master/varComp_lmer.md)
* [Module 2: Variance Components (BGLR)](https://github.com/gdlc/ARVALIS/blob/master/varComp_bglr.md)
* [Module 3: Structure of genotypes and of environments](https://github.com/gdlc/ARVALIS/blob/master/eigen.md)
* [Module 4: Genomic Regressions using BGLR](https://github.com/gdlc/ARVALIS/blob/master/full_data_models.md)
* [Module 5: Training-Testing evaluations](https://github.com/gdlc/ARVALIS/blob/master/training_testing.md)
* [Module 6: Splines]()
* [Moudle 7: Simulations]()
