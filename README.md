# Temperature effects on mass-scaling exponents in colonial animals: a manipulative test

This repository contains code and data needed to reproduce the article:

**Barneche DR, White CR, Marshall DJ**, Temperature effects on mass-scaling exponents in colonial animals: a manipulative test. *Ecology* (accepted 2016-10-04)

[![DOI](https://zenodo.org/badge/46690645.svg)](https://zenodo.org/badge/latestdoi/46690645)

## Instructions

All analyses were done in `R`. To compile the paper, including figures and tables we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies=TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies=TRUE)
```

You also need to install the software `JAGS`, as per these instructions [here](https://sourceforge.net/projects/mcmc-jags/?source=typ_redirect).

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

And then install the package `fontcm`, via `extrafont`. This installs the font `CM Roman` we use in our figures (for more information on see [these instructions](https://cran.r-project.org/web/packages/fontcm/README.html):

```r
extrafont::font_install('fontcm')
```

Then, to generate all figures, analyses, and manuscript (.docx, using Rmarkdown), simply do:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that the Bayesian analysis in this paper will take 20–30 min to run on a regular computer.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename="build.R")
```

### The paper can be reproduced using the following software and associated packages:
```
R version 3.3.1 (2016-06-21)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11.4 (El Capitan)

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] fontcm_1.1         extrafont_0.17     RColorBrewer_1.1-2 MASS_7.3-45        png_0.1-7          LoLinR_0.0.0.9000  rmarkdown_0.9.6    knitr_1.13         R2jags_0.5-7       rjags_4-6         
[11] coda_0.18-1        plyr_1.8.4         lme4_1.1-12        Matrix_1.2-6      

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.5      formatR_1.4      nloptr_1.0.4     tools_3.3.1      boot_1.3-18      digest_0.6.10    remake_0.2.0     evaluate_0.9     nlme_3.1-128     lattice_0.20-33  yaml_2.1.13     
[12] parallel_3.3.1   Rttf2pt1_1.3.4   storr_1.0.1      stringr_1.0.0    lmtest_0.9-34    grid_3.3.1       R6_2.2.0         minqa_1.2.4      extrafontdb_1.0  magrittr_1.5     codetools_0.2-14
[23] R2WinBUGS_2.1-21 htmltools_0.3.5  splines_3.3.1    abind_1.4-3      stringi_1.1.1    crayon_1.3.2     zoo_1.7-13      
```
### Please report if you run into problems or spot a bug on the code:
d13g0 DOT b4rn3ch3 AT m0n4sh DOT 3du (replace the 0 for o, 1 for i, 3 for e, 4 for a)  

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  

### Acknowledgements:  
Many thanks to [Rich FitzJohn](https://github.com/richfitz), [Remko Duursma](https://github.com/RemkoDuursma), and [Daniel Falster](http://danielfalster.com) for providing suggestions and excellent examples on how to implement this work flow using `remake`.
