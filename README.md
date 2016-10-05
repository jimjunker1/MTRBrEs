**Barneche DR, White C, Marshall DJ *(2016)* Temperature effects on mass-scaling exponents in colonial animals: a manipulative test. Ecology, in press**

One can click [here](...) for the article URL  

### Overview  
This repository contains all data, analyses, figures and tables presented at the above-mentioned paper.  

* first I advise that you open your R GUI (R, RStudio, or starting on terminal) from the `MTRBrEs.RData` file...  
* ...this is an empty file that sets the absolute path to your project so everything will work independently on any machine;   
* pay attention to the required packages at the beginning of each `.R` file in the project root directory;  
* in the project root directory, you can reproduce all the analyses and outputs by running `analyses.R`, `numbers.R`, `tables.R` and then `figures.R`;  
* the file called `figures.R` reproduces all the figures exactly as they are shown in the paper **once** you have already reproduced all the outputs by running `analyses.R` above;  
* the file called `numbers.R` contains specific statistics exactly as presented in the paper (need to run `analyses.R` first too);  
* the file called `tables.R` reproduces the numbers presented in Tables 1 and 2, though not in a pretty format (need to run `analyses.R` first too);  
* the figures will be automatically placed in a directory called output (it is going to be automatically created for you);  
* **Importantly**, the alternative figures will only work if you have the packages `extrafont` and `fontcm` installed. Follow the instructions [here](https://cran.r-project.org/web/packages/fontcm/README.html) to install the font `CM Roman`;  
* notice that the Bayesian analysis in `analyses.R` may take a while to run on a regular computer;  

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
 [1] plyr_1.8.4         fontcm_1.1         extrafont_0.17     RColorBrewer_1.1-2 MASS_7.3-45        png_0.1-7         
 [7] LoLinR_0.0.0.9000  R2jags_0.5-7       rjags_4-6          coda_0.18-1        lme4_1.1-12        Matrix_1.2-6      

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.5      R2WinBUGS_2.1-21 lattice_0.20-33  zoo_1.7-13       lmtest_0.9-34    grid_3.3.1       Rttf2pt1_1.3.4  
 [8] nlme_3.1-128     minqa_1.2.4      extrafontdb_1.0  nloptr_1.0.4     boot_1.3-18      splines_3.3.1    abind_1.4-3     
[15] parallel_3.3.1  
```
### Please report if you run into problems or spot a bug on the code:
d13g0 DOT b4rn3ch3 AT m0n4sh DOT 3du (replace the 0 for o, 1 for i, 3 for e, 4 for a)  

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  
