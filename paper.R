library(knitr)
library(bibtex)
library(knitcitations)
library(plyr)

rm(list=ls())
###########
# MAIN TEXT 
###########
knit('text/MS.Rmd', output=file.path(getwd(), 'text/MS.md'), quiet=TRUE, encoding = 'utf-8')
# to word
system('pandoc -o text/MS.docx text/MS.md -s -S --bibliography library.bib --csl ecology.csl')

##############
# COVER LETTER 
##############
knit('text/coverLetter.Rmd', output=file.path(getwd(), 'text/coverLetter.md'), quiet=TRUE, encoding = 'utf-8')
# to word
system('pandoc -o text/coverLetter.docx text/coverLetter.md')
