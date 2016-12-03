## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(concordance=TRUE)

## ----install1, eval=FALSE, results='asis'--------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("interactiveDisplay")

## ----install2, eval=FALSE, results='asis'--------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  useDevel(TRUE)
#  biocLite("interactiveDisplay")

## ----citation, eval=TRUE, results='tex'----------------------------------
citation("interactiveDisplay")

## ----libraries, eval=TRUE, results='asis'-------------------------------------
options(width=80)
options(continue=" ")
suppressMessages(library(ggplot2))
suppressMessages(library(interactiveDisplay))
suppressMessages(library(Biobase))

## ----data_mmgr, eval=TRUE, results='tex'--------------------------------------
data(mmgr)
mmgr

## ----grl, eval=TRUE, results='tex'--------------------------------------------
grl <- GenomicRangesList(list(mmgr,mmgr))

## ----data_mmgrl, eval=TRUE, results='tex'-------------------------------------
data(mmgrl)
mmgrl

## ----display_mmgr, eval=FALSE, results='asis'---------------------------------
#  display(mmgr)

## ----display_mmgrl, eval=FALSE, results='asis'--------------------------------
#  display(mmgrl)

## ----display_mmgr2, eval=FALSE, results='asis'--------------------------------
#  new_mmgr <- display(mmgr)

## ----display_mmgrl2, eval=FALSE, results='asis'-------------------------------
#  new_mmgrl <- display(mmgrl)

## ----data_expr, eval=TRUE, results='tex'--------------------------------------
data(expr)
expr

## ----exprs_expr, eval=TRUE, results='tex'-------------------------------------
exprs(expr)[1:10,1:7]

## ----display_expr, eval=FALSE, results='asis'---------------------------------
#  display(expr)

## ----data_se, eval=TRUE, results='tex'----------------------------------------
data(se)
se

## ----display_se, eval=FALSE, results='asis'-----------------------------------
#  display(se)

## ----display_mtcars, eval=FALSE, results='asis'-------------------------------
#  display(mtcars)

## ----display_mtcars2, eval=FALSE, results='asis'------------------------------
#  new_mtcars <- display(mtcars)

## ----simplenet_mtcars, eval=FALSE, results='asis'-----------------------------
#  simplenet(mtcars)

## ----plot_mtcars, eval=FALSE, results='asis'----------------------------------
#  data(mtcars)
#  qp <- qplot(mpg, data=mtcars, geom="density", fill=factor(cyl), alpha=I(.4))
#  plot(qp)

## ----gridsvgjs_qp, eval=FALSE, results='asis'---------------------------------
#  gridsvgjs(qp)

## ----sessionInfo, eval=TRUE, results='tex'------------------------------------
sessionInfo()

