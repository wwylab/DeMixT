## ----setup, include=FALSE------------------------------------------------
# knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
# devtools::install_github("wwylab/DeMixT")

## ---- results="hide", message=FALSE--------------------------------------
library(DeMixT)
data(test.data1.y)
data(test.data1.comp1)
res <- DeMixT(data.Y = test.data1.y,
    data.comp1 = test.data1.comp1,
    if.filter = FALSE,
    output.more.info = TRUE)

## ------------------------------------------------------------------------
res$pi
head(res$ExprT, 3)
head(res$ExprN1, 3)
head(res$Mu, 3)
head(res$Sigma, 3)
res$pi.iter
res$gene.name

## ---- results="hide", message=FALSE--------------------------------------
# data(test.data2.y)
# data(test.data2.comp1)
# data(test.data2.comp2)
# res <- DeMixT(data.Y = test.data2.y,
#     data.comp1 = test.data2.comp1,
#     data.comp2 = test.data2.comp2,
#     if.filter = FALSE)

## ---- results="hide", message=FALSE--------------------------------------
library(DeMixT)
data <- as.matrix(read.table("input.lcm.txt", header = FALSE))
normal <- data[, 1:25]
adm <- data[, 26:48]
tumor <- data[, 49:73]

nrows <- nrow(adm); ncols <- ncol(adm)
lcm.data.adm <- matrix(2^adm, nrows)
lcm.data.adm <- SummarizedExperiment(assays=list(counts=lcm.data.adm))

nrows <- nrow(tumor); ncols <- ncol(tumor)
lcm.data.tumor <- matrix(2^tumor, nrows)
lcm.data.tumor <- SummarizedExperiment(assays=list(counts=lcm.data.tumor))

nrows <- nrow(normal); ncols <- ncol(normal)
lcm.data.normal <- matrix(2^normal, nrows)
lcm.data.normal <- SummarizedExperiment(assays=list(counts=lcm.data.normal))

testr.TA <- DeMixT(data.Y = lcm.data.adm, data.comp1 = lcm.data.tumor,
    niter = 20, nbin = 60, if.filter = FALSE, tol = 10^-6)
testr.SA <- DeMixT(data.Y = lcm.data.adm, data.comp1 = lcm.data.normal,
    niter = 20, nbin = 60, if.filter = FALSE, tol = 10^-6)

## ------------------------------------------------------------------------
# plot A
dt_purT <- 1- as.numeric(testr.SA$pi)
dt_purS <- 1- as.numeric(testr.TA$pi)
plot(1 - dt_purS, dt_purT,
    col = "blue", pch = 1, xlim = c(0, 1), ylim = c(0, 1),
    xlab = expression(1 - hat(pi)[S]), ylab = expression(hat(pi)[T]))
abline(0, 1, col = "red", lwd = 2)

# Plot - Mean expressions for Tumor
OB_St <- log2(read.table("lcm_normal.txt", header = FALSE))
OB_Tu <- log2(read.table("lcm_tumor.txt", header = FALSE))
DT_Tu_mu <- as.numeric(testr.SA$Mu[, 1])
DT_St_mu <- as.numeric(testr.TA$Mu[, 1])
DT_Tu_sg <- as.numeric(testr.SA$Sigma[, 1])
DT_St_sg <- as.numeric(testr.TA$Sigma[, 1])
OB_St_m <- apply(OB_St, 1, mean)
OB_Tu_m <- apply(OB_Tu, 1, mean)

# filter out genes with large estimated standard deviations
condSt <- (DT_St_sg < 0.99)
condTu <- (DT_Tu_sg < 0.99)
DT_Tu_m <- as.numeric(apply(log2(testr.SA$ExprT), 1, mean))
DT_St_m <- as.numeric(apply(log2(testr.TA$ExprT), 1, mean))
OB_St_m <- OB_St_m[condSt]
OB_Tu_m <- OB_Tu_m[condTu]
DT_St_m <- DT_St_m[condSt]
DT_Tu_m <- DT_Tu_m[condTu]

# Plot - Mean expressions for Tumor
smoothScatter((DT_Tu_m + OB_Tu_m) / 2, DT_Tu_m - OB_Tu_m,
    ylab = "Estimate - Truth", xlab = "(Estimate + Truth)/2",
    xlim = c(2,16), ylim = c(-1.2,1.2),
    main = "Mean expressions for Tumor",
    pch = 1, nrpoints = 0, col = 'yellow',
    colramp=colorRampPalette(c("white","yellow",
    "yellow1","orange","orange1")))

tmp01 <- lowess((DT_Tu_m - OB_Tu_m) ~ ((DT_Tu_m + OB_Tu_m) / 2))
lines(tmp01$x, tmp01$y, col="blue", lwd = 5)
abline(h = 0, col = 'red', lty = 2)

# Plot - Mean expressions for Stroma
smoothScatter((DT_St_m + OB_St_m) / 2, DT_St_m - OB_St_m,
    ylab = "Estimate - Truth", xlab = "(Estimate + Truth)/2",
    xlim = c(2,16), ylim = c(-1.2,1.2),
    main = "Mean expressions for Stroma", pch = 1, nrpoints = 0,
    col = 'yellow',
    colramp=colorRampPalette(c("white","yellow",
    "yellow1","orange","orange1")))

tmp01 <- lowess((DT_St_m - OB_St_m) ~ ((DT_St_m + OB_St_m) / 2))
lines(tmp01$x, tmp01$y, col="blue", lwd = 5)
abline(h = 0, col = 'red', lty = 2)

## ------------------------------------------------------------------------
sessionInfo(package = "DeMixT")

