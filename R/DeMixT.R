require(parallel)
require(SummarizedExperiment)

DeMixT <- function(data.Y, data.comp1, data.comp2 = NULL, 
                   niter = 10, nbin = 50, 
                   if.filter = TRUE, ngene.selected.for.pi = 250, mean.diff.in.CM = 0.25,
                   tol = 10^(-5), output.more.info = FALSE, nthread = detectCores() - 1) {
  
  cat("Step 1: Estimation of Proportions\n")
  res.S1 <- DeMixT.S1(data.Y = data.Y, data.comp1 = data.comp1, data.comp2 = data.comp2, 
                      niter = niter, nbin = nbin, 
                      if.filter = if.filter, ngene.selected.for.pi = ngene.selected.for.pi, mean.diff.in.CM = mean.diff.in.CM,
                      tol = tol, nthread = nthread)
  
  cat("Step 2: Deconvolution of Expressions\n")
  res.S2 <- DeMixT.S2(data.Y = data.Y, data.comp1 = data.comp1, data.comp2 = data.comp2, 
                      givenpi = c(t(res.S1$pi)), nbin = nbin, nthread = nthread)
  
  cat("Deconvolution is finished\n")
  
  if (is.null(data.comp2)) { # two-component
    if (output.more.info) return(list(pi = res.S1$pi, ExprT = res.S2$decovExprT, ExprN1 = res.S2$decovExprN1, Mu = res.S2$decovMu, Sigma = res.S2$decovSigma, pi.iter = res.S1$pi.iter, gene.name = res.S1$gene.name))
    return(list(pi = res.S1$pi, ExprT = res.S2$decovExprT, ExprN1 = res.S2$decovExprN1, Mu = res.S2$decovMu, Sigma = res.S2$decovSigma))
  } else { # three-component
    if (output.more.info) return(list(pi = res.S1$pi, ExprT = res.S2$decovExprT, ExprN1 = res.S2$decovExprN1, ExprN2 = res.S2$decovExprN2, Mu = res.S2$decovMu, Sigma = res.S2$decovSigma, pi.iter = res.S1$pi.iter, gene.name = res.S1$gene.name))
    return(list(pi = res.S1$pi, ExprT = res.S2$decovExprT, ExprN1 = res.S2$decovExprN1, ExprN2 = res.S2$decovExprN2, Mu = res.S2$decovMu, Sigma = res.S2$decovSigma))
  }
}