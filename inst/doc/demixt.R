## ----setup, include=FALSE, message=FALSE--------------------------------------
library(ggplot2)
library(DeMixT)
plot.PCA =function(indata, batch, figure.dir, PCA.fig.title, legend.position = 'bottomleft',
                   label = FALSE, xlimit = NULL, ofile=TRUE, lines = TRUE) {
# data is a data matrix with samples in columns and genes in rows.
# batch is a vector with the order matching the order in indata.
    #batch=as.numeric(batch)
    N.batch = length(unique(batch))
    if(file.exists(paste0('Batch_Effect/', PCA.fig.title, '.RData'))){
      load(paste0('Batch_Effect/', PCA.fig.title, '.RData'))
    }else{
      pca <- SamplePCA(indata, usecor=F, center=T)
      save(pca, file = paste0('Batch_Effect/', PCA.fig.title, '.RData'))
    }
    pct1 <- round (pca@variances[1]/sum(pca@variances), digits=3)*100
    pct2 <- round (pca@variances[2]/sum(pca@variances), digits=3)*100
    xlab.text = paste("First Comp: ", as.character(pct1), "% variance", sep="")
    ylab.text = paste("Second Comp: ", as.character(pct2), "% variance", sep="")    
    
    #jpeg(file=file.path(figure.dir, paste("PCA_", PCA.fig.title, ".jpeg", sep="")), width = 600, height = 600, quality=100, pointsize=16)
    if(ofile) pdf(file=file.path(figure.dir, paste("PCA_", PCA.fig.title, ".pdf", sep="")))
    plot(pca@scores[,1], pca@scores[,2],  cex=1, xlab=xlab.text, ylab=ylab.text, col=rainbow(N.batch)[as.numeric(factor(batch))],
    pch=as.numeric(factor((batch))),lwd=1.5, main=PCA.fig.title)
    if(label == TRUE) {
      library(calibrate)
      textxy(pca@scores[,1], pca@scores[,2],colnames(indata), cex=0.75)
    }
    abline(h=0, v=0, col="brown", lty=2)
    abline(h=0, v=0, col="brown", lty=2)
    center1<-tapply(pca@scores[,1], factor(batch), mean)
    center2<-tapply(pca@scores[,2], factor(batch), mean)
    if(lines){
      for (ii in 1:length(center1)) {
        groupi<-pca@scores[as.numeric(factor(batch))==ii, 1:2]
        if (class(groupi)=="matrix") {
            for (j in (1:nrow(groupi))) {
                segments( groupi[j,1], groupi[j,2], center1[ii], center2[ii], col=rainbow(N.batch)[ii] , lwd=0.3)
            }
        }else {
            segments( groupi[1], groupi[2], center1[ii], center2[ii], col=rainbow(N.batch)[ii] , lwd=0.3)
        }
       }
    }
    points(center1, center2, pch=7, lwd=1.5,col=rainbow(N.batch))
    legend(legend.position,legend=names(table(factor(batch))), text.col=rainbow(N.batch), pch=1:8, col=rainbow(N.batch), lty=1)
   if(ofile) invisible(dev.off())
}

## ----install_DeMixT-----------------------------------------------------------
# devtools::install_github("wwylab/DeMixT")

## ----Algorithm, echo=FALSE, out.width='100%'----------------------------------
knitr::include_graphics(path = paste0("Algorithm.png"))

## ---- sim_2comp_GS, results="hide", message=FALSE-----------------------------
data("test.data.2comp")
# res.GS = DeMixT_GS(data.Y = test.data.2comp$data.Y, 
#                     data.N1 = test.data.2comp$data.N1,
#                     niter = 30, nbin = 50, nspikein = 50,
#                     if.filter = TRUE, ngene.Profile.selected = 150,
#                     mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#                     tol = 10^(-5))
load('Res_2comp/res.GS.RData')

## ----sim_2comp_GS_res---------------------------------------------------------
head(t(res.GS$pi))
head(res.GS$gene.name)

## ---- sim_2comp_S2, results="hide", message=FALSE-----------------------------
data("test.data.2comp")
# res.S2 <- DeMixT_S2(data.Y = test.data.2comp$data.Y, 
#                     data.N1 = test.data.2comp$data.N1,
#                     data.N2 = NULL, 
#                     givenpi = c(t(res.S1$pi[-nrow(res.GS$pi),])), nbin = 50)
load('Res_2comp/res.S2.RData')

## ----sim_2comp_S2_res---------------------------------------------------------
head(res.S2$decovExprT[,1:5],3)
head(res.S2$decovExprN1[,1:5],3)
head(res.S2$decovMu,3)
head(res.S2$decovSigma,3)

## ----read_data_2comp,  warning=FALSE, message=FALSE---------------------------
# ## DeMixT_DE without Spike-in Normal
# res.S1 = DeMixT_DE(data.Y = test.data.2comp$data.Y, 
#                    data.N1 = test.data.2comp$data.N1,
#                    niter = 30, nbin = 50, nspikein = 0,
#                    if.filter = TRUE, 
#                    mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#                    tol = 10^(-5))
# ## DeMixT_DE with Spike-in Normal
# res.S1.SP = DeMixT_DE(data.Y = test.data.2comp$data.Y, 
#                      data.N1 = test.data.2comp$data.N1,
#                      niter = 30, nbin = 50, nspikein = 50,
#                      if.filter = TRUE, 
#                      mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#                      tol = 10^(-5))
# ## DeMixT_GS with Spike-in Normal
# res.GS.SP = DeMixT_GS(data.Y = test.data.2comp$data.Y,
#                      data.N1 = test.data.2comp$data.N1,
#                      niter = 30, nbin = 50, nspikein = 50,
#                      if.filter = TRUE, ngene.Profile.selected = 150,
#                      mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#                      tol = 10^(-5))
load('Res_2comp/res.S1.RData'); load('Res_2comp/res.S1.SP.RData'); 
load('Res_2comp/res.GS.RData'); load('Res_2comp/res.GS.SP.RData'); 

## ----sim_2comp, fig.height = 4, fig.width = 6, fig.align='center', warning=FALSE----
res.2comp = as.data.frame(cbind(round(rep(t(test.data.2comp$pi[2,]),3),2), 
                                round(c(t(res.S1$pi[2,]),t(res.S1.SP$pi[2,]), t(res.GS.SP$pi[2,])),2),
                                rep(c('DE','DE-SP','GS-SP'), each = 100)), num = 1:2)
res.2comp$V1 <- as.numeric(as.character(res.2comp$V1))
res.2comp$V2 <- as.numeric(as.character(res.2comp$V2))
res.2comp$V3 = as.factor(res.2comp$V3)
names(res.2comp) = c('True.Proportion', 'Estimated.Proportion', 'Method')
## Plot
ggplot(res.2comp, aes(x=True.Proportion, y=Estimated.Proportion, group = Method, color=Method, shape=Method)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", lwd = 0.5) +
  xlim(0,1) + ylim(0,1)  +
  scale_shape_manual(values=c(seq(1:3))) +
  labs(x = 'True Proportion', y = 'Estimated Proportion')  

## ----read_data_3comp,  warning=FALSE, message=FALSE---------------------------
data("test.data.3comp")
# res.S1 <- DeMixT_DE(data.Y = test.data.3comp$data.Y, data.N1 = test.data.3comp$data.N1,
#                    data.N2 = test.data.3comp$data.N2, if.filter = TRUE)
load('Res_3comp/res.S1.RData'); 

## ----sim_3comp, fig.height = 4, fig.width = 6, fig.align='center', warning=FALSE----
res.3comp= as.data.frame(cbind(round(t(matrix(t(test.data.3comp$pi), nrow = 1)),2), 
                                round(t(matrix(t(res.S1$pi), nrow = 1)),2), 
                                rep(c('N1','N2','T'), each = 20)))
res.3comp$V1 <- as.numeric(as.character(res.3comp$V1))
res.3comp$V2 <- as.numeric(as.character(res.3comp$V2))
res.3comp$V3 = as.factor(res.3comp$V3)
names(res.3comp) = c('True.Proportion', 'Estimated.Proportion', 'Component')
## Plot
ggplot(res.3comp, aes(x=True.Proportion, y=Estimated.Proportion, group = Component, color=Component, shape=Component)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", lwd = 0.5) +
  xlim(0,1) + ylim(0,1)  +
  scale_shape_manual(values=c(seq(1:3))) +
  labs(x = 'True Proportion', y = 'Estimated Proportion')  

## ----read_data_PRAD,  warning=FALSE, message=FALSE----------------------------
load('res.PRAD.RData'); 

## ----PRAD, fig.height = 4, fig.width = 6, fig.align='center', warning=FALSE----
res.PRAD.df = as.data.frame(cbind(res.PRAD$res.GS.PRAD, res.PRAD$res.GS.SP.PRAD))
res.PRAD.df$V1 <- as.numeric(as.character(res.PRAD.df$V1))
res.PRAD.df$V2 <- as.numeric(as.character(res.PRAD.df$V2))
names(res.PRAD.df) = c('Estimated.Proportion', 'Estimated.Proportion.SP')
## Plot
ggplot(res.PRAD.df, aes(x=Estimated.Proportion, y=Estimated.Proportion.SP)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", lwd = 0.5) +
  xlim(0,1) + ylim(0,1)  +
  scale_shape_manual(values=c(seq(1:3))) +
  labs(x = 'Estimated Proportion', y = 'Estimated Proportion After Spike-in')  

## ----GTEx, echo=FALSE, out.width='65%', fig.align='center'--------------------
knitr::include_graphics(path = paste0("GTEx_normal.png"))

## -----------------------------------------------------------------------------
sessionInfo(package = "DeMixT")

