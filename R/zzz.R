#zzz.R

.onLoad <- function(libname, pkgname) {
  .C("checkopenmp")
  vig_list = tools::vignetteEngine(package = 'knitr')
  vweave <- vig_list[['knitr::knitr']][c('weave')][[1]]
  vtangle <- vig_list[['knitr::knitr']][c('tangle')][[1]]
  tools::vignetteEngine(pkgname, weave = vweave, tangle = vtangle,
                        pattern = "[.]Rmd$", package = pkgname)
}