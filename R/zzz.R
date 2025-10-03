#zzz.R

#' @useDynLib DeMixT, .registration = TRUE
.onAttach <- function(...) {
    result <- .C("checkopenmp", numthread=as.integer(0))
    packageStartupMessage("OpenMP installed.\nNumber of threads = ", 
        result$numthread)
}
