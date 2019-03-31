#zzz.R

.onAttach <- function(...) {
    result <- .C("checkopenmp", numthread=as.integer(0))
    packageStartupMessage("OpenMP installed.\nNumber of threads = ", 
        result$numthread)
}