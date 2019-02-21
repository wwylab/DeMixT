#zzz.R

.onLoad <- function(...) {
  .C("checkopenmp")
}