me1 <- function (Ne, L) {(2 * Ne * L) / log(4 * Ne * L)}
me2 <- function (Ne, L, n) {(2 * Ne * L) / log(2 * Ne * (L / n))}
me3 <- function (Ne, L, n) {(2 * Ne * L) / log(Ne * (L / n))}
me4 <- function (Ne, L) {(2 * Ne * L)}
me5 <- function (Ne, L) {(4 * Ne * L)}


Ne <- 100
L <- 30
me1(Ne, L)
me2(Ne, L, 30)
me3(Ne, L, 30)
me4(Ne, L)
me5(Ne, L)

