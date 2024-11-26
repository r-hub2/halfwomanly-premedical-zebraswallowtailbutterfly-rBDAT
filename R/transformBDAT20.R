#' @title transform BDAT20-matrix
#' @description transforms an intermediate BDAT20-matrix into a data.frame,
#' names it appropriately and removes <zero>-entries.
#' @param assort a matrix, produced by calling Fortran-BDAT20-subroutine and
#' extracting one list element to a matrix
#' @param value character vector indicating return type: either "Vol", "Skl"
#' or "Fix"
#' @details Fortran subroutine BDAT20 returns a list with many entries for each
#' tree; most relevant are 'Vol', 'FixLng' and maybe 'Skl'. These elements are
#' transformed into a matrix and properly whiped into shape.
#' @keywords internal
#' @return a transformed BDAT20-result, usually a data.frame

transformBDAT20 <- function(assort, value) {
  if (identical(value, "Vol")) {
    assort <- as.data.frame(assort)
    cnames <- c("Vfm", "X", "Sth", "Ab", "Ind", "nvDh", "EV")
    colnames(assort) <- cnames
    assort$tree <- 1:nrow(assort)
    assort <- assort[, c("tree", cnames)]
  } else if (identical(value, "Skl")) {
    assort <- as.data.frame(assort)
    cnames <- c("X", "Xab", "Sth", "Sthab", "Ab", "Abab")
    colnames(assort) <- cnames
    assort$tree <- 1:nrow(assort)
    assort <- assort[, c("tree", cnames)]
  } else if (identical(value, "LDSort")) {
    assort <- as.data.frame(matrix(as.numeric(t(assort)), ncol = 4, byrow = T))
    colnames(assort) <- c("height", "length", "midD", "topD")
    assort$tree <- rep(1:(nrow(assort) / 5), each = 5)
    assort$No <- rep(1:5, nrow(assort) / 5)
    assort$Sort <- rep(c("X", "Sth", "Ab", "Ind", "nvDh"), nrow(assort) / 5)
    assort <- assort[, c(
      "tree", "No", "Sort", "height", "length", "midD",
      "topD"
    )]
  } else if (identical(value, "FixLng")) {
    assort <- as.data.frame(matrix(as.numeric(t(assort)), ncol = 6, byrow = T))
    cnames <- c("No", "height", "length", "midD", "topD", "Vol")
    colnames(assort) <- cnames
    assort$tree <- rep(1:(nrow(assort) / 30), each = 30)
    x <- formatC(1:30, width = 2, format = "d", flag = "0")
    assort$Sort <- rep(paste0("Fix", x), nrow(assort) / 30)
    assort <- assort[assort$No != 0, c("tree", "No", "Sort", cnames[-1])]
  }
  return(assort)
}
