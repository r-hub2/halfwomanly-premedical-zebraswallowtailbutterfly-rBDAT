#' clear errors from subroutine BDAT20
#'
#' Returns the NA'd data set
#'
#' @param data result of subroutine BDAT20
#' @return a corrected, i.e. NA'd version of the input
#' @keywords internal
clearError <- function(data){
  errInd <- which(data$Ifeh > 0)

  n <- 6
  idx <- (rep(errInd, each=n)-1) * n + 1:n
  data$Skl[idx] <- NA

  n <- 7
  idx <- (rep(errInd, each=n)-1) * n + 1:n
  data$Vol[idx] <- NA

  n <- 20
  idx <- (rep(errInd, each=n)-1) * n + 1:n
  data$LDSort[idx] <- NA

  data$BHD[errInd] <- NA

  # n <- 180
  # idx <- (rep(errInd, each=n)-1) * n + 1:n
  # data$FixLng[idx] <- NA

  # data$NFixLng[errInd] <- NA
  return(data)
}

