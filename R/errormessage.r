#' Return dedicated error message
#'
#' Returns the correct error message of SUBROUTINE BDAT20 according to help
#' file of BDAT
#'
#' @param code error codes from subroutine BDAT20 (parameter 'Ifeh')
#' @return returns the error message of the respective codes
#' @keywords internal
errormessage <- function(code){
  error <- c("Invalid tree species code",
             "missing tree height: H < 0",
             "lower diameter missing: D1 <= 0",
             "lower diameter height incorrect: H1 > 2.5 m",
             "(not used)",
             "stem height > 0 (i.e. Skz=1 & Hsh>0 & spp>14)",
             "twining below 7m and stem height > 7m in hardwoods: Skz = 2 and Hsh > 7; spp > 14",
             "Deciduous trees without recognizable trunk but with indication of a trunk height:
             Skz = 3 and (Hsh > 3 m or Hsh < 1.3 m); spp < 14",
             "Incorrect trunk height for broken or dead deciduous tree: Skz = 4 and H < Hsh; spp > 14",
             "Top-shafted coniferous tree with positive stem height: Skz = 0|1 and Hsh > 0; spp < 15",
             "Dead or broken coniferous wood without indication of a top or crown break: Skz = 4 and Hkz = 0; spp < 15",
             "no q0.30 estimation possible for given stem dimensions (H, D1)",
             "no volume table values for given stem dimensions (H, D1)",
             "q*0.30 estimation inaccurate, iteration does not pass D7m",
             "Dimensions too small, no volume table values for Dbh < 7cm",
             "unspecified error",
             "no error")
  code <- ifelse(code < 0 | code > 15, 16, code)
  code <- ifelse(code==0 , 17, code)
  return(error[code])
}

