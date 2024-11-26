#' @title replace NAMESPACE from 'rBDATPRO' to 'rBDAT'.
#' @description function replaces calls to package 'rBDATPRO' by calls to
#' package 'rBDAT', simply by replacing the package name.
#' @param inpath path of file to process
#' @param outpath either a file name or a directory, both or NULL, see details.
#' @details This function merely exists to account for the renaming of the
#' package from \code{rBDATPRO} (which was internally used for some period of
#' time) to \code{rBDAT}. Its sole purpose is to update R-scripts which use
#' \code{rBDATPRO} and now should be updated to use \code{rBDAT}. Internally,
#' \code{gsub} is used.
#'
#' \code{outpath} can be (i) a filename, then the newly generated file
#' is stored under that name in the \code{inpath} directory, (ii) a directory,
#' then the \code{inpath}-filename is used with prefixed \code{rbdat_}, (iii) a
#' complete path (directory name and file name), then this is used to store the
#' file, (iv) NULL, then the inpath is used with prefixed \code{rbdat_} to the
#' inpath filename.
#' @return a character holding path and filename of the updated file
#' @examples
#' \dontrun{
#' p <- tempdir()
#' f <- "rbdatpro.r"
#' tx <- c("require(rBDATPRO)", "library(rBDATPRO)",
#'         "rBDATPRO::getDiameter(list(spp=1, D1=30, H=27))")
#' pf <- file.path(p, f)
#' writeLines(tx, con=pf)
#' file.exists(pf)
#' list.files(p)
#' # file.show(pf)
#'
#' ## define different output specs
#' outpath1 <- file.path(tempdir(), "devel/rbdatScript.r")
#' outpath2 <- p
#' outpath3 <- "rbdatScript.r"
#'
#' (updated_file <- updateBdatNamespace(pf, outpath = NULL))
#' list.files(p)
#' # file.show(updated_file)
#'
#' (updated_file <- updateBdatNamespace(pf, outpath1))
#' list.files(file.path(p, "devel"))
#' # file.show(updated_file)
#'
#' (updated_file <- updateBdatNamespace(pf, outpath2))
#' list.files(p)
#' # file.show(updated_file)
#'
#' (updated_file <- updateBdatNamespace(pf, outpath3))
#' list.files(p)
#' # file.show(updated_file)
#' }
#' @export

updateBdatNamespace <- function(inpath=file.choose(), outpath=NULL){

  if(!file.exists(inpath)){
    stop("file not found!")
  }

  ## checking outpath
  if(is.null(outpath)){
    outpath <- file.path(dirname(inpath), paste0("rbdat_", basename(inpath)))
  }
  if(!dir.exists(dirname(outpath))){
    de <- try(dir.create(dirname(outpath)), silent = TRUE) # dir exists
    if(inherits(de, "try-error")){
      stop("'outpath' not available!")
    }
  } else if(!grepl(".r", tolower(basename(outpath)))) {
    ## outpath is a directory name
    outpath <- file.path(outpath, paste0("rbdat_", basename(inpath)))
  } else if(identical(dirname(outpath), ".")){
    ## 'outpath' is a filename
    outpath <- file.path(dirname(inpath), outpath)
  }

  ## replace 'rBDATPRO' by 'rBDAT'
  tx = gsub(pattern = "rBDATPRO", replacement = "rBDAT", x = readLines(inpath))
  writeLines(tx, con=outpath)
  return(outpath)
}
