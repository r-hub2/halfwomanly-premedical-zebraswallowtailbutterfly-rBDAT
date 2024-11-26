#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname){
  packageStartupMessage(paste0("This is rBDAT ",
                               utils::packageVersion("rBDAT")))
}

.onUnload <- function (libpath) {
  library.dynam.unload("rBDAT", libpath)
}
