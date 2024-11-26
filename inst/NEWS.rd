\name{NEWS}
\title{NEWS}

\section{Changes in version 0.10.0}{
\itemize{
\item added information and example about 'inv=0' in getForm to return volume tables
(so-called 'Massetafeln' of Grundner und Schwappach, 1921) equivalent form
\item coerce 'bark'-parameter in getDiameter, getVolume and getHeight into logical
of length one.
\item removed LazyData entry in Description file.
\item updated buildTree(, check="biomass") to calculate D13 and D03 for correct
evaluation of the biomass function. It is now possible to hand over a
percentile of the D03-distribution with D2=0 and H2=]0,100[. Consequently,
the default behaviour has changed in case no D2 information is given: formerly
NFI1-form was assumed but now volume table-equivalence is assumed. This is now
coherent to all other functions in BDAT.
\item The consequence of allowing D1=0 in building the data set for evaluating
biomass for small trees, are wrong estimates for those small trees with
respect to taper curve evaluation and, hence, to other characteristics like
volume, double bark thickness and diameter. Now, those small trees correctly
return zero in those cases.
\item added PACKAGE argument to all .Fortran() calls
}
}

\section{Changes in version 0.9.9}{
\itemize{
\item now the error indication of subroutine BDAT20 is evaluated and respective
trees get treated accordingly, ie. set to NA
\item now allow for D1 = 0 and 0 < H1 < 2.5 so that very small trees can be
evaluated for biomass as indeed implemented.
\item Added citation for biomass functions.
}
}

\section{Changes in version 0.9.8}{
\itemize{
\item added more tolerance for tests on CRAN windows old-rel
}
}

\section{Changes in version 0.9.7}{
\itemize{
\item modifications to avoid warnings in LTO (link time optimization of gcc, option
-flto=10), especially in rBDAT_init.c
\item in Fortran function FnBiomass, use explicit declaration and corrected passing
of variables with required type
}
}

\section{Changes in version 0.9.6}{
\itemize{
\item added doi of citation in Description file
\item corrected par() re-setting in plot.datBDAT examples
}
}

\section{Changes in version 0.9.5}{
\itemize{
\item improved examples in several functions
}
}

\section{Changes in version 0.9.4}{
\itemize{
\item remove BDAT docs from package; now available at
https://gitlab.com/vochr/rbdat/-/blob/master/bdatdocs//
\item small updates in vignette and package info
}
}

\section{Changes in version 0.9.3}{
\itemize{
\item prepare for CRAN
\item update of Fortran Code to remove compile warnings
\item small updates in vignette
}
}

\section{Changes in version 0.9.2}{
\itemize{
\item call of getDiameter via S3-methods dispatch
}
}

\section{Changes in version 0.9.1}{
\itemize{
\item added function to easily update R-scripts: use 'rBDAT' instead of 'rBDATPRO'
\item added startup message
\item added .onUnload to detach .so/.dll
}
}

\section{Changes in version 0.9.0}{
\itemize{
\item This is a clone of https://gitlab.com/vochr/rbdatpro, hence, a rename of
R-Package 'rBDATPRO' v0.8.1.9000. With that, the package is named more
appropriately (as BDATPRO was a specific version of BDAT), shorter and this
is the canonical name of the provided functionality.
}
\subsection{in the devel version of rBDATPRO}{
\itemize{
\item Problems when calling getBiomass with very small trees: internal call to
getDiameter(, Hx=0.3*H, ) returns NA and hence getBiomass stops with error.
As first measure the precalculation of D03 inside getBiomass is stopped and
left to be done via FnBiomasse in Fortran. In consequence, D2 can not be
given as q03-quantile (if so, assumed form is set to H2=50, i.e. mean NFI2
taper form).
}
}
}

\section{Changes in version 0.8.0}{
\itemize{
\item Added original BDAT documentation to directory vignettes/.
(\code{pkgbuild::build()} will delete inst/doc prior to building the tarball.
Just don't use it if you don't have to,
see https://github.com/r-lib/pkgbuild/issues/58 or vignettes/build/issue58.html)
While building the tarball the BDAT documentation (via an entry in \code{vignettes/.install_extras})
will be copied into inst/doc, from where it will be installed into doc/ when
installing the tarball (set \code{build_vignettes = TRUE} if using \code{devtools::install()},
but avoid devtools whenever you can).
\item Fixed stale fortran codes.
\item Infected with packager.
\item Added a \code{NEWS.md} file to track changes to the package.
}
}

\section{Changes in version 0.7.3}{
\itemize{
\item error correction in getVolume: if nothing else than the tree parameter
is given, the function estimates standing coarse wood volume including bark;
if parameter bark was given as FALSE, still the volume in bark was returned,
now correct so that bark parameter is evaluated in that context
\item small correction in buildTree in explaining how parameter Az is
estimated if not given
}
}

\section{Changes in version 0.7.2}{
\itemize{
\item in plotting function the default names using dbh and height are now
rounded to first decimal
\item in getVolume test of inclusion of A and B parameter in parameter tree
now is more general
}
}

\section{Changes in version 0.7.1}{
\itemize{
\item some adjustments for CRAN submission
}
}

\section{Changes in version 0.7.0}{
\itemize{
\item In function BDATVOLDHOR changed the interpretation of the parameter
DHGrz into Dob (Diameter over bark) instead of Dub (Diameter under bark)
\item getVolume now returns total coarse wood volume over bark (Vfm m.R.)
instead of harvested coarse wood volume under bark (Efm o.R.)
\item \code{getBiomass} now returns appropriate results if H2 is given as
quantile or D2 as form quotient
\item consolidated the answers of the get\*-functions when one or several
trees are given and the optional parameter hx / dx is used
}
}

\section{Changes in version 0.6.5}{
\itemize{
\item error in plot function: no assortments for subsetted dat.BDAT (unless
first is selected), because wrong index was used. Corrected and example
added
}
}

\section{Changes in version 0.6.4}{
\itemize{
\item correction in v0.6.3 led to false calculations of stem height of
classical assortments (esp. \code{nvDh}) in certain cases
\item added missing Volume information for \code{X-Holz} in case D1<10
\item added assortment information into plot function (which has partly been
rewritten)
}
}

\section{Changes in version 0.6.3}{
\itemize{
\item in case of tree having D1 < 10 only industrial wood is assorted. Added
missing length-middiameter-topdiameter information for that case (calculated
in separate part of the Fortran code)
}
}

\section{Changes in version 0.6.2}{
\itemize{
\item corrected error in (adapted) Fortran code to test for almost zero
instead of to be equal zero (out-variable \code{LDSort} in BDAT20 subroutine)
}
}

\section{Changes in version 0.6.1}{
\itemize{
\item Corrected error in data checking for class type. It was possible to
call Fortran routines with erroneous data, which led to crashing R
}
}

\section{Changes in version 0.6.0}{
\itemize{
\item Functions to estimate the form factor (\code{getForm()}), i.e. mean q03 are
implemented
}
}

\section{Changes in version 0.5.0}{
\itemize{
\item added an experimental plotting function
}
}

\section{Changes in version 0.4.1}{
\itemize{
\item correct some faulty meta-data (help, examples)
}
}

\section{Changes in version 0.4.0}{
\itemize{
\item Vectorized fortran functions have been added for a quicker calculation
of many trees at once instead of using sapply (at least ten times faster)
}
}

\section{Changes in version 0.2.1}{
\itemize{
\item corrected error in transform BDAT20-output to proper format for output
of \code{getAssortment()}
}
}

\section{Changes in version 0.2.0}{
\itemize{
\item added extra output to \code{getAssortment()}: length and diameter
information about the classical assortments, now comparable to the output of
fix length assortments
}
}

\section{Changes in version 0.1.0}{
\itemize{
\item first version with Fortran scripts being compiled when installing
(instead of using pre-compiled Win-DLL), i.e usage of BDAT is now OS
independent and 64bit available
}
}

