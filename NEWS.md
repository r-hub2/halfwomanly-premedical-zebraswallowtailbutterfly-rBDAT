# rBDAT 1.1.0

* updated English species names in getSpeciesCode() to match to the names in 
  German NFI; the old version can still be accessed.
* bug fix in assortment function for cases where lX (length of unusable wood at 
  stem food) was larger than height of minimum cutting diameter
* adjusted check on lX in buildTree()

# rBDAT 1.0.1

* removing significant warning in r-devel-linux-x86_64-fedora-gcc when compiling
  the fortran code: variable declaration for implied DO loop in data block
* bug fix in error message referring to wrong variable (Skz instead of Sokz)
* added exception for cases where assortment parameter lX is too big, especially
  bigger than 80% of tree height.

# rBDAT 1.0.0

* current state of package considered version 1.0.0
* adjusted/removed comments in Fortran code
* updated helpfiles
* correct check on variable 'inv' in function getForm()

# rBDAT 0.10.0

* added information and example about 'inv=0' in getForm to return volume tables
  (so-called 'Massetafeln' of Grundner und Schwappach, 1921) equivalent form
* coerce 'bark'-parameter in getDiameter, getVolume and getHeight into logical 
  of length one.
* removed LazyData entry in Description file.
* updated buildTree(, check="biomass") to calculate D13 and D03 for correct 
  evaluation of the biomass function. It is now possible to hand over a 
  percentile of the D03-distribution with D2=0 and H2=]0,100[. Consequently, 
  the default behaviour has changed in case no D2 information is given: formerly
  NFI1-form was assumed but now volume table-equivalence is assumed. This is now
  coherent to all other functions in BDAT.
* The consequence of allowing D1=0 in building the data set for evaluating 
  biomass for small trees, are wrong estimates for those small trees with 
  respect to taper curve evaluation and, hence, to other characteristics like
  volume, double bark thickness and diameter. Now, those small trees correctly
  return zero in those cases.
* added PACKAGE argument to all .Fortran() calls

# rBDAT 0.9.9

* now the error indication of subroutine BDAT20 is evaluated and respective 
  trees get treated accordingly, ie. set to NA
* now allow for D1 = 0 and 0 < H1 < 2.5 so that very small trees can be 
  evaluated for biomass as indeed implemented. 
* Added citation for biomass functions.

# rBDAT 0.9.8

* added more tolerance for tests on CRAN windows old-rel

# rBDAT 0.9.7

* modifications to avoid warnings in LTO (link time optimization of gcc, option
  -flto=10), especially in rBDAT_init.c
* in Fortran function FnBiomass, use explicit declaration and corrected passing
  of variables with required type

# rBDAT  0.9.6

* added doi of citation in Description file
* corrected par() re-setting in plot.datBDAT examples

# rBDAT 0.9.5

* improved examples in several functions

# rBDAT 0.9.4

* remove BDAT docs from package; now available at
  https://gitlab.com/vochr/rbdat/-/blob/master/bdatdocs//
* small updates in vignette and package info

# rBDAT 0.9.3

* prepare for CRAN
* update of Fortran Code to remove compile warnings
* small updates in vignette

# rBDAT 0.9.2

* call of getDiameter via S3-methods dispatch

# rBDAT 0.9.1

* added function to easily update R-scripts: use 'rBDAT' instead of 'rBDATPRO'
* added startup message
* added .onUnload to detach .so/.dll

# rBDAT 0.9.0

* This is a clone of https://gitlab.com/vochr/rbdatpro, hence, a rename of 
  R-Package 'rBDATPRO' v0.8.1.9000. With that, the package is named more
  appropriately (as BDATPRO was a specific version of BDAT), shorter and this 
  is the canonical name of the provided functionality.

## in the devel version of rBDATPRO
* Problems when calling getBiomass with very small trees: internal call to 
  getDiameter(, Hx=0.3*H, ) returns NA and hence getBiomass stops with error.
  As first measure the precalculation of D03 inside getBiomass is stopped and
  left to be done via FnBiomasse in Fortran. In consequence, D2 can not be
  given as q03-quantile (if so, assumed form is set to H2=50, i.e. mean NFI2 
  taper form).

# rBDATPRO 0.8.0

* Added original BDAT documentation to directory vignettes/.
  (`pkgbuild::build()` will delete inst/doc prior to building the tarball. 
   Just don't use it if you don't have to, 
   see https://github.com/r-lib/pkgbuild/issues/58 or vignettes/build/issue58.html)
  While building the tarball the BDAT documentation (via an entry in `vignettes/.install_extras`)
  will be copied into inst/doc, from where it will be installed into doc/ when
  installing the tarball (set `build_vignettes = TRUE` if using `devtools::install()`, 
  but avoid devtools whenever you can).
* Fixed stale fortran codes.
* Infected with packager.
* Added a `NEWS.md` file to track changes to the package.



# rBDATPRO 0.7.3
-  error correction in getVolume: if nothing else than the tree parameter
is given, the function estimates standing coarse wood volume including bark;
if parameter bark was given as FALSE, still the volume in bark was returned,
now correct so that bark parameter is evaluated in that context
-  small correction in buildTree in explaining how parameter Az is
estimated if not given

# rBDATPRO 0.7.2
-  in plotting function the default names using dbh and height are now
rounded to first decimal
-  in getVolume test of inclusion of A and B parameter in parameter tree
now is more general

# rBDATPRO 0.7.1
-  some adjustments for CRAN submission

# rBDATPRO 0.7.0
-  In function BDATVOLDHOR changed the interpretation of the parameter
DHGrz into Dob (Diameter over bark) instead of Dub (Diameter under bark)
-  getVolume now returns total coarse wood volume over bark (Vfm m.R.)
instead of harvested coarse wood volume under bark (Efm o.R.)
-  \code{getBiomass} now returns appropriate results if H2 is given as
quantile or D2 as form quotient
-  consolidated the answers of the get\*-functions when one or several
trees are given and the optional parameter hx / dx is used

# rBDATPRO 0.6.5
-  error in plot function: no assortments for subsetted dat.BDAT (unless
first is selected), because wrong index was used. Corrected and example
added

# rBDATPRO 0.6.4
-  correction in v0.6.3 led to false calculations of stem height of
classical assortments (esp. `nvDh`) in certain cases
-  added missing Volume information for `X-Holz` in case D1<10
-  added assortment information into plot function (which has partly been
rewritten)

# rBDATPRO 0.6.3
-  in case of tree having D1 < 10 only industrial wood is assorted. Added
missing length-middiameter-topdiameter information for that case (calculated
in separate part of the Fortran code)

# rBDATPRO 0.6.2
-  corrected error in (adapted) Fortran code to test for almost zero
instead of to be equal zero (out-variable `LDSort` in BDAT20 subroutine)

# rBDATPRO 0.6.1
-  Corrected error in data checking for class type. It was possible to
call Fortran routines with erroneous data, which led to crashing R

# rBDATPRO 0.6.0
-  Functions to estimate the form factor (`getForm()`), i.e. mean q03 are
implemented

# rBDATPRO 0.5.0
-  added an experimental plotting function

# rBDATPRO 0.4.1
-  correct some faulty meta-data (help, examples)

# rBDATPRO 0.4.0
-  Vectorized fortran functions have been added for a quicker calculation
of many trees at once instead of using sapply (at least ten times faster)

# rBDATPRO 0.2.1
-  corrected error in transform BDAT20-output to proper format for output
of `getAssortment()`

# rBDATPRO 0.2.0
-  added extra output to `getAssortment()`: length and diameter
information about the classical assortments, now comparable to the output of
fix length assortments

# rBDATPRO 0.1.0
-  first version with Fortran scripts being compiled when installing
(instead of using pre-compiled Win-DLL), i.e usage of BDAT is now OS
independent and 64bit available
