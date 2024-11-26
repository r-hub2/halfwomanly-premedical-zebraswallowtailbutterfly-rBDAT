---
title: "BDAT Implementation in R"
author: "christian.vonderach@forst.bwl.de"
date: 2020-10-08, 16:43:33
output: rmarkdown::html_vignette
# output: pdf_document
vignette: >
  %\VignetteIndexEntry{BDAT Implementation in R}
  \usepackage[utf8]{inputenc}
  %\VignetteEngine{knitr::knitr}
---



## BDAT: a Fortran library for taper curves

*BDAT* is a Fortran program library to estimate volume, diameter, height, bark
thickness and assortments of a specified tree. Hence, it is a library about
taper functions. It has been developed at FVA-BW on request of BMELV for the 
first german national forest inventory (BWI 1, 1987). Adjustements have been 
made for the BWI2 (2002). *BDATPro* is a version including these latest 
adjustments but also an additionaly GUI. The latest version includes biomass
functions and form parameters for NFI3 (=BWI3). Here the latest version was
implemented.

Colloquially the library is called **BDAT**.

The Fortran library has been extensively used in several application, of which
the most important might be WEHAM (WaldEntwicklungs- and HolzAufkommensModell:
forest developement and timber stock prediction model) but also in other forest
inventories. At that time, R was not emerged, yet, later direct calls to the
Windows-DLL were used. With the present R-package, the usage of BDAT in R was 
highly simplified, ported to other operating systems (32bit and 64bit) and 
better documented. There are some original documents on my gitlab repository
about the methodology and application, but all have been written in German.

The BDAT library is based on a spline representation of the taper function of
different tree species. A methodological improvement has been made using 
B-Splines and mixed-modelling, see Kublin et al. 2013 and the R-package TapeR.
Still, BDAT is in use and keeps being used since the library offers more than
just the taper form itself: estimation of diameters, height, double bark 
thickness, volume of sections defined by height and/or diameter, assortments
given parameters and all already parameterised for a whole bunch of tree species
including deciduous tree species with their complex tree crowns. It is based on
approximately 30.000 measured trees and studies on bark thickness and tree crown
volume.

## further references
Further references like BDAT Documentations, related articles and reports can be
found at https://gitlab.com/vochr/rbdat especially at
https://gitlab.com/vochr/rbdat/-/tree/master/bdatdocs, 
unfortunately these resources are all written in German.

## Functions
The R-package contains all relevant functions from the Fortran-library and uses 
vectorized evaluation. It is recommended to use the get*-functions, but for
convenience wrapper functions using the names of the Fortan-subroutines are 
included so that older scripts can easily be adapted to the use of the 
R-package.

Beside the core-functions (buildTree, getDiameter, getHeight, getVolume, 
getAssortment, getBiomass, getBark, getForm and getSpeciesCode) there is a 
plotting function. 

### preliminary note
The use of the BDAT Fortran functions requires the preparation of the data to
conform with what Fortran is expecting to come. Within this R-Package, this is
implemented in two ways: either one can prepare all necessary variables within
a list or data frame and pass it over to the respective function via the 
`tree` parameter or one can pass a tree definition into parameter `tree` of each
function and use the second parameter (i.e. `Dx`, `Hx`, `AB` or `sort`, here 
called `vars`) to hand over the required - function specific - information. In 
the first case, the function returns one result for each row given. In the 
second case, a cross join / cartesian product between `tree` and
`vars` is calculated. If `vars` is of size one, the results is the same as in 
the first case. If the size of `vars` is bigger than one, the functions return
one value for each given tree (e.g. 3) and element of `vars` (e.g. 4), in the
example this is 12. For functions returning a scalar, a matrix is returned with
`trees` given in rows and `vars` given in columns (e.g. 3x4-matrix). The
assortment function, usually returning a data frame, still returns a data frame,
but now the order (and naming) of trees is different from the given input `tree`
parameter, since internally the tree object is now expanded. Hence, the first
tree is repeated (length-of-vars) number of times, before the second tree is
process and returned and so on.

### Getting started
In BDAT, a tree is specified at its minimum by its species code, dbh (diameter 
in breast height, i.e. 1.3m) and height. Implicitly, this defines a second 
diameter in 30% of three height according to the "Masse-Tafeln" (Volume-Tables) 
from Grundner und Schwappach (1921). For the possibility of more precisely
specifying taper form see further below.


```r
library(rBDAT)
tree <- buildTree(tree = list(spp=1, D1=30, H=27))
str(tree)
```

```
## Classes 'datBDAT' and 'data.frame':	1 obs. of  6 variables:
##  $ spp: num 1
##  $ D1 : num 30
##  $ H  : num 27
##  $ H1 : num 1.3
##  $ D2 : num 0
##  $ H2 : num 0
```

One can visualise the taper curve of a given tree using the plot-function:


```r
plot(tree)
```

![plot of chunk pressure](figure/pressure-1.png)
Here, the taper curve over bark (black) and under bark (grey) is drawn.

### getting BDAT-species code
BDAT has been parameterised for 36 tree species, more or less common in Germany
based on about 30.000 trees. These 36 tree species are index and each posesses 
its own BDAT-species code. This code and the respective species name (short and 
long format), english species name and scientific name can be retrieved by the
getSpeciesCode-function:

```r
getSpeciesCode()
```

```
##    ID kurz            lang short                long            scientific
## 1   1   Fi          Fichte    NS       Norway spruce           Picea abies
## 2   2   SF    Sitka-Fichte    SS        Sitka spruce      Picea sitchensis
## 3   3   Ta           Tanne   ESF European silver fir            Abies alba
## 4   4   KT    Kuestentanne    GF           Grand fir         Abies grandis
## 5   5  Kie          Kiefer    SP         Scotch pine      Pinus sylvestris
## 6   6   SK   Schwarzkiefer   AUP       Austrian pine           Pinus nigra
## 7   7   WK Weymouthskiefer   WEP       Weymouth pine         Pinus strobus
## 8   8   DG       Douglasie    DF         Douglas fir Pseudotsuga menziesii
## 9   9   La         Laerche   XLA               larch            Larix spp.
## 10 10   EL  Europ. Laerche   ELA      European larch         Larix decidua
## 11 11   JL    Jap. Laerche   JLA      Japanese larch       Larix kaempferi
## 12 12   Th           Thuja    RC   Western red cedar         Thuja plicata
## 13 13   Ts           Tsuga    WH     Western hemlock    Tsuga heterophylla
## 14 14   SN             sNB    XC      other conifers     Coniferales trees
## 15 15   Bu           Buche    BE               beech       Fagus sylvatica
## 16 16   HB       Hainbuche   HBM            hornbeam      Carpinus betulus
## 17 17   Ei           Eiche    OK oak (robur/petraea)          Quercus spp.
## 18 18   RE        Roteiche   ROK             Red oak         Quercus rubra
## 19 19   Pa          Pappel   XPO              poplar          Populus spp.
## 20 20   BP    Balsampappel   BPO       Balsam poplar   Populus balsamifera
## 21 21   Es           Esche    AH                 ash    Fraxinus excelsior
## 22 22   Ah           Ahorn   XAH               maple             Acer spp.
## 23 23   BA       Bergahorn    SY            sycamore   Acer pseudoplatanus
## 24 24   SA      Spitzahorn   NOM        Norway maple      Acer platanoides
## 25 25   FA       Feldahorn    FM         Field maple        Acer campestre
## 26 26   Bi           Birke   XBI               birch           Betula spp.
## 27 27   Li           Linde    LI           lime tree            Tilia spp.
## 28 28   Er            Erle    AR               alder            Alnus spp.
## 29 29  Kir         Kirsche   WCH         Wild cherry          Prunus avium
## 30 30   Ul            Ulme    EM                 elm            Ulmus spp.
## 31 31   Ro         Robinie    BL        black locust  Robinia pseudoacacia
## 32 32   El        Elsbeere   WST   Wild service tree     Sorbus torminalis
## 33 33   Ka        Kastanie    SC      Sweet chestnut       Castanea sativa
## 34 34   We           Weide   XWL              willow            Salix spp.
## 35 35   LB             sLB    XB   other broadleaves   Magnoliopsida trees
## 36 36   VB      Vogelbeere   ROW               rowan      Sorbus aucuparia
```
If the function is called without any parameter, a data.frame is returned 
holding the information which can be retrieved. 
One can alternatively specify the type of `input` and `output`. `Input`
must be one of the data entries, `output` must be one of the column names.

```r
getSpeciesCode(1)
```

```
##   ID kurz   lang short          long  scientific
## 1  1   Fi Fichte    NS Norway spruce Picea abies
```

```r
getSpeciesCode("Fi")
```

```
## [1] 1
```

```r
getSpeciesCode("NS") # english abbreviation of Norway spruce
```

```
## [1] 1
```

```r
getSpeciesCode(1, "scientific")
```

```
## [1] "Picea abies"
```

```r
getSpeciesCode(c(1:36), "short")
```

```
##  [1] "NS"  "SS"  "ESF" "GF"  "SP"  "AUP" "WEP" "DF"  "XLA" "ELA" "JLA" "RC" 
## [13] "WH"  "XC"  "BE"  "HBM" "OK"  "ROK" "XPO" "BPO" "AH"  "XAH" "SY"  "NOM"
## [25] "FM"  "XBI" "LI"  "AR"  "WCH" "EM"  "BL"  "WST" "SC"  "XWL" "XB"  "ROW"
```


### getting diameter in given height
The taper functions specify a curve being a function of height within tree and 
returns the respective diameter. Hence, one can evalutate this function at a 
given height and receive the diameter. There are fortran functions to return 
either diameter over bark and diameter under bark, but this package binds 
these functions together and uses a boolean parameter to switch between both.


```r
getDiameter(tree, Hx = 1.3) # return dbh
```

```
## [1] 30
```

As a default, diameter over bark is returned, but this can easily be changed:


```r
getDiameter(tree, Hx = 1.3, bark = FALSE) # return d.u.b. at 1.3m
```

```
## [1] 28.40971
```

It is possible to request diameter in several heights at once:

```r
getDiameter(tree, Hx = c(1:10))
```

```
##  [1] 31.14922 28.44761 27.57622 26.92106 26.32229 25.76090 25.21790 24.67427
##  [9] 24.11234 23.52139
```

Using this kind of call to the functions, `Hx` is evaluated for each given tree 
separately:

```r
tree2 <- buildTree(list(spp=1, D1=c(30, 35), H=27))
getDiameter(tree2, Hx = c(1, 2))
```

```
##          [,1]     [,2]
## [1,] 31.14922 28.44761
## [2,] 36.54927 32.92451
```
Here, a matrix is returned with one row per tree and one column per requested 
Hx. 

An alternative way is specifying `Hx` directly inside the tree-object:

```r
tree2 <- buildTree(list(spp=1, D1=c(30, 35), H=27, Hx=c(1, 2)))
tree2
```

```
##   spp D1  H Hx  H1 D2 H2
## 1   1 30 27  1 1.3  0  0
## 2   1 35 27  2 1.3  0  0
```

```r
getDiameter(tree2)
```

```
## [1] 31.14922 32.92451
```
Here, tree2 consists of 2 rows, because the give list is internally transformed
into a data frame (common rules for building data.frames apply). Finally, the
object is evaluated row-wise. Because `Hx` is already given, the parameter `Hx`
can be left empty.


### getting double bark thickness
Double bark thickness is that part of a diameter over bark, which is considered 
to consist of bark tissue. The relation between wood and bark with respect to 
diameter can be expressed as {double bark thickness} + {diameter under bark} =
{diameter over bark}. The implemented functions originate from the works of 
Altherr et al. (1974, 1975, 1976, 1978 and 1979).

```r
dub <- getDiameter(tree, Hx=1.3, bark = FALSE)
dob <- getDiameter(tree, Hx=1.3, bark = TRUE)
dbt <- getBark(tree, Hx=1.3)
dub + dbt == dob
```

```
## [1] TRUE
```

Again, it is possible to either include parameter `Hx` into the tree-object or 
pass it separately. In the second case, an matrix is returned if Hx is longer 
than one.

```r
getBark(tree2, Hx = 1:5)
```

```
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 1.636124 1.527308 1.491416 1.464177 1.439089
## [2,] 1.842470 1.705602 1.661812 1.629562 1.600471
```


### getting height for given diameter
The diameter-height-relation from above, where we evaluated the function for a
given height-value, can also be evaluated for a given diameter. The internal
function is an iterative procedure to determine height. Diameter can be
specified over and under bark. Default is bark = TRUE. 

```r
getHeight(tree, Dx=30) # height of diameter over bark
```

```
## [1] 1.3
```

```r
getHeight(tree, Dx=30, bark = FALSE) # height of diameter under bark == 30
```

```
## [1] 0.8916153
```

Again, it is possible to vary passing of the `Dx`-parameter: it can take one or 
several values, but can also be included to the tree parameter. In the first
case, a cross join / cartesian product between `tree` and `Dx` is created, 
in the second case the `tree` object is processed *as is*.

```r
getHeight(tree2, Dx=c(30, 20, 10)) # returns value in meters
```

```
##          [,1]     [,2]     [,3]
## [1,] 1.300000 14.82080 22.65465
## [2,] 5.370549 17.27376 23.05898
```

```r
tree2$Dx = c(30, 20)
getHeight(tree2)
```

```
## [1]  1.30000 17.27376
```
As one can see, in the first case a matrix with one row per tree and one column 
per `Dx` is returned and in the second call, a vector with one element for each
row in `tree2`.


### getting volume
The maybe most interesting function is the one returning volume. The function 
includes a switch to return volume with (the default) or without bark volume as
well. Volume is calculated using middiameter-formula (in Germany called 
"Huber'sche Formel") from 2m-sections (default), but section lengths can be 
varied. Additionally, it is necessary to specify the section for which volume
is required. This can be done either using diameters or heights, or a mixture of
both. If parameter `AB` is not given, i.e. NULL, then the function assumes 
coarse wood volume over bark is required (i.e. from forest floor up to diameter
over bark of 7cm):

```r
getVolume(tree) # get coarse wood, which is the same as next line
```

```
## [1] 0.9168622
```

```r
getVolume(tree, AB = list(A=0, B=7), iAB=c("H", "Dob"), bark = T) 
```

```
## [1] 0.9168622
```

One can precisely specify the section for which volume under or over bark is
required:

```r
# volume including bark between height 1m and 2m
getVolume(tree2, AB=list(A=1, B=2)) 
```

```
## [1] 0.06184182 0.08197092
```

```r
# volume excluding bark between 30 and 7cm in diameter over bark
getVolume(tree2, AB=list(A=30, B=7), iAB="dob", bark = F) 
```

```
## [1] 0.7153726 0.6583801
```

The section length for which middiameter-formula is applied is 2m as a default.
It is possible to change that behaviour by setting parameter `sl` inside the 
`AB`-argument:

```r
# again: coarse wood volume
getVolume(tree) 
```

```
## [1] 0.9168622
```

```r
# identical
getVolume(tree, AB=list(A=0.27, B=7, sl=2), iAB=c("H", "Dob"), bark=F) 
```

```
## [1] 0.7892534
```

```r
# using sl=0.1, that is section length for volume calculation set to be 0.1m
getVolume(tree, AB=list(A=0.27, B=7, sl=0.1), iAB=c("H", "Dob"), bark=F) 
```

```
## [1] 0.7971188
```

If one wants to get the volume for several sections of a tree, one could either
build a suited data.frame on its own, repeating the tree attributes and
adjusting the `A`and `B` values as needed. The more elegant way is to let the
functions do that work for you. In the easiest case, we saw that already in the
example above, specifying one tree and one section, where the function returns
exactly one volume. Internally, the `AB` data is merged to the tree data by a
cross join (or cartesion product), hence we can make use of that behaviour by
defining several sections at once:

```r
getVolume(tree, AB=list(A=0:9, B=1:10))
```

```
##  [1] 0.09056835 0.06184182 0.06128590 0.05816483 0.05563769 0.05319688
##  [7] 0.05102357 0.04886994 0.04674566 0.04458112
```

```r
getVolume(tree = tree2, AB=list(A=0:9, B=1:10))
```

```
##            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
## [1,] 0.09056835 0.06184182 0.06128590 0.05816483 0.05563769 0.05319688
## [2,] 0.12786366 0.08197092 0.08171326 0.07714033 0.07358742 0.07018656
##            [,7]       [,8]       [,9]      [,10]
## [1,] 0.05102357 0.04886994 0.04674566 0.04458112
## [2,] 0.06725162 0.06435400 0.06151795 0.05862010
```


### getting biomass
Beside getting volume section information, there is a function to get total
aboveground biomass. These biomass-functions were fitted independently of the
parameterisation of BDAT and were developed in preparation of the german NFI3.
These functions are based on 983 analysed trees of a subset of species only.
The other species are either fit at a synthetical data set or in worst case
subsumed to other species. These functions are the official ones used during
reporting of results of the NFI3.

The call is identical to the already shown pattern:

```r
getBiomass(tree)
```

```
## [1] 406.1523
```

```r
getBiomass(list(spp=1, D1=30, H=27, D2=c(23, 24, 25))) 
```

```
## [1] 350.5038 375.2445 400.5624
```




### getting assortments
One very nice feature of the BDAT program library is its ability to use the
presented functions to simulate roundwoods from given trees. For that purpose
a separate function was written, which is called via `getAssortment()`. 
Similarly to other functions, it requires data of one or more trees
and optionally parameters to control the sorting process. There are quite a few
parameters to be specified. It's easiest to show that using some examples:

```r
getAssortment(tree2) # using standard assortment parameters
```

```
##    tree No Sort height    length      midD      topD        Vol
## 1     1  1    X   0.00  0.000000  0.000000  0.000000 0.00000000
## 2     1  2  Sth   0.27 20.000000 21.250875 12.429831 0.70937115
## 3     1  3   Ab   0.00  0.000000  0.000000  0.000000 0.00000000
## 4     1  4  Ind  20.47  1.000000 11.460602 10.725959 0.01031584
## 5     1  5 nvDh  21.47  2.717499  8.529173  5.998856 0.01552648
## 6     2  1    X   0.00  0.000000  0.000000  0.000000 0.00000000
## 7     2  2  Sth   0.27 20.000000 24.499969 13.956373 0.94286811
## 8     2  3   Ab   0.00  0.000000  0.000000  0.000000 0.00000000
## 9     2  4  Ind  20.47  1.000000 12.778620 11.893189 0.01282501
## 10    2  5 nvDh  21.47  2.904999  9.106683  5.990857 0.01892153
```
By default, a data.frame is returned with one row for each tree and roundwood 
piece. It keeps information about the roundwood foot position inside the stem,
its length (without add-on), mid-diameter under bark (midD), top-diameter under
bark (topD) and the respective volume under bark. The standard assortments 
comprise stem-wood (Sth), second length stem wood (Ab, after optionally 
specified cutting diameter for Sth or after transportation-cut, i.e. 20m),
industrial roundwood (Ind) and residual coarse wood (nvDh, between cutting 
diameter for Ind and 7cm over bark). Assortment X is unusuable wood at the stem
foot.

If one is interested in raw BDAT return, one can specify `value='raw'`:

```r
getAssortment(tree2, value = "raw") # usually of little interest...
```

```
## $n
## [1] 2
## 
## $BDATArtNr
## [1] 1 1
## 
## $D1
## [1] 30 35
## attr(,"Csingle")
## [1] TRUE
## 
## $H1
## [1] 1.3 1.3
## attr(,"Csingle")
## [1] TRUE
## 
## $D2
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $H2
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $H
## [1] 27 27
## attr(,"Csingle")
## [1] TRUE
## 
## $lX
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $Hkz
## [1] 0 0
## 
## $Skz
## [1] 0 0
## 
## $Az
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $Hsh
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $Zsh
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $Zab
## [1] 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $Sokz
## [1] 1 1
## 
## $Skl
##  [1] 0 0 2 0 0 0 0 0 2 0 0 0
## 
## $Vol
##  [1] 0.91686225 0.00000000 0.70937115 0.00000000 0.01031584 0.01552648
##  [7] 0.18164879 1.21056700 0.00000000 0.94286811 0.00000000 0.01282501
## [13] 0.01892153 0.23595232
## attr(,"Csingle")
## [1] TRUE
## 
## $LDSort
##  [1]  0.000000  0.000000  0.000000  0.000000  0.270000 20.000000 21.250875
##  [8] 12.429831  0.000000  0.000000  0.000000  0.000000 20.470001  1.000000
## [15] 11.460602 10.725959 21.470001  2.717499  8.529173  5.998856  0.000000
## [22]  0.000000  0.000000  0.000000  0.270000 20.000000 24.499969 13.956373
## [29]  0.000000  0.000000  0.000000  0.000000 20.470001  1.000000 12.778620
## [36] 11.893189 21.470001  2.904999  9.106683  5.990857
## attr(,"Csingle")
## [1] TRUE
## 
## $BHD
## [1] 30 35
## attr(,"Csingle")
## [1] TRUE
## 
## $Ifeh
## [1] 0 0
## 
## $FixLngDef
## [1] 0 0 0 0 0 0 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $NMaxFixLng
## [1] 0 0
## 
## $FixLng
##   [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
##  [38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
##  [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [112] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [149] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [186] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [223] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [260] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [297] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [334] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## attr(,"Csingle")
## [1] TRUE
## 
## $NFixLng
## [1] 0 0
```
Other options are also available like 'Vol' (Volume), 'Skl' (Stärkeklassen, i.e.
diameter classes), 'Fix' (fixed length assortments, if specified) and 'LDSort'
(an added feature w.r.t. the original BDAT Fortran code, keeping length and 
diameter of the assortments, which was missing in the original BDAT library,
hence only available in rBDAT!), which return a subset of the raw BDAT 
return value. The default is 'merge' which produces an aggregated data.frame of
relevant information about the roundwoods produced for each tree.

As said before, assortment rules can be specified using several parameters, 
which in detail are given in the help file. Some examples follow:

```r
getAssortment(tree, sort=list(Az=15)) # minimum diameter o.b. for assortments
```

```
##   tree No Sort height    length     midD      topD        Vol
## 1    1  1    X  0.000  0.000000  0.00000  0.000000 0.00000000
## 2    1  2  Sth  0.270 19.100000 21.51180 13.584801 0.69418710
## 3    1  3   Ab  0.000  0.000000  0.00000  0.000000 0.00000000
## 4    1  4  Ind  0.000  0.000000  0.00000  0.000000 0.00000000
## 5    1  5 nvDh 19.561  4.626499 10.10419  5.998856 0.03709757
```

```r
getAssortment(tree, sort=list(Az=15, Hsh=10)) # Hsh= max. height of sawlog quality
```

```
## Warning in getAssortment.datBDAT(tree, sort = list(Az = 15, Hsh = 10)): error indication in subroutine BDAT20:
## tree 1: Top-shafted coniferous tree with positive stem height: Skz = 0|1 and Hsh > 0; spp < 15
```

```
##   tree No Sort height length midD topD Vol
## 1    1  1    X     NA     NA   NA   NA  NA
## 2    1  2  Sth     NA     NA   NA   NA  NA
## 3    1  3   Ab     NA     NA   NA   NA  NA
## 4    1  4  Ind     NA     NA   NA   NA  NA
## 5    1  5 nvDh     NA     NA   NA   NA  NA
```
Additionally one can specify a fixed length assortment at stem foot, which is
cut before sawlogs.

```r
## fixN = number of roundwood pieces
## fixL = length of roundwood pieces in m
## fixA = absolute length addition (good for measure) in cm
## fixR = relative length addition in %
## fixZ = minimum cutting diameter for this assortment (Z=Zopf) in cm
getAssortment(tree, sort=list(Az=15, fixN=2, fixL=4, fixA=10, fixZ=20)) 
```

```
##   tree No  Sort    height  length     midD      topD       Vol
## 1    1  1     X  0.000000  0.0000  0.00000  0.000000 0.0000000
## 2    1  2   Sth  8.470001 10.9000 19.04127 13.584801 0.3103904
## 3    1  3    Ab  0.000000  0.0000  0.00000  0.000000 0.0000000
## 4    1  4   Ind  0.000000  0.0000  0.00000  0.000000 0.0000000
## 5    1  5  nvDh 19.479000  4.7085 10.16843  5.998856 0.0382367
## 6    1  6 Fix01  0.270000  4.0000 25.86922 24.547409 0.2102405
## 7    1  7 Fix02  4.370000  4.0000 23.40202 22.309429 0.1720507
```
In case the considered tree exhibits rotten/decaying wood at stem foot, this
can also be specified (in german termed X-Holz):

```r
getAssortment(tree, sort=list(lX=1.4)) # remove 1.4m from stump upwards
```

```
##   tree No Sort height    length      midD      topD        Vol
## 1    1  1    X   0.27  1.400000 28.892862 26.715826 0.09179077
## 2    1  2  Sth   1.67 20.000000 20.387375 10.421552 0.65289372
## 3    1  3   Ab   0.00  0.000000  0.000000  0.000000 0.00000000
## 4    1  4  Ind   0.00  0.000000  0.000000  0.000000 0.00000000
## 5    1  5 nvDh  21.87  2.317499  8.178885  5.998856 0.01217580
```

Additionally, one can plot the taper curve and include the assortments:

```r
assort <- getAssortment(tree, sort=list(Az=15, fixN=2, fixL=4, fixA=10, fixZ=20))
plot(tree, assort=assort)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)

#### using more than one assortment specification 
`buildTree` uses the parameters `tree` and `vars` which are merged
by a cross join / cartesian product. Hence, if assortment specifications are
extended to a length bigger one, the function should return an estimated 
assortment for each tree and specified assortment. Let's check that:

```r
getAssortment(tree, sort = list(Az=c(10, 7)))
```

```
##    tree No Sort height     length      midD      topD         Vol
## 1     1  1    X   0.00  0.0000000  0.000000  0.000000 0.000000000
## 2     1  2  Sth   0.27 20.0000000 21.250875 12.429831 0.709371150
## 3     1  3   Ab   0.00  0.0000000  0.000000  0.000000 0.000000000
## 4     1  4  Ind  20.47  2.0000000 10.725959  9.139494 0.018071413
## 5     1  5 nvDh  22.47  1.7174988  7.639552  5.998856 0.007872671
## 6     2  1    X   0.00  0.0000000  0.000000  0.000000 0.000000000
## 7     2  2  Sth   0.27 20.0000000 21.250875 12.429831 0.709371150
## 8     2  3   Ab   0.00  0.0000000  0.000000  0.000000 0.000000000
## 9     2  4  Ind  20.47  3.0000000  9.953113  7.379699 0.023341509
## 10    2  5 nvDh  23.47  0.7174988  6.702212  5.998856 0.002531322
```

```r
getAssortment(tree2, sort = list(Az=c(10, 7)))
```

```
##    tree No Sort height     length      midD      topD         Vol
## 1     1  1    X   0.00  0.0000000  0.000000  0.000000 0.000000000
## 2     1  2  Sth   0.27 20.0000000 21.250875 12.429831 0.709371150
## 3     1  3   Ab   0.00  0.0000000  0.000000  0.000000 0.000000000
## 4     1  4  Ind  20.47  2.0000000 10.725959  9.139494 0.018071413
## 5     1  5 nvDh  22.47  1.7174988  7.639552  5.998856 0.007872671
## 6     2  1    X   0.00  0.0000000  0.000000  0.000000 0.000000000
## 7     2  2  Sth   0.27 20.0000000 24.499969 13.956373 0.942868114
## 8     2  3   Ab   0.00  0.0000000  0.000000  0.000000 0.000000000
## 9     2  4  Ind  20.47  2.0000000 11.893189 10.009459 0.022218592
## 10    2  5 nvDh  22.47  1.9049988  8.071947  5.990857 0.009748576
## 11    3  1    X   0.00  0.0000000  0.000000  0.000000 0.000000000
## 12    3  2  Sth   0.27 20.0000000 21.250875 12.429831 0.709371150
## 13    3  3   Ab   0.00  0.0000000  0.000000  0.000000 0.000000000
## 14    3  4  Ind  20.47  3.0000000  9.953113  7.379699 0.023341509
## 15    3  5 nvDh  23.47  0.7174988  6.702212  5.998856 0.002531322
## 16    4  1    X   0.00  0.0000000  0.000000  0.000000 0.000000000
## 17    4  2  Sth   0.27 20.0000000 24.499969 13.956373 0.942868114
## 18    4  3   Ab   0.00  0.0000000  0.000000  0.000000 0.000000000
## 19    4  4  Ind  20.47  3.0000000 10.970320  7.971595 0.028356308
## 20    4  5 nvDh  23.47  0.9049988  6.997648  5.990857 0.003480503
```
Here, the first call returns the estimated assortments for both assortment
specifications, the tree order inside the resulting object is clear: tree 1 uses
Az=10 and tree 2 uses Az=7. In the second call, the order is the same, although
this is not as clear as in the first example.

### more advanced specification of trees
As shown, a tree is specified by its species code, its dbh and height. If done
so, this assumes a predefined diameter in 30% of tree height, which in turn 
defines the taper form or the relation between diameter in 30% and 5% of the 
tree height. But easily, this can be changed: 
the BWI-equivalent taper form can be specified using H2=50. 
The D2 and H2 parameter can take a very flexible parameterisation being diameter
and heights but also quantiles of the q03-distribution (H2) or the q03-parameter
itself (D2). q03 is the quotient of diameter in 30% of three height and in 5% of
tree height. Hence, quite a variety of taper forms can be represented. All 
BDAT-functions respond sensitively to this fourth parameter.

```r
tree <- buildTree(tree = list(spp=1, D1=30, H=27, H2=c(30, 50, 99, 0), 
                              D2=c(0, 0, 0, -0.8), Hx=0.3*27))
str(tree)
```

```
## Classes 'datBDAT' and 'data.frame':	4 obs. of  7 variables:
##  $ spp: num  1 1 1 1
##  $ D1 : num  30 30 30 30
##  $ H  : num  27 27 27 27
##  $ H2 : num  30 50 99 0
##  $ D2 : num  0 0 0 -0.8
##  $ Hx : num  8.1 8.1 8.1 8.1
##  $ H1 : num  1.3 1.3 1.3 1.3
```

```r
getDiameter(tree) # Hx specified beforehand being 30% of tree height
```

```
## [1] 23.19664 24.00090 27.59044 23.84774
```

With the getForm-function, one can specify trees mimicking the mean sample trees
from different inventories. This function returns the expected q03 of a given
diameter-height-class, here of dbh=30 and h=27:

```r
getForm(tree[1,], inv=c(1, 2, 3, 4))
```

```
## [1] 0.8049681 0.7996023 0.7913954 0.8107125
```

