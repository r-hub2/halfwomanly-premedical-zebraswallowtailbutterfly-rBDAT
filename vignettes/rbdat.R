## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(rBDAT)
tree <- buildTree(tree = list(spp=1, D1=30, H=27))
str(tree)

## ----pressure, fig.dim=c(7, 4), echo=TRUE-------------------------------------
plot(tree)

## -----------------------------------------------------------------------------
getSpeciesCode()

## -----------------------------------------------------------------------------
getSpeciesCode(1)
getSpeciesCode("Fi")
getSpeciesCode("NS") # english abbreviation of Norway spruce
getSpeciesCode(1, "scientific")
getSpeciesCode(c(1:36), "short")

## -----------------------------------------------------------------------------
getDiameter(tree, Hx = 1.3) # return dbh

## -----------------------------------------------------------------------------
getDiameter(tree, Hx = 1.3, bark = FALSE) # return d.u.b. at 1.3m

## -----------------------------------------------------------------------------
getDiameter(tree, Hx = c(1:10))

## -----------------------------------------------------------------------------
tree2 <- buildTree(list(spp=1, D1=c(30, 35), H=27))
getDiameter(tree2, Hx = c(1, 2))

## -----------------------------------------------------------------------------
tree2 <- buildTree(list(spp=1, D1=c(30, 35), H=27, Hx=c(1, 2)))
tree2
getDiameter(tree2)

## -----------------------------------------------------------------------------
dub <- getDiameter(tree, Hx=1.3, bark = FALSE)
dob <- getDiameter(tree, Hx=1.3, bark = TRUE)
dbt <- getBark(tree, Hx=1.3)
dub + dbt == dob

## -----------------------------------------------------------------------------
getBark(tree2, Hx = 1:5)

## -----------------------------------------------------------------------------
getHeight(tree, Dx=30) # height of diameter over bark
getHeight(tree, Dx=30, bark = FALSE) # height of diameter under bark == 30

## -----------------------------------------------------------------------------
getHeight(tree2, Dx=c(30, 20, 10)) # returns value in meters
tree2$Dx = c(30, 20)
getHeight(tree2)

## -----------------------------------------------------------------------------
getVolume(tree) # get coarse wood, which is the same as next line
getVolume(tree, AB = list(A=0, B=7), iAB=c("H", "Dob"), bark = T) 

## -----------------------------------------------------------------------------
# volume including bark between height 1m and 2m
getVolume(tree2, AB=list(A=1, B=2)) 
# volume excluding bark between 30 and 7cm in diameter over bark
getVolume(tree2, AB=list(A=30, B=7), iAB="dob", bark = F) 

## -----------------------------------------------------------------------------
# again: coarse wood volume
getVolume(tree) 
# identical
getVolume(tree, AB=list(A=0.27, B=7, sl=2), iAB=c("H", "Dob"), bark=F) 
# using sl=0.1, that is section length for volume calculation set to be 0.1m
getVolume(tree, AB=list(A=0.27, B=7, sl=0.1), iAB=c("H", "Dob"), bark=F) 

## -----------------------------------------------------------------------------
getVolume(tree, AB=list(A=0:9, B=1:10))
getVolume(tree = tree2, AB=list(A=0:9, B=1:10))

## -----------------------------------------------------------------------------
getBiomass(tree)
getBiomass(list(spp=1, D1=30, H=27, D2=c(23, 24, 25))) 

## -----------------------------------------------------------------------------
getAssortment(tree2) # using standard assortment parameters

## -----------------------------------------------------------------------------
getAssortment(tree2, value = "raw") # usually of little interest...

## -----------------------------------------------------------------------------
getAssortment(tree, sort=list(Az=15)) # minimum diameter o.b. for assortments
getAssortment(tree, sort=list(Az=15, Hsh=10)) # Hsh= max. height of sawlog quality

## -----------------------------------------------------------------------------
## fixN = number of roundwood pieces
## fixL = length of roundwood pieces in m
## fixA = absolute length addition (good for measure) in cm
## fixR = relative length addition in %
## fixZ = minimum cutting diameter for this assortment (Z=Zopf) in cm
getAssortment(tree, sort=list(Az=15, fixN=2, fixL=4, fixA=10, fixZ=20)) 

## -----------------------------------------------------------------------------
getAssortment(tree, sort=list(lX=1.4)) # remove 1.4m from stump upwards

## ---- fig.dim=c(7, 4), echo=TRUE----------------------------------------------
assort <- getAssortment(tree, sort=list(Az=15, fixN=2, fixL=4, fixA=10, fixZ=20))
plot(tree, assort=assort)

## -----------------------------------------------------------------------------
getAssortment(tree, sort = list(Az=c(10, 7)))
getAssortment(tree2, sort = list(Az=c(10, 7)))

## -----------------------------------------------------------------------------
tree <- buildTree(tree = list(spp=1, D1=30, H=27, H2=c(30, 50, 99, 0), 
                              D2=c(0, 0, 0, -0.8), Hx=0.3*27))
str(tree)
getDiameter(tree) # Hx specified beforehand being 30% of tree height

## -----------------------------------------------------------------------------
getForm(tree[1,], inv=c(1, 2, 3, 4))

