

!*************************************************************************
      Subroutine vBDATHxDx (n,vBDATBArtNr,vD1,vH1,vD2,vH2,vH,
     1    vDx,vHx,vIErr)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vH(n)
      real vDx(n)
      real vHx(n)
      INTEGER vIErr(n)

      real FNBDATHxDx

      Do 10 i=1, n
          vHx(i) = FNBDATHxDx(vBDATBArtNr(i), vD1(i), vH1(i), vD2(i),
     1        vH2(i), vH(i), vHx(i), vDx(i), vIErr(i))
  10  continue

      END Subroutine vBDATHxDx


!*************************************************************************
      Subroutine vBDATHxDxoR (n,vBDATBArtNr,vD1,vH1,vD2,vH2,vH,
     1    vDx,vHx,vIErr)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vH(n)
      real vDx(n)
      real vHx(n)
      INTEGER vIErr(n)

      real FNBDATHxDxoR

      Do 10 i=1, n
          vHx(i) = FNBDATHxDxoR(vBDATBArtNr(i), vD1(i), vH1(i),
     1        vD2(i), vH2(i), vH(i), vHx(i), vDx(i), vIErr(i))
  10  continue

      END Subroutine vBDATHxDxoR


!***********************************************************************
       Subroutine vBiomasse (n,vBdatBart,vD13,vD2,vH2,vH,vBiom)
!***********************************************************************


       implicit none
       INTEGER i
       INTEGER n
       INTEGER vBdatBart(n)
       real vD13(n)
       real vD2(n)
       real vH2(n)
       real vH(n)
       real vBiom(n)

       real FNBiomasse

       Do 10 i=1, n
          vBiom(i) = FNBiomasse(vBdatBart(i), vD13(i), vD2(i),
     1        vH2(i), vH(i))
  10   continue

      END subroutine vBiomasse

!***********************************************************************
      Subroutine vBDAT20(n, vBDATBArtNr, vD1, vH1, vD2, vH2, vH, vHx,
     1    vHkz, vSkz, vAz, vHsh, vZsh, vZab, vSokz, vSkl,
     2    vVol, vLDSort, vBhd, vIfeh, vFixLngDef,
     3    vNMaxFixLng, vFixLng, vNFixLng)
!***********************************************************************


      implicit none
      INTEGER n
      INTEGER i, j, ijSkl, ijVol, ijLDSort, ijFLD, ijFL
      INTEGER vSkl(n*6)
      real vVol(n*7)
      real vLDSort(n*20)
      real vBhd(n)
      INTEGER vIfeh(n)
      real vFixLngDef(n*4)
      real vFixLng(n*180)
      INTEGER vNFixLng(n)

      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vH(n)
      real vHx(n)
      INTEGER vHkz(n)
      INTEGER vSkz(n)
      real vAz(n)
      real vHsh(n)
      real vZsh(n)
      real vZab(n)
      INTEGER vSokz(n)
      INTEGER aSkl(1:6)
      DATA  aSkl /6*0/ !AUS
      real aVol(1:7)
      DATA  aVol  /7*0/ !AUS
      real aLDSort(1:20)
      DATA  aLDSort  /20*0/ !AUS
      real aBhd
      DATA  aBhd /0/ !AUS
      INTEGER aIfeh
      DATA  aIfeh /0/ !AUS
      real aFixLngDef(1:4)
      DATA  aFixLngDef /10, 5, 10, 1/
      INTEGER vNMaxFixLng(n)
      real aFixLng(1:180)
      DATA  aFixLng /180*0/ !AUS
      INTEGER aNFixLng
      DATA  aNFixLng /0/ !AUS

!=======================================================================
      ijSkl=0
      ijVol=0
      ijLDSort=0
      ijFLD=0
      ijFL=0

      Do 10 i=1, n
          Do 30, j=1, 6
              aSkl(j) = 0
  30      continue

          Do 31, j=1, 20
              aLDSort(j) = 0
  31      continue

          Do 32, j=1, 180
              aFixLng(j) = 0
  32      continue

          Do 33, j=1, 7
              aVol(j) = 1
  33      continue

          Do 34, j=1, 4
              ijFLD=ijFLD+1
              aFixLngDef(j) = vFixLngDef(ijFLD)
  34      continue

          call BDAT20(vBDATBArtNr(i), vD1(i), vH1(i), vD2(i), vH2(i),
     1        vH(i), vHx(i), vHkz(i), vSkz(i), vAz(i),
     2        vHsh(i), vZsh(i), vZab(i), vSokz(i),
     3        aSkl, aVol, aLDSort, aBhd, aIfeh, aFixLngDef,
     4        vNMaxFixLng(i), aFixLng, aNFixLng)

          Do j=1, 6
              ijSkl = ijSkl+1
              vSkl(ijSkl) = aSkl(j)
          end do
          Do j=1, 7
              ijVol = ijVol+1
              vVol(ijVol) = aVol(j)
          end do
          Do j=1, 20
              ijLDSort = ijLDSort+1
              vLDSort(ijLDSort) = aLDSort(j)
          end do
          vBhd(i) = aBhd
          vIfeh(i) = aIfeh
          Do j=1, 180
              ijFL = ijFL+1
              vFixLng(ijFL) = aFixLng(j)
          end do
          vNFixLng(i)=aNFixLng

  10  continue

      end subroutine vBDAT20

!*************************************************************************
      Subroutine vBDATRINDE2Hx (n, vBDATBArtNr,vD1,vH1,vD2,vH2,vHges,
     1    vHx,vIErr,vRinde2Hx)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vHges(n)
      real vHx(n)
      INTEGER vIErr(n)
      real vRinde2Hx(n)

      real xFNBDATRinde2Hx

      Do 10 i=1, n
          vRinde2Hx(i) = xFNBDATRinde2Hx(vBDATBArtNr(i), vD1(i),
     1        vH1(i), vD2(i),vH2(i), vHges(i), vHx(i),
     2        vIErr(i), vRinde2Hx(i))
  10  continue

      END Subroutine vBDATRINDE2Hx

!*************************************************************************
      Subroutine vBDATDmRHx (n, vBDATBArtNr,vD1,vH1,vD2,vH2,vHges,
     1    vHx,vIErr,vDmRHx)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vHges(n)
      real vHx(n)
      INTEGER vIErr(n)
      real vDmRHx(n)

      real xFNBDATDmRHx

      Do 10 i=1, n
          vDmRHx(i) = xFNBDATDmRHx(vBDATBArtNr(i), vD1(i),
     1        vH1(i), vD2(i),vH2(i), vHges(i), vHx(i),
     2        vIErr(i), vDmRHx(i))
  10  continue

      END Subroutine vBDATDmRHx

!*************************************************************************
      Subroutine vBDATDoRHx (n, vBDATBArtNr,vD1,vH1,vD2,vH2,vHges,
     1    vHx,vIErr,vDoRHx)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vHges(n)
      real vHx(n)
      INTEGER vIErr(n)
      real vDoRHx(n)

      real xFNBDATDoRHx

      Do 10 i=1, n
          vDoRHx(i) = xFNBDATDoRHx(vBDATBArtNr(i), vD1(i),
     1        vH1(i), vD2(i),vH2(i), vHges(i), vHx(i),
     2        vIErr(i), vDoRHx(i))
  10  continue

      END Subroutine vBDATDoRHx

!*************************************************************************
      Subroutine vBDATVOLABmR (n, vBDATBArtNr,vD1,vH1,vD2,vH2,vHges,
     1    vA,vB,vSekLng,vIErr,vVolABmR)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vHges(n)
      real vA(n)
      real vB(n)
      real vSekLng(n)
      INTEGER vIErr(n)
      real vVolABmR(n)

      real xFNBDATVolABmR

      Do 10 i=1, n
          vVolABmR(i) = xFNBDATVolABmR(vBDATBArtNr(i),vD1(i),vH1(i),
     1        vD2(i), vH2(i), vHges(i), vA(i),
     2        vB(i), vSekLng(i), vIErr(i),
     3        vVolABmR(i))
  10  continue

      END Subroutine vBDATVOLABmR

!*************************************************************************
      Subroutine vBDATVOLABoR (n, vBDATBArtNr,vD1,vH1,vD2,vH2,vHges,
     1    vA,vB,vSekLng,vIErr,vVolABoR)
!*************************************************************************

      implicit none
      INTEGER i
      INTEGER n
      INTEGER vBDATBArtNr(n)
      real vD1(n)
      real vH1(n)
      real vD2(n)
      real vH2(n)
      real vHges(n)
      real vA(n)
      real vB(n)
      real vSekLng(n)
      INTEGER vIErr(n)
      real vVolABoR(n)

      real xFNBDATVolABoR

      Do 10 i=1, n
          vVolABoR(i) = xFNBDATVolABoR(vBDATBArtNr(i),vD1(i),vH1(i),
     1        vD2(i), vH2(i), vHges(i), vA(i),
     2        vB(i), vSekLng(i), vIErr(i),
     3        vVolABoR(i))
  10  continue

      END Subroutine vBDATVOLABoR

!***********************************************************************
       Subroutine vBDATFormTarif(n,vTarif,vBDATBArtNr,vD,vH,vMwQ03BWI)
!***********************************************************************


       implicit none
       INTEGER i
       INTEGER n
       INTEGER vTarif(n)
       INTEGER vBDATBArtNr(n)
       real vD(n)
       real vH(n)
       real vMwQ03BWI(n)

       real FNBDATFormTarif

       Do 10 i=1, n
          vMwQ03BWI(i) = FNBDATFormTarif(vTarif(i), vBDATBArtNr(i),
     1        vD(i), vH(i), vMwQ03BWI(i))
  10   continue

      END subroutine vBDATFormTarif


!***********************************************************************
      Subroutine test (n, a, b, c, d, e, f)
!***********************************************************************


      implicit none
      INTEGER i
      INTEGER n
      INTEGER a(n)
      real b(n)
      real c(n)
      real d(n)
      real e(n)
      real f(n)

      Do 10 i=1, n
            call Biomasse(a(i),b(i),c(i),d(i),e(i),f(i))
  10  continue

      END subroutine test
