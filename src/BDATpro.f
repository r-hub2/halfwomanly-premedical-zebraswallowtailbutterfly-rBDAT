
!     ##################################################################
      SUBROUTINE BDAT20(BDATBArtNr, D1, H1, D2, H2, H,
     1    Hx, Hkz, Skz, Az, Hsh, Zsh, Zab, Sokz,
     2    Skl, Vol, LDSort, Bhd, Ifeh,
     3    FixLngDef,NMaxFixLng,FixLng,NFixLng)
!     ##################################################################

!  Aenderung <<22.09.03>> :--------------------------------------------------------------

       parameter (NParFixLngDef=4, NParFixLng=6, MMaxFixLng=30)
       parameter (MFixLng=NParFixLng*MMaxFixLng)
       parameter (StammFussPrz =1)

!    implicit logical (a-z)

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H

!  Aenderung <<03.11.05>> :--------------------------------------------------------------
       REAL   Hx
       INTEGER Hkz
       INTEGER Skz
       REAL   Az
       REAL   Hsh
       REAL   Zsh
       REAL   Zab
       INTEGER Sokz
       INTEGER Skl(1:6)
       REAL   Vol(1:7)
       REAL  LDSort(1:20) ! Ergaenzung Christian Vonderach 24.07.2018 / for Output
       REAL   Bhd
       INTEGER Ifeh

       REAL   FixLngDef(1:NParFixLngDef)
       INTEGER NMaxFixLng
       REAL   FixLng(1:NParFixLng*MMaxFixLng)
       INTEGER NFixLng

! ----------------------------------------------------------------------------------------

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wH

       REAL   wHx
       INTEGER wHkz
       INTEGER wSkz
       REAL   wAz
       REAL   wHsh
       REAL   wZsh,wwZsh

!  Aenderung <<03.11.05>> :--------------------------------------------------------------

       REAL   wZab
       INTEGER wSokz

       INTEGER wSkl(1:6)
       REAL   wVol(1:7)
       DATA  wSkl/6*0/, wVol/7*0/

       REAL   wBhd
       INTEGER wIfeh

       INTEGER wNMaxFixLng
       REAL   wFixLng(1:MMaxFixLng, 1:NParFixLng)

! ----------------------------------------------------------------------------------------

       INTEGER i,j,ij, IErr
       REAL   Hges, H0FixLng, HMaxFixLng
       REAL   wH0, wL0
       REAL   DMoR,DZoR,wDx,VoloR,SuVoloR

       REAL  FixLngZ
       REAL  FixLngM
       REAL  FixLngZugPrz
       REAL  FixLngZugCm
       DATA FixLngZ/7/, FixLngM/2/, FixLngZugPrz/0/, FixLngZugCm/0/
       REAL  Pi
       DATA Pi/3.14159E0/
! ---------------------------------------------------------------------------------------

       REAL   wHAz
       DATA    wHAz    /0/ !   Hoehe Aufarbeitungszopf
       REAL   wDHGrz
       DATA    wDHGrz    /7/ ! DerholzGrenze in cm
       REAL   wHDHGrz
       DATA    wHDHGrz    /0/ ! Hoehe der DerbHolzGrenze (aus)
       REAL   wSekLng
       DATA    wSekLng    /2/ !   SektionsL?ngen VolBerechnung

! ---------------------------------------------------------------------------------------

!  Aenderung <<22.09.03>> :--------------------------------------------------------------

       REAL   A,B,SekLng,VolABmR,VolDHmR
       REAL   wH0FixLng

! ---------------------------------------------------------------------------------------

! ...<<18.09.03>> : Aenderung :.......................................................

       REAL   HStockEnde
       REAL   HStHAnfang
       REAL   LngStH
       REAL   HStHLzEnde
       REAL   HBDATGes
       REAL   glLSort(1:5), glDSort(1:5) ! christian vonderach, 23.07.2018
       REAL  tmp ! christian vonderach 24.07.2018

       COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
     1     , HBDATGes
       COMMON /glLDSort/ glLSort, glDSort ! christian vonderach, 23.07.2018

! ...<<26.10.05>> : Aenderung :.......................................................

       wHAz = H
       if (Az.gt.0) then
        wHAz = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,
     2  wHAz,Az,IErr)
       end if

       wHDHGrz = H

       if (wDHGrz .gt. 0) then
        wHDHGrz = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,
     2  wHDHGrz,wDHGrz,IErr)
       else
        wDHGrz=7
        wHDHGrz = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,
     2  wHDHGrz,wDHGrz,IErr)
       end if

! ...<<03.11.05>> : Aenderung :.......................................................

       If (Hsh.gt.0) then
        wHsh=MIN(Hsh,H,wHAz,wHDHGrz)
       else
        wHsh=MIN(H,wHAz,wHDHGrz)
       end if
	   
	   ! cv 05.06.2025
	   If((Hx+0.01*H).gt.MIN(wHDHGrz, wHAz)) then ! für den Fall falscher X-Holz Laenge
	    Hx = MIN(wHDHGrz, wHAz) - 0.01*H
		IF(Hx.gt.3) Hx=int(Hx*10)/10
	   end if

       wZsh= xFNBDATDoRHx (BDATBArtNr,D1, H1, D2, H2,
     1        H, wHsh, IErr, wZsh)

!  Aenderung <<05.11.05>> :--------------------------------------------------------------
       wZsh=xFNBDATDxFoRu(wZsh)

!  Kontrolle DmoR
!  Hx=xFNBDATHxDxoRFoRu(BDATBArtNr,D1,H1,D2,H2,H,
! 2        Hx,wZsh,IErr)

       if(Zsh.gt.0) then
        wwZsh=MAX(Zsh,wZsh)
       else
        wwZsh=wZsh
       end if


!   Zsh=MAX(Az,Zsh,wDHGrz)

! ...<<18.09.03>> : Aenderung :.......................................................

       wBDATBArtNr=BDATBArtNr
       wD1=D1
       wH1=H1
       wD2=D2
       wH2=H2
       wH=H
       wHx=Hx
       wHkz=Hkz
       wSkz=Skz
       wAz=Az
       wHsh=Hsh
       wZsh=wwZsh
       wZab=Zab
       wSokz=Sokz

       FixLngZ=FixLngDef(1)
       FixLngM=FixLngDef(2)
       FixLngZugCm=FixLngDef(3)
       FixLngZugPrz=FixLngDef(4)

       if (FixLngM<2) FixLngM=2

!  if (FixLngZugPrz*0.01*FixLngM*0.01>FixLngZugCm)
! 1  FixLngZugCm=FixLngZugPrz*0.01*FixLngM*0.01

! ...<<10.06.05> : Aenderung :------------------------------------------------------------

       if (FixLngZugPrz*FixLngM>FixLngZugCm)
     1  FixLngZugCm=FixLngZugPrz*FixLngM

       do i=1,6
        wSkl(i)=0
       end do
       do i=1,7
        wVol(i)=-1
       end do
! ...Initialisierung glDSort(1:5), glLSort(1:5) / christian Vonderach 23.07.2018
       do i=1,5
        glDSort(i) = 0
        glLSort(i) = 0
       end do
! ...Initialisierung LDSort(1:20) / christian vonderach, 24.07.2018
       do i=1,20
        LDSort(i)=0
       end do

       Pi = 3.14159

       wBhd=0
       wIfeh=0

! ...<<17.01.03>> : Aenderung :------------------------------------------------------------

! ...StammHoehe Hges = H(Hkz) :........................................................

       if (Hkz.eq.1) then
        Hges=H+2
       else if (Hkz.eq.2) then
        if (30.lt.D1) then
         Hges=30+(D1-30)*0.3
        else
         Hges=D1
        end if
        if (H.gt.Hges-3) Hges=H+4
       else
        Hges=H
       end if

       HBDATGes=Hges
       HStockEnde=Hges*StammFussPrz*0.01
       HStockEnde=MIN(HStockEnde,H)

!  BDAT 1.0 - X-Holz:------------------------------------------------------------------

! ...<<17.01.03>> : Aenderung :------------------------------------------------------------

       call xBDATD2H2Trans (wBDATBArtNr,wD1,wH1,wD2,wH2,Hges)
!  call xBDATD2H2Trans (wBDATBArtNr,wD1,wH1,wD2,wH2,wH)

!  write (*,*)
!   write (*,*) " BDAT2.0 --> BDAT (X-Holz):"
!   write (*,*)

       Call BDAT(wBDATBArtNr,wD1,wH1,wD2,wH2,wH,wHx,wHkz,
     1    wSkz,wAz,wHsh,wZsh,wZab,wSokz,
     2    wSkl(1),wVol(1),wBhd,wIfeh)

!  X-HolzVolumen (Efm oR) :.............................................................

       Vol(2) = wVol(2)
       Skl(1) = wSkl(1)
       Skl(2) = wSkl(2)
! ...X-Holz am Stammfuß, Christian Vonderach, 24.01.2019
       LDSort(1) = 0.01*wH ! Fussposition im Stamm
       LDSort(2) = Hx ! Laenge des Sortiments [m]: Eingabe
       if(Hx.gt.0.001) then
         LDSort(3) = xFNBDATDoRHx(BDATBArtNr, D1, H1, D2, H2, H,
     1     LDSort(1) + LDSort(2)/2, iErr, tmp) ! Mittendurchmesser [cm]
         LDSort(3) = Rund(LDSort(3))
         LDSort(4) = xFNBDATDoRHx(BDATBArtNr, D1, H1, D2, H2, H,
     1     LDSort(1) + LDSort(2), iErr, tmp) ! Zopfdurchmesser
         LDSort(4) = Rund(LDSort(4))
       end if
! ...ende x-holz Ergaenzung

       BHD  = wBHD

       Ifeh = wIfeh

       SuVoloR = Vol(2)

!  Aenderung <<16.08.02>> :--------------------------------------------------------------

!  DerbHolzVolumen (Vfm mR) :...........................................................

!  Vol(1) = wVol(1)
!  vol(1)  = xFNBDATVolDHmR(wBDATBArtNr,wD1,wH1,wD2,wH2,wH,
!       wDHGrz,wHDHGrz,wSekLng,wIfeh,wVolDHmR)

!  Aenderung <<22.09.03>> :--------------------------------------------------------------

       wHDHGrz = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,Hges,
     2 wHDHGrz,wDHGrz,IErr)
       wHDHGrz = MIN(wHDHGrz,H)

       A=0
       B=wHDHGrz
       SekLng=2

       VolDHmR = xFNBDATVolABmR(BDATBArtNr,D1,H1,D2,H2,Hges,
     2       A,B,SekLng,IErr,VolABmR)

!  Aenderung <<22.09.03>> :--------------------------------------------------------------


!   vol(1)  = xFNBDATVolDHmR(wBDATBArtNr,wD1,wH1,wD2,wH2,Hges,
! 2       wDHGrz,wHDHGrz,wSekLng,wIfeh,wVolDHmR)

       vol(1) = VolDHmR

!  Aenderung <<22.09.03>> :--------------------------------------------------------------

       if((0<Ifeh).and.(Ifeh<5)) return

!  BDAT 2.0 - FixL?ngen:---------------------------------------------------------------

       wBDATBArtNr=BDATBArtNr
       wD1=D1
       wH1=H1
       wD2=D2
       wH2=H2
       wH=H
       wHx=Hx
       wHkz=Hkz
       wSkz=Skz
       wAz=Az
       wHsh=Hsh

!  Aenderung <<05.11.05>> :--------------------------------------------------------------
       wZsh=wwZsh
       wZab=Zab
       wSokz=Sokz

       do i=1,6
        wSkl(i)=0
       end do
       do i=1,7
        wVol(i)=-1
       end do

       wBhd=0
       wIfeh=0

       do i=1,MMaxFixLng
        do j=1,NParFixLng
         wFixLng(i,j)=0
        end do
       end do

!  Aenderung <<22.09.03>> :--------------------------------------------------------------

!  ...Aufarbeitungsgrenze HFixLngMax  = H( Hges, Skz) ! Hges bei BruchHoehen falsch:....

       if  (Skz.eq.0) then
              HMaxFixLng=H
       elseif (Skz.eq.1) then

! ...<<10.05.03>> : Aenderung :............................................................

        if (BDATBArtNr.gt.14) then
              HMaxFixLng = H*0.7
        else           ! SKZ Nadelholz
              Skz=0
              HMaxFixLng=H
        end if

       else if (Skz.eq.2) then
              HMaxFixLng = 5
       else if (Skz.eq.3) then
              HMaxFixLng = 0
       else if (Skz.eq.4) then
              HMaxFixLng = H
       else
              HMaxFixLng = 0
       end if

!  Aenderung <<22.09.03>> :--------------------------------------------------------------

!  ...Fixl?ngenaushaltung nur im Stammholzbereich :...............................

       if ((wHsh>0).and.(HMaxFixLng>wHsh)) then
        HMaxFixLng=wHsh
       end if

       wNMaxFixLng=NMaxFixLng

       if (wNMaxFixLng<0) wNMaxFixLng=0
       if (MMaxFixLng<wNMaxFixLng) wNMaxFixLng=MMaxFixLng

       NFixLng=0
       H0FixLng = 0
       wH0FixLng = Hges*StammFussPrz*0.01 + Hx

       Do
       if ((wNMaxFixLng<=NFixLng).or.(HMaxFixLng<=wH0FixLng))
     1 goto 100

        wHx = wH0FixLng
     1   + FixLngM
     2   + FixLngZugCm*0.01


!   wHx = wH0FixLng
! 1   + FixLngM
! 2   + FixLngZugCm
!

       if (HMaxFixLng<wHx)
     1 goto 100

        wHx = wH0FixLng+FixLngM

!   ...Zopfdurchmesser o.R. (stat. Erwartungswert gerundet) pr?fen :................

        DZoR=xFNBDATDoRHx(wBDATBArtNr,wD1,wH1,wD2,wH2,Hges,
     1      wHx, iErr,DZoR)

        DZoR= Rund(DZoR) ! christian vonderach 24.01.2019: Werte>20: -0.75

       if (DZoR<FixLngZ)
     1 goto 100

!  ...Fixl?ngenaushaltung nur im Stammholzbereich :...............................

       if ((wZsh>0).and.(DZoR<wZsh))
     1 goto 100

        NFixLng = NFixLng + 1

!   ...MittenDurchmesser o.R. Fixl?nge (stat. Erwartungswert gerundet) :............

        wHx = wH0FixLng+FixLngM*0.5

        DMoR=xFNBDATDoRHx(wBDATBArtNr,wD1,wH1,wD2,wH2,Hges,
     1      wHx, iErr,DMoR)

        if (DMoR<20) then
         DMoR=DMoR-0.5
        else
         DMoR=DMoR-0.75
        end if

!   ...Volumen o.R. Fixl?nge (stat. Erwartungswert gerundet) :......................

        wDx = DMoR * 0.01
        VoloR = Pi * 0.25 * wDx * wDx * FixLngM

!   --------------------------------------------------------------------------------
        wFixLng(NFixLng,1) = NFixLng
        wFixLng(NFixLng,2) = wH0FixLng
        wFixLng(NFixLng,3) = FixLngM
        wFixLng(NFixLng,4) = DMoR
        wFixLng(NFixLng,5) = DZoR
        wFixLng(NFixLng,6) = VoloR
!   --------------------------------------------------------------------------------

        wH0FixLng = wH0FixLng
     1     + FixLngM
     2     + FixLngZugCm*0.01
        H0FixLng = H0FixLng
     1     + FixLngM
     2     + FixLngZugCm*0.01

        SuVoloR    = SuVoloR + VoloR

       end do
100    continue

       ij=0
       do i=1,MMaxFixLng
        do j=1,NParFixLng
         ij=ij+1
         FixLng(ij)=wFixLng(i,j)
        end do
       end do

!  BDAT 1.0 - Sortierung StammHolz IndustrieHolz :-------------------------------------

       wBDATBArtNr=BDATBArtNr
       wD1=D1
       wH1=H1
       wD2=D2
       wH2=H2
       wH=H
       wHx=Hx + H0FixLng
       wHkz=Hkz
       wSkz=Skz
       wAz=Az
       wHsh=Hsh
       wZsh=wwZsh
       wZab=Zab
       wSokz=Sokz

       do i=1,6
        wSkl(i)=0
       end do
       do i=1,7
        wVol(1:7)=0
       end do

       wBhd=0
       wIfeh=0

       call xBDATD2H2Trans (wBDATBArtNr,wD1,wH1,wD2,wH2,Hges)

       wIFeh=0

!  Aenderung <<17.01.03>> :...schwaches Stangenholz Du < 10 cm --------------------------

!   ...MittenDurchmesser o.R. Fixl?nge (stat. Erwartungswert gerundet) :............

       if (wD1 < 10) then

        wH0 = Hges*StammFussPrz*0.01 + wHx

        if(wH0>H) wH0=H

!   <<18.09.03>> : Aenderung :.......................................................

!   FixLng(NParFixLng*MMaxFixLng)= HStHLzEnde

!   Aenderung <<22.09.03>> :--------------------------------------------------------------

!    wHDHGrz = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,Hges,
!  2  wHDHGrz,wDHGrz,IErr)

!   Aenderung <<01.04.04>> :----------------------------------------------------

        wHDHGrz = xFNBDATHxDx(BDATBArtNr,wD1,wH1,wD2,wH2,Hges,
     2  wHDHGrz,wDHGrz,IErr)

        wHDHGrz = MIN(wHDHGrz,H)

!   Aenderung <<22.09.03>> :----------------------------------------------------

!    <<20.01.03>> : Aenderung :..................................................

        wL0 = (wHDHGrz-wH0)

        if(wL0>wHDHGrz) wL0=wHDHGrz
        if(wL0<0) wL0=0

        wHx = wH0+wL0*0.5

        if (wHx>H) wHx=H

        DMoR=xFNBDATDoRHx(wBDATBArtNr,wD1,wH1,wD2,wH2,Hges,
     1       wHx, iErr,DMoR)

        if (DMoR<20) then
         DMoR=DMoR-0.5
        else
         DMoR=DMoR-0.75
        end if

        wDx = DMoR * 0.01
        VoloR = Pi * 0.25 * wDx * wDx * wL0

!     Ergaenzung Christian Vonderach 22.01.2019
!     MDM, Laenge + zopf IndustrieHolz ins output schreiben
      LDSort(13) = LDSort(1) + LDSort(2) + H0FixLng ! Fußpunkt Ind
      LDSort(14) = wL0 ! Länge Ind
      LDSort(15) = DMoR ! MDM Ind
      LDSort(16) = xFNBDATDoRHx(wBDATBArtNr,wD1,wH1,wD2,wH2,Hges,
     1       wH0+wL0, iErr,DMoR) ! Zopf Ind
	  if(LDSort(16).lt.20) then ! forstlich abrunden
	   LDSort(16) = LDSort(16) - 0.5
	  else
	   LDSort(16) = LDSort(16) - 0.75
	  end if
	  
      if(LDSort(2).le.0.0001) then
       LDSort(1)=0 ! kein X-holz=>Nullsetzen
      else
       Vol(2)= Pi * 0.25 * LDSort(3)/100 * LDSort(3)/100 * LDSort(2)
       SuVoloR = SuVoloR + Vol(2)
      end if

        do i=3,6
         Skl(i)=0
        end do

        do i=3,7
         Vol(i)=0
        end do

        if (VoloR>(Vol(1)-SuVoloR)) VoloR = Vol(1)-SuVoloR

        if (wSokz > 0) then

         SuVoloR = SuVoloR+VoloR
         vol(5)=VoloR

        else
         vol(5)=0
        end if

        vol(7) = Vol(1) - SuVoloR

        RETURN

       end if

!  Aenderung <<17.01.03>> :--------------------------------------------------------

!  ------------------------------------------------
!   write (*,*)
!  write (*,*) " VOL(5): ", VOL(5)
!   write (*,*) " BDAT2.0 --> BDAT (Sortierung):"
!   write (*,*)
!  ------------------------------------------------


       Call BDAT(wBDATBArtNr,wD1,wH1,wD2,wH2,wH,wHx,wHkz,
     1    wSkz,wAz,wHsh,wZsh,wZab,wSokz,
     2    wSkl(1),wVol(1),wBhd,wIfeh)

!  ------------------------------------------------
!   write (*,*) " BDAT --> BDAT2.0 (Sortierung):"
!  write (*,*) " VOL(5): ", wVOL(5)
!  ------------------------------------------------

       do i=3,6
        Skl(i)=wSkl(i)
        Vol(i)=wVol(i)
        SuVoloR = SuVoloR + Vol(i)
       end do

! Aenderung  <<20.09.02>> :----------------------------------------------------------------

       if ((SuVoloR>0).and.(SuVoloR<Vol(1))) then
        Vol(7) = Vol(1)-SuVoloR
        if (Vol(7)<0) Vol(7)=0
       else
        Vol(7) = 0
       end if

       Ifeh = wIfeh

! ...Zuweisung COMMON /glLDSort/ zu LDSort(1:20) christian vonderach 24.07.2018
! ...aktualisierung X-holz cv 06.06.2025
      LDSort(2) = glLSort(1) ! Laenge des Sortiments [m]
      LDSort(3) = glDSort(1) ! Mittendurchmesser [cm]
! ...Stammholz, zzgl Fixlaengen am Stammfuss
      LDSort(5) = LDSort(1) + LDSort(2) + H0FixLng! Fussposition im Stamm
      LDSort(6) = glLSort(2) ! Laenge des Sortiments [m]
      LDSort(7) = glDSort(2) ! Mittendurchmesser [cm]
! ...Abschnitt
      LDSort(9) = LDSort(5) + LDSort(6)*1.01 ! Fussposition im Stamm
      LDSort(10) = glLSort(3) ! Laenge des Sortiments [m]
      LDSort(11) = glDSort(3) ! Mittendurchmesser [cm]
! ...Industrieholz
      LDSort(13) = LDSort(9) + LDSort(10)*1.01 ! Fussposition im Stamm
      LDSort(14) = glLSort(4) ! Laenge des Sortiments [m]
      LDSort(15) = glDSort(4) ! Mittendurchmesser [cm]
! ...nvD-Holz (Obs: Ih ohne Zugabe)
      LDSort(17) = LDSort(13) + LDSort(14) ! Fussposition im Stamm
      LDSort(18) = glLSort(5) ! Laenge des Sortiments [m]
      LDSort(19) = glDSort(5) ! Mittendurchmesser [cm]
      do i=4, 20, 4 ! Berechnung Zopfdurchmesser für Sortimente
       if (LDSort(i-2).gt.0.0001) then ! Abfang Rundungsfehler 10.12.2018 cv
        tmp = xFNBDATDoRHx(BDATBArtNr, D1, H1, D2, H2, H,
     1     LDSort(i-3) + LDSort(i-2), iErr, tmp)
        tmp = Rund(tmp)
       else
        tmp=0 ! Zopfdurchmesser 0, falls Laenge=0
        LDSort(i-3)=0 ! Fussposition auf 0 setzen, falls Laenge=0
        LDSort(i-2)=0 ! Laenge auf exakt 0 setzen
        LDSort(i-1)=0 ! Mittendurchmesser auf 0 setzen, falls Laenge=0
       end if
       LDSort(i) = tmp ! Zopfdurchmesser [cm] der Sortimente
      end do
! ...ende ergaenzung christian vonderach 24.07.2018...

! ...<<18.09.03>> : Aenderung :.......................................................

! FixLng(NParFixLng*MMaxFixLng)= HStHLzEnde

! ...<<18.09.03>> : Aenderung :.......................................................

      END SUBROUTINE BDAT20


! ****************************************************************************************
      subroutine BDATVolDHmR(BDATBArtNr, D1, H1, D2, H2, Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHmR)
! ****************************************************************************************

! Derbholzvolumen zu gegebenem DerbholzGrenzDurchmesser :.................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   DHGrz     ! DerholzGrenze in cm
       REAL   HDHGrz     ! Hoehe der DerbHolzGrenze (aus)
       REAL   SekLng
       INTEGER IErr
       REAL   VolDHmR     ! VolumenDerbHolz  (aus)

! -----------------------------------------------------------------------------------------
       VolDHmR = xFNBDATVolDHmR(BDATBArtNr, D1, H1, D2, H2,Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHmR)
! -----------------------------------------------------------------------------------------

      end subroutine BDATVolDHmR


! ****************************************************************************************
      REAL function FNBDATVolDHmR(BDATBArtNr, D1, H1, D2, H2, Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHmR)
! ****************************************************************************************

! Derbholzvolumen zu gegebenem DerbholzGrenzDurchmesser :.................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   DHGrz     ! DerholzGrenze in cm
       REAL   HDHGrz     ! Hoehe der DerbHolzGrenze (aus)
       REAL   SekLng
       INTEGER IErr
       REAL   VolDHmR     ! VolumenDerbHolz  (aus)

! ----------------------------------------------------------------------------------------
       VolDHmR = xFNBDATVolDHmR(BDATBArtNr, D1, H1, D2, H2,Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHmR)
! ----------------------------------------------------------------------------------------

       FNBDATVolDHmR=VolDHmR

      end function FNBDATVolDHmR


! ****************************************************************************************
      REAL function xFNBDATVolDHmR(BDATBArtNr, D1, H1, D2, H2, Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHmR)
! ****************************************************************************************

! Derbholzvolumen zu gegebenem DerbholzGrenzDurchmesser :.................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   DHGrz     ! DerholzGrenze in cm
       REAL   HDHGrz     ! Hoehe der DerbHolzGrenze (aus)
       REAL   SekLng
       INTEGER IErr
       REAL   VolDHmR     ! VolumenDerbHolz  (aus)

       REAL   Dx,Hx,A,B,VolABmR

! -----------------------------------------------------------------------------------------

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

       Dx=DHGrz
       HDHGrz= xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,Hges,Hx,Dx,IErr)

       A=0
       B=HDHGrz
       SekLng=2

       VolDHmR = xFNBDATVolABmR(BDATBArtNr,D1,H1,D2,H2,Hges,
     2       A,B,SekLng,IErr,VolABmR)

       xFNBDATVolDHmR=VolDHmR

      end function xFNBDATVolDHmR

! *** insertion by cv 10.07.2018
! ****************************************************************************************
      subroutine BDATVolDHoR(BDATBArtNr, D1, H1, D2, H2, Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHoR)
! ****************************************************************************************

! Derbholzvolumen zu gegebenem DerbholzGrenzDurchmesser :.................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   DHGrz     ! DerholzGrenze in cm
       REAL   HDHGrz     ! Hoehe der DerbHolzGrenze (aus)
       REAL   SekLng
       INTEGER IErr
       REAL   VolDHoR     ! VolumenDerbHolz ohnr Rinde (aus)

! -----------------------------------------------------------------------------------------
       VolDHoR = xFNBDATVolDHoR(BDATBArtNr, D1, H1, D2, H2,Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHoR)
! -----------------------------------------------------------------------------------------

      end subroutine BDATVolDHoR


! ****************************************************************************************
      REAL function FNBDATVolDHoR(BDATBArtNr, D1, H1, D2, H2, Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHoR)
! ****************************************************************************************

! Derbholzvolumen zu gegebenem DerbholzGrenzDurchmesser :.................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   DHGrz     ! DerholzGrenze in cm
       REAL   HDHGrz     ! Hoehe der DerbHolzGrenze (aus)
       REAL   SekLng
       INTEGER IErr
       REAL   VolDHoR     ! VolumenDerbHolz  (aus)

! ----------------------------------------------------------------------------------------
       VolDHoR = xFNBDATVolDHoR(BDATBArtNr, D1, H1, D2, H2,Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHoR)
! ----------------------------------------------------------------------------------------

       FNBDATVolDHoR=VolDHoR

      end function FNBDATVolDHoR


! ****************************************************************************************
      REAL function xFNBDATVolDHoR(BDATBArtNr, D1, H1, D2, H2, Hges,
     2      DHGrz, HDHGrz, SekLng, IErr,VolDHoR)
! ****************************************************************************************

! Derbholzvolumen zu gegebenem DerbholzGrenzDurchmesser :.................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   DHGrz     ! DerholzGrenze in cm
       REAL   HDHGrz     ! Hoehe der DerbHolzGrenze (aus)
       REAL   SekLng
       INTEGER IErr
       REAL   VolDHoR     ! VolumenDerbHolz  (aus)

       REAL   Dx,Hx,A,B,VolABoR

! -----------------------------------------------------------------------------------------

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

       Dx=DHGrz
       HDHGrz= xFNBDATHxDxoR(BDATBArtNr,D1,H1,D2,H2,Hges,Hx,Dx,IErr)

       A=0
       B=HDHGrz
       SekLng=2

       VolDHoR = xFNBDATVolABoR(BDATBArtNr,D1,H1,D2,H2,Hges,
     2       A,B,SekLng,IErr,VolABoR)

       xFNBDATVolDHoR=VolDHoR

      end function xFNBDATVolDHoR

! ** end insertion by cv 10.07.2018

! ****************************************************************************************
      subroutine BDATVolABmR(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wA, wB, wSekLng, wIErr,
     2      wVolABmR)
! ****************************************************************************************

! ... BaumVolumen [m?] für den Abschnitt StammHoehe A - B durch sektionsweise Kubierung
!  in Form von Walzen der L?nge SekLng.

!  INTEGER BDATBArtNr - BDATBaumArtSchl?ssel
!  REAL   D1    - Durchmesser (cm) in Hoehe wH1 (wH1 = 0 : H1 = 1.3 m)
!  REAL   H1    - Hoehe (m) Durchmesser 1 (0: wH1 = 1.3)
!  REAL   D2    - Durchmesser(cm) in Hoehe wH2 (vgl. BDAT)
!  REAL   H2    - Hoehe (m) Durchmesser 2 (0: wH2 = 7.0)
!  REAL   Hges   - StammHoehe (m)
!  REAL   A    - Abschnittsanfang (m) für die Voluminierung
!  REAL   B    - Abschnittsende (m) für die Voluminierung
!  REAL   SekLng   - Sektionsl?nge (m)
!  INTEGER IErr   - Fehlerindikator
!  REAL   VolABmR  - Abschnittsvolumen (m?) in Rinde

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wA
       REAL   wB
       REAL   wSekLng
       INTEGER wIErr
       REAL   wVolABmR

! ----------------------------------------------------------------------------------------
       wVolABmR = xFNBDATVolABmR(wBDATBArtNr,wD1,wH1,wD2,wH2,wHges,
     2        wA, wB, wSekLng,wIErr, wVolABmR)
! ----------------------------------------------------------------------------------------

      End Subroutine BDATVolABmR


! ****************************************************************************************
      REAL function FNBDATVolABmR(wBDATBArtNr, wD1, wH1, wD2,
     1        wH2, wHges, wA, wB, wSekLng,
     2        wIErr, wVolABmR)
! ****************************************************************************************

! ... BaumVolumen [m?] für den Abschnitt StammHoehe A - B durch sektionsweise Kubierung in
!  Form von Walzen der L?nge SekLng.

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wA
       REAL   wB
       REAL   wSekLng
       INTEGER wIErr
       REAL   wVolABmR

! ----------------------------------------------------------------------------------------
       FNBDATVolABmR = xFNBDATVolABmR(wBDATBArtNr,wD1,wH1,wD2,wH2,
     2      wHges,wA, wB, wSekLng,wIErr, wVolABmR)
! ----------------------------------------------------------------------------------------

      End Function FNBDATVolABmR


! ****************************************************************************************
      REAL function xFNBDATVolABmR(wBDATBArtNr, wD1, wH1, wD2,
     1        wH2, wHges, wA, wB, wSekLng,
     2        wIErr, wVolABmR)
! ****************************************************************************************

! ... BaumVolumen [m?] für den Abschnitt StammHoehe A - B durch sektionsweise Kubierung in
!  Form von Walzen der L?nge SekLng.

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wA
       REAL   wB
       REAL   wSekLng
       INTEGER wIErr
       REAL   wVolABmR

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   A
       REAL   B
       REAL   SekLng
       INTEGER IErr

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   volum(7)
       DATA    volum /7*0.0/

       REAL   Pi
       DATA    Pi   /3.14159E0/

! ........................................................................................

       REAL   H
       REAL   Hx
       REAL   HxM
       REAL   VolAx
       REAL   VolBx
       REAL   Ax
       REAL   Bx
       REAL   ABx

! ----------------------------------------------------------------------------------------

       BDATBArtNr = wBDATBArtNr

       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges
       H=Hges

! ----------------------------------------------------------------------------------------

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

       Call BDAT(BDATBArtNr, D1, H1, D2, H2, Hges, Stxu,
     1    Hkz,Sk, Azop, sthh, Zost, Zoab, sok,
     2    klasse(1), volum(1), BHDz,ifeh)

! ----------------------------------------------------------------------------------------

       A = wA
       B = wB
       SekLng = wSekLng
       IErr = wIErr

! ----------------------------------------------------------------------------------------

       VolBx = 0
       VolAx = 0

       Ax = A

       If (Ax > H) Then
        Ax = H
       End If

       Bx = B
       If (Bx <= Ax) Then
        Bx = Ax
        xFNBDATVolABmR = 0
        wVolABmR = xFNBDATVolABmR
        return
       End If

       Hx = 0

       ABx = SekLng

       If (Ax > 0) Then

        LoopExit = 0
        Hx = 0

        Do While (LoopExit .eq. 0)

         If (Hx + ABx <= Ax) Then

          HxM = Hx + ABx * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ----------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
!     ----------------------------------

          Dx = tmp

          VolAx = VolAx + Pi * 0.25 * Dx * 0.01 * Dx *
     1     0.01 * ABx
          Hx = Hx + ABx

         Else

          HxM = (Ax + Hx) * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ----------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
!     ----------------------------------

          Dx = tmp

          VolAx = VolAx + Pi * 0.25 * Dx * 0.01 * Dx * 0.01
     1     * (Ax - Hx)
          LoopExit = 1

         End If

        end do

       Else
        VolAx = 0
       End If

       If (Bx > 0) Then

        LoopExit = 0
        Hx = 0

        Do While (LoopExit .eq. 0)

         If (Hx + ABx <= Bx) Then

          HxM = Hx + ABx * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ----------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
!     ----------------------------------

          Dx = tmp

          VolBx = VolBx + Pi * 0.25 * Dx * 0.01 * Dx * 0.01
     1     * ABx
          Hx = Hx + ABx

         Else

          HxM = (Bx + Hx) * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ----------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
!     ----------------------------------

          Dx = tmp

          VolBx = VolBx + Pi * 0.25 * Dx * 0.01 * Dx * 0.01
     1     * (Bx - Hx)

          LoopExit = 1

         End If

        end do
       Else
        VolBx = 0
       End If

       If (VolBx > VolAx) Then
        xFNBDATVolABmR = VolBx - VolAx
       Else
        xFNBDATVolABmR = 0
       End If

        wVolABmR = xFNBDATVolABmR

      End Function xFNBDATVolABmR


! ****************************************************************************************
      subroutine BDATVolABoR(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wA, wB, wSekLng, wIErr,
     2      wVolABoR)
! ****************************************************************************************

! ... BaumVolumen [m?] für den Abschnitt StammHoehe A - B durch sektionsweise Kubierung
!  in Form von Walzen der L?nge SekLng.

!  INTEGER BDATBArtNr - BDATBaumArtSchl?ssel
!  REAL   D1    - Durchmesser (cm) in Hoehe wH1 (wH1 = 0 : H1 = 1.3 m)
!  REAL   H1    - Hoehe (m) Durchmesser 1 (0: wH1 = 1.3)
!  REAL   D2    - Durchmesser(cm) in Hoehe wH2 (vgl. BDAT)
!  REAL   H2    - Hoehe (m) Durchmesser 2 (0: wH2 = 7.0)
!  REAL   Hges   - StammHoehe (m)
!  REAL   A    - Abschnittsanfang (m) für die Voluminierung
!  REAL   B    - Abschnittsende (m) für die Voluminierung
!  REAL   SekLng   - Sektionsl?nge (m)
!  INTEGER IErr   - Fehlerindikator
!  REAL   VolABoR  - Abschnittsvolumen (m?) in Rinde
!
! ----------------------------------------------------------------------------------------
       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wA
       REAL   wB
       REAL   wSekLng
       INTEGER wIErr
       REAL   wVolABoR

! ----------------------------------------------------------------------------------------
          wVolABoR=xFNBDATVolABoR(wBDATBArtNr, wD1, wH1, wD2,
     1        wH2, wHges, wA, wB, wSekLng,
     2        wIErr, wVolABoR)
! ----------------------------------------------------------------------------------------

      End Subroutine BDATVolABoR


! ****************************************************************************************
      REAL function FNBDATVolABoR(wBDATBArtNr, wD1, wH1, wD2,
     1        wH2, wHges, wA, wB, wSekLng,
     2        wIErr, wVolABoR)
! ****************************************************************************************

! ... BaumVolumen [m?] für den Abschnitt StammHoehe A - B durch sektionsweise Kubierung in
!  Form von Walzen der L?nge SekLng.

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wA
       REAL   wB
       REAL   wSekLng
       INTEGER wIErr
       REAL   wVolABoR


! ----------------------------------------------------------------------------------------
       FNBDATVolABoR = xFNBDATVolABoR(wBDATBArtNr, wD1, wH1, wD2,
     1        wH2, wHges, wA, wB, wSekLng,
     2        wIErr, wVolABoR)
! ----------------------------------------------------------------------------------------


      End Function FNBDATVolABoR

! ****************************************************************************************
      REAL function xFNBDATVolABoR(wBDATBArtNr, wD1, wH1, wD2,
     1        wH2, wHges, wA, wB, wSekLng,
     2        wIErr, wVolABoR)
! ****************************************************************************************

! ... BaumVolumen [m?] für den Abschnitt StammHoehe A - B durch sektionsweise Kubierung in
!  Form von Walzen der L?nge SekLng.

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wA
       REAL   wB
       REAL   wSekLng
       INTEGER wIErr
       REAL   wVolABoR

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   A
       REAL   B
       REAL   SekLng
       INTEGER IErr

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   volum(7)
       DATA    volum /7*0.0/

       REAL   Pi
       DATA    Pi   /3.14159E0/

! ........................................................................................

       REAL   H
       REAL   Hx
       REAL   HxM
       REAL   VolAx
       REAL   VolBx
       REAL   Ax
       REAL   Bx
       REAL   ABx

! ----------------------------------------------------------------------------------------

       BDATBArtNr = wBDATBArtNr

       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges
       H=Hges

! ----------------------------------------------------------------------------------------

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

       Call BDAT(BDATBArtNr, D1, H1, D2, H2, Hges, Stxu,
     1    Hkz,Sk, Azop, sthh, Zost, Zoab, sok,
     2    klasse(1), volum(1), BHDz,ifeh)

! ----------------------------------------------------------------------------------------

       A = wA
       B = wB
       SekLng = wSekLng
       IErr = wIErr

! ----------------------------------------------------------------------------------------

       VolBx = 0
       VolAx = 0

       Ax = A

       If (Ax > H) Then
        Ax = H
       End If

       Bx = B
       If (Bx <= Ax) Then
        Bx = Ax
        xFNBDATVolABoR = 0
        wVolABoR = xFNBDATVolABoR
        return
       End If

       Hx = 0

       ABx = SekLng

       If (Ax > 0) Then

        LoopExit = 0
        Hx = 0

        Do While (LoopExit .eq. 0)

         If (Hx + ABx <= Ax) Then

          HxM = Hx + ABx * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ----------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
          Call RINDE (1 - HxM / Hges, tmp, R2, 0, 0)
!     -----------------------------------------

          If (tmp < 0) Then
           tmp = 0.0
          End If

          Dx = tmp

          VolAx = VolAx + Pi * 0.25 * Dx * 0.01 * Dx *
     1     0.01 * ABx
          Hx = Hx + ABx

         Else

          HxM = (Ax + Hx) * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ----------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
          Call RINDE (1 - HxM / Hges, tmp, R2, 0, 0)
!     -----------------------------------------

          If (tmp < 0) Then
           tmp = 0.0
          End If

          Dx = tmp

          VolAx = VolAx + Pi * 0.25 * Dx * 0.01 * Dx * 0.01
     1     * (Ax - Hx)
          LoopExit = 1

         End If

        end do

       Else
        VolAx = 0
       End If

       If (Bx > 0) Then

        LoopExit = 0
        Hx = 0

        Do While (LoopExit .eq. 0)

         If (Hx + ABx <= Bx) Then

          HxM = Hx + ABx * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ------------------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
          Call RINDE (1 - HxM / Hges, tmp, R2, 0, 0)
!     ------------------------------------------

          If (tmp < 0) Then
           tmp = 0.0
          End If

          Dx = tmp

          VolBx = VolBx + Pi * 0.25 * Dx * 0.01 * Dx * 0.01
     1     * ABx
          Hx = Hx + ABx

         Else

          HxM = (Bx + Hx) * 0.5

          If (HxM > Hges) Then
           HxM = Hges
          EndIf

!     ------------------------------------------
          Call KUWERT(1 - HxM / Hges, tmp)
          Call RINDE (1 - HxM / Hges, tmp, R2, 0, 0)
!     ------------------------------------------

          If (tmp < 0) Then
           tmp = 0.0
          End If

          Dx = tmp

          VolBx = VolBx + Pi * 0.25 * Dx * 0.01 * Dx * 0.01
     1     * (Bx - Hx)

          LoopExit = 1

         End If

        end do
       Else
        VolBx = 0
       End If

       If (VolBx > VolAx) Then
        xFNBDATVolABoR = VolBx - VolAx
       Else
        xFNBDATVolABoR = 0
       End If

        wVolABoR = xFNBDATVolABoR

      End Function xFNBDATVolABoR


! ****************************************************************************************
      subroutine BDATDmRHx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wHx, wIErr, wDmRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDmRHx

! ----------------------------------------------------------------------------------------
       wDmRHx = xFNBDATDmRHx(wBDATBArtNr, wD1, wH1, wD2, wH2
     1       ,wHges,wHx, wIErr, wDmRHx)
! ----------------------------------------------------------------------------------------

      End subroutine BDATDmRHx


! ****************************************************************************************
      REAL function FNBDATDmRHx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1        wHges, wHx, wIErr, wDmRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDmRHx

! ----------------------------------------------------------------------------------------
       FNBDATDmRHx = xFNBDATDmRHx(wBDATBArtNr, wD1, wH1, wD2, wH2
     1        ,wHges,wHx, wIErr, wDmRHx)
! ----------------------------------------------------------------------------------------

      End Function FNBDATDmRHx


! ****************************************************************************************
      REAL function xFNBDATDmRHx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wHx, wIErr, wDmRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDmRHx

! ----------------------------------------------------------------------------------------
       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   Hx
       INTEGER IErr
       REAL   DmRHx

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   Pi
       DATA    Pi   /3.14159E0/
       REAL   volum(7)
       DATA    volum /7*0.0/

       REAL   H

! ----------------------------------------------------------------------------------------

       BDATBArtNr = wBDATBArtNr

       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges
       Hx = wHx
       IErr = wIErr
       DmRHx = wDmRHx

       H=Hges

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

! ....Liefert den Durchmesser an der Stelle Hx im Schaft als Funktion der Baumart
!  Durchmesser(Hoehe 1/2) (D1,H1) (D2,H2)und der GesamtHoehe H

! ----------------------------------------------------------------------------------------
       Call BDAT(BDATBArtNr, D1, H1, D2, H2, Hges, Stxu, Hkz,
     1 Sk, Azop, sthh, Zost, Zoab, sok, klasse(1), volum(1), BHDz,
     2 ifeh)
! ----------------------------------------------------------------------------------------

       wIErr = ifeh

       If (Hx .gt. Hges) Then
        Hx = Hges
       EndIf

!  ----------------------------------
       Call KUWERT(1 - Hx / Hges, tmp)
!  ----------------------------------

       If (tmp .lt. 0) Then
        tmp = 0
       endif

       wDmRHx = tmp
       xFNBDATDmRHx = wDmRHx

      End Function xFNBDATDmRHx

! ****************************************************************************************
      REAL function yFNBDATDmRHx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wHx, wIErr, wDmRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDmRHx

! ----------------------------------------------------------------------------------------
       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   Hx
       INTEGER IErr
       REAL   DmRHx

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   Pi
       DATA    Pi   /3.14159E0/
       REAL   volum(7)
       DATA    volum /7*0.0/

       REAL   H

! ----------------------------------------------------------------------------------------

       BDATBArtNr = wBDATBArtNr

       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges
       Hx = wHx
       IErr = wIErr
       DmRHx = wDmRHx

       H=Hges

!  call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

! ....Liefert den Durchmesser an der Stelle Hx im Schaft als Funktion der Baumart
!  Durchmesser(Hoehe 1/2) (D1,H1) (D2,H2)und der GesamtHoehe H

! ----------------------------------------------------------------------------------------
       Call BDAT(BDATBArtNr, D1, H1, D2, H2, Hges, Stxu, Hkz,
     1 Sk, Azop, sthh, Zost, Zoab, sok, klasse(1), volum(1), BHDz,
     2 ifeh)
! ----------------------------------------------------------------------------------------

       wIErr = ifeh

       If (Hx .gt. Hges) Then
        Hx = Hges
       EndIf

!  ----------------------------------
       Call KUWERT(1 - Hx / Hges, tmp)
!  ----------------------------------

       If (tmp .lt. 0) Then
        tmp = 0
       endif

       wDmRHx = tmp
       yFNBDATDmRHx = wDmRHx

      End Function yFNBDATDmRHx


! ****************************************************************************************
      Subroutine BDATDoRHx (wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wHx, wIErr, wDoRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDoRHx

! ----------------------------------------------------------------------------------------
       wDoRHx=xFNBDATDoRHx (wBDATBArtNr, wD1, wH1, wD2, wH2,
     1      wHges, wHx, wIErr, wDoRHx)
! ----------------------------------------------------------------------------------------

      End Subroutine BDATDoRHx


! ****************************************************************************************
      REAL function FNBDATDoRHx (wBDATBArtNr, wD1, wH1, wD2, wH2,
     1        wHges, wHx, wIErr, wDoRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDoRHx

! ----------------------------------------------------------------------------------------
       FNBDATDoRHx=xFNBDATDoRHx (wBDATBArtNr, wD1, wH1, wD2, wH2,
     1        wHges, wHx, wIErr, wDoRHx)
! ----------------------------------------------------------------------------------------

      End Function FNBDATDoRHx


! ****************************************************************************************
      REAL function xFNBDATDoRHx (wBDATBArtNr, wD1, wH1, wD2, wH2,
     1        wHges, wHx, wIErr, wDoRHx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wDoRHx

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   Hx
       INTEGER IErr
       REAL   DoRHx

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   Pi
       DATA    Pi   /3.14159E0/
       REAL   volum(7)
       DATA    volum /7*0.0/

       REAL   kw
       DATA    kw   /0.0/

       REAL   H

! ----------------------------------------------------------------------------------------

       BDATBArtNr = wBDATBArtNr
       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges
       Hx = wHx
       IErr = wIErr
       DoRHx = wDoRHx

       H = Hges

! ----------------------------------------------------------------------------------------

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)


       Call BDAT(BDATBArtNr, D1, H1, D2, H2, Hges, Stxu,
     1    Hkz,Sk, Azop, sthh, Zost, Zoab, sok,
     2    klasse(1), volum(1), BHDz,ifeh)

! ----------------------------------------------------------------------------------------

!  -----------------------------------------
       Call KUWERT(1 - Hx / Hges, kw)
       Call RINDE (1 - Hx / Hges, Kw, R2, 0, 0)
!  -----------------------------------------

       If (Kw < 0) Then
        Kw = 0.0
       End If

       wDoRHx = Kw
       xFNBDATDoRHx = Kw

      End Function xFNBDATDoRHx


! ****************************************************************************************
      subroutine BDATRinde2Hx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1       wHges, wHx, wIErr, wRinde2Hx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wRinde2Hx

! ----------------------------------------------------------------------------------------
       wRinde2Hx = xFNBDATRinde2Hx(wBDATBArtNr, wD1, wH1, wD2,
     1       wH2,wHges, wHx, wIErr, wRinde2Hx)
! ----------------------------------------------------------------------------------------

      End subroutine BDATRinde2Hx


! ****************************************************************************************
      REAL function FNBDATRinde2Hx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1        wHges, wHx, wIErr, wRinde2Hx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wRinde2Hx

! ----------------------------------------------------------------------------------------
       FNBDATRinde2Hx = xFNBDATRinde2Hx(wBDATBArtNr, wD1, wH1, wD2,
     1       wH2,wHges, wHx, wIErr, wRinde2Hx)
! ----------------------------------------------------------------------------------------

      End Function FNBDATRinde2Hx


! ****************************************************************************************
      REAL function xFNBDATRinde2Hx(wBDATBArtNr, wD1, wH1, wD2, wH2,
     1        wHges, wHx, wIErr, wRinde2Hx)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
       REAL   wHx
       INTEGER wIErr
       REAL   wRinde2Hx

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges
       REAL   Hx
       INTEGER IErr
       REAL   Rinde2Hx

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   Pi
       DATA    Pi   /3.14159E0/
       REAL   volum(7)
       DATA    volum /7*0.0/

       REAL   kw
       DATA    kw   /0.0/

       REAL   H

! ----------------------------------------------------------------------------------------

       BDATBArtNr = wBDATBArtNr

       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges
       Hx = wHx
       IErr = wIErr
       Rinde2Hx = wRinde2Hx

       H = Hges

       call xBDATD2H2Trans (BDATBArtNr,D1,H1,D2,H2,Hges)

! ----------------------------------------------------------------------------------------
       Call BDAT(BDATBArtNr, D1, H1, D2, H2, Hges, Stxu, Hkz,
     1 Sk, Azop, sthh, Zost, Zoab, sok, klasse(1), volum(1), BHDz
     2 , ifeh)
! ----------------------------------------------------------------------------------------

       IErr = ifeh

       If (Hx .gt. Hges) Then
        Hx = Hges
       EndIf

!  -----------------------------------------
       Call KUWERT(1 - Hx / Hges, kw)
       Call RINDE (1 - Hx / Hges, kw, tmp, 0, 0)
!  -----------------------------------------

       wRinde2Hx = tmp
       xFNBDATRinde2Hx = tmp


      End Function xFNBDATRinde2Hx


! ****************************************************************************************
      subroutine BDATD2H2Trans (wBDATBArtNr, wD1, wH1, wD2, wH2,wHges)
! ****************************************************************************************

! ----------------------------------------------------------------------------------------
       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges
! ----------------------------------------------------------------------------------------
       call xBDATD2H2Trans (wBDATBArtNr, wD1, wH1, wD2, wH2,wHges)
! ----------------------------------------------------------------------------------------

      End subroutine BDATD2H2Trans


! ****************************************************************************************
      subroutine xBDATD2H2Trans (wBDATBArtNr, wD1, wH1, wD2, wH2,wHges)
! ****************************************************************************************

       INTEGER wBDATBArtNr
       REAL   wD1
       REAL   wH1
       REAL   wD2
       REAL   wH2
       REAL   wHges

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   Hges

       INTEGER Hkz
       DATA  Hkz   /0/
       INTEGER sok
       DATA  sok   /0/
       INTEGER Sk
       DATA  Sk   /0/
       INTEGER ifeh
       DATA  ifeh  /0/
       INTEGER iz
       DATA  iz   /0/
       INTEGER klasse(6)
       DATA  klasse /6*0/

       REAL   Du
       DATA    Du   /0.0/
       REAL   Hu
       DATA    Hu   /0.0/
       REAL   ddo
       DATA    ddo   /0.0/
       REAL   Ho
       DATA    Ho   /0.0/
       REAL   BHDz
       DATA    BHDz  /0.0/
       REAL   tmp
       DATA    tmp   /0.0/
       REAL   Stxu
       DATA    Stxu  /0.0/
       REAL   Azop
       DATA    Azop  /0.0/
       REAL   sthh
       DATA    sthh  /0.0/
       REAL   Zost
       DATA    Zost  /0.0/
       REAL   Zoab
       DATA    Zoab  /0.0/
       REAL   Pi
       DATA    Pi   /3.14159E0/
       REAL   volum(7)
       DATA    volum /7*0.0/
       REAL   kw
       DATA    kw   /0.0/

! ----------------------------------------------------------------------------------------

       REAL   Hx, Dx, Q03Pct,MwQ03BWI,StDevQ03BWI,MwQ03BWIPct
       INTEGER iERR
       DATA  iERR  /0/
       REAL   H,D2t,H2t,D2u,D2o
! ----------------------------------------------------------------------------------------

       if (wH1<=0) then
        wH1=1.3
       end if

       BDATBArtNr = wBDATBArtNr

       D1 = wD1
       H1 = wH1
       D2 = wD2
       H2 = wH2
       Hges = wHges

! ----------------------------------------------------------------------------------------

       if  (D2 > 0) then

!   D2 (ein) ~ Durchmesserwert

        if (H2>0) then

!    D2 (aus) = D2 (ein)
!    H2 (aus) = H2 (ein)

        else

!    H2 (aus) = 7 m (Standardwert BWI)

         H2 = 7

        end if

!  ...<<12.11.02>> : Aenderung :........................................................

        H=Hges
        D2t=-0.400
        H2t=0.3*H
        Hx=H2

        D2u=yFNBDATDmRHx(BDATBArtNr,D1,H1,D2t,H2t,H,Hx,iERR,Dx)

        H=Hges
        D2t=-0.9500
        H2t=0.3*H
        Hx=H2

        D2o=yFNBDATDmRHx(BDATBArtNr,D1,H1,D2t,H2t,H,Hx,iERR,Dx)

        H=Hges
        D2t=D2
        H2t=H2
        Hx=H2

        Dx=yFNBDATDmRHx(BDATBArtNr,D1,H1,D2t,H2t,H,Hx,iERR,Dx)
        D2t=Dx

        if (ABS(D2o-D2)>ABS(D2u-D2)) then
         if (ABS(D2t-D2)>ABS(D2u-D2)) then
          D2=-0.4000
          H2=0.3*H
         end if
        else
         if (ABS(D2t-D2)>ABS(D2o-D2)) then
          D2=-0.9500
          H2=0.3*H
         end if
        end if

       else if ((-1 < D2).and.(D2 < 0)) then

!   Formigkeit: q0.3:= D0.30/D095

!   D2 (ein) ~ - q0.3;

!   D2 (aus) = D2 (ein)
!   H2 (aus) = 0.3*Hges

        H2 = 0.3*Hges

       else if (D2 <= -1)     then

!    Formigkeit: BWI-?quivalenz (MEDIAN q0.3)

!    D2 (aus) = EST MED [q0.3 |BHD,Hges;BWI]
!    H2 (aus) = 0.3*Hges

         Q03Pct=0.50

         call xBDATMwQ03BWI (BDATBArtNr,D1,Hges,Q03Pct,
     1      MwQ03BWI,StDevQ03BWI,MwQ03BWIPct)

         D2 = -MwQ03BWIPct
         H2 = 0.3*Hges

       else

!   D2 (ein) = 0 : D2 ~ ?ber <<H2>> festlegen

        if (H2 <= 0) then

!    Formigkeit: MassenTafel-?quivalenz (BDAT 1.0)

         D2=0
         H2=0

        else if ((0<H2).and.(H2<100)) then

!    Formigkeit: BWI-?quivalenz (PCTL q0.3)

!    H2 (ein) ~ PercentilWert * 100
!
!    D2 (aus) = EST PCTL [q0.3 | BHD,Hges,H2(ein);BWI]
!    H2 (aus) = 0.3*Hges

         Q03Pct=wH2*0.01

         call xBDATMwQ03BWI (BDATBArtNr,D1,Hges,Q03Pct,
     1      MwQ03BWI,StDevQ03BWI,MwQ03BWIPct)

         D2 = -MwQ03BWIPct
         H2 = 0.3*Hges

        else

!    Formigkeit: BWI-?quivalenz (MEDIAN q0.3)
!
!    D2 (aus) = EST MED [q0.3 |BHD,Hges;BWI]
!    H2 (aus) = 0.3*Hges

         Q03Pct=0.50

         call xBDATMwQ03BWI (BDATBArtNr,D1,Hges,Q03Pct,
     1      MwQ03BWI,StDevQ03BWI,MwQ03BWIPct)

         D2 = -MwQ03BWIPct
         H2 = 0.3*Hges

        end if

       end if

!     AusgabeParameter :------------------------------------------------------------------

       wD1=D1
       wH1=H1

       wD2=D2
       wH2=H2

      End subroutine xBDATD2H2Trans


! ****************************************************************************************
      subroutine BDATMwQ03BWI(BDATBArtNr,D,H,Q03Pct,
     1      MwQ03BWI,StDevQ03BWI,MwQ03BWIPct)
! ****************************************************************************************

!  q0.3 - Percentile = MwQ03BWIPct = F(D,H,PctlWert) / MW / STD = F(D,H)
!
!  BDATBArtNr INT*2 EIN - BDAT BaumArt Nummer 1-36
!   D   REAL  EIN - BHD [cm]
!  H   REAL  EIN - BaumHoehe [m]
!  Q03Pct  REAL  EIN - PercentilWert (eps,1-eps) für MwQ03BWIPct Schätzung
!  MwQ03BWI      REAL  AUS - Mittelwert q0.3 (gesch?tzt)
!  StDevQ03BWI      REAL  AUS - Standardabweichung Mwq0.3 - Verteilung
!  MwQ03BWIPct      REAL  AUS - Mittlerer q0.3 - Percentilwert (gesch?tzt)

! ****************************************************************************************

! ----------------------------------------------------------------------------------------
       INTEGER BDATBArtNr
       REAL   D, H, Q03Pct
       REAL   MwQ03BWI,StDevQ03BWI, MwQ03BWIPct
! ----------------------------------------------------------------------------------------
       call xBDATMwQ03BWI(BDATBArtNr,D,H,Q03Pct,
     1      MwQ03BWI,StDevQ03BWI,MwQ03BWIPct)
! ----------------------------------------------------------------------------------------

      end subroutine BDATMwQ03BWI


! ****************************************************************************************
      subroutine xBDATMwQ03BWI(BDATBArtNr,D,H,Q03Pct,
     1      MwQ03BWI,StDevQ03BWI,MwQ03BWIPct)
! ****************************************************************************************

!  q0.3 - Percentile = MwQ03BWIPct = F(D,H,PctlWert) / MW / STD = F(D,H)
!
!  BDATBArtNr INT*2 EIN - BDAT BaumArt Nummer 1-36
!   D   REAL  EIN - BHD [cm]
!  H   REAL  EIN - BaumHoehe [m]
!  Q03Pct  REAL  EIN - PercentilWert (eps,1-eps) für MwQ03BWIPct Schätzung
!  MwQ03BWI      REAL  AUS - Mittelwert q0.3 (gesch?tzt)
!  StDevQ03BWI      REAL  AUS - Standardabweichung Mwq0.3 - Verteilung
!  MwQ03BWIPct      REAL  AUS - Mittlerer q0.3 - Percentilwert (gesch?tzt)

! ****************************************************************************************

       INTEGER BDATBArtNr
       REAL   D, H, Q03Pct
       REAL   MwQ03BWI,StDevQ03BWI, MwQ03BWIPct

       REAL   EQP(1:8,1:2,1:7)
       REAL   SQP(1:8,1:6)

       INTEGER BDATSKNrList(1:36)

       INTEGER BDATSKNr

       REAL   Q
       REAL   Q1, Q2, Q3, sQ1, sQ2, sQ3
       REAL   EQ03,StDevQ03

       REAL   EQ03uG
       DATA    EQ03uG /0.4000/
       REAL   EQ03oG
       DATA    EQ03oG /0.9800/

       REAL   a11, a12, a13, h11, h12, H13, D1
       REAL   a21, a22, a23, h21, h22, H23, D2

       REAL   Phi

       REAL   Z1, Z2

       REAL   CDFx, x

       REAL   eps
       DATA    eps /0.001/

! ----------------------------------------------------------------------------------------

       data BDATSKNrList       /
     1 1,  1,  2,  2,  4,  4,  4,  3,  5,  5,
     2 5,  1,  1,  1,  6,  6,  7,  8,  8,  8,
     3 6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
     4 6,  6,  7,  6,  6,  6      /

! -------------------------------- Stand 23.08.01 ----------------------------------------

       data (((EQP(i,j,k),k=1,7),j=1,2),i=1,8)  /
     1  20, 10, 50, 0.650, 0.875, 0.850, 0.250,
     1 55, 10, 50, 0.525, 0.860, 0.775, 0.000,
     2 20, 10, 50, 0.750, 0.950, 0.860, 0.500,
     2 70, 10, 50, 0.670, 0.875, 0.780, 0.000,
     3 20, 10, 50, 0.650, 0.925, 0.860, 0.300,
     3 70, 10, 50, 0.500, 0.825, 0.700, 0.000,
     4 20, 10, 50, 0.700, 0.800, 0.770, 0.250,
     4 50, 10, 50, 0.700, 0.790, 0.770, 0.000,
     5 20, 10, 50, 0.725, 0.950, 0.875, 0.500,
     5 60, 10, 50, 0.630, 0.875, 0.750, 0.000,
     6 20, 10, 50, 0.700, 0.900, 0.830, 0.750,
     6 60, 10, 50, 0.650, 0.870, 0.820, 0.000,
     7 20, 10, 50, 0.700, 0.850, 0.840, 0.750,
     7 60, 10, 50, 0.675, 0.840, 0.825, 0.000,
     8 20, 10, 50, 0.775, 0.850, 0.810, 1.000,
     8 60, 10, 50, 0.725, 0.800, 0.760, 0.000 /

       data ((SQP(i,j),j=1,6),i=1,8)      /
     1  0.50, 0.75, 1.00, 0.2500, 0.0700, 0.0000,
     2 0.50, 0.75, 1.00, 0.2500, 0.0800, 0.0000,
     3 0.50, 0.75, 1.00, 0.2500, 0.0600, 0.0000,
     4 0.50, 0.75, 1.00, 0.3000, 0.0550, 0.0000,
     5 0.50, 0.75, 1.00, 0.2000, 0.0600, 0.0000,
     6 0.50, 0.75, 1.00, 0.3000, 0.0900, 0.0000,
     7 0.50, 0.80, 1.00, 0.2500, 0.0700, 0.0000,
     8 0.50, 0.80, 1.00, 0.0300, 0.0300, 0.0300  /

! ----------------------------------------------------------------------------------------

       BDATSKNr = BDATSKNrList(BDATBArtNr)

       a11 = EQP(BDATSKNr, 1, 4)
       a12 = EQP(BDATSKNr, 1, 5)
       a13 = EQP(BDATSKNr, 1, 6)

       h11 = EQP(BDATSKNr, 1, 2)
       h12 = EQP(BDATSKNr, 1, 3)
       H13 = (h12 + h11) * 0.5

       D1 = EQP(BDATSKNr, 1, 1)

       a21 = EQP(BDATSKNr, 2, 4)
       a22 = EQP(BDATSKNr, 2, 5)
       a23 = EQP(BDATSKNr, 2, 6)

       h21 = EQP(BDATSKNr, 2, 2)
       h22 = EQP(BDATSKNr, 2, 3)
       H23 = (h22 + h21) * 0.5

       D2 = EQP(BDATSKNr, 2, 1)

       Phi = EQP(BDATSKNr, 1, 7)

!     ****************************************************************************************
!     * Z(H|D(i)) = MW [Q0.3| H | a(H(i,j)|D(i))]; j=1,2,3; i=1,2                            *
! * Ratkowsky, D.A. (1990) (4.3.9), S97                                                  *
!     ****************************************************************************************

       Q1 = 2 * (H - h11) / (h12 - h11)
       Z1 = a11 + (a12 - a11) * (1 - ((a12 - a13)/(a13 - a11))** Q1)
     1 / (1 - ((a12 - a13) / (a13 - a11)) ** 2)

       Q2 = 2 * (H - h21) / (h22 - h21)
       Z2 = a21 + (a22 - a21) * (1 - ((a22 - a23)/(a23 - a21)) ** Q2)
     1 / (1 - ((a22 - a23) / (a23 - a21)) ** 2)

!     ****************************************************************************************
!     * EQ0.3(D,H) =  E [Q0.3| D, Z(H|D(i)); i=1,2] * Ratkowsky, D.A. (1990) (4.3.23), S104  *
!     ****************************************************************************************

       EQ03 = Z1 * Z2 * (D2 ** Phi - D1 ** Phi)
     1 / (Z2 * (D2 ** Phi - D ** Phi) + Z1 * (D ** Phi - D1 ** Phi))

       If (EQ03 < EQ03uG) Then
        EQ03 = EQ03uG
       End If

       If (EQ03 > EQ03oG) Then
        EQ03 = EQ03oG
       End If

!  ***************
       MwQ03BWI = EQ03
!  ***************

!     ****************************************************************************************
!     * sQ0.3(D,H) =  s [Q0.3| D,sQ(Q|i)); i=1,3] * Ratkowsky, D.A. (1990) (4.3.29), S106    *
!     ****************************************************************************************


!  Call BDATStDevQ03(BDATBArtNr, MwQ03, BDATSKNr, StDevQ03)

       Q1 = SQP(BDATSKNr, 1)
       Q2 = SQP(BDATSKNr, 2)
       Q3 = SQP(BDATSKNr, 3)

       sQ1 = SQP(BDATSKNr, 4)
       sQ2 = SQP(BDATSKNr, 5)
       sQ3 = SQP(BDATSKNr, 6)

       Q = EQ03

       if (ABS(sQ3-sQ1) < eps) then
        StDevQ03 = sQ3
       else
        StDevQ03 = (Q-Q3)*(Q2-Q1)*sQ1*sQ2 + (Q-Q2)*(Q1-Q3)*sQ1*sQ3
     1  + (Q-Q1)*(Q3-Q2)*sQ2*sQ3

        StDevQ03 = StDevQ03/((Q-Q1)*(Q2-Q3)*sQ1 + (Q - Q2)*(Q3-Q1)
     1  * sQ2 + (Q-Q3)*(Q1-Q2)*sQ3)
       end if

!  ********************
       StDevQ03BWI=StDevQ03
!  ********************

       CDFx = Q03Pct

       if (CDFx < eps) then
        CDFx=0.5
       end if

       if (CDFx > 1 - eps) then
        CDFx=0.5
       end if

!  *****************************************
       Call CDFNORMInv(EQ03, StDevQ03, CDFx, x)
!  *****************************************


!  ***************
       MwQ03BWIPct = x
!  ***************

       if (MwQ03BWIPct>1) then
        MwQ03BWIPct=1
       end if

       if (MwQ03BWIPct<0) then
        MwQ03BWIPct=0
       end if

      end subroutine xBDATMwQ03BWI


! ****************************************************************************************
      subroutine BDATPctQ03BWI(BDATBArtNr,D,H,Q03,
     1      MwQ03BWI,StDevQ03BWI,PctQ03BWI)
! ****************************************************************************************

!  q0.3 - PercentilWert = PctQ03BWI = F(D,H,q0.3) / Mittelwert / Standardabweichung = F(D,H)

!  BDATBArtNr INT*2 EIN - BDAT BaumArt Nummer 1-36
!   D   REAL  EIN - BHD [cm]
!  H   REAL  EIN - BaumHoehe [m]
!  Q03     REAL  EIN - Formquotient für die Percentilbestimmung
!  MwQ03BWI      REAL  AUS - Mittelwert q0.3 (gesch?tzt)
!  StDevQ03BWI      REAL  AUS - Standardabweichung Mwq0.3 - Verteilung
!  PctQ03BWI      REAL  AUS - q0.3 - Percentilwert (gesch?tzt) zu q0.3, D, H

! ****************************************************************************************

       INTEGER BDATBArtNr
       REAL   D, H, Q03
       REAL   MwQ03BWI, StDevQ03BWI, PctQ03BWI

       REAL   EQP(1:8,1:2,1:7)
       REAL   SQP(1:8,1:6)

       INTEGER BDATSKNrList(1:36)

       INTEGER BDATSKNr

       REAL   Q
       REAL   Q1, Q2, Q3, sQ1, sQ2, sQ3

       REAL   EQ03,StDevQ03
       REAL   EQ03uG
       DATA    EQ03uG /0.4000/
       REAL   EQ03oG
       DATA    EQ03oG /0.9800/

       REAL   a11, a12, a13, h11, h12, H13, D1
       REAL   a21, a22, a23, h21, h22, H23, D2

       REAL   Phi

       REAL   Z1, Z2

       REAL   CDFx, x

       REAL   eps
       DATA    eps /0.001/

! ----------------------------------------------------------------------------------------

       data BDATSKNrList      /
     1 1,  1,  2,  2,  4,  4,  4,  3,  5,  5,
     2 5,  1,  1,  1,  6,  6,  7,  8,  8,  8,
     3 6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
     4 6,  6,  7,  6,  6,  6     /

! -------------------------------- Stand 23.08.01 ----------------------------------------

       data (((EQP(i,j,k),k=1,7),j=1,2),i=1,8)  /
     1  20, 10, 50, 0.650, 0.875, 0.850, 0.250,
     1 55, 10, 50, 0.525, 0.860, 0.775, 0.000,
     2 20, 10, 50, 0.750, 0.950, 0.860, 0.500,
     2 70, 10, 50, 0.670, 0.875, 0.780, 0.000,
     3 20, 10, 50, 0.650, 0.925, 0.860, 0.300,
     3 70, 10, 50, 0.500, 0.825, 0.700, 0.000,
     4 20, 10, 50, 0.700, 0.800, 0.770, 0.250,
     4 50, 10, 50, 0.700, 0.790, 0.770, 0.000,
     5 20, 10, 50, 0.725, 0.950, 0.875, 0.500,
     5 60, 10, 50, 0.630, 0.875, 0.750, 0.000,
     6 20, 10, 50, 0.700, 0.900, 0.830, 0.750,
     6 60, 10, 50, 0.650, 0.870, 0.820, 0.000,
     7 20, 10, 50, 0.700, 0.850, 0.840, 0.750,
     7 60, 10, 50, 0.675, 0.840, 0.825, 0.000,
     8 20, 10, 50, 0.775, 0.850, 0.810, 1.000,
     8 60, 10, 50, 0.725, 0.800, 0.760, 0.000  /

       data ((SQP(i,j),j=1,6),i=1,8)     /
     1 0.50, 0.75, 1.00, 0.2500, 0.0700, 0.0000,
     2 0.50, 0.75, 1.00, 0.2500, 0.0800, 0.0000,
     3 0.50, 0.75, 1.00, 0.2500, 0.0600, 0.0000,
     4 0.50, 0.75, 1.00, 0.3000, 0.0550, 0.0000,
     5 0.50, 0.75, 1.00, 0.2000, 0.0600, 0.0000,
     6 0.50, 0.75, 1.00, 0.3000, 0.0900, 0.0000,
     7 0.50, 0.80, 1.00, 0.2500, 0.0700, 0.0000,
     8 0.50, 0.80, 1.00, 0.0300, 0.0300, 0.0300  /

! -------------------------------- Stand 23.08.01 ----------------------------------------

       BDATSKNr = BDATSKNrList(BDATBArtNr)

       a11 = EQP(BDATSKNr, 1, 4)
       a12 = EQP(BDATSKNr, 1, 5)
       a13 = EQP(BDATSKNr, 1, 6)

       h11 = EQP(BDATSKNr, 1, 2)
       h12 = EQP(BDATSKNr, 1, 3)
       H13 = (h12 + h11) * 0.5

       D1 = EQP(BDATSKNr, 1, 1)

       a21 = EQP(BDATSKNr, 2, 4)
       a22 = EQP(BDATSKNr, 2, 5)
       a23 = EQP(BDATSKNr, 2, 6)

       h21 = EQP(BDATSKNr, 2, 2)
       h22 = EQP(BDATSKNr, 2, 3)
       H23 = (h22 + h21) * 0.5

       D2 = EQP(BDATSKNr, 2, 1)

       Phi = EQP(BDATSKNr, 1, 7)

!     ****************************************************************************************
!     * Z(H|D(i)) = MW [Q0.3| H | a(H(i,j)|D(i))]; j=1,2,3; i=1,2                            *
! * Ratkowsky, D.A. (1990) (4.3.9), S97                                                  *
!     ****************************************************************************************

       Q1 = 2 * (H - h11) / (h12 - h11)
       Z1 = a11 + (a12 - a11) * (1 - ((a12 - a13) / (a13 - a11))**Q1)
     1 / (1 - ((a12 - a13) / (a13 - a11)) ** 2)

       Q2 = 2 * (H - h21) / (h22 - h21)
       Z2 = a21 + (a22 - a21) * (1 - ((a22 - a23) / (a23 - a21))**Q2)
     1 / (1 - ((a22 - a23) / (a23 - a21)) ** 2)

!     ****************************************************************************************
!     * EQ0.3(D,H) =  E [Q0.3| D, Z(H|D(i)); i=1,2] * Ratkowsky, D.A. (1990) (4.3.23), S104  *
!     ****************************************************************************************

       EQ03 = Z1 * Z2 * (D2 ** Phi - D1 ** Phi)
     1 / (Z2 * (D2 ** Phi - D ** Phi) + Z1 * (D ** Phi - D1 ** Phi))


       If (EQ03 < EQ03uG) Then
        EQ03 = EQ03uG
       End If

       If (EQ03 > EQ03oG) Then
        EQ03 = EQ03oG
       End If


!  ***************
       MwQ03BWI = EQ03
!  ***************

!  Call BDATStDevQ03(BDATBArtNr, MwQ03, BDATSKNr, StDevQ03)

       Q1 = SQP(BDATSKNr, 1)
       Q2 = SQP(BDATSKNr, 2)
       Q3 = SQP(BDATSKNr, 3)

       sQ1 = SQP(BDATSKNr, 4)
       sQ2 = SQP(BDATSKNr, 5)
       sQ3 = SQP(BDATSKNr, 6)

       Q = EQ03

       if (ABS(sQ3-sQ1) < eps) then
        StDevQ03 = sQ3
       else
        StDevQ03 = (Q-Q3)*(Q2-Q1)*sQ1*sQ2 + (Q-Q2)*(Q1-Q3)*sQ1*sQ3
     1  + (Q-Q1)*(Q3-Q2)*sQ2*sQ3

        StDevQ03 = StDevQ03/((Q-Q1)*(Q2-Q3)*sQ1 + (Q - Q2)*(Q3-Q1)
     1  * sQ2 + (Q-Q3)*(Q1-Q2)*sQ3)
       end if

!  ********************
       StDevQ03BWI=StDevQ03
!  ********************

       if (Q03 < eps) then
        PctQ03BWI=0
       else if (Q03 > 1 - eps) then
        PctQ03BWI=1
       else

!   **************************************
        x=Q03
        Call CDFNORM(EQ03, StDevQ03, x, CDFx)
        PctQ03BWI = CDFx
!   **************************************

       end if


      end subroutine BDATPctQ03BWI


! ****************************************************************************************
      subroutine EQ03ParIni(WEQP,WSQP)
! ****************************************************************************************

       REAL  WEQP(1:8,1:2,1:7)
       REAL  WSQP(1:8,1:6)

       REAL  EQP(1:8,1:2,1:7)
       REAL  SQP(1:8,1:6)


!  COMMON / EQ03 / EQP, SQP

! -------------------------------- Stand 23.08.01 ----------------------------------------

       data (((EQP(i,j,k),k=1,7),j=1,2),i=1,8)  /
     1 20, 10, 50, 0.650, 0.875, 0.850, 0.250,
     1 55, 10, 50, 0.525, 0.860, 0.775, 0.000,
     2 20, 10, 50, 0.750, 0.950, 0.860, 0.500,
     2 70, 10, 50, 0.670, 0.875, 0.780, 0.000,
     3 20, 10, 50, 0.650, 0.925, 0.860, 0.300,
     3 70, 10, 50, 0.500, 0.825, 0.700, 0.000,
     4 20, 10, 50, 0.700, 0.800, 0.770, 0.250,
     4 50, 10, 50, 0.700, 0.790, 0.770, 0.000,
     5 20, 10, 50, 0.725, 0.950, 0.875, 0.500,
     5 60, 10, 50, 0.630, 0.875, 0.750, 0.000,
     6 20, 10, 50, 0.700, 0.900, 0.830, 0.750,
     6 60, 10, 50, 0.650, 0.870, 0.820, 0.000,
     7 20, 10, 50, 0.700, 0.850, 0.840, 0.750,
     7 60, 10, 50, 0.675, 0.840, 0.825, 0.000,
     8 20, 10, 50, 0.775, 0.850, 0.810, 1.000,
     8 60, 10, 50, 0.725, 0.800, 0.760, 0.000 /

      data ((SQP(i,j),j=1,6),i=1,8)      /
     1 0.50, 0.75, 1.00, 0.2500, 0.0700, 0.0000,
     2 0.50, 0.75, 1.00, 0.2500, 0.0800, 0.0000,
     3 0.50, 0.75, 1.00, 0.2500, 0.0600, 0.0000,
     4 0.50, 0.75, 1.00, 0.3000, 0.0550, 0.0000,
     5 0.50, 0.75, 1.00, 0.2000, 0.0600, 0.0000,
     6 0.50, 0.75, 1.00, 0.3000, 0.0900, 0.0000,
     7 0.50, 0.80, 1.00, 0.2500, 0.0700, 0.0000,
     8 0.50, 0.80, 1.00, 0.0300, 0.0300, 0.0300  /

! -------------------------------------------------------------------------------------------

       do i= 1,8,1
        do j= 1,2,1
         do k = 1,7,1
          WEQP(i,j,k)=EQP(i,j,k)
         end do
        end do
       end do

       do i= 1,8,1
        do j= 1,6,1
         WSQP(i,j)=SQP(i,j)
        end do
       end do

      end subroutine EQ03ParIni


! q03 - Fortschreibung :_____________________________________________________________________


! *******************************************************************************************
      REAL function FNBDATEstQ032(BDATBArtNr, BHD1, D71, H1, BHD2, H2,
     1       Estq032, EstD72, iErr)
! *******************************************************************************************

       INTEGER BDATBArtNr
       REAL   BHD1
       REAL   D71
       REAL   H1
       REAL   BHD2
       REAL   H2
       REAL   Estq032
       REAL   EstD72
       INTEGER iErr

       REAL   D1,HD1,D2,HD2,H,Hx
       REAL   q03Pct,MwQ03BWI,StDevQ03BWI,MwQ03BWIPct
       REAL   q031,MWq031,MWq032

       REAL   D005,D03

! --------------------------------------------------------------------------------------------

       iErr = 0

       If ((BDATBArtNr < 1).Or.(BDATBArtNr > 36)) Then
            BDATBArtNr = 1
            iErr = 1
       End If

!        '...q03 = F(BHD,D7=0,H,INV 2) :...............................................

       D1 = BHD2
       HD1 = 1.3
       D2 = 0
       HD2 = 50
       H = H2

       q03Pct = 0.5

       Call xBDATMwQ03BWI(BDATBArtNr, BHD2, H2, q03Pct, MwQ03BWI,
     1      StDevQ03BWI, MwQ03BWIPct)
       MWq032 = MwQ03BWI

       If (D71 <= 0) Then

        Estq032 = MWq032

        D1 = BHD2
        HD1 = 1.3
        D2 = -Estq032
        HD2 = 0.3 * H2
        H = H2

        Hx = 7
        Dx = 0

        Dx=xFNBDATDmRHx(BDATBArtNr, D1, HD1, D2, HD2, H, Hx, iErr
     1  , Dx)

        EstD72 = Dx

       Else

!  '...q03 = F(BHD,D7=0,H,INV 1) :..............................................

        D1 = BHD1
        HD1 = 1.3
        D2 = 0
        HD2 = 50
        H = H1

        q03Pct = 0.5

        Call xBDATMwQ03BWI(BDATBArtNr, BHD1, H1, q03Pct, MwQ03BWI,
     1       StDevQ03BWI, MwQ03BWIPct)
        MWq031 = MwQ03BWI

!  '...q03 = F(BHD,D7,H,INV 1) :................................................

        D1 = BHD1
        HD1 = 1.3
        D2 = D71
        HD2 = 7
        H = H1

        Hx = 0.05 * H
        Dx = xFNBDATDmRHx(BDATBArtNr, D1, HD1, D2, HD2, H, Hx,
     1      iErr, D005)

        Hx = 0.3 * H
        Dx = xFNBDATDmRHx(BDATBArtNr, D1, HD1, D2, HD2, H, Hx,
     1      iErr, D03)

        q031 = D03 / D005

!  '...Fortschreibung :........................................................

        Estq032 = q031 + (MWq032 - MWq031)

       End If

       D1 = BHD2
       HD1 = 1.3
       D2 = -Estq032
       HD2 = 0.3 * H2
       H = H2

       Hx = 7
       Dx = 0

       Dx = xFNBDATDmRHx(BDATBArtNr, D1, HD1, D2, HD2, H, Hx, iErr, Dx)

!  ******************************
       If (Dx > D71) Then
        EstD72 = Dx
       Else
        iErr = iErr + 1
        EstD72 = D71 * BHD2 / BHD1
       End If
!  ******************************


!  ******************************
       FNBdatEstq032 = Estq032
!  ******************************

      End Function FNBDATEstQ032

! subroutine ergaenzt für BDATHXDX, 10.06.2018 cv
! ****************************************************************************************
      Subroutine BDATHxDx (BDATBArtNr, D1, H1, D2, H2,
     1      H, Dx, Hx, IErr)
! ****************************************************************************************

       INTEGER BDATBArtNr
       REAL  D1
       REAL  H1
       REAL  D2
       REAL  H2
       REAL  H
       REAL  Hx
       REAL  Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------
       Hx=FNBDATHxDx( BDATBArtNr, D1, H1, D2, H2,
     1      H, Hx, Dx, IErr)
! ----------------------------------------------------------------------------------------

      End Subroutine BDATHxDx


! ***************************************************************************************
      REAL function FNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,
     2        Hx,Dx,IErr)
! ***************************************************************************************

!     ... Hoehe (Hx) zu gegebenem Durchmesser (Dx) :...........................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H
       REAL   Hx
       REAL   Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------
       FNBDATHxDx = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,Hx,Dx,IErr)
! ----------------------------------------------------------------------------------------

      END Function FNBDATHxDx


! ***************************************************************************************
      REAL function xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,
     2        Hx,Dx,IErr)
! ***************************************************************************************

!     ... Hoehe (Hx) zu gegebenem Durchmesser (Dx) :...........................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H
       REAL   Hx
       REAL   Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------

       INTEGER NSFktNr
       DATA  NSFktNr   /1/
       REAL   NSFktPar
       DATA    NSFktPar  /0/

       REAL   A
       DATA    A    /0/
       REAL   B
       REAL   NSFktAbsErr
       DATA    NSFktAbsErr  /1e-3/
       REAL   NSFktXAbsErr
       DATA    NSFktXAbsErr /0/
       REAL   NSFktXRelErr
       DATA    NSFktXRelErr /1e-3/

       INTEGER MIt
       DATA  MIt    /100/
       REAL   x1
       REAL   x2

! ----------------------------------------------------------------------------------------

       A = 0
       B = H

       call BDATNullStellenSuche
     1  (BDATBArtNr, D1, H1, D2, H2, H,
     2   NSFktNr, Dx,A, B,
     3   NSFktAbsErr, NSFktXAbsErr, NSFktXRelErr,
     4   MIt, x1, x2, Hx, IErr
     5  )

       xFNBDATHxDx = Hx

      END Function xFNBDATHxDx


! ...<<07.03.03>> : Aenderung :............................................................

! subroutine ergaenzt für BDATHXDX, 14.06.2018 cv
! ****************************************************************************************
      Subroutine BDATHxDxoR (BDATBArtNr, D1, H1, D2, H2,
     1      H, Dx, Hx, IErr)
! ****************************************************************************************

       INTEGER BDATBArtNr
       REAL  D1
       REAL  H1
       REAL  D2
       REAL  H2
       REAL  H
       REAL  Hx
       REAL  Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------
       Hx=FNBDATHxDxoR( BDATBArtNr, D1, H1, D2, H2,
     1      H, Hx, Dx, IErr)
! ----------------------------------------------------------------------------------------

      End Subroutine BDATHxDxoR

! ***************************************************************************************
      REAL function FNBDATHxDxoR(BDATBArtNr,D1,H1,D2,H2,H,
     2        Hx,Dx,IErr)
! ***************************************************************************************

!     ... Hoehe (Hx) zu gegebenem Durchmesser oR (Dx) :.......................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H
       REAL   Hx
       REAL   Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------
       FNBDATHxDxoR = xFNBDATHxDxoR(BDATBArtNr,D1,H1,D2,H2,H,Hx,Dx,
     1 IErr)
! ----------------------------------------------------------------------------------------

      END Function FNBDATHxDxoR


! ***************************************************************************************
      REAL function xFNBDATHxDxoR(BDATBArtNr,D1,H1,D2,H2,H,
     2        Hx,Dx,IErr)
! ***************************************************************************************

!     ... Hoehe (Hx) zu gegebenem Durchmesser oR (Dx) :.......................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H
       REAL   Hx
       REAL   Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------

       INTEGER NSFktNr
       DATA  NSFktNr   /2/
       REAL   NSFktPar
       DATA    NSFktPar  /0/

       REAL   A
       DATA    A    /0/
       REAL   B
       REAL   NSFktAbsErr
       DATA    NSFktAbsErr  /1e-3/
       REAL   NSFktXAbsErr
       DATA    NSFktXAbsErr /0/
       REAL   NSFktXRelErr
       DATA    NSFktXRelErr /1e-3/

       INTEGER MIt
       DATA  MIt    /100/
       REAL   x1
       REAL   x2

! ----------------------------------------------------------------------------------------

       A = 0
       B = H

       call BDATNullStellenSuche
     1  (BDATBArtNr, D1, H1, D2, H2, H,
     2   NSFktNr, Dx,A, B,
     3   NSFktAbsErr, NSFktXAbsErr, NSFktXRelErr,
     4   MIt, x1, x2, Hx, IErr
     5  )

       xFNBDATHxDxoR = Hx

      END Function xFNBDATHxDxoR


! ...<<07.03.03>> : Aenderung :............................................................


! ***************************************************************************************
      REAL function FNBDATHxDxoRFoRu(BDATBArtNr,D1,H1,D2,H2,H,
     2        Hx,Dx,IErr)
! ***************************************************************************************

!     ... Hoehe (Hx) zu gegebenem Durchmesser oR mit forstlicher Rundung [Dx]:................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H
       REAL   Hx
       REAL   Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------
       FNBDATHxDxoRFoRu = xFNBDATHxDxoRFoRu(BDATBArtNr,D1,H1,D2,H2,H
     1 ,Hx,Dx,IErr)
! ----------------------------------------------------------------------------------------

      END Function FNBDATHxDxoRFoRu


! ***************************************************************************************
      REAL function xFNBDATHxDxoRFoRu(BDATBArtNr,D1,H1,D2,H2,H,
     2        Hx,Dx,IErr)
! ***************************************************************************************

!     ... Hoehe (Hx) zu gegebenem Durchmesser (Dx) :...........................................

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H
       REAL   Hx
       REAL   Dx
       INTEGER IErr

! ----------------------------------------------------------------------------------------

       INTEGER NSFktNr
       DATA  NSFktNr   /3/
       REAL   NSFktPar
       DATA    NSFktPar  /0/

       REAL   A
       DATA    A    /0/
       REAL   B
       REAL   NSFktAbsErr
       DATA    NSFktAbsErr  /1e-3/
       REAL   NSFktXAbsErr
       DATA    NSFktXAbsErr /0/
       REAL   NSFktXRelErr
       DATA    NSFktXRelErr /1e-3/

       INTEGER MIt
       DATA  MIt    /100/
       REAL   x1
       REAL   x2

! ----------------------------------------------------------------------------------------

       A = 0
       B = H

       call BDATNullStellenSuche
     1  (BDATBArtNr, D1, H1, D2, H2, H,
     2   NSFktNr, Dx,A, B,
     3   NSFktAbsErr, NSFktXAbsErr, NSFktXRelErr,
     4   MIt, x1, x2, Hx, IErr
     5  )

       xFNBDATHxDxoRFoRu = Hx

      END Function xFNBDATHxDxoRFoRu


! ...<<07.03.03>> : Aenderung :............................................................


! ****************************************************************************************
      REAL function FNBDATQ03VHDx
     1  (BDATBArtNr, D1, H1, H, Dx, VolHDx,
     2   MIt, q031, q032, q03x, IErr
     3  )
! ****************************************************************************************

! ---------------------------------------------------------------------------------------
!     Formquotient q0.3x zu VolHDx (Volumen bis zum GrenzDurchmesser Dx)
! ---------------------------------------------------------------------------------------

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   H
       REAL   Dx
       REAL   VolHDx

       INTEGER MIt
       REAL   q031
       REAL   q032
       REAL   q03x      ! Nullstelle

       INTEGER IErr

! ----------------------------------------------------------------------------------------

       REAL   D2
       REAL   H2

       REAL   NSFktAbsErr
       DATA    NSFktAbsErr    /1E-5/
       REAL   NSFktXAbsErr
       DATA    NSFktXAbsErr   /0.00/
       REAL   NSFktXRelErr
       DATA    NSFktXRelErr   /1E-3/

       REAL   x1
       REAL   x2
       REAL   XNs

       REAL   SekLng
       DATA    SekLng     /2/
       REAL   VolABmR
       DATA    VolABmR     /0/

       REAL   HDx1
       REAL   HDx2
       DATA    HDx2     /0/
       REAL   HDx3
       DATA    HDx3     /0/

       REAL   A
       DATA    A      /0/
       REAL   B
       DATA    B      /0/

       REAL   VolHDx1
       DATA    VolHDx1     /0/
       REAL   VolHDx2
       DATA    VolHDx2     /0/
       REAL   VolHDx3
       DATA    VolHDx3     /0/

! ----------------------------------------------------------------------------------------

       REAL   X3
       REAL   f1
       REAL   f2
       REAL   f3
       REAL   s12
       INTEGER It

! ----------------------------------------------------------------------------------------

       It = 0
       A=0
       x1 = q031
       x2 = q032

! GrenzHoehe/Volumen zu Dx und q031 :......................................................

       H2=0
       D2=-x1

       HDx1 = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,Hx,Dx,IErr)

       A=0
       B=HDx1

       VolHDx1 = xFNBDATVolABmR(BDATBArtNr,D1,H1,D2,H2,H,A,HDx1
     2      ,SekLng,IErr,VolABmR)

       f1= VolHDx1 - VolHDx

! GrenzHoehe/Volumen zu Dx und q032 :......................................................

       H2=0
       D2=-x2
       HDx2 = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,Hx,Dx,IErr)

       A=0
       B=HDx2
       VolHDx2 = xFNBDATVolABmR(BDATBArtNr,D1,H1,D2,H2,H,A,B,SekLng,
     2      IErr,VolABmR)

       f2= VolHDx2 - VolHDx

       If (f1 * f2.lt.0) Then
       ElseIf (f1 * f2.gt.0) Then
        IErr = 1 ! falsche Startwerte, beide kleiner oder groesser als Zielvariable
        q03x = 0 ! Rueckgabe von Null

        FNBDATQ03VHDx = q03x

        return
       Else
        IErr = 0
        If (Abs(f1).lt.Abs(f2)) Then
         XNs = x1
        Else
         XNs = x2
        End If
        q03x = XNs

        FNBDATQ03VHDx = q03x

        return
       End If

       It = 0

       Do While (It.lt.MIt)
        It = It + 1
        If (Abs(f2).lt.NSFktAbsErr) Then
         XNs = x2
         IErr = 0
         q03x = XNs
         FNBDATQ03VHDx = q03x
         return
        ElseIf
     1   (Abs(x2 - x1).le.(Abs(x2) * NSFktXRelErr+NSFktXAbsErr))
     2  Then
         XNs = x2
         If (Abs(f1).lt.Abs(f2)) Then
          XNs = x1
         End If

! ...<<11.12.02>> : Aenderung :............................................................

         iErr=0
         q03x = XNs
         FNBDATQ03VHDx = q03x
         return
        Else
         s12 = (f2 - f1) / (x2 - x1)
         X3 = x2 - f2 / s12

!   GrenzHoehe/Volumen zu Dx und q032 :..............................................

         H2=0
         D2=-x3
         HDx3 = xFNBDATHxDx(BDATBArtNr,D1,H1,D2,H2,H,Hx,Dx
     1         ,IErr)

         A=0
         B=HDx3
         VolHDx3 = xFNBDATVolABmR(BDATBArtNr,D1,H1,D2,H2,H,A,B,
     2        SekLng,IErr,VolABmR)

         f3= VolHDx3 - VolHDx

         If (f2 * f3.le.0) Then
          x1 = x2
          f1 = f2
         Else
          f1 = f1 * f2 / (f2 + f3)
         End If
         x2 = X3
         f2 = f3
        End If
       end do

! ...Zul?ssige Anzahl von Iterationen erreicht ohne Konvergenz

       IErr = 2 ! keine Konvergenz
       q03x = 0 ! Rueckgabe von Null
       FNBDATQ03VHDx = q03x
      End Function FNBDATQ03VHDx


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++         lokale Prozeduren          +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! lokal: CDFNORM (NUM)  CDFNORMInv (NUM)
!
! ****************************************************************************************
      subroutine CDFNORM (mw,StDev,x,CDFx)
! ****************************************************************************************

! Summenfunktion (CDF) der Normalverteilung :.............................................

!  USE numerical_libraries
       REAL  mw
       REAL  StDev
       REAL  x
       REAL  CDFx

       REAL  N01

       if (StDev <=0) then
        CDFx = -1
        return
       end if

       N01 = (x-mw)/StDev

       CDFx=ANORDF(N01)
!  CDFx=0.5

      end subroutine CDFNORM


! ****************************************************************************************
      subroutine CDFNORMInv (mw,StDev,CDFx,x)
! ****************************************************************************************

! Inverse Summenfunktion (Percentile) der Normalverteilung :..............................

!  USE numerical_libraries

       REAL  mw
       REAL  StDev
       REAL  CDFx
       REAL  x

       REAL  N01

       if (CDFx <= 0.0001) then
        x=-999999
        return
       else if (CDFx>0.9999) then
        X=999999
        return
       else
        N01=dinvnorm(CDFx)
!   N01=1
        x=mw+N01*StDev
       end if

      end subroutine CDFNORMInv


! ****************************************************************************************
      Subroutine BDATNullStellenSuche
     1  (BDATBArtNr, D1, H1, D2, H2, H,
     2   NSFktNr, NSFktPar,A, B,
     3   NSFktAbsErr, NSFktXAbsErr, NSFktXRelErr,
     4   MIt, x1, x2, XNs, IErr
     5  )
! ****************************************************************************************

! ---------------------------------------------------------------------------------------
!     NullStelle einer stetigen Funktion nach dem PegasusVerfahren vgl. Engeln-M?llges P2.8.2
! ---------------------------------------------------------------------------------------

       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H

       INTEGER NSFktNr
       REAL   NSFktPar

       REAL   A
       REAL   B
       REAL   NSFktAbsErr
       REAL   NSFktXAbsErr
       REAL   NSFktXRelErr

       INTEGER MIt
       REAL   x1
       REAL   x2
       REAL   XNs      ! Nullstelle
       INTEGER IErr

! ----------------------------------------------------------------------------------------

       REAL   X3
       REAL   f1
       REAL   f2
       REAL   f3
       REAL   s12
       INTEGER It

! ----------------------------------------------------------------------------------------

       It = 0
       x1 = A
       x2 = B

       Call BDATNullStellenFkt(BDATBArtNr, D1, H1, D2, H2, H,
     1       NSFktNr, NSFktPar, x1, f1)

       Call BDATNullStellenFkt(BDATBArtNr, D1, H1, D2, H2, H,
     1       NSFktNr, NSFktPar, x2, f2)

       If (f1 * f2.lt.0) Then
       ElseIf (f1 * f2.gt.0) Then
        IErr = 1
        return
       Else
        IErr = 0
        If (Abs(f1).lt.Abs(f2)) Then
         XNs = x1
        Else
         XNs = x2
        End If
        return
       End If

       It = 0

       Do While (It.lt.MIt)
        It = It + 1
        If (Abs(f2).lt.NSFktAbsErr) Then
         XNs = x2
         IErr = 0
         return
        ElseIf
     1   (Abs(x2 - x1).le.(Abs(x2) * NSFktXRelErr+NSFktXAbsErr))
     2  Then
         XNs = x2
         If (Abs(f1).lt.Abs(f2)) Then
          XNs = x1
         End If

! ...<<11.12.02>> : Aenderung :............................................................

         iErr=0
         return
        Else
         s12 = (f2 - f1) / (x2 - x1)
         X3 = x2 - f2 / s12

         Call BDATNullStellenFkt(BDATBArtNr, D1, H1, D2, H2, H,
     1        NSFktNr, NSFktPar, x3, f3)

         If (f2 * f3.le.0) Then
          x1 = x2
          f1 = f2
         Else
          f1 = f1 * f2 / (f2 + f3)
         End If
         x2 = X3
         f2 = f3
        End If
       end do

! ...Zul?ssige Anzahl von Iterationen erreicht ohne Konvergenz

       IErr = 3

      End Subroutine BDATNullStellenSuche


! ***************************************************************************************
      Subroutine BDATNullStellenFkt (BDATBArtNr,D1,H1,D2,H2,H,
     1        NSFktNr,NSFktPar,X,Fx)
! ***************************************************************************************


!     ...NullStellen Funktion für das PegasusVerfahren :.....................................


       INTEGER BDATBArtNr
       REAL   D1
       REAL   H1
       REAL   D2
       REAL   H2
       REAL   H

       INTEGER NSFktNr
       REAL   NSFktPar
       REAL   X
       REAL   Fx

! ---------------------------------------------------------------------------------------

       INTEGER IErr
       REAL   Hx
       REAL   Dx
       REAL   DmRHx,DoRHx

! ----------------------------------------------------------------------------------------



       If (NSFktNr .eq. 1) Then

!   ...X: Hoehe zu gegebenem Durchmesser: Hx = H (Dx) :...............................

        Dx = NSFktPar
        Hx = X

        DmRHx = xFNBDATDmRHx(BDATBArtNr, D1, H1, D2, H2, H, Hx
     1       ,IErr,DmRHx)
        Fx = Dx - DmRHx

       Else If (NSFktNr .eq. 2) Then


!   ...X: Hoehe zu gegebenem Durchmesser oR: Hx = H (Dx) :............................


        Dx = NSFktPar
        Hx = X

        DoRHx = xFNBDATDoRHx(BDATBArtNr, D1, H1, D2, H2, H, Hx
     1       ,IErr,DoRHx)
        Fx = Dx - DoRHx

       Else If (NSFktNr .eq. 3) Then


!   ...X: Hoehe zu gegebenem Durchmesser oR mit forstlicher Rundung Hx = H (Dx) :.....


        Dx = NSFktPar
        Hx = X

        DoRHx = xFNBDATDoRHx(BDATBArtNr, D1, H1, D2, H2, H, Hx
     1       ,IErr,DoRHx)

        Fx = Dx - xFNBDATDxFoRu(DoRHx)

       Else

        Dx = NSFktPar
        Hx = X

        DmRHx = xFNBDATDmRHx(BDATBArtNr, D1, H1, D2, H2, H, Hx
     1       ,IErr,DmRHx)

        Fx = Dx - DmRHx

       End If

      END Subroutine BDATNullStellenFkt


! ...<<07.03.03>> : Aenderung :............................................................


! ***************************************************************************************
      REAL function FNBDATDxFoRu(Dx)
! ***************************************************************************************

!     ... Durchmesser forstlich gerundet [Dx] :..............................................

       REAL   Dx

        FNBDATDxFoRu=xFNBDATDxFoRu(Dx)

      END Function FNBDATDxFoRu


! ***************************************************************************************
      REAL function xFNBDATDxFoRu(Dx)
! ***************************************************************************************

!     ... Durchmesser forstlich gerundet [Dx] :..............................................

       REAL   Dx

       if (20.le.Dx) then
        xFNBDATDxFoRu=Dx-0.75
       else
        xFNBDATDxFoRu=Dx-0.50
       end if

      END Function xFNBDATDxFoRu

! ...<<18.09.03>> : Aenderung :.......................................................

!c ***************************************************************************************
! REAL function FNBDATHStockEnde
!c ***************************************************************************************
!
!  REAL   HStockEnde
!  REAL   HStHAnfang
!  REAL   LngStH
!  REAL   HStHLzEnde
!  REAL   HBDATGes
!
!  COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
! 1     , HBDATGes
!
!  FNBDATHStockEnde=HStockEnde
!
! END Function FNBDATHStockEnde

!c ***************************************************************************************
! REAL function FNBDATHStHAnfang
!c ***************************************************************************************
!
!  REAL   HStockEnde
!  REAL   HStHAnfang
!  REAL   LngStH
!  REAL   HStHLzEnde
!  REAL   HBDATGes
!
!  COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
! 1     , HBDATGes
!
!  FNBDATHStHAnfang=HStHAnfang
!
! END Function FNBDATHStHAnfang

!c ***************************************************************************************
! REAL function FNBDATLngStH
!c ***************************************************************************************
!
!  REAL   HStockEnde
!  REAL   HStHAnfang
!  REAL   LngStH
!  REAL   HStHLzEnde
!  REAL   HBDATGes
!
!  COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
! 1     , HBDATGes
!
!  FNBDATLngStH=LngStH
!
! END Function FNBDATLngStH

!c ***************************************************************************************
! REAL function FNBDATHStHLzEnde
!c ***************************************************************************************
!
!  REAL   HStockEnde
!  REAL   HStHAnfang
!  REAL   LngStH
!  REAL   HStHLzEnde
!  REAL   HBDATGes
!
!  COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
! 1     , HBDATGes
!
!  FNBDATHStHLzEnde=HStHLzEnde
!
! END Function FNBDATHStHLzEnde

!c ***************************************************************************************
! REAL function FNBDATHBDATGes
!c ***************************************************************************************
!
!  REAL   HStockEnde
!  REAL   HStHAnfang
!  REAL   LngStH
!  REAL   HStHLzEnde
!  REAL   HBDATGes
!
!  COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
! 1     , HBDATGes
!
!  FNBDATHBDATGes=HBDATGes
!
! END Function FNBDATHBDATGes

! #############################  BWI BDAT 1.0  ###########################################

! Aenderungen: <<19.01.04>> :--------------------------------------------------------------
! Aenderungen: <<05.11.02>> :--------------------------------------------------------------
! Aenderungen: <<01.07.02>> :--------------------------------------------------------------

!....                   Version:17.10.01
! Aenderungen: <<19.01.04>> :--------------------------------------------------------------

!      REAL wH095

! ...<<17.10.01>> :...(Rindenfunktion)  ..................................................

!.... Aenderung 23.10.00 Kublin :.........................................

!     M?glichkeit (Ddo,Hho) abweichend vom D7 oder q0.30 vorzugeben

!.... Aenderung 23.10.00 Kublin :.........................................

!     *******************************************************************
!     * Voluminierung & Sortierung  ***  VOL = INTEGRAL(SPLINE(Du,Do,H) *
!     *******************************************************************
!     * FILE: BDAT_B.FOR * STAND: 13.12.90 * AUTOREN: SCHARNAGL/KUBLIN  *
!     *******************************************************************


!     ###################################################################
      SUBROUTINE BDAT(ibba,Ddu,Hhu,Ddo,Hho,Hoe,Stx,Hkz,Sk,Azopf,
     *                     ksthh,zzost,Zzoab,Sok,Klasse,Volum,Bhd,ifeh)
!     ###################################################################


!                                               VERSION: 13.12.90 Kublin
!  Eingabe(E)/ Ausgabe(A) - Variable
!
!        Name(bei Eingabe-Ausgabe)
!             Name (im Programm) bei m?glicher WertAenderung
!                        Variablentyp  (i=Integer,R=Real)
!                               Beschreibung
!                               ============
!  E     Ibba(Iba)       I      Baumart nach Baumartenschluessel
!                               s. Block data
!                               iba bleibt konstant !
!  E     Ddu (Du)        R      Unterer Durchmesser (cm)
!  E     Hhu (Hu)        R      Hoehe des  unt.Durchmessers am Stamm (cm)
!  E     Ddo (Do)        R      Oberer Durchmesser (cm)
!                               Ddo > 0 und Hho=0 und ->                     Hho = 7
!                               Ddo < 0               -> q0.3 = ABS(Ddo) und Hho = 0.3*Hoe
!            (D7m)       R                  "

!.....Aenderung 23.10.00 Kublin :..........................................

!  E     HHo                    Hoehe oberer Durchmesser 0 und Ddo>0 -> Hho = 7m
!                                                             Ddo<0 -> Q0.3 = ABS (Ddo) Hho = 0.3*Hoe

!.....Aenderung 23.10.00 Kublin :..........................................

!  E     Hoe (H)         R      Hoehe des Baumes (m)
!            (Hbr)       R      Hilfsv. fuer H,H wird durch Hkz
!                               ge?ndert
!  E     Stx (Stxu)      R      L?nge X-Holz am Stammfuss
!  E     Hkz             I      Hoehenkennziffer 1 = Wipfelbruch
!                                              2 = Gipfelbruch
!  E     Sk              I      Stammkennziffer  s. Dokum. S.83
!  E     Azopf(Azop)     R      Aufarbeitungszopf(cm)
!  E     Ksthh(sthh)     R      Hoehe bis zu der Stammholz ausgehalten
!                                     werden darf
!  E     Zzost(zost)     R      Zopf bis zu dem Stammholz ausgehalten
!                               werden darf
!  E     Zzoab(zoab)     R      Zopf bis zu dem Abschnitt ausgehalten
!                               werden darf
!  E     Sok(Sokz)       I      Sortierkennziffer 0=keine Sortierung
!                                    1=Mittenst. 2=Heilbronner
!  A     Klasse          I      St?rkeklassen (Feld 1-6)
!         (1...Kl )      I      Klasse     Stammholz          a
!         (2...Ukl)      I      Unterklasse    "     werden fuer Ausg.
!         (3...Klx)      I      Klasse       X-Holz         in Klasse
!         (4...Uklx)     I      Unterklasse    "          uebertragen
!         (5...Kla )     I      Klasse    Abschnitt
!         (6...Ukla)     I      Unterklasse    "
!
!
!  A     Volum           R      Volumenfeld  (Feld 1-7)
!         ( 1...Vol )           Hilsgr??e für Indh.ber.Laubholz
!                                (Stammholz o.R. incl. X-Holz)
!         ( 2...Volx)           X-Holz Efm. o. Ri.
!         ( 3...Vols)           Stammholz   "
!         ( 4...Vola)           Stammteil(-abschnitt)  "
!         ( 5...Voli)           Industrieholz          "
!         ( 6...Volu)           nicht verwertetes Derbholz  "
!         ( 7...Volr)           Ernteverlust                "
!
!  A     Bhd             R      BHD aus Schaftk. wenn Du <>1,30
!  A     Ifeh            I      Fehler der Eingabegr?ssen
!                               nach Fehlertabelle s.Doc.  S.84 ff.
!
!==================================================================
!   Variable im Programm:
!        Hazop           R      Hoehe (Lage) des Aufarbeitungszopfes
!        Hdgren          R       "            "  Derbholzgrenze
!        Hsthzop         R       "            "  Stammholzzopfes
!        Habzop          R       "            "  Abschnittzopfes
!        Ba              I      zugeordneter Baumarten-Index fuer
!                               best.Funktion
!
!   Schaft
!
!==================================================================
!        Konstante :
!        ===========
!  K     Dgr (Dgrenz)    R      Derbholzgrenze  = 7 cm
!
! =================================================================
!
!  K     Volk(6,30,3)    I      Volumen nach Krenn (ba,mm-St,h?-St)
!                               s.BD.
!  K     Hoehr(6,6)      R      Hoehenrahmen fur Krenn-Vol.(h?-st)
!  K     Azo(7,3)        R      Durchschn. Aufarbeitungszopf der
!                               Baumarten
!  K     Rin(28,4,3)     R      Koeffizienten fuer Rindenfunktionen
!                               Rin(Ba,Hoehen-St im Baum,Koeffz-Nr)
!  K     Rinh(3,5,3)     R      wie Rin jedoch fuer H-Sortierung
!  K     Ban (36,7)      I      Baumartenindex verschiedene
!                               Funktionen nach Zuordnungstabelle
!
!  Schaftkurvenvariable:                 s.Dok.
!  K     Nnp,Np,Nxk,B,Xk
!  K     nxk95a,nxk95b,xk95a,xk95b,a95,b95
!  K     Add07,Md07,Snxkn7,Sxk07,Sdo7
!  K     Add07(14,5:45,3)I      Addresse der D07-Werte im Feld
!                               Md07(20333)   Add07(Ba,BHD,  )
!  K     Md07(20333)     I      Feld der Massentafel-D07 (*1000)
!
!
!        Yy              R      Hilfsvariable
!        Durel           R
!        Hurel           R      relative Lage des unteren Durchmessers
!        Horel           R      relative Lage des oberen Durchmessers
!  K     Ho              R      Hoehe (Lage) des oberen Durchmessers
!                                auf 7m festgesetzt
!        D095            R      Gesch?tzter Durchmesser in 5% der Hoehe
!        Dnorm           R      Normungsdurchmesser (=Du)
!        D07             R      Relativdurchmesser in 30% der Hoehe
!                                 iteriert oder aus Tabelle (SUB D07tab)
!        Xsi             R      Hilfsgr?sse bei D07-Iteration
!        Ifehl           I      Fehler-Zeiger zeigt an, ob Iteration
!                                 erfolgreich war
!
!
!        D07lu           R      untere Iterationsgrenze fuer D07
!        D07lo           R      obere  Iterationsgrenze fuer D07
!
!
      implicit logical (a-z)


! ...<<11.11.02>> : Aenderung :.......................................................
!
! REAL DuTmp,HuTmp,DoTmp,HoTmp,HTmp,DoLoTmp,DoLuTmp,DTmp

!............................................Aenderung: 17.01.91 Kublin

      REAL ABS,AMAX1,AMIN1
!     REAL hhxx

!.....................................................................

      INTEGER Nnp,Np,Nxk,Sk,Sokz,Sok,Hkz,I
      REAL B(1:80,1:8),Xk(1:6) ,Azopf,Stx,Hoe,Ddu,Ddo,Hhu,X,D7m,Ksthh

!.....Aenderung 23.10.00 Kublin :..........................................

      REAL Hho, wDdu

!.....Aenderung 23.10.00 Kublin :..........................................

      REAL Rin(1:28,1:4,1:3)
      REAL Rinh(1:3,1:5,1:3)
      REAL a95(1:24,1:8),b95(1:24,1:8)
      REAL xk95a(7,8),xk95b(7,8)
      INTEGER nxk95a(8),nxk95b(8)
      REAL Yy,Durel,Xsi
      INTEGER Nbaum,ibba,ifeh
      REAL D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
      REAL Spline,hhrel,kw
      REAL Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zzost,Zost
      REAL Zzoab,Zoab,Hbr

! Aenderungen: <<19.01.04>> :--------------------------------------------------------------
      REAL wH095
! Aenderungen: <<19.01.04>> :--------------------------------------------------------------

!
      INTEGER Ifehl,Ianz
      REAL F1,F2,D07lo,D07lu,Bb,Aa
      REAL Azo(7,0:3)
      INTEGER Kl,Ukl,Klx,Uklx,Kla,Ukla
      REAL Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      INTEGER Ba,Iba,Ban(36,7)
      REAL sthh,Stxu
      INTEGER Volk(6,30,3),Klasse(6)
      REAL Volum(7),bhd ,dgr,hoehr(6,6)
      INTEGER md07(20333)
      INTEGER Add07(14,5:45,3)
      REAL Sxk07(7,8),Sd07(24,8)
      INTEGER Snxkn7(8)

! ...<<18.09.03>> : Aenderung :.......................................................

      REAL   HStockEnde
      REAL   HStHAnfang
      REAL   LngStH
      REAL   HStHLzEnde
      REAL   HBDATGes
!   REAL glLSort(1:5), glDSort(1:5) ! christian vonderach, 24.07.2018

      COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
     1     , HBDATGes

! ...<<18.09.03>> : Aenderung :.......................................................

      COMMON /Volk/Volk,hoehr
      COMMON /Wert1/ sthh,Sokz,Stxu
      COMMON /Klasse/ Kl,Ukl,Klx,Uklx,Kla,Ukla
      COMMON /Volum/ Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      COMMON /Duazop/ Azo
      COMMON /It/Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      COMMON /Rind/Rin,Rinh
      COMMON /Baum/ Ba,Iba,Ban
      COMMON /Schaft/  Nnp,Np,Nxk,B,Xk
      COMMON /d95/nxk95a,nxk95b,xk95a,xk95b,a95,b95
      COMMON /D07/ Add07,Md07,Snxkn7,Sxk07,Sd07
      COMMON /Sk/ Yy,Durel
      COMMON /Baum0/ Nbaum
      COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
!   COMMON /glLDSort/ glLSort, glDSort ! christian vonderach, 24.07.2018

! ----------------------------------------------------------------------------------------

!     return
!     write(*,*) '***** Hallo *****'
!     return
!
      do 10 i=1,6
        Klasse(i)=0
        Volum(i)=0
   10 continue

      hstzop=0
      habzop=0
      hazop=0
      ifeh=0
      Volum(7)=0

      Nbaum=1
!                                        ....Aenderung 31.01.91 Kublin
      Iba=ibba
      if(iba.lt.1.or.iba.gt.36)ifeh=1
      Stxu=Stx
      Azop=Azopf
      dgr=7
      Dgrenz=dgr

      Du=Ddu
      Bhd=Ddu

!.....Aenderung 23.10.00 Kublin :..........................................

      wDdu=Ddu
      Do=Ddo

!.....Aenderung 23.10.00 Kublin :..........................................

!.....                   oberer Durchmesser darf nicht > unterer sein

      if(do.gt.du*.99)do=du*.99

      Hu=Hhu
      H=Hoe
      if(h.le.0.and.du.ge.10)ifeh=2

!  Aenderung <<15.12.03>> :--------------------------------------------------------------

      if(du.le.0)ifeh=3
!     if(hu.gt.2.5)ifeh=4
      if(ifeh.gt.0)return
      sthh=Ksthh
      Zost=Zzost
      Zoab=Zzoab
      Hbr=h
      Sokz=Sok

!.....                                     Hoehenkorrektur bei Hkz > 0

      if(hkz.eq.1)H=h+2
      if(hkz.eq.2)THEN
        IF(Du.gt.30)THEN
          H=30+(DU-30)*.3
        ELSE
          H=Du
        END IF

!.....                              Mindestzuschlag bei Wipfelbruch 4 m

        IF(Hbr.gt.H-3)H=Hbr+4
      END IF
!     write(10,*) 'h,hbr',h,hbr
      IF (Hu.le.0)  Hu=1.3


! DuTmp=DDu
! HuTmp=Hu
! DoTmp=DDo
! HoTmp=HHo
! HTmp = H


!.....................................................................

      IF(h.ge.3)then

        Ho=Hho
        Horel=1.-Ho/H
        Hurel=1.-Hu/H

        D7m=Do
        dnorm=du

!.....  Aenderung 23.10.00 Kublin :....................................

!       IF(Sk.eq.0.AND.Du.lt.20) D7m=0
!
!.....  Aenderung 23.10.00 Kublin :....................................

!.....  ............D095-Schätzung aus Du u H  (Du wird =D1.3 gesetzt)
!

! Aenderungen: <<19.01.04>> :--------------------------------------------------------------
!     REAL wH095
! Aenderungen: <<19.01.04>> :--------------------------------------------------------------

        if (h>60) then
       wH095=60
        else
          wH095=H
        end if

        Ba=Ban(iba,1)
        Aa=Spline(wH095,Ba,Nxk95a,Xk95a,A95)
        Bb=Spline(wH095,Ba,Nxk95b,Xk95b,B95)
        D095=Aa+Bb*Du

! Aenderungen: <<19.01.04>> :--------------------------------------------------------------

!
        IF (D7m.GT.0.)  THEN

!.....    d07-Startwert :..............................................

!         Do=D7m                            (Aenderung: 12.12.90 Kublin)

          D07lu=.4
          D07lo=.99

          Ifehl=-1

          CALL Fkt(D07lu,F1)
          CALL Fkt(D07lo,F2)

!  F=Fdrel(D095,X,H,1-Ho/H)/Fdrel(D095,X,H,1-Hu/H)-Do/Du

          Ianz=0
          Xsi=0

          CALL Pegasu(Ifehl,Ianz,F1,F2,Xsi,d07lu,d07lo)

!         write(10,*)'ifehl,d07lu,xsi,d07lo',ifehl,d07lu,xsi,d07lo

          IF (Ifehl.LT.3.AND.Ifehl.GE.0) THEN
            D07=AMAX1(D07lu,AMIN1(Xsi,D07lo))
          ELSE

            IF(ifehl.eq.-1)THEN

!   ...<<11.11.02>> : Aenderung :...............................................

!   d07=D07lu
!   call Kuwert(HHo,DoLuTmp)

!   REAL DuTmp,HuTmp,DoTmp,HoTmp,HTmp,DoLoTmp,DoLuTmp


!   DoLuTmp = xFNBDATDmRHx(Iba, 135.0, 1.3,65.0,
! 1      7.0,HTmp,7,ifehl, DTmp)
!   DoTmp=-0.90
!   HoTmp=0.0
!   DuTmp=135
!   HuTmp=0.0
!   HTmp=60
!   x=7.0
!   Iba=3
!   ifehl=0
!    DoLuTmp = xFNBDATDmRHx(Iba,DuTmp, HuTmp, DoTmp,
!  1       HoTmp,HTmp, x, ifehl, DTmp)

! --------------------------------------------------------------------
!   DoTmp=-.40
!   HoTmp=0.0
!   DuTmp=135
!   HuTmp=0.0
!   HTmp=60
!   x=7
!   Iba=3
!   ifehl=0
!   DTmp=0.0

!    DoLoTmp = xFNBDATDmRHx(Iba, DuTmp, HuTmp, DoTmp,
!  1       HoTmp,HTmp,x, ifehl, DTmp)

!     if (DoLoTmp>65) then
!     d07=.9333
!     else
!        d07=.40
!     endif

! --------------------------------------------------------------------

        d07=D07lo

!
! ...<<11.12.02>> : Aenderung :............................................................

!    call Kuwert(HHo,DoLoTmp)


!.....                                       (Aenderung: 12.12.90 Kublin)
                 ifeh=12

                 CALL D07tab(ifehl)
                 if (ifehl.eq.0) then
                   ifeh=13
                 end if
!.....
            ELSE
               ifeh=14
               d07=Xsi
            END IF
          END IF

        ELSE

!.....    Aenderung 23.10.00 Kublin :....................................

          IF (D7m.LT.0.)        THEN

        Horel=.7
        D07 = ABS(D7m)
        Ho=H*.3

!.....  Aenderung 29.08.01 (Kublin) :....................................................

        if (1<=D07) then
         CALL D07tab(ifehl)
         iFEH = 15
        end if

!   Aenderung 03.09.01 :.......................................
!
!    D07=AMAX1(D07lu,AMIN1(D07,D07lo))

          ELSE

!.....      Aenderung 23.10.00 Kublin :..................................

            Horel=.7
            CALL D07tab(ifehl)
!.....                                       (Aenderung: 12.12.90 Kublin)
            if (ifehl.eq.0) then
              ifeh=13
            end if
!.....
            Ho=H*.3
          END IF
        END IF

!........................Unterer Durchmesser nicht bei 1.30m gemessen :

        IF(ABS(hu-1.3).gt.0.005) then
!                                        ....Aenderung 31.01.91 Kublin
          hhrel=1.-1.3/h
          call Kuwert(Hhrel,Kw)
          bhd=kw
        END IF

      ELSE
!       ...................................(Aenderung: 12.12.90 Kublin)

!       H < 3 m : Anpassen einer Schaftkurve für Testzwecke. Wird bei
!       bei der BWI nicht gebraucht, da für diese B?ume kein Volumen
!       berechnet wird

!       ..............................................................

!       Ba=Ban(iba,1)
!       Aa=Spline(H,Ba,Nxk95a,Xk95a,A95)
!       Bb=Spline(H,Ba,Nxk95b,Xk95b,B95)
!       D095=Aa+Bb*Du
!
!       if(hu.le.0) then
!         hu=1.3
!       end if
!       hurel=1.0-hu/h
!
!       Horel=.7
!       CALL D07tab(ifehl)
!
!       if (ifehl.eq.0) then
!         ifeh=13
!       end if
!
!       Ho=H*.3
!
!       dnorm=du
!.....  .............................................................

      END IF

      IF(D07.gt..95) D07=.95
      if(d07.lt..40) D07=.40

!............................................Aenderung 11.12.90 Kublin

      nbaum = 1

!............................................Aenderung 17.01.91 Kublin

!
!     Ausgabe der Baumparameter und der Durchmesser bei 0.0(0.1)1.0
!     für Testzwecke. Ausgabe auf FILE # 10.
!
!.....................................................................
!
!     write(10,'(3(/),5x,a5,f10.3  )') 'H:   ',H
!     write(10,'(     5x,a5,f10.3  )') 'Hu:  ',hu
!     write(10,'(     5x,a5,f10.3  )') 'Du:  ',du
!     write(10,'(     5x,a5,f10.3,/)') 'BHD: ',bhd
!
!     write(10,'(     5x,a5,f10.3  )') 'D095:',d095
!     write(10,'(     5x,a5,f10.3,/)') 'D070:',d07
!
!     if (du.gt.0) then
!
!       hhxx  = 1.-1./h
!       call kuwert(hhxx,kw)
!       write(10,'(5x,a5,3f10.3)')     'D1 : ',1.0-hhxx,kw/du,kw
!
!       hhxx=1.-1.3/h
!       call kuwert(hhxx,kw)
!       write(10,'(5x,a5,3f10.3)')     'BHD:  ',1.0-hhxx,kw/du,kw
!
!       if (h.gt.7) then
!         hhxx=1.-7./h
!         call kuwert(hhxx,kw)
!         write(10,'(5x,a5,3f10.3)')    'D7 : ',1.0-hhxx,kw/du,kw
!       else
!         write(10,*)
!       end if
!
!       write(10,'(/)')
!
!       do 12300 hhxx=0.00,1.0001,.1
!         call kuwert(1.-hhxx,kw)
!         write(10,'(10x,3f10.3)') hhxx,kw/du,kw
!12300  continue
!
!     end if
!
!......................................................................

!     Liegt der Du bei 1.30m Hoehe und ist < 10, wird eine Durchmesser-
!     abrundung durch + 0.5 ausgeglichen :

!.....Aenderung 23.10.00 Kublin :..........................................

      if(Bhd.eq.Ddu)wDdu=Ddu+.5

!...................................................Volumen nach Krenn
!..............................................Aenderung 11.03.94 St?hr

      IF(wDdu.lt.10)then
       IF(wDdu.lt.7) then
           ifeh=15
           RETURN
         ENDIF
         ba=ban(iba,6)
       Volr=Volk(ba,INT((wDdu-7)*10.0001)+1,2)*0.0001
         Volum(1)=Volr
         Volum(7)=volr*.3
         Volum(5)=volr-volum(7)

!.....Aenderung 23.10.00 Kublin :..........................................

         RETURN
      ELSE
        IF(H.lt.3)then
        ifeh=15
          RETURN
        ENDIF
      ENDIF

!        Schätzung des Aufarbeitungszopfes mit Funktionen falls nicht
!          eingegeben:

! Aenderungen: <<05.11.02>> :--------------------------------------------------------------

!     IF(Azop.eq.0) THEN

      IF(Azop.lt.0.001) THEN
        Ba=Ban(Iba,3)
        X=Bhd
        Azop=EXP(Azo(Ba,0)+Azo(Ba,1)*LOG(X)+Azo(Ba,2)*(LOG(X)**2)+
     *           Azo(Ba,3)*X)
      END IF

      IF(Sk.gt.4) THEN
       RETURN
      END IF

!     Ueberpruefung der Sortier-Eingabeparameter bei Laubholz:.......

      IF(Iba.gt.14) THEN
        if(sokz.gt.0) sokz=3
        IF(Sk.eq.1) then
          IF(sthh.gt.0) ifeh=6
          sthh=hbr*.7
        end if
!.....                                    (Aenderung: 12.12.90 Kublin)

        IF(Sk.eq.2) then
          if (sthh.gt.0) then
            IF(sthh.gt.7.or.sthh.lt.1.3) then
              ifeh=7
              sthh=5
            endif
          else
            sthh=5.0
          end if
        endif
!.....                                     (Aenderung: 12.12.90 Kublin)
        IF(Sk.eq.3.)then
          if (sthh.gt.0) then
            if (sthh.gt.3.or.sthh.lt.1.3) then
              ifeh=8
              sthh=.1
            end if
          else
            sthh=0.1
          end if
        end if
!.....                                   (Aenderung: 12.12.90 Kublin)
        IF(Sk.eq.4) then
          if (sthh.gt.0) then
            IF(sthh.gt.Hbr*.66) then
              ifeh = 9
              sthh = hbr*.66
              hazop= hbr*.70
            else
              hazop  = hbr*.7
            end if
          else
            sthh   = hbr*.66
            hazop  = hbr*.7
          end if
        end if

      ELSE

!       .................................................Nadelholz :

        IF(Sokz.gt.2) sokz=1

!       write(10,*)'sokz,sthh,zost'
!       write(10,*) sokz,sthh,zost

        IF(Sokz.eq.1.and.Sthh.le.0.and.zost.le.0) Sokz=4

!       write(10,*)'sokz'
!       write(10,*) sokz

        if(sk.eq.0.or.sk.eq.1)then
          IF(sthh.gt.0) THEN
             ifeh=10
          ELSE
             Sthh=Hbr
          END IF
        END IF
!
        if(sk.eq.2)then
          sthh=5
          habzop=5
          hstzop=5
        end if
        if(sk.eq.3)then
          sthh=.1
          habzop=.1
          hstzop=.1
        end if

!.....                                     (Aenderung: 12.12.90 Kublin)

        if(sk.eq.4)then
          if (hkz.eq.0) then
            ifeh=11
          end if
          if (sthh.gt.0) then
            if(sthh.gt.hbr*.66)then
              ifeh=9
              sthh = hbr*.66
              hazop= hbr*.70
            else
              hazop  = hbr*.7
            end if
          else
            sthh   = hbr*.66
            hazop  = hbr*.7
          end if
        end if
      END IF
!
!.................................................Volumen Vfm mit Rinde
!
         volx=0
         volu=0
         vols=0
         voli=0
         vola=0
         kl=0
         ukl=0
         kla=0
         ukla=0
         klx=0
         uklx=0

!.....Berechnung Vorrats, Lage Grenzdm und -Hoehen fuer die Sortierung...

        CALL Iter
        Volum(1)=Volr

! .......................................Aufruf eines Sortierprogrammes:

      IF(Sokz.gt.0) then
          IF(iba.gt.14)sokz=3
          IF(Sokz.eq.1)then
             CALL Sortmi
          ELSE
             IF(sokz.eq.2.or.sokz.eq.4)then
                CALL sorthl
             ELSE
                CALL sortlb
             ENDIF
          ENDIF
!
!...................................Klassen und Unterklassen der Sorten:
!
!    X-Holz
         Klasse(1)=Klx
         Klasse(2)=Uklx
!    Stammholz
         Klasse(3)=Kl
         Klasse(4)=Ukl
!    Stammholz-Abschnitt
         Klasse(5)=Kla
         Klasse(6)=Ukla
      do 20  i=1,6,2
      if(klasse(i).gt.6) then
        klasse(i)=6
        klasse(i+1)=1
      endif
  20  continue

!  write (*,*) " BDAT ! VOLi: ", VOLi
!
! Volumina
!    X-holz
         Volum(2)=Volx
!    Stammholz
         Volum(3)=Vols
!    Stammholz-Abschnitt
         Volum(4)=Vola
!    Industrieholz
         Volum(5)=Voli
!    nicht verwertetes Derbholz
         Volum(6)=Volu
!    Ernteverlust
         Volum(7)=Volum(1)-volx-vols-vola-voli-volu
      ENDIF
      RETURN
!     debug subchk
      END SUBROUTINE


!     ################## BDAT - Unterprogramme #########################


!     ###################################################################
      SUBROUTINE D07tab(mt)
!     ###################################################################

      implicit logical (a-z)
!.....                                       (Aenderung: 12.12.90 Kublin)
      INTEGER mt
!.....
      INTEGER md07(20333)
      INTEGER Add07(14,5:45,3)
      INTEGER Ba,Iba,Ban(36,7)
      INTEGER dd
      REAL Spline
      INTEGER hh,Ad1
      REAL Sxk07(7,8),Sd07(24,8)
      INTEGER Snxkn7(8)
      REAL D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
!
      COMMON /D07/ Add07,Md07,Snxkn7,Sxk07,Sd07
      COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
      COMMON /Baum/ Ba,Iba,Ban
!
      Ba=Ban(Iba,7)

!.....Rundung auf ganze Hoehen und Durchmesser fuer den Tabellenzugriff:

      Hh=int(H+.5)
      Dd=int(Du+.5)

!.....                 Begrenzung auf Minimal und Maximal-Hoehen-Werte :

      IF(Hh.gt.45)Hh=45
      IF(Hh.lt.5)Hh=5

!.....                        Durchmesser innerhalb/ausserhalb  Tabelle:

      IF(Dd.gt.Add07(Ba,Hh,3).or.Dd.lt.Add07(Ba,HH,2).or.
     *   Add07(Ba,Hh,1).le.0)Then
         Ba=Ban(iba,1)

!.....                           Schätzung ueber Spline-Schaetzfunktion:


         D07=Spline(H,Ba,Snxkn7,Sxk07,Sd07)
!        write(10,*)'d07spline  ',d07,h,ba

!.....                                       (Aenderung:12.12.90 Kublin)
         mt = 0
!.....

      ELSE

         Ad1=Add07(Ba,HH,1)+Dd-Add07(Ba,Hh,2)
         D07=Md07(Ad1)*.001
!        write(10,*)'md07  ', d07,ad1,ba,add07(ba,hh,1),add07(ba,hh,2)

!.....                                       (Aenderung:12.12.90 Kublin)
         mt = 1
!.....

      END IF
      END SUBROUTINE


!     ###################################################################
      SUBROUTINE KUWERT(Hhrel,Kw)
!     ###################################################################

!     Das Unterprogramm gibt an der relativen Lage (Hhrel) des Schaftes
!     den Durchmesser aus.
!
      implicit logical (a-z)
        REAL D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
        INTEGER Nbaum
        INTEGER Ba,Iba,Ban(36,7)
        REAL Yy,Durel,D1rel
! REAL Ddo,Ddu,Fdrel,Hhrel
      REAL Fdrel,Hhrel
        REAL Kw ,eins
        COMMON /Baum0/ Nbaum
        COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
        COMMON /Sk/ Yy,Durel
        COMMON /Baum/ Ba,Iba,Ban
!
        ba=ban(iba,1)
        eins=1.
        IF (Nbaum.GT.0) THEN

!.....         Durel,Yy,Ddo,Ddu werden nur beim ersten Aufruf berechnet
!                                                         (wenn Nbaum>0

          Nbaum=0
          Durel=Fdrel(D095,D07,H,Hurel)
          Yy=1./Durel

!.....    ...................................(Aenderung: 12.12.90 Kublin)
!         Ddo=Fdrel(D095,D07,H,Horel)*Yy*Dnorm
!         Ddu=Fdrel(D095,D07,H,Hurel)*Yy*Dnorm
!.....    ..............................................................

        END IF
!
        IF (Hhrel.GT.Hurel) THEN
            D1rel=Fdrel(D095,D07,H,eins)
            IF(D1rel.lt.Durel)THEN
               Kw=Dnorm
               RETURN
            ELSE
               Kw=Fdrel(d095,d07,H,Hhrel)*Yy*Dnorm
            ENDIF
        ELSE
          Kw=Fdrel(d095,d07,H,Hhrel)*Yy*Dnorm
        END IF
      RETURN
!     debug subchk
      END SUBROUTINE


!     ###################################################################
      SUBROUTINE Pegasu(Ifehl,Ianz,F1,F2,Xsi,A,Bb)
!     ###################################################################

!     Iterationsprogramm für die D07 Berechnung :

!      delt             Genauigkeitsminimum
!      Nmax             Maximale Anzahl der Iterationen
!
      implicit logical (a-z)
        REAL S12,F1,F2,X1,X2,Xsi,A,BB,Delt,X3,F3
        INTEGER Ifehl,Ianz,Nmax,I
        data delt/.01/,Nmax/20/
!
        Ianz=0
        X1=A
        X2=Bb
        IF (F1*F2.GT.0) THEN
          Ifehl=-1
          RETURN
        END IF
        IF (F1*F2.EQ.0) THEN
          Ifehl=0
          RETURN
        END IF
!----------------------------------------------------------
        DO 30   I= 1, Nmax
          Ianz=I
          IF (F2.EQ.0) THEN
            Xsi=X2
            Ifehl=1
            RETURN
          END IF
          IF (ABS(X2-X1).LE.Delt) THEN
            Xsi=X2
            IF (ABS(F1).LT.ABS(F2))  Xsi=X1
            Ifehl=2
            RETURN
          ELSE
            S12=(F2-F1)/(X2-X1)
            X3=X2-F2/S12
            CALL Fkt(X3,F3)
            IF (F2*F3.LE.0)  THEN
              X1=X2
              F1=F2
            ELSE
              F1=F1*F2/(F2+F3)
            END IF
            X2=X3
            F2=F3
            Xsi=X3
          END IF
 30     CONTINUE
        Ifehl=3
      RETURN
!     debug subchk
      END SUBROUTINE


!     ###################################################################
      SUBROUTINE Fkt(X,F)
!     ###################################################################

      implicit logical (a-z)
        REAL D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
        COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
        REAL X,F,Fdrel
        F=Fdrel(D095,X,H,1-Ho/H)/Fdrel(D095,X,H,1-Hu/H)-Do/Du
!
      RETURN
      END SUBROUTINE


!     ###################################################################
      REAL FUNCTION Fdrel(D095,D07,H,Hhrel)
!     ###################################################################

        implicit logical (a-z)
        INTEGER Nnp,Np,Nxk,Ftint
        REAL B(1:80,1:8),Xk(1:6),D095,D07,h,Hhrel
        INTEGER ba,ban(1:36,1:7),iba
        REAL X(1:4)
        COMMON /Schaft/  Nnp,Np,Nxk,B,Xk
        COMMON /Baum/ Ba,Iba,Ban
        INTEGER I,Ii,Ix
        REAL A,Bb,C,D,T,U
        REAL Klein
!
!       write (10,*) ba,ban(iba,1)
        Klein=1.E-4
        fdrel=0.
        X(1)=1.0
        X(2)=H
        X(3)=D095
        X(4)=D07
        Ix=Ftint(Hhrel)
        DO 40  I= 1,4
          Ii=(Ix-1)*4+(I-1)*Nnp
          T=(Hhrel-Xk(Ix))/(Xk(Ix+1)-Xk(Ix))
          U=1.-T
!....                                              Splinekoeffizienten:
          A=B(Ii+1,ba)
          Bb=B(Ii+2,ba)
          C=B(Ii+3,ba)
          D=B(Ii+4,ba)
          Fdrel=Fdrel+X(I)*(A*U+Bb*T+C*U*U*U+D*T*T*T)
 40     CONTINUE
        IF (Hhrel.LT.Klein) THEN
          Fdrel=0
        END IF
        RETURN
      END FUNCTION


!     ###################################################################
      INTEGER FUNCTION Ftint(Hhrel)
!     ###################################################################

        implicit logical (a-z)
        INTEGER Nnp,Np,Nxk
        REAL B(1:80,1:8),Xk(1:6),Hhrel
        COMMON /Schaft/  Nnp,Np,Nxk,B,Xk
        INTEGER I
        DO 50  I=1, nxk
          IF (Hhrel.LT.Xk(I)) THEN
            Ftint=I-1
            GOTO 60
          END IF
 50     CONTINUE
          Ftint=Nxk-1
          Hhrel=Xk(Nxk)
 60     IF (FTINT.EQ.0) THEN
          Ftint=1
          Hhrel=0
        END IF
        RETURN
      END FUNCTION


!     ###################################################################
      REAL FUNCTION   Spline(x,ba,nk,Xk,B)
!     ###################################################################

!
!  Berechnung (d/dx)**a [SPLINE(x)] ; a=0,1,2
!=======================================================================
!
! X         - EIN: Abszissenwert
! Nk        - EIN: # Knotenpunkte der Zerlegung des Def.-Ber. <a,b>
! Xk(nk)    - EIN: Knotenpunkte der Zerlegung: a=xk(1)<...<xk(nk)=b
! B(4(nk-1))- EIN: Splinekoeffizienten
!
!=======================================================================
!
        implicit logical (a-z)
       INTEGER I,Iu,Io,Fehler,nk(8),ba
       REAL Eps,U,T,x,xk(7,8),b(24,8)
       Eps=1.E-6
       Fehler=0
!
!=======================================================================
!
       Spline=0.
!
       if (x.lt.Xk(1,ba)-Eps)then
         Fehler=1
       ELSE IF(x.lt.Xk(1,ba)+Eps)then
         T=0.
         I=1
       ELSE IF(x.lt.Xk(Nk(ba),ba)-Eps)then
         Iu=1
         Io=Nk(ba)
         CALL Bisekt(X,ba,Xk,Iu,Io,I)
         T=(X-Xk(I,ba))/(Xk(I+1,ba)-Xk(I,ba))
       ELSE IF(x.lt.Xk(Nk(ba),ba)+Eps)then
         T=1.
         I=Nk(ba)-1
       ELSE
         Fehler=1
       ENDIF
       IF (Fehler.eq.0) THEN
         U=1.-T
         Iu=(I-1)*4
!...y=Spline(x)
      Spline=B(Iu+1,ba)*U+B(Iu+2,ba)*T+B(Iu+3,ba)*U*U*U+B(Iu+4,ba)*T*T*T
       ENDIF
       RETURN
!
!     debug subchk
      END FUNCTION


!     ###################################################################
      SUBROUTINE Bisekt(X,ba,Xk,Iu,Io,Iku)
!     ###################################################################


!     ==================================================================
!     Iku: xk(iku) <= x < xk(iku+1) / Bin?res Suchen des Teilintervalls
!     ==================================================================
!
        implicit logical (a-z)
       REAL X,Xk(7,8)
       INTEGER Ba,Iko,Mitte,Iu,Io,Iku
       Iku=Iu
       Iko=Io
!
  10   IF((Iko-Iku).lt.2)return
         Mitte=INT((Iko+Iku)/2)
         IF (X.lt.Xk(Mitte,ba))THEN
           Iko=Mitte
         ELSE
           Iku=Mitte
         END IF
       go to 10
!
!     debug subchk
      END SUBROUTINE
!

!     ###################################################################
      SUBroutine RINDE(Hhrel,Kw,Ri,Hsga,Zo)
!     ###################################################################

!     Ausgabe der Rindenst?rken ..................
!     Hsg=1,2,3 ...nach Heilbronner Sortierung ...
!                  Hsg gibt Funktions-Nr an (min.Aush-L?nge=1,max.Kl-L?nge=2
!                                         max.Draufholz=3)
!                  Baumart mu? stimmen (Fi,Ta,Dgl)  sonst Hsg=0 !
!                  Zo=1 Rindenfunktion für Zopfdurchmesser
!                  Zo=0    "            "  Mittendurchmesser
!     Hsg=0 ..."normale Rindenfunktion"
!               Algorithmus für Hoehenlagenbest.
!
!
      implicit logical (a-z)
        INTEGER Hsg,Zo,hsga
        INTEGER Ba,Iba,Ban(36,7)
        REAL Hhrel,Ri
        REAL Kw
        REAL Rin(1:28,1:4,1:3),Rinh(1:3,1:5,1:3)
        COMMON /Baum/ Ba,Iba,Ban
        COMMON /Rind/Rin,Rinh
        INTEGER Bi,I
        hsg=hsga
!
        IF(zo.gt.0.and.ban(iba,1).gt.3) hsg=0
!
        IF (Hsg.eq.0) THEN

! ...<<17.10.01>> :...........................................................................

!          I=3-INT((.8999-hhrel)/.3)
!          IF(I.eq.0) I=1
!          IF(I.gt.3) I=3

       if (hhrel<=0.4) then
        i=3
       elseif(hhrel<=0.7) then
        i=2
       else
        i=1
       end if

! ...<<17.10.01>> :...........................................................................

          Bi=Ban(Iba,2)
          Ri=(Rin(Bi,i,1)+Rin(Bi,i,2)*Kw+Rin(Bi,i,3)*Kw*Kw)*.1
        ELSE
          IF(Zo.gt.0) THEN
            I=Hsg
            Bi=Ban(Iba,1)
            Ri=(Rinh(Bi,I,1)+Rinh(Bi,I,2)*Kw+Rinh(Bi,I,3)*Kw*Kw)*.1
          ELSE
            Bi=Ban(Iba,2)
            I=Hsg
            Ri=(Rin(Bi,I,1)+Rin(Bi,I,2)*Kw+Rin(Bi,I,3)*Kw*Kw)*.1
          END IF
        END IF
        Kw=Kw-Ri
      RETURN
      END SUBroutine


!     ###################################################################
      REAL FUNCTION Rund(kw)
!     ###################################################################

      REAL Kw
      IF (Kw.ge.20) THEN
         Rund=Kw-0.75
      ELSE
         Rund=Kw-0.5
      END IF
	  IF (Rund.lt.0.0) then !CV 05.06.2025
		 Rund = 0
	  END IF
      RETURN
!     debug subchk
      END FUNCTION


!     ###################################################################
      SUBroutine Sorthl
!     ###################################################################

      implicit logical (a-z)
      REAL D095,D07,Du,Do,Dnorm,H,Hu,
     *hurel,Ho,Horel
      REAL Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      REAL Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      REAL sthh,Stxu,Hakt,Rund,Lxu,Ri,Met
! add integer litmp variable for DO CONTINUE, 10.6.2018 cv
      INTEGER Rif,Ikl,Lind,Sokz, Litmp
      REAL Hhrel,Kw,Li,Hsth,Lsth,Zdrh,Hdrh,Labsch,Luv
      REAL Dhs(1:6),Hhs(1:6)
      INTEGER Kl,Ukl,Klx,Uklx,Kla,Ukla
        REAL glLSort(1:5), glDSort(1:5) ! christian vonderach, 24.07.2018

! ...<<18.09.03>> : Aenderung :.......................................................

      REAL   HStockEnde
      REAL   HStHAnfang
      REAL   LngStH
      REAL   HStHLzEnde
      REAL   HBDATGes

      COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
     1     , HBDATGes
      COMMON /glLDSort/ glLSort, glDSort ! christian vonderach, 24.07.2018

! ...<<18.09.03>> : Aenderung :.......................................................
!
      COMMON /Klasse/ Kl,Ukl,Klx,Uklx,Kla,Ukla
      COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
      COMMON/It/Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      COMMON /Volum/ Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      COMMON /Wert1/ sthh,Sokz,Stxu

!     Erlaeuterung der lokalen Variablen :
!
!        Hakt     R.... Hoehe der aktuellen Sortierungsstelle
!        Met      R.... Rundungsdivisor
!        Lxu      R.... Laenge X-Holz
!        Dhs(6)   R.....Durchmesser der H-Sortierklassen
!        Hhs(6)   R.....Hoehe der H-Sortierklassen
!        Rif      I.....Rindenfunktion Nr
!        Lind     I.....Laenge Industrie-Holz
!        Ikl      I.....Klasse
!
!        Hsthh    R.....Hoehe Stammholz
!        Lsthh    R.....Laenge    "
!        Zdrh     R.....Zaehler Draufholz (Meterschritte)
!        Hdrh     R.....Hoehe       "
!        Hstzop   R.....Hoehe Stammholzzopf
!        Habzop   R.....Hoehe Abschnittzopf
!        Luv      R.....Laenge unverwertetes Derbholz
!
!
!
!
!..................Durchm. HS:
      DATA Dhs/10.,12.,14.,17.,22.,30./
!..................L?ngen HS:
      DATA Hhs/8.,10.,14.,16.,18.,18./
!       Klasse 1   2   3   4   5   6
!
! ........... Met  1= Sortierung auf 0.10 m
!                10      "         1 m
!
      Met=1
      IF(Sokz.eq.4)Met=10
!
!
      sthh=hazop
!
!................................Klasse u Vol.des X-Holzes-unten
      IF((stxu+h*.01).gt.hazop)stxu=hazop-h*.01 !cv 06.06.2025
      Hakt=Stxu+H*.01
      Lxu=Stxu
      IF(Lxu.gt.3) Lxu=int(Lxu*10.0001)*.1
      Rif=1
      IF(Lxu.gt.0)then
         Hhrel=1.-(Lxu*.5+H*.01)/H
         CALL Kuwert(Hhrel,Kw)
         CALL Rinde(Hhrel,Kw,Ri,0,0)
         Kw=Rund(Kw)
         Klx=INT(Kw*.1)
         Uklx=INT((Kw-Klx*10)*.2)
         hakt=h*.01+lxu
         Volx=Lxu*Kw*Kw*PI
        glLSort(1)=Lxu ! christian vonderach 24.07.2018
        glDSort(1)=Kw ! christian vonderach 24.07.2018
         Rif=2
        ELSE
        glLSort(1)=0 ! christian vonderach 24.07.2018
        glDSort(1)=0 ! christian vonderach 24.07.2018
      END IF
!     write(10,*)'Hl-  volx,kw,lxu,hakt',volx,kw,lxu,hakt

!...................................... H - Sortierung

      DO 100 Ikl=6 ,1,-1
       if(sthh.ge.hhs(ikl)+hakt)then
        Hhrel=1.-(Hhs(Ikl)+hakt)/H
        if(hhrel.gt.0) then
          CALL Kuwert(Hhrel,Kw)
          CALL Rinde(Hhrel,Kw,Ri,5,1)
          kw=rund(kw)
!     write(10,*)'kw,dhs(ikl),hhs(ikl),ikl,sthh,hhrel'
!     write(10,*)kw,dhs(ikl),hhs(ikl),ikl,sthh,hhrel
          if(Kw.ge.Dhs(Ikl))then
            Kl=Ikl
            GOTO 180
          END IF
        END IF
       END IF
 100  CONTINUE
!..............<<< Holz zu schwach ... kein Stammholz
      GOTO 250
!     ......max. Klassenlaenge
 180  CONTINUE
! ** DO-CONTINUE mit ganzzahligem Ausdruck Litmp implementiert. 10.06.2018 cv
!      DO 120 Li=Hhs(Kl)+hakt+1,sthh,1
      Li=INT(Hhs(Kl)+hakt+1)
      DO 120 Litmp=INT(Hhs(Kl)+hakt+1),INT(sthh),1
        Li=Litmp
        Hhrel=1.-Li/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,5,1)
        Kw=rund(kw)
        IF(Kw.lt.Dhs(Kl))goto 190
 120  CONTINUE
 190  CONTINUE
      Hsth=Li-1
      IF(Hsth.gt.Hhs(kl)+Hakt+.1)Rif=2
      IF(Sokz.eq.4)THEN
! ** DO-CONTINUE mit ganzzahligem Ausdruck Litmp implementiert. 10.06.2018 cv
!         DO 121 Li=Hsth+.1, sthh,.1
         DO 121 Litmp=INT((Hsth+.1)*10), INT(sthh*10),1
           Li=Litmp/10
           Hhrel=1.-Li/H
           CALL Kuwert(Hhrel,Kw)
           CALL Rinde(Hhrel,Kw,Ri,5,1)
           Kw=rund(kw)
           IF(Kw.lt.Dhs(Kl))goto 191
 121     CONTINUE
 191     hsth=li-.1
      END IF
!.......................................max.Draufholz
      Zdrh=0
      IF(Kl.gt.1)THEN
! ** DO-CONTINUE mit ganzzahligem Ausdruck Litmp implementiert. 10.06.2018 cv
        DO 130 Litmp=INT(Hsth+1) , INT(sthh), 1
          Li=Litmp
          Hhrel=1.-Li/H
          CALL Kuwert(Hhrel,Kw)
          CALL Rinde(Hhrel,Kw,Ri,5,1)
          Kw=Rund(Kw)
          IF(Kw.lt.Dhs(kl-1))go to 200
          Zdrh=Zdrh+1
  130   CONTINUE
  200   CONTINUE
        Hdrh=Li-1
        IF(Hdrh.gt.Hsth+.1)Rif=3
        IF(Kl.gt.4.AND.Zdrh.ge.4) THEN
           Hdrh=Hsth+4
        ELSE
           IF(Sokz.eq.4)THEN
! ** DO-CONTINUE mit ganzzahligem Ausdruck Litmp implementiert. 10.06.2018 cv
!              DO 131 Li=Hdrh+.1,sthh,.1
              DO 131 Litmp=INT((Hdrh+.1)*10) , INT(sthh*10),1
                Li=Litmp/10
                Hhrel=1.-Li/H
                CALL Kuwert(Hhrel,Kw)
                CALL Rinde(Hhrel,Kw,Ri,5,1)
                Kw=Rund(Kw)
                IF(Kw.lt.Dhs(kl-1))go to 201
  131         CONTINUE
  201         Hdrh=Li-.1
           END IF
        END IF
      ELSE
        Hdrh=Hsth
      END IF
      Lsth=Hdrh-hakt

! ...<<18.09.03>> : Aenderung :.......................................................

      LngStH = LStH
      HStHAnfang = Hakt

! ...<<18.09.03>> : Aenderung :.......................................................

!
!     write(10,*)'hakt lsth',hakt,lsth
!....................... Stammholzvolumen:
      if(lsth.gt.0)then
        Hhrel=1.-(hakt+Lsth*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,Rif,1)
!       write(10,*) 'rif' ,rif ,kw
        Kw=Rund(Kw)
!       write (10,*) kw
        Hakt=Hakt+Lsth*1.01
        IF(Sokz.eq.4)THEN
          Kl=INT(Kw*.1)
          Ukl=INT((Kw-Kl*10)*.2)
        END IF
        Vols=Kw*Kw*Lsth*PI
       glLSort(2)=Lsth ! christian vonderach 24.07.2018
       glDSort(2)=Kw ! christian vonderach 24.07.2018
!     write(10,*)'Hl-  vols,kw,lsth,hakt',vols,kw,lsth,hakt
!     write(10,*)'hhrel,ri',hhrel,ri,rif
        else
          glLSort(2)=0 ! christian vonderach 24.07.2018
       glDSort(2)=0 ! christian vonderach 24.07.2018
      end if
!
!     ...................zus?tzlich Kronenst?ck nach Mittenst?rkens.:
 250  continue

! ...<<18.09.03>> : Aenderung :.......................................................

      HStHLzEnde = MIN(HAkt,Hbr)

! ...<<16.06.05> : Aenderung :------------------------------------------------------------

      if (sthh.gt.0) then
       HStHLzEnde = MIN(HStHLzEnde,sthh)
      else
      end if

! ...<<18.09.03>> : Aenderung :.......................................................

      Labsch=INT((habzop+.0001-Hakt)*Met)/Met
      IF(Labsch.ge.3) THEN
        Hhrel=1.-(Hakt+Labsch*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(Kw)
        Kla=INT(Kw*.1)
        Ukla=INT((Kw-Kla*10)*.2)

!       Volumenabsch.

        Vola=Labsch*Kw*Kw*PI
       glLSort(3)=Labsch ! christian vonderach 24.07.2018
       glDSort(3)=Kw ! christian vonderach 24.07.2018
        Hakt=Hakt+Labsch*1.01
!       write(10,*)'Hl-  vola,kw,labsch,hakt',vola,kw,labsch,hakt
!       write (10,*)'ri hhrel',ri,hhrel
        ELSE
          glLSort(3)=0 ! christian vonderach 24.07.2018
       glDSort(3)=0 ! christian vonderach 24.07.2018
      END IF
!..................................Industrieholz:
      Voli=0
      IF(Hazop.gt.Hakt+.999)THEN
        Lind=int(Hazop+.0001-Hakt)
        Hhrel=1.-(Hakt+Lind*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(kw)
        Voli=Lind*Kw*Kw*PI
       glLSort(4)=Lind ! christian vonderach 24.07.2018
       glDSort(4)=Kw ! christian vonderach 24.07.2018
        Hakt=Lind+Hakt
!     write(10,*)'Hl-  voli,kw,lind,hakt',voli,kw,lind,hakt
        ELSE
          glLSort(4)=0 ! christian vonderach 24.07.2018
       glDSort(4)=0 ! christian vonderach 24.07.2018
      END IF
!..................................... UVDerbholz:
      Luv=Hdgren-Hakt
      IF(Luv.gt.0)THEN
        Hhrel=1.-(Hakt+Luv*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(Kw)
        Volu=Luv*Kw*Kw*PI
       glLSort(5)=Luv ! christian vonderach 24.07.2018
       glDSort(5)=Kw ! christian vonderach 24.07.2018
!     write(10,*)'Hl-  volu,kw,luv,hakt',volu,kw,luv,hakt
        ELSE
          glLSort(5)=0 ! christian vonderach 24.07.2018
       glDSort(5)=0 ! christian vonderach 24.07.2018
      END IF
      RETURN
!     debug subchk
      END SUBroutine


!     ###################################################################
      SUBroutine Sortmi
!     ###################################################################

      implicit logical (a-z)
      REAL D095,D07,Du,Do,Dnorm,H,
     *hu,hurel,Ho,Horel
      REAL Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      REAL Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      REAL sthh,Stxu,Luvd,Rund,Hhrel
      REAL Kw,Lxu,Lsth,Hakt,Labsch,Ri,Met
      INTEGER Ba,Iba,Ban(36,7),Lind,Sokz
      INTEGER Kl,Ukl,Klx,Uklx,Kla,Ukla
        REAL glLSort(1:5), glDSort(1:5) ! christian vonderach, 24.07.2018

! ...<<18.09.03>> : Aenderung :.......................................................

      REAL   HStockEnde
      REAL   HStHAnfang
      REAL   LngStH
      REAL   HStHLzEnde
      REAL   HBDATGes

      COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
     1     , HBDATGes
      COMMON /glLDSort/ glLSort, glDSort ! christian vonderach, 24.07.2018

! ...<<18.09.03>> : Aenderung :.......................................................

      COMMON /Klasse/ Kl,Ukl,Klx,Uklx,Kla,Ukla
!
      COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
      COMMON/It/Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      COMMON /Volum/ Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      COMMON /Wert1/ sthh,Sokz,Stxu
      COMMON /Baum/ Ba,Iba,Ban
!

!       Erlaeuterung der lokalen Variablen :
!
!        Hakt     R.... Hoehe der aktuellen Sortierungsstelle
!        Met      R.... Rundungsdivisor
!        Lxu      R.... L?nge X-Holz
!        Rif      I.....Rindenfunktion Nr
!        Lind     I.....L?nge Industrie-Holz
!        Ikl      I.....Klasse
!
!        Hsthh    R.....Hoehe Stammholz
!        Lsthh    R.....L?nge    "
!        Hdrh     R.....Hoehe       "
!        Hstzop   R.....Hoehe Stammholzzopf
!        Habzop   R.....Hoehe Abschnittzopf
!        Labsch   R.....L?nge Abschnitt
!        Luv      R.....L?nge unverwertetes Derbholz
!
!
!
!     .....Rundung auf 10cm
      Met=10
!................................Klasse u Vol.des X-Holzes-unten
      IF((stxu+h*.01).gt.hazop)stxu=hazop-h*.01 !cv 06.06.2025
      Hakt=Stxu+H*.01
      Lxu=Stxu
      IF(Lxu.gt.3) Lxu=int(Lxu*10.0001)*.1
      IF(Lxu.gt.0)THEN
         Hhrel=1.-(Lxu*.5+H*.01)/H
         CALL Kuwert(Hhrel,Kw)
         CALL Rinde(Hhrel,Kw,Ri,0,0)
         Kw=Rund(Kw)
         Klx=INT(Kw*.1)
         Uklx=INT((Kw-Klx*10)*.2)
         hakt=h*.01+lxu
         Volx=Kw*Kw*Lxu*PI
        glLSort(1)=Lxu ! christian vonderach 24.07.2018
        glDSort(1)=Kw ! christian vonderach 24.07.2018
        ELSE
        glLSort(1)=0 ! christian vonderach 24.07.2018
        glDSort(1)=0 ! christian vonderach 24.07.2018
      END IF
!     write(10,*)'Mi-  volx,kw,lxu,hakt',volx,kw,lxu,hakt
! ....Stammholz.................................................



      Lsth=INT((hstzop+.0001-Hakt)*Met)/Met

! Aenderung <<16.06.04>> :------------------------------------------------------------------

! (Trennschnitt 20m )

      if (Lsth>20) then
       Lsth=20
      end if

      Hstzop=Lsth+hakt

! ...<<18.09.03>> : Aenderung :.......................................................

      LngStH = LStH
      HStHAnfang = Hakt

! ...<<18.09.03>> : Aenderung :.......................................................


      IF(Lsth.ge.3)THEN
        Hhrel=1.-(hakt+Lsth*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(Kw)
        Kl=INT(Kw*.1)
        Ukl=INT((Kw-Kl*10)*.2)
        Hakt=Hakt+Lsth*1.01
! Stammholzvolumen
        Vols=Lsth*Kw*Kw*PI
       glLSort(2)=Lsth ! christian vonderach 24.07.2018
       glDSort(2)=Kw ! christian vonderach 24.07.2018
!     write(10,*)'Mi-  vols,kw,lsth,hakt',vols,kw,lsth,hakt
        ELSE
          glLSort(2)=0 ! christian vonderach 24.07.2018
       glDSort(2)=0 ! christian vonderach 24.07.2018
      END IF

! ...<<18.09.03>> : Aenderung :.......................................................

      HStHLzEnde = MIN(HAkt,Hbr)

! ...<<16.06.05> : Aenderung :------------------------------------------------------------

      if (sthh.gt.0) then
       HStHLzEnde = MIN(HStHLzEnde,sthh)
      else
      end if

! ...<<18.09.03>> : Aenderung :.......................................................

! zus?tzliches  Kronenstammholzst?ck :
! Abschnitte:.................................................


      Labsch=INT((Habzop+.0001-Hakt)*Met)/Met

! Aenderung <<16.06.04>> :------------------------------------------------------------------

! (Trennschnitt 20 m)

      if (Labsch>20) then
       Labsch=20
      end if

      Habzop=Hakt+Labsch
      IF(Labsch.ge.3) THEN
        Hhrel=1.-(Hakt+Labsch*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(Kw)
        Kla=INT(Kw*.1)
        Ukla=INT((Kw-Kla*10)*.2)
!   Volumenabsch:
        Vola=Labsch*Kw*Kw*PI
        Hakt=Habzop+Labsch*.01
       glLSort(3)=Labsch ! christian vonderach 24.07.2018
       glDSort(3)=Kw ! christian vonderach 24.07.2018
!     write(10,*)'Mi-  vola,kw,labsch,hakt',vola,kw,labsch,hakt
        ELSE
          glLSort(3)=0 ! christian vonderach 24.07.2018
       glDSort(3)=0 ! christian vonderach 24.07.2018
      END IF
! Industrieholz ...........................................
      Voli=0
      Lind=INT(Hazop+.0001-Hakt)
      IF(Lind.gt.0.999)THEN
        Hhrel=1.-(Hakt+(Lind*.5))/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(Kw)
        Voli=Lind*Kw*Kw*PI
        Hakt=Hakt+Lind
       glLSort(4)=Lind ! christian vonderach 24.07.2018
       glDSort(4)=Kw ! christian vonderach 24.07.2018
!     write(10,*)'Mi-  voli,kw,lind,hakt',voli,kw,lind,hakt
        ELSE
          glLSort(4)=0 ! christian vonderach 24.07.2018
       glDSort(4)=0 ! christian vonderach 24.07.2018
      END IF
! Unverwertbares Derbholz .................................
      Volu=0
      IF(Hakt.lt.Hdgren) THEN
        Luvd=Hdgren-Hakt
        Hhrel=1.-(Hakt+Luvd*.5)/H
        CALL Kuwert(Hhrel,Kw)
        CALL Rinde(Hhrel,Kw,Ri,0,0)
        Kw=Rund(Kw)
        Volu=(Hdgren-Hakt)*Kw*Kw*PI
       glLSort(5)=Luvd ! christian vonderach 24.07.2018
       glDSort(5)=Kw ! christian vonderach 24.07.2018
!     write(10,*)'Mi-  volu,kw,hdgren-hakt',volu,kw,hdgren-hakt
        ELSE
          glLSort(5)=0 ! christian vonderach 24.07.2018
       glDSort(5)=0 ! christian vonderach 24.07.2018
      END IF
      RETURN
!     debug subchk
      END SUBroutine


!     ###################################################################
      SUBroutine Iter
!     ###################################################################

      implicit logical (a-z)
      REAL D095,D07,Du,Do,Dnorm,H,
     *hu,hurel,Ho,Horel
      REAL Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      REAL sthh,Stxu,Hhrel,stx
      REAL Kw,Ri,rest
      REAL Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      INTEGER Ba,Iba,Ban(36,7)
      INTEGER I,J,Sokz,Vorz,Anzit,anfit
      REAL Wit(1:4),Hwit(1:4),Dis,Po,rund,bis
      INTEGER Ix
!
      COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
      COMMON/It/Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      COMMON /Wert1/ sthh,Sokz,Stxu
      COMMON /Volum/ Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      COMMON /Baum/ Ba,Iba,Ban

!
!       Erlaeuterung der lokalen Variablen :
!         Hwit(4)       R.....Hoehe des iterierten Durchmessers
!                            1 Hoehe der Derbholzgrenze
!                            2  "       Aufarbeitungszopfs
!                            3  "       Abschnittszopfs
!                            4  "       Stammholzzopfs
!         Wit (4)       R     entsprechende Durchmesser
!         Anfit         I     Anfangsindex der Iteration
!         anzit         I     Anzahl der Interationen (step=1)
!         Po            R     aktuelle IterationsHoehe (Position)
!         Dis           R     Iterationsschrittweite
!         Vorz          I     (Vorzeichen) Iterationsrichtung
!         Rest          R     Hilfsgroesse
!         Stx           R     Hilfsgroesse für Stx+H*0.01(=Stock)
!

      do 1 i=1,4
         hwit(i)=0
  1   continue
!
!..BHD - Zopfklasse
      IF(Sokz.eq.0) then
         anfit=1
         anzit=1
      end if


      IF(Sokz.eq.3) then

         anfit=1

! ...<<17.09.03>> : Aenderung :............................................................

         anzit=4
!        anzit=2

      end if
      IF(Sokz.eq.1) then
         anfit=1
         anzit=4
      end if
      IF(Sokz.eq.2.or.Sokz.eq.4) then
         anfit=1
         anzit=3
      end if
      Wit(1)=Dgrenz
      Wit(2)=Azop
      if(Zoab.eq.0)then
        Wit(3)=14.
      else
        Wit(3)=Zoab
      end if
! Stammholzzopf
      Wit(4)=Zost
!
      DO 160 J=anfit,Anzit
        DO 140 I=INT(Hbr) ,1, -1
          Hhrel=1.-I/H
          CALL Kuwert(Hhrel,Kw)
          IF(j.gt.2) then

!  Durchmesser ohne Rinde forstlich gerundet :...................................

            CALL Rinde(Hhrel,Kw,Ri,0,0)
            kw=rund(kw)
          END IF
          IF(Kw.gt.Wit(J))goto  220
 140    CONTINUE
        Hwit(J)=1
 220    CONTINUE
!
        IF(Hwit(J).lt.1) THEN
          Po=I+.5
          Dis=.5
          DO 150 I=1 ,18
            Hhrel=1.-Po/H
            CALL Kuwert(Hhrel,Kw)
            if(j.gt.2) then
               CALL Rinde(Hhrel,Kw,Ri,0,0)
               kw=rund(kw)
            END IF
            IF(Kw.lt.Wit(J)+.1.AND.Kw.gt.Wit(J)-.1) goto 230
            Vorz=1
            IF(Kw.lt.Wit(J)) Vorz=(-1)
            Dis=ABS(Dis)*.5*Vorz
            Po=Po+Dis
 150      CONTINUE
 230   CONTINUE
       if(hbr.lt.po)then
         hwit(j)=hbr
       else
         Hwit(j)=po
       END IF
!
      END IF

!     write (10,*)'hwit(j),wit(j),azop,j'
!     write (10,*) hwit(j),wit(j),azop,j

 160  CONTINUE
!
!     write (10,*)'hdgren,hazop,habzop,hstzop'
!     write (10,*) hdgren,hazop,habzop,hstzop

      if(hazop.le.0.or.hazop.gt.hwit(2)) Hazop=Hwit(2)
      Hdgren=Hwit(1)
      if (habzop.le.0.or.habzop.gt.hwit(3))Habzop=Hwit(3)
      IF(Habzop.gt.Hazop) Habzop=Hazop
      if(hstzop.le.0.or.hstzop.gt.hwit(4))Hstzop=Hwit(4)
      if(hstzop.le.0.or.hstzop.gt.sthh)hstzop=sthh
      IF(hstzop.gt.hazop)hstzop=hazop
!........................................AbschnittzopfHoehe (1.2.94)
      if(habzop.le.0.or.habzop.gt.sthh)habzop=sthh
      if(habzop.le.0)habzop=hstzop
      IF(Hbr.lt.Hdgren)Hdgren=Hbr
      IF((stxu+h*.01).gt.hazop)stxu=hazop-h*.01

!     write (10,*)'hdgren,hazop,habzop,hstzop,stxu'
!     write (10,*) hdgren,hazop,habzop,hstzop,stxu
!
      Volr=0
      Vol=0
!..............................Volumenberechnung:   Vol Vfm.m.R.
      if (hdgren.le.2)then
        Hhrel=1.-Hdgren*.5/H
        CALL Kuwert(Hhrel,Kw)
        Volr= Kw*Kw*PI*Hdgren
      ELSE
        DO 240 Ix=1,INT(Hdgren)-1,2
          Hhrel=1.-Ix/H
          CALL Kuwert(Hhrel,Kw)
          Volr=Volr+Kw*Kw*PI*2

!         write (10,*)'volr,kw,hhrel'
!         write (10,*) volr,kw,hhrel

 240    CONTINUE
        rest=hdgren-(ix-1)
        Hhrel=1.-(Hdgren-rest*.5)/H
        CALL Kuwert(Hhrel,Kw)
        Volr=Volr+Kw*Kw*PI*rest
      END IF
!
      IF(Iba.gt.14)THEN
          stx=stxu+h*.01
          bis=stx
          IF(hstzop-stx.ge.2) bis=(int(hstzop*10.0001)*.1)*1.01
          if(stx.ge.bis)then
             hstzop=bis*.9999
             bis=stx
             IF(bis.ge.3)bis=(int(bis*10.0001)*.1)*1.01
          END IF
          vol=0

!          write(10,*)'hstzop,hazop ',hstzop,hazop

          DO 250 Ix=1, INT(bis)-1 , 2
            Hhrel=1.-Ix/H
            CALL Kuwert(Hhrel,Kw)
            Vol=Vol+Kw*Kw*PI*2
 250      CONTINUE
          rest=bis-(ix-1)

!  write(10,*)'bis,rest,ix'
!  write(10,*) bis,rest,ix

          Hhrel=1.-(bis-rest*.5)/H

          CALL Kuwert(Hhrel,Kw)

          Vol=Vol+Kw*Kw*Pi*rest

!  Aenderungen: <<01.07.02>> :--------------------------------------------------------------
!
!  SUBroutine Iter: IF(H.lt.22) Volr=Volr*.98
!         IF(H.lt.22)Volr=Volr*.98

      END IF

!     write(10,*)'bis,rest,kw,vol'
!     write(10,*) bis,rest,kw,vol
!     debug subchk


      END SUBroutine


!     ###################################################################
      SUBroutine Sortlb
!     ###################################################################

      implicit logical (a-z)

           REAL  D095,D07,Du,Do,Dnorm,H,Volir,hu,hurel,ho,horel
           REAL  zopf
           REAL  Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost
           REAL  Zoab,Hbr
           REAL  Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
           REAL  Volrer,Volur
          REAL glLSort(1:5), glDSort(1:5) ! christian vonderach, 24.07.2018

! <<1>> ################################### 01.04.04 #################################
! REAL Volrer,Volur
! ######################################### 01.04.04 #################################

           REAL  sthh,Stxu,Hhrel,Lxu,Kw,Hakt,Lsth,W1,W0,Vl
! ####cv: function Rund returns REAL NOT double precision==double precision ############
      REAL Rund

!     INTEGER Ba,Iba,Ban(36,7),Ib,Ia,Sokz,ix

      INTEGER Ba,Iba,Ban(36,7),Ib,Ia,Sokz
           REAL   Unvd(2,27,33),Unvpr,Ddur,Riproz,Ri,Dduru,Dduri
           REAL  D1,fakt
      INTEGER Kl,Ukl,Klx,Uklx,Kla,Ukla

! ...<<18.09.03>> : Aenderung :.......................................................

      REAL   HStockEnde
      REAL   HStHAnfang
      REAL   LngStH
      REAL   HStHLzEnde
      REAL   HBDATGes

      COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
     1     , HBDATGes
      COMMON /glLDSort/ glLSort, glDSort ! christian vonderach, 24.07.2018

! ...<<18.09.03>> : Aenderung :.......................................................

      COMMON /Klasse/ Kl,Ukl,Klx,Uklx,Kla,Ukla
!
      COMMON /Baum1/ D095,D07,Du,Do,Dnorm,H,Hu,Hurel,Ho,Horel
      COMMON/It/Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      COMMON /Volum/ Vol,Volx,Vols,Vola,Voli,Volu,Volr,Volsu
      COMMON /Wert1/ sthh,Sokz,Stxu
      COMMON /Baum/ Ba,Iba,Ban
      COMMON /Uvd/ Unvd

! ---------------------------------------------
!
!       Erlaeuterung der lokalen Variablen :
!
!        Hakt     R.... Hoehe der aktuellen Sortierungsstelle
!        Lxu      R.... Laenge X-Holz
!        Lind     I.....Laenge Industrie-Holz
!        Ikl      I.....Klasse
!
!        Hsthh    R.....Hoehe Stammholz
!        Lsthh    R.....L?nge    "
!        Zdrh     R.....Z?hler Draufholz (Meterschritte)
!        Hdrh     R.....Hoehe       "
!        Hstzop   R.....Hoehe Stammholzzopf
!        Luv      R.....L?nge unverwertetes Derbholz

!
!...................................Klasse u Vol.des X-Holzes-unten:

      VolRer=0
      VoluR=0
      IF((stxu+h*.01).gt.hazop)stxu=hazop-h*.01 !cv 06.06.2025
      Hakt=Stxu+H*.01
      Lxu=Stxu

      IF(Lxu.ge.3) lxu=int(lxu*10.0001)*.1

      IF(Lxu.gt.0)THEN

         Hhrel=1.-(Lxu*.5+H*.01)/H
         CALL Kuwert(Hhrel,Kw)
         CALL Rinde(Hhrel,Kw,Ri,0,0)
         Kw=Rund(Kw)
         Klx=INT(Kw*.1)
         Uklx=INT((Kw-Klx*10)*.2)
         hakt=h*.01+ lxu
         Volx=Kw*Kw*Lxu*PI
        glLSort(1)=Lxu ! christian vonderach 24.07.2018
        glDSort(1)=Kw ! christian vonderach 24.07.2018
        ELSE
           glLSort(1)=0 ! christian vonderach 24.07.2018
        glDSort(1)=0 ! christian vonderach 24.07.2018
      END IF

! ...<<18.09.03>> : Aenderung :.......................................................

      HStHAnfang = Hakt

! ...<<18.09.03>> : Aenderung :.......................................................

!     write(10,*)'b- volr,hazop,hstzop',volr,hazop,hstzop
!     write(10,*)'Lb-  volx,kw,lxu,hakt',volx,kw,lxu,hakt

!.............................................Stammholz:

      if(hstzop.le.0) hstzop=stxu+h*.01
      Vols=0

      IF(hstzop.gt.0) THEN

         lsth=int((hstzop+.0001-hakt)*10)*.1

! Aenderung <<16.06.04>> :------------------------------------------------------------------

       if (Lsth>20) then
        Lsth=20
       end if

         IF(Lsth.ge.2) THEN

            Hhrel=1.-(hakt+Lsth*.5)/H

            CALL Kuwert(Hhrel,Kw)
            CALL Rinde(Hhrel,Kw,Ri,0,0)

            kw=Rund(Kw)
            Kl=INT(Kw*.1)
            Ukl=INT((Kw-Kl*10)*.2)

!......................... Stammholzvolumen

            Vols=Lsth*Kw*Kw*PI
        glLSort(2)=Lsth ! christian vonderach 24.07.2018
        glDSort(2)=Kw ! christian vonderach 24.07.2018

!............Stammholzzopf....................................Aenderung 1.2.94

            Hstzop=hakt+Lsth
            Hakt=Hakt+Lsth*1.01

         ELSE
            Hstzop=hakt
        glLSort(2)=0 ! christian vonderach 24.07.2018
        glDSort(2)=0 ! christian vonderach 24.07.2018
         END IF
        LngStH = LStH
      END IF

!     write(10,*)'Lb-  vols,kw,lsth,hakt',vols,kw,lsth,hakt
      Vola=0

! ...<<18.09.03>> : Aenderung :.......................................................

      HStHLzEnde = MIN(HAkt,Hbr)

! ...<<16.06.05> : Aenderung :------------------------------------------------------------

      if (sthh.gt.0) then
       HStHLzEnde = MIN(HStHLzEnde,sthh)
      else
      end if

! ...<<18.09.03>> : Aenderung :.......................................................

!.............................. Industrieholz:

! Rest-Holz:

! ---------------------------------------------
!   write (*,*)
! write (*,*) " BDAT ! SORTLB Block3 "
!  write (*,*) " BDAT ! SORTLB VOL   : ", VOL
! write (*,*) " BDAT ! SORTLB VOLx  : ", VOLx
! write (*,*) " BDAT ! SORTLB VOLs  : ", VOLs
! write (*,*) " BDAT ! SORTLB VOLa  : ", VOLa
! write (*,*) " BDAT ! SORTLB VOLi  : ", VOLi
! write (*,*) " BDAT ! SORTLB VOLu  : ", VOLu
! write (*,*) " BDAT ! SORTLB VOLr  : ", VOLr
! write (*,*) " BDAT ! SORTLB VOLuR : ", VOLuR
! write (*,*) " BDAT ! SORTLB VOLsu : ", VOLsu
! write (*,*) " BDAT ! SORTLB VOLrer: ", VOLrer
!   write (*,*)
! ---------------------------------------------

      volrer=volr-vol

! ---------------------------------------------

! <<2>> #################################### 01.04.04 ################################
!   write (*,*)
! ######################################### 01.04.04 #################################

!     write(10,*)'Lb- volr,vol,volrer',volr,vol,volrer
! Unverwertbares Derbholz mit Rinde

      volu=0
      volur=0

      IF(azop.gt.8)THEN

        Zopf=azop
        if(zopf.gt.40)zopf=40

!  geschaetzt aus  Aufarbeitungsgrenze und BHD

        Hhrel=1.-1.3/H
        CALL Kuwert(Hhrel,Kw)
        IF(Kw.gt.60)Kw=60
        Ia=INT(Zopf+.0001-7)
        Ib=INT((Kw-8)*.5)+1
        D1=Kw-(Ib*2+6)
        Ba=Ban(Iba,4)
        W1=Unvd(Ba,Ib,Ia)
        IF(Ib.gt.26)then
          W0=W1
        ELSE
          W0=Unvd(Ba,Ib+1,Ia)
        END IF

        Unvpr=(W1-W0)*(2-D1)/2+W0

!       write(10,*)'Unvpr=(W1-W0)*(2-D1)/2+W0,ib,ia,ba'
!       write(10,*)Unvpr,W1,W0,2,D1,2,W0,ib,ia,ba

        fakt=volr/volrer*unvpr*.01

        if(fakt.gt.1) fakt=1.0

        volur=volrer*fakt

! ######################################### 01.04.04 #################################
!     volur=volrer
! ######################################### 01.04.04 #################################

!  write(10,*)'Lb-  volu,volr,unvpr',volu,volr,unvpr,

      END IF

!     Mittlere ?stedurchmesser mit Rinde

      vl=hstzop
      if(vl.lt.h*.1)vl=h*.15
      IF(Ban(Iba,5).eq.1) THEN
        Ddur=8.866-.2721*vl+.2035*Du+.3546*H/vl+1.4139*vl/Du
      ELSE
        Ddur=6.473-.119*vl+.1578*Du+1.4236*H/vl
      END IF

!.....Industrieholz Volumen mit/ohne Rinde.....................Aenderung 1.2.94

      voli=0
      volir=0

! ######################################### 01.04.04 #################################
!    volur=volrer
! ######################################### 01.04.04 #################################

      volir=volrer-volur
      Zopf=0
      hhrel=1.-(hstzop/h)

      CALL kuwert(hhrel,zopf)


! Aenderung <<01.04.04>> :-------------------------------------------------------------

      IF(ABS(volir).gt.1e-8) THEN

! #################### 01.04.04 ############## FALSE bei <<1>> oder <<2>> ############
!   IF(volir.gt.0) THEN
! ######################################### 01.04.04 #################################

        IF(hazop.lt.hdgren) THEN
          dduri=ddur+((volrer-volir)/volrer)*(zopf-ddur)
        ELSE
          Dduri=Ddur
        END IF

        CALL Rinde(Hhrel,Dduri,Ri,0,0)

        Riproz=(Dduri)**2/(Dduri+ri)**2

        Voli=Volir*Riproz

! ######################################### 01.04.04 #################################
!    Voli=-99999
! ######################################### 01.04.04 #################################
!
      END IF

!   write (*,*)
!     write (*,*) '========================='
!   write (*,*) 'BDAT!SortLb:'c
!  write (*,*) 'voli   = ',voli
!  write (*,*) 'volir  = ',volir
!  write (*,*) 'riproz = ', riproz
!  write (*,*) '========================='
!  write (*,*)

!
      Hakt=Hazop

! Unverwertbares Derbholz ohne Rinde...........................Aenderung 1.2.94

      Volu=0
      IF(hazop.lt.hdgren) THEN
         dduru=ddur-((volrer-volur)/volrer)*(ddur-7)
         CALL Rinde(Hhrel,Dduru,Ri,0,0)
         Riproz=Dduru**2/(Dduru+ri)**2
         Volu=Volur*Riproz
!        write(10,*)'Lb-  volu,volur,riproz',volu,volur,riproz
      END IF
!
      RETURN

      END SUBroutine

!***********************************************************************
!  von Bernhard boesch, ergaenzt durch christian vonderach am 09.06.2018
      real FUNCTION ANORDF( x )
          ! This function evaluates the standard normal (Gaussian) cumulative distribution function.
          ! -> http://www.roguewave.com/Portals/0/products/imsl-numerical-libraries/fortran-library/docs/6.0/stat/default.htm?turl=anorin.htm
          !REAL ANORDF
          REAL x ! Argument for which the normal cumulative distribution function is to be evaluated.   (Input)
          double precision xtmp
          !0 CALL ShowNYIMessage('ANORDF')
          xtmp=x/SQRT(2.d0)
          ANORDF= real((1.d0+ERF(xtmp) )/2.d0)
      END FUNCTION ANORDF

      real FUNCTION ANORIN( prob )
      ! This function evaluates the inverse of the standard normal (Gaussian)
      ! cumulative distribution function.
      ! -> http://www.roguewave.com/Portals/0/products/imsl-numerical-libraries/fortran-library/docs/6.0/stat/default.htm?turl=anorin.htm
      !REAL ANORIN
      double precision prob, xtmp


      !0 CALL ShowNYIMessage('ANORIN')
        !xtmp=-ERFCInv( prob * 2. )
      !  ANORIN1 = -ERFCInv( prob * 2. ) * SQRT(2.d0)
        xtmp=1.0d0 - prob*2.d0
        ANORIN = real(-SQRT(2.d0)*ERFInv(xtmp))
        !xtmp1=ERFInv(xtmp)
      END FUNCTION ANORIN



      FUNCTION ERFInv( x )
      ! inverse Error Function
      REAL ERFInv
      double precision x    ! nur fuer x /= 1
      !
      ! mittels Approximation (gemaess email vom 22.12.2011, von Dr. Boesch,
      ! Forstliche Versuchs- und Forschungsanstalt, Wonnhalde 4, 79100 Freiburg;
      ! Tel. 0761 4018 193   Fax: 0761 4018 333
      ! Email: bernhard.boesch@forst.bwl.de<mailto:bernhard.boesch@forst.bwl.de>)
      !
      ! vgl. http://en.wikipedia.org/wiki/Error_function
      ! und
      ! qt18\...\mma7\Approximation_InverseErrorFunction.nb
      !
      double precision ln1x, c2PiA, x1, a
      double precision PI, A0

      PARAMETER  (A0 = 0.14001228868666660600424949138612d0)
      PARAMETER  (A1 = 0.147)
      PARAMETER (PI = 3.1415926535897932384626433832795d0)


      ! This is very accurate in a neighborhood of 0 and a neighborhood of infinity,
      ! and the error is less than 0.00035 for all x.
      !

      ! Using the alternate value a ~ 0.147 reduces the maximum error to about 0.00012.

      ln1x = LOG(1.D0 - x*x)  ! Abk.
      IF ( x < 0.602469 ) THEN   ! siehe: qt18\...\mma7\Approximation_InverseErrorFunction.nb
         a = A1   !Aenderung Boesch vorerst Symmetrie
      ELSE
         a = A1
      END IF
      c2PiA = 2.D0/(PI*a)
      x1 = SQRT( (c2PiA + ln1x/2.D0)**2 - ln1x/a )
      ERFInv = real(SIGN(1.D0,x) * SQRT( x1 - (c2PiA + ln1x/2) ))
      END FUNCTION ERFInv


      FUNCTION ERFCInv( x )
      ! inverse komplementaere Error Function
      REAL ERFCInv
      double precision x    ! nur fuer x < 0
      !
      ! erfcInv(1-z) = erfInv(z)
      ! mit x = 1 - z
      ! -> z = 1 - x
      ! -> erfcInv(x) = erfInv(1 - x)
      ! vgl. http://en.wikipedia.org/wiki/Error_function
      ERFCInv = ERFInv(1. - x)
      END FUNCTION ERFCInv

!
!######################################################################
!
!
      real function dinvnorm(pi)
!      normal inverse translate from a routine written by john herrero
      real pi
      double precision p,p_low,p_high
      double precision a1,a2,a3,a4,a5,a6
      double precision b1,b2,b3,b4,b5
      double precision c1,c2,c3,c4,c5,c6
      double precision d1,d2,d3,d4
      double precision z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low

      p=pi

      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/
     &(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=real(z)
      return
      end
!


!**********************************************************************
!*** FILE: BDAT_C.FOR *** Tabellen BLOCK DATA *** STAND: 13.12.90   ***
!**********************************************************************

      BLOCK DATA

!======================================================================
      REAL Rin(1:28,1:4,1:3)
      REAL Rinh(1:3,1:5,1:3)
      REAL a95(1:24,1:8),b95(1:24,1:8)
      REAL xk95a(7,8),xk95b(7,8)
      INTEGER nxk95a(8),nxk95b(8)
      INTEGER Nnp,Np,Nxk
      REAL B(1:80,1:8),Xk(1:6)
      REAL Yy
      REAL Azo(7,0:3)
      REAL Unvd(2,27,33)
      REAL Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
      INTEGER Ba,Iba,Ban(36,7)
      INTEGER Volk(6,30,3)
      INTEGER md07(20333)
      INTEGER Add07(14,5:45,3)
      REAL hoehr(6,6)
      REAL Sxk07(7,8),Sd07(24,8)
      INTEGER Snxkn7(8)
      INTEGER i, j ! christian vonderach, 08.05.2025
      REAL glLSort(1:5), glDSort(1:5) ! christian vonderach, 23.07.2018
      REAL   HStockEnde
      DATA    HStockEnde /0/
      REAL   HStHAnfang
      DATA    HStHAnfang /0/
      REAL   LngStH
      DATA    LngStH  /0/
      REAL   HStHLzEnde
      DATA    HStHLzEnde /0/
      REAL   HBDATGes
      DATA    HBDATGes /0/
!
      COMMON /XtrComPar/ HStockEnde, HStHAnfang, LngStH, HStHLzEnde
     1     , HBDATGes
      COMMON /D07/ Add07,Md07,Snxkn7,Sxk07,Sd07
      COMMON /Volk/Volk,hoehr
      COMMON/It/Azop,Dgrenz,Hazop,Hdgren,Hstzop,Habzop,PI,Zost,Zoab,Hbr
        COMMON /glLDSort/ glLSort, glDSort ! christian vonderach, 23.07.2018
      COMMON /Uvd/ Unvd
      COMMON /Duazop/ Azo
      COMMON /Baum/ Ba,Iba,Ban
      COMMON /Rind/ Rin,Rinh
      COMMON /d95/nxk95a,nxk95b,xk95a,xk95b,a95,b95
      COMMON /Schaft/  Nnp,Np,Nxk,B,Xk
      COMMON /Sk/ Yy,Durel
      COMMON /itr/ i, j ! christian vonderach, 08.05.2025
!
!.....Splinekoeffizienten für die Schaftkurven..........................
!
        DATA NP/ 80/
        DATA  NXK/6/,XK/0.,.3,.5,.7,.9,1./ , NnP/20/
!
!.....Fichte............................................................
!
      data (B(i,1),i=1,80) /
     *-.663092E-01,-.503199E-01, .663092E-01,-.221933E-01,-.626495E-01,
     *-.125967    ,-.986368E-02, .247511E-01,-.125967    ,-.407789E-01,
     * .247511E-01, .407789E-01,-.407789E-01, .289083    , .407789E-01,
     * .215608    , .450789    , 1.10084    , .539020E-01, 1.73959    ,
     * .114868E-01, .461916E-03,-.114868E-01, .191777E-02, .152734E-02,
     * .569999E-03, .852343E-03,-.741014E-05, .569999E-03,-.431807E-03,
     *-.741014E-05, .431807E-03,-.431807E-03, .115723E-02, .431807E-03,
     *-.144520E-02, .733330E-04,-.238384E-02,-.361299E-03, .960332E-02,
     *-.885276E-02, .101540E-03, .885276E-02,-.157632E-02,-.774196E-03,
     *-.590717E-04,-.700588E-03,-.765393E-04,-.590717E-04, .196816E-03,
     *-.765393E-04,-.196816E-03, .196816E-03,-.728193E-03,-.196816E-03,
     * .270672E-03,-.525189E-03,-.378682E-03, .676680E-04, .354782E-02,
     * .152552    , .726242    ,-.152552    ,-.262736E-01, .711645    ,
     * 1.00653    ,-.116772E-01,-.401196E-01, 1.00653    , 1.06069    ,
     *-.401196E-01,-.606896E-01, 1.06069    , .750715    ,-.606896E-01,
     *-.207830    , .594842    ,-.277629E-01,-.519576E-01,-2.21636/
!
!.....Tanne....................................................................
!
      data (B(i,2),i=1,80) /
     *-.272913    ,-.535819E-01, .272913    ,-.554505E-01,-.843877E-01,
     *-.123002    ,-.246447E-01, .150921E-01,-.123002    ,-.710641E-01,
     * .150921E-01, .710641E-01,-.710641E-01, .407259    , .710641E-01,
     * .215010    , .568516    , 1.29145    , .537524E-01, .506389    ,
     * .129034E-01, .714102E-02,-.129034E-01,-.146077E-03, .705986E-02,
     * .273133E-02,-.649231E-04, .129084E-03, .273133E-02,-.822691E-03,
     * .129084E-03, .822691E-03,-.822691E-03, .559431E-03, .822691E-03,
     *-.449196E-03, .222534E-03,-.970964E-04,-.112299E-03,-.389450E-03,
     *-.597541E-02,-.245982E-02, .597541E-02,-.154806E-03,-.254582E-02,
     *-.718117E-03,-.688026E-04,-.118285E-03,-.718117E-03, .399879E-03,
     *-.118285E-03,-.399879E-03, .399879E-03,-.881400E-03,-.399879E-03,
     * .126996E-03,-.786153E-03,-.114105E-02, .317490E-04, .767707E-02,
     * .363365    , .691837    ,-.363365    , .207776E-01, .703381    ,
     * .991621    , .923447E-02,-.312238E-01, .991621    , 1.09252    ,
     *-.312238E-01,-.925185E-01, 1.09252    , .638305    ,-.925185E-01,
     *-.229105    , .466476    ,-.276117    ,-.572764E-01,-.704160/
!
!.....Douglasie.........................................................
!
      data (B(i,3),i=1,80) /
     * .916011E-01, .267276E-01,-.916011E-01,-.245526E-01, .130873E-01,
     *-.112004    ,-.109123E-01, .248836E-01,-.112004    ,-.877935E-01,
     * .248836E-01, .877935E-01,-.877935E-01, .463178    , .877935E-01,
     * .159385    , .582717    , 1.21682    , .398463E-01, .762010    ,
     * .899184E-02, .449597E-02,-.899184E-02,-.369141E-03, .429089E-02,
     * .631694E-04,-.164063E-03, .693708E-03, .631694E-04,-.230363E-05,
     * .693708E-03, .230363E-05,-.230363E-05,-.539549E-04, .230363E-05,
     * .233716E-04,-.364262E-04,-.966559E-05, .584291E-05, .178525E-03,
     *-.582458E-02,-.264850E-02, .582458E-02, .186921E-03,-.254465E-02,
     * .195809E-03, .830762E-04,-.499295E-03, .195809E-03,-.595029E-04,
     *-.499295E-03, .595029E-04,-.595029E-04, .422024E-04, .595029E-04,
     * .206368E-04, .576800E-04, .154965E-03, .515921E-05,-.855741E-03,
     *-.187996    , .547549    , .187996    ,-.136054E-01, .539991    ,
     * .985003    ,-.604687E-02,-.529516E-01, .985003    , 1.11231    ,
     *-.529516E-01,-.112306    , 1.11231    , .565774    ,-.112306    ,
     *-.183433    , .428199    ,-.257790    ,-.458582E-01,-.635779/
!
!....Kiefer.............................................................
!
      data (B(i,4),i=1,80) /
     *-.104056    , .278676E-01, .104056    ,-.561985E-01,-.335377E-02,
     *-.102733    ,-.249771E-01, .271358E-01,-.102733    ,-.392980E-01,
     * .271358E-01, .392980E-01,-.392980E-01, .259925    , .392980E-01,
     * .317489    , .498042    , 1.36200    , .793722E-01, .480447    ,
     * .818093E-02, .343902E-03,-.818093E-02, .132793E-02, .108164E-02,
     * .283377E-03, .590190E-03, .686929E-04, .283377E-03,-.102728E-03,
     * .686929E-04, .102728E-03,-.102728E-03, .127535E-03, .102728E-03,
     *-.176758E-02,-.119815E-02,-.506009E-02,-.441896E-03, .254749E-01,
     *-.403881E-02,-.332653E-03, .403881E-02,-.603638E-03,-.668007E-03,
     *-.209360E-03,-.268283E-03,-.463799E-04,-.209360E-03,-.289915E-04,
     *-.463799E-04, .289915E-04,-.289915E-04, .325326E-03, .289915E-04,
     *-.117108E-04, .316543E-03, .467352E-03,-.292769E-05,-.313265E-02,
     * .188379    , .673190    ,-.188379    , .435361E-03, .673432    ,
     * .998091    , .193494E-03,-.454820E-01, .998091    , 1.04986    ,
     *-.454820E-01,-.498577E-01, 1.04986    , .802479    ,-.498577E-01,
     *-.320644    , .561996    ,-.283143    ,-.801610E-01,-1.03525/
!
!.....L?rche............................................................
!
      data (B(i,5),i=1,80) /
     *-.142477    , .130680    , .142477    ,-.994457E-01, .754323E-01,
     *-.739485E-01,-.441981E-01, .272001E-01,-.739485E-01,-.601290E-01,
     * .272001E-01, .601290E-01,-.601290E-01, .314464    , .601290E-01,
     * .254681    , .505475    , 1.26580    , .636702E-01, .851216    ,
     * .951228E-02, .212846E-02,-.951228E-02, .963436E-03, .266370E-02,
     * .952612E-03, .428194E-03, .540570E-04, .952612E-03,-.434139E-03,
     * .540570E-04, .434139E-03,-.434139E-03, .783942E-03, .434139E-03,
     *-.111104E-02,-.493394E-04,-.194014E-02,-.277761E-03, .823570E-02,
     *-.556938E-02,-.121968E-02, .556938E-02,-.539264E-03,-.151927E-02,
     *-.417015E-03,-.239673E-03,-.751449E-04,-.417015E-03, .234369E-03,
     *-.751449E-04,-.234369E-03, .234369E-03,-.520462E-03,-.234369E-03,
     *-.758713E-04,-.577366E-03,-.112549E-02,-.189678E-04, .683040E-02,
     * .200393    , .491092    ,-.200393    , .650944E-01, .527256    ,
     * .938037    , .289308E-01,-.445546E-01, .938037    , 1.08149    ,
     *-.445546E-01,-.814898E-01, 1.08149    , .736004    ,-.814898E-01,
     *-.258429    , .542182    ,-.212027    ,-.646074E-01,-1.25601/
!
!.....Buche.............................................................
!
      data (B(i,6),i=1,80) /
     *-.188531    ,-.574869E-01, .188531    ,-.221598E-01,-.697980E-01,
     *-.563012E-01,-.984882E-02,-.811468E-02,-.563012E-01,-.914926E-01,
     *-.811468E-02, .914926E-01,-.914926E-01, .422272    , .914926E-01,
     * .211683    , .581034    , 1.31420    , .529207E-01, .366134    ,
     * .149121E-02, .837835E-03,-.149121E-02,-.556769E-04, .806903E-03,
     * .185733E-03,-.247453E-04, .212077E-04, .185733E-03,-.308190E-03,
     * .212077E-04, .308190E-03,-.308190E-03, .104703E-02, .308190E-03,
     *-.885825E-03, .382660E-03,-.932838E-03,-.221456E-03, .242217E-02,
     * .233504E-02, .386838E-02,-.233504E-02,-.907703E-03, .336410E-02,
     * .136064E-02,-.403424E-03, .903357E-04, .136064E-02,-.100794E-03,
     * .903357E-04, .100794E-03,-.100794E-03,-.957468E-03, .100794E-03,
     * .405613E-03,-.653258E-03,-.168966E-03, .101403E-03, .318749E-02,
     * .207604    , .617870    ,-.207604    , .627061E-03, .618218    ,
     * .893818    , .278694E-03,-.752490E-02, .893818    , 1.12427    ,
     *-.752490E-02,-.124270    , 1.12427    , .609103    ,-.124270    ,
     *-.220643    , .443621    ,-.310410    ,-.551609E-01,-.477682/
!
!.....Eiche.............................................................
!
      data (B(i,7),i=1,80) /
     *-.172467    ,-.510311E-02, .172467    ,-.372934E-01,-.258217E-01,
     *-.385571E-01,-.165749E-01,-.715138E-02,-.385571E-01,-.942007E-01,
     *-.715138E-02, .942007E-01,-.942007E-01, .415360    , .942006E-01,
     * .229863    , .587757    , 1.35973    , .574658E-01, .152587    ,
     * .565412E-03, .108904E-02,-.565412E-03,-.270753E-03, .938617E-03,
     * .385190E-03,-.120335E-03,-.179621E-04, .385190E-03,-.276010E-03,
     *-.179621E-04, .276010E-03,-.276010E-03, .718850E-03, .276010E-03,
     *-.964489E-03,-.451681E-05,-.167719E-02,-.241122E-03, .696794E-02,
     * .320367E-02, .317107E-02,-.320367E-02,-.600605E-03, .283740E-02,
     * .813647E-03,-.266936E-03, .199659E-03, .813647E-03,-.121489E-04,
     * .199659E-03, .121489E-04,-.121489E-04,-.765051E-03, .121489E-04,
     * .252746E-03,-.575492E-03,-.383263E-03, .631866E-04, .377183E-02,
     * .216229    , .592582    ,-.216229    , .122861E-01, .599407    ,
     * .891263    , .546049E-02,-.995866E-02, .891263    , 1.12337    ,
     *-.995866E-02,-.123367    , 1.12337    , .615269    ,-.123367    ,
     *-.233533    , .440120    ,-.339378    ,-.583832E-01,-.344582/
!
!.....Roteiche..........................................................
!
      data (B(i,8),i=1,80) /
     *-.129150    , .260865E-02, .129150    ,-.317103E-01,-.150082E-01,
     *-.328701E-01,-.140935E-01,-.977679E-02,-.328701E-01,-.109393    ,
     *-.977679E-02, .109393    ,-.109393    , .470441    , .109393    ,
     * .185881    , .609852    , 1.31800    , .464703E-01, .242118    ,
     *-.589486E-03,-.193826E-02, .589486E-03, .647802E-03,-.157837E-02,
     *-.318212E-03, .287912E-03,-.178726E-03,-.318212E-03,-.130412E-03,
     *-.178726E-03, .130413E-03,-.130412E-03, .839862E-03, .130413E-03,
     *-.789992E-03, .247368E-03,-.104498E-02,-.197498E-03, .338793E-02,
     * .144288E-02, .260640E-02,-.144288E-02,-.665381E-03, .223674E-02,
     * .794479E-03,-.295725E-03, .140658E-03, .794479E-03, .196169E-03,
     * .140658E-03,-.196169E-03, .196169E-03,-.157916E-02,-.196169E-03,
     * .681781E-03,-.106782E-02,-.421477E-03, .170445E-03, .578675E-02,
     * .202464    , .631965    ,-.202464    ,-.644796E-02, .628383    ,
     * .893223    ,-.286576E-02,-.363381E-02, .893223    , 1.13626    ,
     *-.363381E-02,-.136261    , 1.13626    , .561734    ,-.136261    ,
     *-.191540    , .418080    ,-.300148    ,-.478849E-01,-.423841/
!
!....Regressionskoeffizienten fuer die D0.95- / d0.7-Schätzung:.....
!
      data nxk95a/5,5,2,2,5,6,2,5/
      data nxk95b/5,5,5,5,5,5,5,5/
!
      data xk95a/
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 60.0000,     0. ,     0. ,     0. ,     0. ,     0. ,
     *0.0000, 60.0000,     0. ,     0. ,     0. ,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 15.0000, 25.0000, 40.0000, 60.0000,     0. ,
     *0.0000, 60.0000,     0. ,     0. ,     0. ,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. /
!
      data xk95b/
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 30.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. ,
     *0.0000, 10.0000, 25.0000, 40.0000, 60.0000,     0. ,     0. /
!
!....D0.95 = a(H) + b(H)*D1.3 m / a(H)=Spline(H)........................
!
!.....FICHTE
!
      data (a95(i,1),i=1,16)/
     *    -.2500 ,  -.2500 ,  0.0000 ,  0.0000 ,
     *    -.2500 ,  -.2500 ,  0.0000 ,   .2989 ,
     *    -.2500 ,  1.5432 ,   .2989 ,  2.2087 ,
     *    -.1747 , 22.8312 ,  3.9267 ,  0.0000 /
!
!.....TANNE
!
      data (a95(i,2),i=1,16)/
     *     .1056 ,  -.1785 ,  0.0000 ,  0.0000 ,
     *    -.1785 ,  -.6047 ,  0.0000 ,   .5564 ,
     *    -.6047 ,  2.3076 ,   .5564 ,  0.0000 ,
     *    2.3076 ,  6.1906 ,  0.0000 ,  0.0000 /
!
!.....DOUGLASIE
!
      data (a95(i,3),i=1,4)/
     *    -.4938 ,   .9066 ,  0.0000 ,  0.0000 /
!
!.....KIEFER
!
      data (a95(i,4),i=1,4)/
     *    -.6966 ,  1.3798 ,  0.0000 ,  0.0000 /
!
!.....L?RCHE
!
      data (a95(i,5),i=1,16)/
     *    -.5000 ,  -.5000 ,  0.0000 ,  0.0000 ,
     *    -.5000 ,  -.5000 ,  0.0000 ,   .5225 ,
     *    -.5000 ,  2.6350 ,   .5225 ,  0.0000 ,
     *    2.6350 ,  6.8149 ,  0.0000 ,  0.0000 /
!
!.....BUCHE
!
      data (a95(i,6),i=1,20)/
     *    -.1097 ,  -.0903 ,  0.0000 , 0.0000    ,
     *    -.0903 ,  -.0806 ,  0.0000 ,-4.8637E-3 ,
     *    -.0660 ,  -.1342 ,  -.0195 ,  .1207    ,
     *    -.2851 ,   .9706 ,   .2716 , 0.0000    ,
     *     .9706 ,  2.6450 ,  0.0000 , 0.0000    /
!
!.....EICHE
!
      data (a95(i,7),i=1,4)/
     *    -.2299 ,   .4612 ,  0.0000 ,  0.0000 /
!
!.....ROTEICHE
!
      data (a95(i,8),i=1,16)/
     *    -.4213 ,  -.3170 ,  0.0000 ,  0.0000 ,
     *    -.3170 ,  -.1605 ,  0.0000 ,   .2115 ,
     *    -.1605 ,  1.2647 ,   .2115 ,  0.0000 ,
     *    1.2647 ,  3.1650 ,  0.0000 ,  0.0000 /
!
!.....D0.95 = a(H) + b(H)*D1.3 m / b(H)=Spline(H)
!
!.....FICHTE
!
      data (b95(i,1),i=1,16)/
     *    1.0967 ,  1.0666 ,  0.0000 ,  0.0000 ,
     *    1.0666 ,  1.0215 ,  0.0000 ,  -.0187 ,
     *    1.0215 ,   .8644 ,  -.0187 ,  0.0000 ,
     *     .8644 ,   .6550 ,  0.0000 ,  0.0000 /
!
!.....TANNE
!
      data (b95(i,2),i=1,16)/
     *    1.0806 ,  1.0569 ,  0.0000 ,  0.0000 ,
     *    1.0569 ,  1.0213 ,  0.0000 ,  -.0163 ,
     *    1.0213 ,   .8879 ,  -.0163 ,  0.0000 ,
     *     .8879 ,   .7101 ,  0.0000 ,  0.0000 /
!
!.....DOUGLASIE
!
      data (b95(i,3),i=1,16)/
     *    1.1958 ,  1.1105 ,  0.0000 ,  0.0000 ,
     *    1.1105 ,   .9826 ,  0.0000 ,   .0160 ,
     *     .9826 ,   .9509 ,   .0160 ,  0.0000 ,
     *     .9509 ,   .9087 ,  0.0000 ,  0.0000 /
!
!.....KIEFER
!
      data (b95(i,4),i=1,16)/
     *    1.2743 ,  1.1536 ,  0.0000 ,  0.0000 ,
     *    1.1536 ,   .9122 ,  0.0000 ,   .0533 ,
     *     .9522 ,   .9515 ,   .0133 ,  0.0000 ,
     *     .9515 ,   .9500 ,  0.0000 ,  0.0000 /
!
!.....L?RCHE
!
      data (b95(i,5),i=1,16)/
     *    1.1716 ,  1.1086 ,  0.0000    ,  0.0000    ,
     *    1.1086 ,  1.0141 ,  0.0000    , -9.5323E-3 ,
     *    1.0141 ,   .8624 , -9.5323E-3 ,  0.0000    ,
     *     .8624 ,   .6601 ,  0.0000    ,  0.0000    /
!
!.....BUCHE
!
      data (b95(i,6),i=1,16)/
     *    1.1185 ,  1.0708 ,  0.0000    ,  0.0000    ,
     *    1.0708 ,   .9992 ,  0.0000    ,  3.32990E-3,
     *     .9992 ,   .9476 ,  3.32990E-3,  0.0000    ,
     *     .9476 ,   .8788 ,  0.0000    ,  0.0000    /
!
!.....EICHE
!
      data (b95(i,7),i=1,16)/
     *    1.1246 ,  1.0723 ,  0.0000    ,  0.0000    ,
     *    1.0723 ,   .9940 ,  0.0000    ,  6.22952E-3,
     *     .9940 ,   .9530 ,  6.22952E-3,  0.0000    ,
     *     .9530 ,   .8984 ,  0.0000    ,  0.0000    /
!
!.....ROTEICHE
!
      data (b95(i,8),i=1,16)/
     *    1.1860 ,  1.1065 ,  0.0000 ,  0.0000 ,
     *    1.1065 ,   .9872 ,  0.0000 ,   .0122 ,
     *     .9872 ,   .9414 ,   .0122 ,  0.0000 ,
     *     .9414 ,   .8803 ,  0.0000 ,  0.0000 /
!
!.....Schätzung des d0.70 = a(H) / a(H)=Spline(H).......................
!

      Data Snxkn7/5,5,5,5,5,5,5,5/

      Data SXK07 /
     * 0.,10.,25.,40.,60.,0.,0.,0.,10.,25.,40.,60.,0.,0.,
     * 0.,10.,25.,40.,60.,0.,0.,0.,10.,25.,40.,60.,0.,0.,
     * 0.,10.,25.,40.,60.,0.,0.,0.,10.,25.,40.,60.,0.,0.,
     * 0.,10.,25.,40.,60.,0.,0.,0.,10.,25.,40.,60.,0.,0./

      Data (sd07(i,1), i=1,16)/
     *     .8623,     .8554,-2.3156E-3,-2.6624E-3,
     *     .8587,     .8183,-5.9904E-3,3.52532E-3,
     *     .8183,     .7991,3.52532E-3,3.30230E-3,
     *     .7965,     .8017,5.87075E-3,-1.7316E-3/

      data (sd07(i,2), i=1,16)/
     *     .9050,     .8901,-4.9698E-3,-3.2377E-3,
     *     .8941,     .8353,-7.2848E-3,1.58498E-3,
     *     .8353,     .7860,1.58498E-3,6.93605E-3,
     *     .7806,     .7797,     .0123,3.22789E-4/

      data (sd07(i,3), i=1,16)/
     *     .8604,     .8593,-3.5844E-4,-8.3832E-3,
     *     .8698,     .7738,    -.0189,     .0169,
     *     .7738,     .7793,     .0169,-7.4833E-4,
     *     .7798,     .7801,-1.3304E-3,-8.0434E-5/

      data (sd07(i,4), i=1,16)/
     *     .6999,     .7002,7.87642E-5,5.53208E-3,
     *     .6932,     .7558,     .0124,-4.4964E-3,
     *     .7558,     .7914,-4.4964E-3,-4.2429E-3,
     *     .7947,     .8026,-7.5429E-3,-2.6275E-3/

      data (sd07(i,5), i=1,16)/
     *     .8406,     .8388,-6.1619E-4,-5.0611E-3,
     *     .8451,     .7854,    -.0114,4.41328E-3,
     *     .7854,     .7522,4.41328E-3,5.02966E-3,
     *     .7482,     .7509,8.94162E-3,-8.7892E-4/

      data (sd07(i,6), i=1,16)/
     *     .7932,     .8136,6.80110E-3,-6.7888E-4,
     *     .8145,     .8374,-1.5275E-3,-3.0239E-3,
     *     .8374,     .8422,-3.0239E-3,-1.1995E-3,
     *     .8432,     .8384,-2.1325E-3,1.58773E-3/

      data (sd07(i,7), i=1,16)/
     *     .7863,     .7974,3.70835E-3,2.97515E-3,
     *     .7937,     .8439,6.69409E-3,-9.1629E-3,
     *     .8439,     .8390,-9.1629E-3,9.63937E-4,
     *     .8383,     .8409,1.71367E-3,-8.5683E-4/

       data (sd07(i,8), i=1,16)/
     *     .7999,     .8002,7.75448E-5,-6.7517E-4,
     *     .8010,     .7938,-1.5191E-3,-9.6160E-4,
     *     .7938,     .7807,-9.6160E-4,1.98990E-3,
     *     .7792,     .7804,3.53761E-3,-4.0580E-4/

!......Rindenabz?ge.....................................................

       data ((Rin(1,j,i),i=1,3),j=1,4)/
     * 1.55540 ,  .55475 ,  -0.00255   ,
     *  .82652 ,  .59424 ,  -0.00212   ,
     *  .17440 ,  .67905 ,  -0.00247   ,
     *  .85149 ,  .60934 ,  -0.00228   /

!      Ba(2)/'Tanne Mittenst.S.'/

       data ((Rin(2,j,i),i=1,3),j=1,4)/
     * 1.67703 , .56074 ,    0.     ,
     *  .82802 , .62504 ,    0.     ,
     *  .67058 , .68492 ,    0.     ,
     * 1.76896 , .59175 ,    0.     /

!      Ba(3)/'Douglasie Mittenst.S. BRD'/

       data ((Rin(3,j,i),i=1,3),j=1,4)/
     *   .26442, .84153, -.00233       ,
     *  -.78130, .80713, -.00290       ,
     * -1.28853, .82882, -.00354       ,
     * -2.13785, .91597, -.00375       /

!      Ba(4)/'Kiefer'/

       data ((Rin(4,j,i),i=1,3),j=1,4) /
     * 5.43367 ,  .62571 ,   0.     ,
     *  .05652 ,  .56149 ,   0.     ,
     * 4.17891 ,  .22292 ,   0.     ,
     * 1.59099 ,  .50146 ,   0.     /

!      Ba(5)/'Sitka-Fichte'/

       data ((Rin(5,j,i),i=1,3),j=1,4) /
     * 1.63833 , .79330 , -.01024      ,
     * 1.87531 , .66048 , -.00701      ,
     *-3.23264 , .99143 , -.01183      ,
     * 0.05167 , .81782 , -.00968      /

!      Ba(6)/'Europ.-L?rche'/

       data ((Rin(6,j,i),i=1,3),j=1,4) /
     * 4.17426 ,  .99857 ,   0.     ,
     * 1.21230 , 1.08527 ,   0.     ,
     * 1.24887 , 1.20858 ,   0.     ,
     * 3.58012 , 1.03147 ,   0.     /

!      Ba(7)/'Jap.L?rche'/

       data ((Rin(7,j,i),i=1,3),j=1,4) /
     * -8.04088 , 1.38585 ,   0.    ,
     * -3.47512 , 1.23849 ,   0.    ,
     * -15.26157, 1.94338 ,   0.    ,
     * -3.77073 , 1.29960 ,   0.    /

!      Ba(8)/'Schwarzkiefer'/

       data ((Rin(8,j,i),i=1,3),j=1,4) /
     * 10.72482 , 1.02208 ,   0.    ,
     * 11.00729 ,  .80471 ,   0.    ,
     * 10.30825 ,  .72255 ,   0.    ,
     *  5.27169 , 1.12602 ,   0.    /

!      Ba(9)/'Weymouthskiefer'/

       data ((Rin(9,j,i),i=1,3),j=1,4) /
     *  4.67218 ,  .52074 ,   0.    ,
     *  4.83026 ,  .44885 ,   0.    ,
     *  4.06786 ,  .46533 ,   0.    ,
     *  3.63666 ,  .50782 ,   0.    /

!      Ba(10)/'Buche'/

       data ((Rin(10,j,i),i=1,3),j=1,4) /
     * 1.97733 ,  .28119 ,   0.     ,
     * 2.25734 ,  .29724 ,   0.     ,
     * 2.69794 ,  .31096 ,   0.     ,
     * 2.61029 ,  .28522 ,   0.     /

!      Ba(11)/'Stiel-Eiche'/

       data ((Rin(11,j,i),i=1,3),j=1,4) /
     *  9.10974 ,  .66266 , 0.      ,
     *  8.94454 ,  .71505 , 0.      ,
     *  9.88377 ,  .75877 , 0.      ,
     * 10.18342 ,  .68997 , 0.      /

!      Ba(12)/'Trauben-Eiche (Bauland)'/

       data ((Rin(12,j,i),i=1,3),j=1,4) /
     *  11.90259 ,  .74264 , 0.     ,
     *  12.80219 ,  .75623 , 0.     ,
     *  13.54592 ,  .78294 , 0.     ,
     *  14.31589 ,  .72699 , 0.     /

!      Ba(13)/'Trauben-Eiche (Vorbergzone)'/

       data ((Rin(13,j,i),i=1,3),j=1,4) /
     *   6.53655 ,  .58931 , 0.     ,
     *   7.41153 ,  .61954 , 0.     ,
     *   9.30269 ,  .64363 , 0.     ,
     *   9.88855 ,  .56734 , 0.     /

!      Ba(14)/'Rot-Eiche'/

       data ((Rin(14,j,i),i=1,3),j=1,4) /
     *  -1.84576 ,  .86920 , -.00523   ,
     *  -3.98295 ,  .95987 , -.00616   ,
     *  -3.44335 ,  .97132 , -.00630   ,
     *  -3.19581 ,  .93891 , -.00596   /

!      Ba(15)/'Bergahorn'/

       data ((Rin(15,j,i),i=1,3),j=1,4) /
     *  -0.60951 ,  .64014 , -.00329   ,
     *  -3.53373 ,  .91611 , -.00707   ,
     *  -4.57300 , 1.06506 , -.00929   ,
     *  -0.62466 ,  .73312 , -.00482   /

!      Ba(16)/'Linde'/

       data ((Rin(16,j,i),i=1,3),j=1,4) /
     *   2.26031 ,  .60113 ,    0.  ,
     *    .38145 ,  .68665 ,    0.  ,
     *    .32491 ,  .72416 ,    0.  ,
     *   1.39565 ,  .65024 ,    0.  /

!      Ba(17)/'Esche'/

       data ((Rin(17,j,i),i=1,3),j=1,4) /
     * -1.14181 ,  .96466 , -.00432    ,
     * -8.40201 , 1.41083 , -.00964    ,
     * -3.62803 , 1.21051 , -.00777    ,
     * -7.97623 , 1.40182 , -.01011    /

!      Ba(18)/'Hainbuche'/

       data ((Rin(18,j,i),i=1,3),j=1,4) /
     *  5.03406 ,  .28277 ,   0.    ,
     *  9.12572 ,  .16002 ,   0.    ,
     *  6.96653 ,  .24105 ,   0.    ,
     *  7.47159 ,  .20957 ,   0.    /

!      Ba(19)/'Robinie'/

       data ((Rin(19,j,i),i=1,3),j=1,4) /
     *  -0.27505 , 1.58820 , 0.     ,
     *  -3.25513 , 1.70546 , 0.     ,
     *  -3.85748 , 1.81024 , 0.     ,
     *  -2.57381 , 1.69622 , 0.     /

!      Ba(20)/'Bergulme'/

       data ((Rin(20,j,i),i=1,3),j=1,4) /
     *   8.18895 ,  .84696 , 0.     ,
     *   8.61699 ,  .88780 , 0.     ,
     *   5.61794 , 1.06130 , 0.     ,
     *   8.26574 ,  .89505 , 0.     /

!      Ba(21)/'Birke'/

       data ((Rin(21,j,i),i=1,3),j=1,4) /
     *  -2.11471 , 0.94927 , 0.     ,
     *   1.59351 ,  .77363 , 0.     ,
     *   3.77295 ,  .69716 , 0.     ,
     *   1.63138 ,  .78958 , 0.     /

!      Ba(22)/'Maril.Pappel'/

       data ((Rin(22,j,i),i=1,3),j=1,4) /
     *  14.26657 ,  .77462 , 0.     ,
     *   6.75576 ,  .90682 , 0.     ,
     *   8.66612 ,  .90762 , 0.     ,
     *  12.74678 ,  .79624 , 0.     /

!      Ba(23)/'Robusta-Pappel'/

       data ((Rin(23,j,i),i=1,3),j=1,4) /
     *    .59848 , 1.20446 , -.00622   ,
     *    .66922 , 1.09135 , -.00480   ,
     *  -5.09109 , 1.35474 , -.00771   ,
     *  -1.26961 , 1.21661 , -.00624   /

!      Ba(24)/'Neupotz-Pappel'/

       data ((Rin(24,j,i),i=1,3),j=1,4) /
     *  -3.63150 ,  1.20021  , -0.00481,
     *  -4.36578 ,  1.13876  , -0.00394,
     *  -8.28526 ,  1.42238  , -0.00720,
     *  -3.68200 ,  1.21090  , -0.00514/

!      Ba(25)/'Regenerata-Pappel'/

       data ((Rin(25,j,i),i=1,3),j=1,4) /
     *  10.09257 ,  0.65481  ,       0.,
     *   7.88366 ,  0.63141  ,       0.,
     *   3.63461 ,  0.72700  ,       0.,
     *   7.20165 ,  0.67110  ,       0./

!      Ba(26)/'Kirsche'/

       data ((Rin(26,j,i),i=1,3),j=1,4) /
     *   4.54768 ,  0.52967  ,       0.,
     *   2.02813 ,  0.66068  ,       0.,
     *   1.54289 ,  0.73426  ,       0.,
     *   4.05603 ,  0.58080  ,       0./

!      Ba(27)/'Spitzahorn'/

       data ((Rin(27,j,i),i=1,3),j=1,4) /
     *   6.68053 ,  0.42792  ,       0.,
     *   4.81679 ,  0.50288  ,       0.,
     *   6.39906 ,  0.50162  ,       0.,
     *   7.43957 ,  0.43244  ,       0./

!      Ba(28)/'Roterle'/

       data ((Rin(28,j,i),i=1,3),j=1,4) /
     *  7.33868 ,  .82437 ,   0.    ,
     *  8.02746 ,  .83030 ,   0.    ,
     *  8.35663 ,  .90946 ,   0.    ,
     *  8.05239 ,  .84910 ,   0.    /

!      Ba(1)/'Fichte  HS'/

       data ((Rinh(1,j,i),i=1,3),j=1,5) /
     *  .44837 , .62253 , -.00288      ,
     *  .86715 , .59173 , -.00207      ,
     *  .08917 , .65407 , -.00286      ,
     *  .10719 , .64706 , -.00299      ,
     *-2.78642 ,1.03228 , -.00928      /

!      Ba(2)/'Tanne'/

       data ((Rinh(2,j,i),i=1,3),j=1,5) /
     * 1.19602 , .5998  ,    0.     ,
     *  .87088 , .62366 ,    0.     ,
     *  .71386 , .63868 ,    0.     ,
     * -.00716 , .68394 , -.00098   ,
     *-1.80354 , .98094 , -.00553   /

!  Douglasie *******  Ba-Wue ******************************

!    * 2.80936 , .61135 ,    0.     ,
!    * 2.73074 , .59815 ,    0.     ,
!    * 2.77691 , .58918 ,    0.     ,
!    * 2.72169 , .60313 ,    0.     ,
!    * 1.12142 , .65961 , -.00103   /

!      Ba(3)/'Douglasie BRD'/

       data ((Rinh(3,j,i),i=1,3),j=1,5) /
     *  .15897 , .75493 , -.00175      ,
     * -.74639 , .80389 , -.00282      ,
     * -.09716 , .74927 , -.00206      ,
     *  .30396 , .73153 , -.00154      ,
     *-1.78499 , .84379 , -.00377      /

!      Ba(4)/'Sitka-Fichte'/

!       data ((Rinh(4,j,i),i=1,3),j=1,5) /
!     * 1.53580 , .71693 , -.00877      ,
!     * 1.53580 , .71693 , -.00877      ,
!     * 1.53580 , .71693 , -.00877      ,
!     * 1.53580 , .71693 , -.00877      ,
!     *-3.00126 ,1.10681 , -.01385      /


!......Zuordnungstabelle................................................


      data Ban /

! 1 Schaftform

     * 1, 1, 2, 2, 4, 4, 4, 3, 5, 5, 5, 1, 1, 1, 6, 6, 7, 8, 8, 8, 6,
     * 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 6, 6, 6,

! 2 Rinde

     * 1, 1, 2, 2, 4, 8, 9, 3, 6, 6, 7, 1, 1, 2,10,18,11,14,25,25,17,
     *15,15,15,21,21,16,28,26,20,19,23,11,12,12,26,

! 3 Duchschni. Aufarbeitungsgrenze (EST)

     * 1, 1, 3, 3, 4, 4, 4, 1, 5, 5, 5, 1, 1, 1, 6, 6, 7, 7, 6, 6, 6,
     * 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7,

! 4 Unverwb. Derbholz.

     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2,
     * 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2,

! 5 durchschn. Astdurchmesser

     * 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2,
     * 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1,

! 6 Massentarife nach Krenn fuer schwache B?ume

     * 1, 1, 2, 2, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 4, 4, 5, 6, 6, 6, 6,
     * 6, 4, 4, 4, 6, 4, 6, 5, 6, 5, 4, 5, 6, 6, 4,

! 7 Massentafeln

     * 1, 1, 2, 2, 4, 4, 5, 3, 6, 6, 7, 1, 1, 1, 8,13, 9,10,14,14,11,
     *11,11,11, 9,13, 8,12, 9,11,13, 9, 9,11,11, 9/
!
!      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
!     Fi SF Ta KT Ki SK WK DG La EL JL Th Ts SN Bu HB Ei RE Pa BP Es
!     22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
!     Ah BA SA FA Bi Li Er Ki Ul Ro El Ka We LB VB
!
!     'FICHTE  S?ddeutschland'

      DATA (Azo(1,j),j=0,3)/
     *    .639936E+0   ,.663158E+0 ,-.380240E-1  , 0./

!     'FICHTE  Norddeutschland'

      DATA (Azo(2,j),j=0,3)/
     *   -.162417E+1   ,.230909E+1 ,-.349112E+0  , .913825E-002/

!     'TANNE'
      DATA (Azo(3,j),j=0,3)/
     *    .923207E+0   ,.615834E+0 ,-.688764E-1  , .781537E-2/

!     'KIEFER'

      DATA (Azo(4,j),j=0,3)/
     *   -.170313E+1   ,.238869E+1 ,-.396125     ,.241501E-1/

!     'L?RCHE'

      DATA (Azo(5,j),j=0,3)/
     *    .113353E+1   ,.167188    , .113906     ,-.161586E-1/

!     'BUCHE'

      DATA (Azo(6,j),j=0,3)/
     *    .241433E+1  ,-.330161E+0 , .129946E+0  ,-.112017E-1/

!     'EICHE'

      DATA (Azo(7,j),j=0,3)/
     *   -.358315E+0   ,.113226E+1 ,-.208216E-1  ,-.150484E-1/
!
!      BUCHE
!    ==========

      DATA (unvd(1,1,i),i=1,33)/
     *   36.00, 80.50,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,2,i),i=1,33)/
     *   15.79, 36.93, 58.90, 85.50,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,3,i),i=1,33)/
     *    7.77, 18.82, 32.68, 50.63, 66.33,
     *  85.82,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,4,i),i=1,33)/
     *    3.43,  8.65, 16.81, 28.18, 39.42,
     *  56.06, 71.01, 83.48,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,5,i),i=1,33)/
     *    2.32,  5.49, 10.24, 17.02, 24.70,
     *  35.46, 47.49, 59.66, 72.84, 84.49,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,6,i),i=1,33)/
     *    1.78,  3.87,  6.59, 10.40, 15.49,
     *  21.82, 30.56, 41.11, 53.39, 66.87, 77.54,
     *   89.87,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,7,i),i=1,33)/
     *    1.36,  2.87,  4.78,  7.20, 10.70,
     *  14.77, 20.22, 27.83, 36.98, 48.02, 58.74,
     *   70.50, 79.38, 88.40,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,8,i),i=1,33)/
     *    1.11,  2.44,  3.98,  5.85,  8.38,
     *  11.40, 14.35, 19.02, 24.62, 31.78, 40.76,
     *   50.54, 61.17, 72.32, 82.02, 92.71,
     *  99.94,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,9,i),i=1,33)/
     *    1.08,  2.15,  3.36,  4.81,  6.61,
     *   8.83, 10.83, 13.86, 17.28, 21.97, 28.61,
     *   35.99, 46.03, 56.67, 66.71, 75.05,
     *  83.40, 93.67, 99.85,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,10,i),i=1,33)/
     *    1.22,  1.95,  2.92,  4.07,  5.38,
     *   7.03,  9.39, 11.93, 14.41, 17.88, 21.81,
     *   26.57, 33.25, 41.19, 51.22, 59.64,
     *  68.32, 75.63, 82.77, 91.59,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,11,i),i=1,33)/
     *    1.27,  1.83,  2.66,  3.64,  4.69,
     *   6.02,  8.43, 10.65, 12.59, 15.30, 17.51,
     *   20.51, 24.60, 30.22, 37.80, 46.48,
     *  54.69, 61.63, 69.11, 78.18, 86.16, 92.58,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,12,i),i=1,33)/
     *    1.18,  1.79,  2.58,  3.51,  4.55,
     *   5.80,  7.69,  9.57, 11.26, 13.50, 15.23,
     *   17.52, 20.39, 24.47, 28.74, 35.58,
     *  42.53, 49.67, 56.86, 65.57, 72.83, 79.91,
     *   87.22, 93.61,100.76,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,13,i),i=1,33)/
     *    1.04,  1.76,  2.58,  3.53,  4.66,
     *   5.92,  7.16,  8.69, 10.27, 12.22, 14.05,
     *   16.15, 18.52, 21.59, 24.01, 27.52,
     *  32.80, 39.75, 46.02, 53.78, 60.27, 67.34,
     *   73.89, 80.07, 87.04, 92.38,100.75,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,14,i),i=1,33)/
     *     .94,  1.74,  2.59,  3.54,  4.70,
     *   5.95,  6.85,  8.01,  9.48, 11.18, 13.04,
     *   14.95, 16.93, 19.19, 21.64, 23.89,
     *  27.49, 31.86, 36.60, 42.79, 48.48, 54.86,
     *   61.35, 68.16, 74.86, 80.93, 87.70,
     *  94.60,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(1,15,i),i=1,33)/
     *     .87,  1.73,  2.59,  3.53,  4.69,
     *   5.89,  6.72,  7.52,  8.90, 10.40, 12.20,
     *   13.93, 15.60, 17.27, 19.60, 21.24,
     *  23.25, 26.01, 28.80, 32.99, 37.85, 43.97,
     *   50.59, 56.88, 64.23, 70.52, 77.31,
     *  82.04, 86.52, 90.71, 94.32,100.00,100.00/

      DATA (unvd(1,16,i),i=1,33)/
     *     .84,  1.73,  2.60,  3.51,  4.62,
     *   5.74,  6.59,  7.23,  8.52,  9.87, 11.52,
     *   13.08, 14.55, 15.84, 17.92, 19.94,
     *  21.07, 22.20, 23.79, 26.57, 30.66, 35.69,
     *   41.62, 47.22, 54.15, 60.15, 66.56,
     *  70.86, 74.58, 79.66, 82.91, 87.62, 91.80/

      DATA (unvd(1,17,i),i=1,33)/
     *     .85,  1.73,  2.60,  3.47,  4.49,
     *   5.50,  6.43,  7.14,  8.35,  9.58, 11.02,
     *   12.40, 13.76, 14.90, 16.58, 18.54,
     *  19.63, 20.43, 21.78, 23.92, 27.31, 30.51,
     *   34.43, 39.19, 44.62, 49.81, 55.47,
     *  61.06, 65.62, 71.30, 75.94, 79.38, 84.24/

      DATA (unvd(1,18,i),i=1,33)/
     *     .88,  1.74,  2.60,  3.44,  4.35,
     *   5.24,  6.26,  7.12,  8.25,  9.41, 10.61,
     *   11.83, 13.13, 14.22, 15.50, 17.20,
     *  18.60, 19.69, 21.25, 23.14, 25.90, 28.46,
     *   30.25, 33.16, 36.89, 41.05, 45.82,
     *  52.62, 57.62, 62.62, 66.40, 70.54, 76.75/

      DATA (unvd(1,19,i),i=1,33)/
     *     .92,  1.77,  2.61,  3.43,  4.25,
     *   5.05,  6.05,  7.07,  8.13,  9.21, 10.22,
     *   11.31, 12.52, 13.58, 14.58, 16.07,
     *  17.67, 18.97, 20.68, 22.37, 24.55, 25.55,
     *   27.31, 29.51, 32.20, 35.40, 39.39,
     *  45.57, 50.59, 54.63, 59.30, 64.09, 70.31/

      DATA (unvd(1,20,i),i=1,33)/
     *     .97,  1.80,  2.63,  3.42,  4.19,
     *   4.93,  5.83,  6.99,  7.97,  8.99,  9.87,
     *   10.83, 11.95, 12.99, 13.82, 15.15,
     *  16.84, 18.28, 20.06, 21.60, 23.26, 24.79,
     *   26.32, 27.94, 30.18, 32.50, 35.83,
     *  39.88, 44.53, 48.32, 52.64, 57.03, 62.15/

      DATA (unvd(1,21,i),i=1,33)/
     *    1.03,  1.86,  2.64,  3.44,  4.17,
     *   4.88,  5.68,  6.87,  7.77,  8.75,  9.54,
     *   10.40, 11.41, 12.45, 13.23, 14.45,
     *  16.12, 17.62, 19.40, 20.83, 22.04, 23.17,
     *   25.55, 26.67, 28.57, 30.24, 32.95,
     *  35.58, 39.44, 42.70, 46.41, 50.36, 54.62/

      DATA (unvd(1,22,i),i=1,33)/
     *    1.09,  1.92,  2.66,  3.46,  4.19,
     *   4.90,  5.61,  6.72,  7.54,  8.48,  9.25,
     *   10.02, 10.91, 11.95, 12.80, 13.95,
     *  15.50, 16.98, 18.70, 20.07, 20.87, 22.70,
     *   23.70, 25.38, 27.01, 28.25, 30.39,
     *  32.24, 35.19, 38.76, 40.81, 44.09, 47.92/

      DATA (unvd(1,23,i),i=1,33)/
     *    1.15,  1.98,  2.70,  3.50,  4.23,
     *   4.95,  5.61,  6.57,  7.32,  8.23,  8.99,
     *    9.69, 10.47, 11.50, 12.45, 13.56,
     *  14.93, 16.34, 17.95, 19.28, 19.77, 21.38,
     *   22.77, 24.12, 25.49, 26.55, 28.15,
     *  29.46, 31.65, 34.64, 36.01, 39.20, 42.06/

      DATA (unvd(1,24,i),i=1,33)/
     *    1.19,  2.03,  2.75,  3.54,  4.27,
     *   5.00,  5.70,  6.46,  7.15,  8.02,  8.77,
     *    9.43, 10.13, 11.08, 12.09, 13.13,
     *  14.35, 15.68, 17.16, 18.43, 18.73, 20.21,
     *   21.77, 22.84, 24.01, 25.12, 26.22,
     *  27.26, 28.83, 30.47, 32.02, 34.71, 37.05/

      DATA (unvd(1,25,i),i=1,33)/
     *    1.21,  2.07,  2.82,  3.57,  4.32,
     *   5.04,  5.84,  6.39,  7.03,  7.86,  8.60,
     *    9.24,  9.89, 10.71, 11.72, 12.68,
     *  13.75, 14.99, 16.33, 17.51, 17.75, 19.18,
     *   20.70, 21.57, 22.57, 23.91, 24.59,
     *  25.55, 26.66, 27.56, 28.81, 30.67, 32.88/

      DATA (unvd(1,26,i),i=1,33)/
     *    1.20,  2.10,  2.90,  3.61,  4.37,
     *   5.09,  5.95,  6.35,  6.97,  7.75,  8.48,
     *    9.11,  9.76, 10.37, 11.33, 12.21,
     *  13.13, 14.27, 15.45, 16.54, 16.83, 18.30,
     *   19.55, 20.29, 21.18, 22.63, 23.06,
     *  24.01, 24.82, 25.61, 26.28, 27.51, 29.58/

      DATA (unvd(1,27,i),i=1,33)/
     *    1.18,  2.11,  2.99,  3.65,  4.42,
     *   5.14,  5.99,  6.35,  6.96,  7.69,  8.40,
     *    9.05,  9.73, 10.08, 10.94, 11.70,
     *   12.51, 13.52, 14.54, 15.5, 15.97, 17.56,
     *  18.33, 19., 19.82, 21.22, 21.59, 22.57,
     *   23.23, 24., 24.4, 25.28, 27.15/
!
!  EICHE
!===========
      DATA (unvd(2, 1,i),i=1,33)/
     *   25.00, 47.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 2,i),i=1,33)/
     *   15.00, 29.00, 63.42, 83.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 3,i),i=1,33)/
     *    8.03, 16.61, 31.00, 49.26, 70.81,
     *  91.33,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 4,i),i=1,33)/
     *    3.56,  8.30, 12.99, 21.55, 31.03,
     *  46.42, 66.69, 83.55, 99.47,102.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 5,i),i=1,33)/
     *    2.42,  5.24,  9.27, 14.41, 21.06,
     *  29.84, 43.72, 57.78, 71.80, 84.97,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 6,i),i=1,33)/
     *    2.28,  4.14,  6.71, 10.40, 14.72,
     *  20.25, 27.79, 38.41, 50.64, 64.84, 76.42,
     *   89.34,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 7,i),i=1,33)/
     *    2.15,  3.51,  5.31,  8.38, 11.13,
     *  14.59, 18.88, 25.45, 33.98, 44.85, 55.04,
     *   65.71, 77.23, 87.59,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 8,i),i=1,33)/
     *    2.00,  3.12,  4.54,  6.63,  8.70,
     *  11.27, 14.22, 17.39, 22.13, 28.45, 38.38,
     *   46.85, 56.96, 66.35, 78.08, 87.67,
     * 100.00,100.00,100.00,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2, 9,i),i=1,33)/
     *    1.83,  2.77,  3.90,  5.27,  6.81,
     *   8.69, 10.98, 12.68, 15.39, 19.11, 26.46,
     *   32.76, 40.24, 47.29, 58.20, 67.22,
     *  78.07, 88.57,100.24,100.00,100.00,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2,10,i),i=1,33)/
     *    1.63,  2.45,  3.39,  4.28,  5.47,
     *   6.87,  9.06, 10.94, 13.15, 15.97, 19.00,
     *   23.25, 28.07, 34.44, 43.01, 51.04,
     *  61.36, 71.11, 79.11, 87.31, 95.05,100.00,
     *  100.00,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2,11,i),i=1,33)/
     *    1.40,  2.17,  3.00,  3.68,  4.67,
     *   5.80,  7.68,  9.79, 11.74, 14.00, 16.45,
     *   18.29, 21.43, 25.79, 32.52, 39.13,
     *  47.70, 55.75, 65.18, 72.88, 80.36, 88.63,
     *   95.76,100.00,100.00,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2,12,i),i=1,33)/
     *    1.15,  1.92,  2.73,  3.45,  4.43,
     *   5.48,  6.71,  8.83, 10.55, 12.35, 13.55,
     *   15.70, 18.29, 21.34, 26.31, 31.49,
     *  38.07, 44.47, 52.45, 60.69, 67.57, 76.17,
     *   81.80, 89.74, 95.32,100.00,100.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2,13,i),i=1,33)/
     *     .92,  1.71,  2.54,  3.41,  4.43,
     *   5.51,  6.15,  8.06,  9.57, 11.01, 12.01,
     *   13.85, 16.38, 19.16, 21.95, 26.48,
     *  31.17, 36.08, 41.77, 49.73, 56.68, 62.80,
     *   69.38, 75.64, 82.40, 88.98, 96.00,
     * 100.00,100.00,100.00,100.00,100.00,100.00/

      DATA (unvd(2,14,i),i=1,33)/
     *     .75,  1.55,  2.39,  3.36,  4.40,
     *   5.49,  6.01,  7.47,  8.81, 10.00, 11.54,
     *   13.11, 15.39, 17.31, 19.02, 22.47,
     *  25.70, 29.41, 34.00, 40.00, 45.67, 51.54,
     *   57.49, 64.55, 71.61, 78.84, 85.03,
     *  91.54, 96.56, 99.95,100.52,100.00,100.00/

      DATA (unvd(2,15,i),i=1,33)/
     *     .64,  1.43,  2.28,  3.28,  4.33,
     *   5.42,  6.21,  7.08,  8.26,  9.30, 11.14,
     *   12.50, 14.35, 15.78, 17.39, 19.46,
     *  21.67, 24.46, 27.89, 31.73, 35.87, 41.37,
     *   47.14, 54.48, 60.93, 66.45, 73.33,
     *  80.11, 85.13, 89.63, 95.07,100.00,100.00/

      DATA (unvd(2,16,i),i=1,33)/
     *     .60,  1.36,  2.20,  3.19,  4.22,
     *   5.29,  6.37,  6.88,  7.93,  8.92, 10.81,
     *   12.00, 13.46, 14.58, 16.19, 17.45,
     *  19.06, 21.21, 23.99, 26.16, 29.15, 34.29,
     *   39.33, 45.42, 51.37, 55.83, 62.21,
     *  68.30, 75.71, 80.61, 86.40, 92.00, 95.84/

      DATA (unvd(2,17,i),i=1,33)/
     *     .61,  1.34,  2.15,  3.09,  4.08,
     *   5.10,  6.42,  6.87,  7.82,  8.86, 10.55,
     *   11.63, 12.75, 13.71, 15.29, 16.44,
     *  17.88, 19.68, 21.05, 23.53, 25.82, 29.32,
     *   33.06, 37.37, 41.93, 46.96, 52.69,
     *  58.11, 64.29, 70.90, 76.50, 81.41, 85.45/

      DATA (unvd(2,18,i),i=1,33)/
     *     .66,  1.34,  2.13,  3.00,  3.95,
     *   4.93,  6.37,  6.95,  7.83,  8.95, 10.36,
     *   11.37, 12.21, 13.16, 14.68, 15.97,
     *  17.45, 19.06, 20.09, 22.45, 24.33, 26.08,
     *   28.53, 31.20, 35.24, 39.85, 44.76,
     *  49.54, 55.87, 61.50, 67.39, 71.36, 75.19/

      DATA (unvd(2,19,i),i=1,33)/
     *     .70,  1.36,  2.14,  2.96,  3.87,
     *   4.82,  6.20,  7.01,  7.85,  9.03, 10.24,
     *   11.23, 11.85, 12.93, 14.37, 15.62,
     *  17.10, 18.55, 19.09, 21.57, 23.13, 24.22,
     *   25.96, 27.75, 30.93, 34.50, 38.41,
     *  42.59, 47.46, 52.41, 58.05, 61.85, 66.08/

      DATA (unvd(2,20,i),i=1,33)/
     *     .74,  1.39,  2.16,  2.96,  3.85,
     *   4.79,  5.95,  7.05,  7.90,  9.08, 10.19,
     *   11.19, 11.65, 12.99, 14.31, 15.36,
     *  16.82, 18.15, 18.93, 20.87, 22.20, 23.56,
     *   25.11, 26.70, 28.59, 30.90, 33.67,
     *  36.26, 39.48, 44.50, 49.36, 53.89, 59.12/

      DATA (unvd(2,21,i),i=1,33)/
     *     .77,  1.43,  2.21,  3.01,  3.89,
     *   4.82,  5.78,  7.08,  7.96,  9.11, 10.13,
     *   11.16, 11.64, 13.07, 14.28, 15.21,
     *  16.62, 17.87, 18.88, 20.37, 21.57, 23.10,
     *   24.53, 26.03, 27.72, 29.07, 30.51,
     *  32.55, 34.36, 37.63, 41.18, 47.48, 53.30/

      DATA (unvd(2,22,i),i=1,33)/
     *     .80,  1.49,  2.28,  3.11,  3.99,
     *   4.92,  5.69,  7.10,  8.05,  9.13, 10.07,
     *   11.12, 11.79, 13.13, 14.23, 15.17,
     *  16.50, 17.70, 18.82, 20.06, 21.22, 22.66,
     *   23.99, 25.41, 26.92, 28.15, 28.94,
     *  31.47, 33.09, 34.81, 38.50, 43.05, 48.28/

      DATA (unvd(2,23,i),i=1,33)/
     *     .82,  1.55,  2.36,  3.22,  4.10,
     *   5.04,  5.70,  7.10,  8.12,  9.13, 10.00,
     *   11.06, 12.01, 13.16, 14.17, 15.16,
     *  16.41, 17.58, 18.74, 19.85, 21.01, 22.25,
     *   23.49, 24.83, 26.19, 27.30, 28.20,
     *  30.25, 31.28, 33.56, 36.34, 40.03, 43.73/

      DATA (unvd(2,24,i),i=1,33)/
     *     .84,  1.59,  2.42,  3.31,  4.19,
     *   5.14,  5.81,  7.08,  8.17,  9.11,  9.92,
     *   11.00, 12.18, 13.16, 14.10, 15.14,
     *  16.30, 17.45, 18.64, 19.64, 20.80, 21.87,
     *   23.03, 24.30, 25.51, 26.53, 27.52,
     *  29.14, 30.53, 32.42, 34.46, 37.42, 40.64/

      DATA (unvd(2,25,i),i=1,33)/
     *     .86,  1.63,  2.47,  3.37,  4.25,
     *   5.21,  5.99,  7.05,  8.20,  9.08,  9.84,
     *   10.92, 12.30, 13.15, 14.01, 15.09,
     *  16.18, 17.32, 18.53, 19.42, 20.59, 21.52,
     *   22.61, 23.82, 24.90, 25.84, 26.91,
     *  28.15, 29.84, 31.38, 32.65, 35.22, 37.00/

      DATA (unvd(2,26,i),i=1,33)/
     *     .87,  1.66,  2.50,  3.42,  4.29,
     *   5.24,  6.15,  7.01,  8.19,  9.02,  9.75,
     *   10.84, 12.37, 13.10, 13.91, 15.02,
     *  16.05, 17.17, 18.40, 19.21, 20.37, 21.19,
     *   22.23, 23.38, 24.35, 25.22, 26.35,
     *  27.27, 29.21, 30.44, 31.90, 32.44, 33.84/

      DATA (unvd(2,27,i),i=1,33)/
     *     .88,  1.68,  2.52,  3.44,  4.31,
     *   5.25,  6.25,  6.94,  8.16,  8.95,  9.66,
     *   10.74, 12.40, 13.04, 13.79, 14.92,
     *   15.91,17.02, 18.26, 19.  , 20.16, 20.89,
     *   21.88, 22.99, 23.86, 24.67, 25.86, 26.5,
     *   27.64, 28.61 ,30.23,31.07,32.13/

      DATA PI/7.8539816E-5/

!     Hoehenrahmen FICHTE

      data (hoehr(1,i),i=1,6)/
     *6.5,7.5,8.5,8.,9.,10./

!     Hoehenrahmen TANNE

      data (hoehr(2,i),i=1,6)/
     *7.,7.5,8.5,8.,9.,10./

!     Hoehenrahmen KIEFER

      data (hoehr(3,i),i=1,6)/
     *7.,8.,9.,9.,10.,11./

!     Hoehenrahmen BUCHE

      data (hoehr(4,i),i=1,6)/
     *9.,10.,11.,10.5,12.,13./

!     Hoehenrahmen EICHE

      data (hoehr(5,i),i=1,6)/
     *8.,9.,10.,9.5,10.5,11.5/

!     Hoehenrahmen ESCHE

      data (hoehr(6,i),i=1,6)/
     *9.5,10.5,11.5,11.,12.,13./

!    * FICHTE-U

      data (Volk(1,I,1),i=1,30)/
     *110,114,118,121,125,129,133,137,140,144,
     *148,155,161,168,174,181,188,194,201,207,
     *214,224,235,245,256,266,276,287,297,308/

!    * FICHTE-M

      data (Volk(1,I,2),i=1,30)/
     *125,129,134,138,143,147,151,155,159,164,
     *168,176,183,191,198,206,213,221,228,236,
     *243,255,267,278,290,302,314,326,337,349/

!    * FICHTE-O

      data (Volk(1,I,3),i=1,30)/
     *140,145,150,154,159,164,169,174,178,183,
     *188,196,205,213,222,230,238,247,255,264,
     *272,285,298,312,325,338,351,364,378,391/

!    * TANNE-U

      data (Volk(2,I,1),i=1,30)/
     * 95,101,106,112,118,124,129,135,141,146,
     *152,160,167,175,182,190,198,205,213,220,
     *228,238,247,257,267,277,286,296,306,315/

!    * TANNE-M

      data (Volk(2,I,2),i=1,30)/
     *108,115,121,128,134,141,147,154,160,167,
     *173,182,190,199,207,216,225,233,242,250,
     *259,270,281,292,303,314,325,336,347,358/

!    * TANNE-O

      data (Volk(2,I,3),i=1,30)/
     *121,128,136,143,150,158,165,172,179,187,
     *194,204,213,223,232,242,252,261,271,280,
     *290,302,315,327,339,352,364,376,388,401/

!    * KIEFER-U

      data (Volk(3,I,1),i=1,30)/
     * 91, 96,101,106,111,116,121,126,131,136,
     *141,148,154,161,167,174,181,187,194,200,
     *207,215,224,232,240,249,257,265,273,282/

!    * KIEFER-M

      data (Volk(3,I,2),i=1,30)/
     *104,110,115,121,126,132,138,143,149,154,
     *160,168,175,183,190,198,205,213,220,228,
     *235,245,254,264,273,283,292,302,311,321/

!    * KIEFER-O

      data (Volk(3,I,3),i=1,30)/
     *116,122,129,135,141,148,154,160,166,173,
     *179,187,196,204,213,221,229,238,246,255,
     *263,274,284,295,306,317,327,338,349,359/

!    *Buche-U

      data (Volk(4,I,1),i=1,30)/
     *34, 40, 45, 51, 56, 62, 67, 73, 78, 84,
     *89, 98,107,115,124,133,142,151,159,168,
     *177,187,196,206,216,226,235,245,255,264/

!    * BUCHE-M

      data (Volk(4,I,2),i=1,30)/
     *39, 45, 51, 58, 64, 70, 76, 82, 87, 94,
     *101,111,121,131,141,151,161,171,181,191,
     *201,212,223,234,245,256,267,278,289,300/

!    * BUCHE-O

      data (Volk(4,I,3),i=1,30)/
     *44, 51, 58, 65, 72, 79, 85, 92, 99,106,
     *113,124,135,147,158,169,180,191,203,214,
     *225,237,250,262,274,287,299,311,323,336/

!    * EICHE-U

      data (Volk(5,I,1),i=1,30)/
     *123,129,135,142,148,154,160,166,173,179,
     *185,192,198,205,211,218,225,231,238,244,
     *251,260,270,279,289,298,307,317,326,336/

!    * EICHE-M

      data (Volk(5,I,2),i=1,30)/
     *140,147,154,161,168,175,182,189,196,203,
     *210,218,225,233,240,248,255,263,270,278,
     *285,296,306,317,328,339,349,360,371,381/

!    * EICHE-O

      data (Volk(5,I,3),i=1,30)/
     *157,164,173,180,188,196,202,212,219,227,
     *235,243,252,260,269,277,285,294,302,311,
     *319,331,343,355,367,379,391,403,415,427/

!    * ESCHE-U

      data (Volk(6,I,1),i=1,30)/
     *26, 34, 42, 50, 58, 67, 75, 83, 91, 99,
     *107,117,126,136,145,155,164,174,183,193,
     *202,213,223,234,244,255,266,276,287,297/

!    *Esche-M

      data (Volk(6,I,2),i=1,30)/
     * 29, 38, 47, 57, 66, 75, 84, 93,103,112,
     *121,132,143,153,164,175,186,197,207,218,
     *229,241,253,265,277,290,302,314,326,338/

!    * ESCHE-O

      data (Volk(6,I,3),i=1,30)/
     *33, 43, 54, 64, 74, 85, 95,105,115,126,
     *136,148,160,172,184,197,209,221,233,245,
     *257,271,284,298,311,325,338,352,365,379/


!.....Massentafelaequivalente d07-Werte..................................


!     Fichte............................................................

      DATA (Md07(i),i=1,9) /580,576,557,522,473,457,456,411,416/
      DATA (Md07(i),i=10,20) /720,671,656,628,595,581,605,592,572,
     *561,544/
      DATA (Md07(i),i=21,35) /604,753,745,716,707,689,663,672,687,
     *673,652,653,646,633,625/
      DATA (Md07(i),i=36,52) /425,820,779,769,761,747,724,728,738,
     *722,714,709,699,692,681,673,669/
      DATA (Md07(i),i=53,70) /855,814,788,782,770,764,753,771,755,
     *745,737,733,731,724,713,710,704,695/
      DATA (Md07(i),i=71,91) /885,820,797,794,783,790,781,793,777,
     *773,764,757,752,748,740,734,730,723,717,712,706/
      DATA (Md07(i),i=92,115) /901,827,808,803,806,811,801,809,801,
     *794,783,779,770,764,759,754,747,741,737,732,726,721,714,708/
      DATA (Md07(i),i=116,141) /894,833,814,808,811,816,814,819,810,
     *803,796,790,785,781,773,770,764,759,751,748,742,735,730,725,
     *720,714/
      DATA (Md07(i),i=142,167) /903,836,817,813,815,818,815,827,817,
     *815,806,803,796,789,784,778,774,769,763,759,754,747,743,737,
     *732,726/
      DATA (Md07(i),i=168,197) /903,838,826,818,817,829,824,834,823,
     *818,814,809,804,800,792,788,784,778,773,769,764,758,754,748,
     *743,739,733,728,723,718/
      DATA (Md07(i),i=198,230) /905,848,826,820,822,830,825,838,832,
     *826,820,815,812,806,802,797,791,786,783,778,774,768,764,758,
     *753,749,744,739,734,730,726,721,715/
      DATA (Md07(i),i=231,265) /912,862,836,832,831,830,831,841,835,
     *828,826,822,818,813,809,804,800,795,791,786,782,777,773,768,
     *764,759,755,750,746,741,737,732,727,723,719/
      DATA (Md07(i),i=266,301) /870,845,828,831,831,831,844,837,834,
     *831,826,823,818,814,811,805,801,797,793,789,784,781,777,772,
     *768,763,759,755,750,747,743,738,734,730,726,721/
      DATA (Md07(i),i=302,338) /846,831,833,832,832,846,843,838,834,
     *832,828,824,819,815,812,808,803,800,796,791,788,784,780,776,
     *772,768,764,760,756,752,748,745,740,737,733,729,725/
      DATA (Md07(i),i=339,378) /834,831,834,833,849,844,839,837,834,
     *832,827,823,820,817,813,810,806,802,799,794,791,787,783,780,
     *776,772,769,765,761,758,754,750,746,742,739,735,732,728,725,
     *721/
      DATA (Md07(i),i=379,421) /833,831,831,851,846,843,841,839,835,
     *832,828,824,821,818,815,811,808,805,801,797,794,790,787,784,
     *780,776,773,769,766,762,759,755,752,749,745,742,738,735,732,
     *728,725,721,718/
      DATA (Md07(i),i=422,468) /834,827,831,852,850,847,843,840,839,
     *835,833,829,826,823,819,817,813,810,807,803,800,797,793,790,
     *787,784,780,777,774,771,767,764,760,757,754,751,748,744,741,
     *738,735,732,729,726,723,719,716/
      DATA (Md07(i),i=469,526) /825,835,854,851,848,846,844,842,838,
     *836,833,830,826,823,821,818,815,812,809,805,803,800,796,793,
     *790,787,784,781,778,775,771,769,766,763,759,756,754,750,747,
     *744,742,739,736,733,730,727,724,721,718,715,713,710,707,704,
     *701,699,696,693/
      DATA (Md07(i),i=527,612) /839,858,852,851,848,846,844,842,839,
     *836,833,831,828,825,823,819,816,814,811,808,805,802,799,796,
     *794,790,788,785,782,779,776,773,770,767,765,762,759,756,753,
     *751,748,745,742,739,737,734,731,729,726,723,721,718,715,713,
     *710,708,705,702,700,697,695,692,690,687,685,682,680,677,675,
     *672,670,667,665,662,660,658,655,653,650,648,646,643,641,639,
     *636,634/
      DATA (Md07(i),i=613,698) /859,855,853,850,849,847,844,842,840,
     *837,834,832,829,827,824,821,818,816,813,810,807,805,802,799,
     *796,794,791,788,786,783,780,777,775,772,769,767,764,761,759,
     *756,753,751,748,746,743,741,738,736,733,730,728,726,723,721,
     *718,716,713,711,708,706,704,701,699,697,694,692,690,687,685,
     *683,680,678,676,673,671,669,667,664,662,660,658,655,653,651,
     *649,646/
      DATA (Md07(i),i=699,783) /856,854,853,851,849,847,845,843,841,
     *838,836,833,830,828,825,823,820,818,815,812,810,807,805,802,
     *799,797,794,792,789,786,784,781,779,776,774,771,769,766,764,
     *761,759,757,754,752,749,747,744,742,740,737,735,733,730,728,
     *726,723,721,719,716,714,712,710,707,705,703,701,699,696,694,
     *692,690,688,686,683,681,679,677,675,673,671,669,666,664,662,
     *660/
      DATA (Md07(i),i=784,867) /856,855,854,852,850,848,846,843,842,
     *839,836,834,831,830,827,824,822,819,817,815,812,810,807,805,
     *802,800,797,795,793,790,788,785,783,780,778,776,773,771,769,
     *766,764,762,759,757,755,753,750,748,746,744,741,739,737,735,
     *733,730,728,726,724,722,720,718,716,713,711,709,707,705,703,
     *701,699,697,695,693,691,689,687,685,683,681,679,677,675,673/
      DATA (Md07(i),i=868,950) /857,855,854,852,850,849,846,844,842,
     *840,838,835,833,831,828,826,824,821,819,817,814,812,810,807,
     *805,803,801,798,796,793,791,789,787,784,782,780,778,775,773,
     *771,769,767,765,762,760,758,756,754,752,750,747,745,743,741,
     *739,737,735,733,731,729,727,725,723,721,719,717,715,713,711,
     *709,707,705,703,702,700,698,696,694,692,690,688,686,685/
      DATA (Md07(i),i=951,1032) /857,856,855,853,851,849,847,845,843,
     *841,839,837,834,832,830,828,826,823,821,819,817,814,812,810,
     *808,806,803,801,799,797,795,792,790,788,786,784,782,780,778,
     *775,773,771,769,767,765,763,761,759,757,755,753,751,749,747,
     *745,743,741,740,738,736,734,732,730,728,726,724,723,721,719,
     *717,715,713,711,710,708,706,704,702,701,699,697,695/
      DATA (Md07(i),i=1033,1113) /858,857,855,854,852,850,848,846,
     *844,842,840,838,836,834,832,829,827,825,823,821,819,817,814,
     *813,810,808,806,804,802,800,798,796,794,792,790,788,786,784,
     *782,780,778,776,774,772,770,768,766,764,762,760,759,757,755,
     *753,751,749,747,746,744,742,740,738,737,735,733,731,729,728,
     *726,724,722,721,719,717,716,714,712,710,709,707,705/
      DATA (Md07(i),i=1114,1193) /859,857,856,855,853,851,849,847,
     *845,843,841,839,837,835,833,831,829,827,825,823,821,819,817,
     *815,813,811,809,807,805,803,801,799,797,795,793,791,789,787,
     *786,784,782,780,778,776,774,773,771,769,767,765,764,762,760,
     *758,756,755,753,751,750,748,746,744,743,741,739,738,736,734,
     *733,731,729,728,726,724,723,721,719,718,716,715/
      DATA (Md07(i),i=1194,1272) /860,859,857,856,854,852,850,849,
     *846,844,843,841,839,837,835,833,831,829,827,825,823,821,819,
     *817,815,814,812,810,808,806,804,802,800,798,797,795,793,791,
     *789,787,786,784,782,780,779,777,775,773,772,770,768,767,765,
     *763,762,760,758,757,755,753,752,750,748,747,745,744,742,740,
     *739,737,736,734,733,731,729,728,726,725,723/
      DATA (Md07(i),i=1273,1350) /861,859,858,856,855,853,851,850,
     *848,846,844,842,840,838,836,835,833,831,829,827,825,823,821,
     *820,818,816,814,812,811,809,807,805,803,802,800,798,796,795,
     *793,791,789,788,786,784,783,781,779,778,776,774,773,771,770,
     *768,766,765,763,762,760,758,757,755,754,752,751,749,748,746,
     *745,743,742,740,739,737,736,734,733,731/
      DATA (Md07(i),i=1351,1427) /862,861,859,858,856,854,852,851,
     *849,847,845,843,842,840,838,836,834,833,831,829,827,825,824,
     *822,820,818,817,815,813,811,810,808,806,805,803,801,800,798,
     *796,795,793,791,790,788,787,785,783,782,780,779,777,776,774,
     *772,771,769,768,766,765,763,762,760,759,757,756,754,753,752,
     *750,749,747,746,744,743,742,740,739/
      DATA (Md07(i),i=1428,1502) /862,860,859,857,855,854,852,850,
     *848,847,845,843,841,840,838,836,834,833,831,829,828,826,824,
     *822,821,819,817,816,814,812,811,809,807,806,804,803,801,799,
     *798,796,795,793,792,790,789,787,786,784,783,781,780,778,777,
     *775,774,772,771,769,768,766,765,764,762,761,759,758,757,755,
     *754,753,751,750,748,747,746/
      DATA (Md07(i),i=1503,1576) /863,862,860,858,856,855,853,851,
     *850,848,846,845,843,841,840,838,836,835,833,831,830,828,826,
     *825,823,821,820,818,817,815,813,812,810,809,807,806,804,803,
     *801,800,798,797,795,794,792,791,789,788,786,785,783,782,781,
     *779,778,776,775,774,772,771,770,768,767,766,764,763,762,760,
     *759,758,756,755,754,752/
      DATA (Md07(i),i=1577,1648) /862,861,859,858,856,854,853,851,
     *850,848,846,845,843,841,840,838,836,835,833,832,830,828,827,
     *825,824,822,821,819,818,816,815,813,812,810,809,807,806,804,
     *803,801,800,798,797,796,794,793,791,790,789,787,786,785,783,
     *782,780,779,778,776,775,774,773,771,770,769,767,766,765,764,
     *762,761,760,759/
      DATA (Md07(i),i=1649,1718) /862,861,859,858,856,854,853,851,
     *849,848,846,845,843,842,840,838,837,835,834,832,831,829,828,
     *826,825,823,822,820,819,817,816,814,813,811,810,809,807,806,
     *804,803,802,800,799,798,796,795,794,792,791,790,788,787,786,
     *784,783,782,780,779,778,777,775,774,773,772,770,769,768,767,
     *766,764/
      DATA (Md07(i),i=1719,1785) /860,859,857,856,854,853,851,850,
     *848,846,845,843,842,840,839,837,836,834,833,831,830,828,827,
     *826,824,823,821,820,818,817,816,814,813,811,810,809,807,806,
     *805,803,802,801,799,798,797,796,794,793,792,790,789,788,787,
     *786,784,783,782,781,779,778,777,776,775,773,772,771,770/
      DATA (Md07(i),i=1786,1851) /862,860,859,857,856,854,853,851,
     *850,848,847,845,844,842,841,839,838,837,835,834,832,831,829,
     *828,827,825,824,822,821,820,818,817,816,814,813,812,810,809,
     *808,806,805,804,803,801,800,799,798,796,795,794,793,792,790,
     *789,788,787,786,784,783,782,781,780,779,778,776,775/
      DATA (Md07(i),i=1852,1913) /859,858,856,854,853,852,850,849,
     *847,846,844,843,841,840,839,837,836,834,833,832,830,829,828,
     *826,825,824,822,821,820,818,817,816,815,813,812,811,810,808,
     *807,806,805,803,802,801,800,799,797,796,795,794,793,792,790,
     *789,788,787,786,785,784,783,781,780/
      DATA (Md07(i),i=1914,1972) /858,856,855,853,852,851,849,848,
     *846,845,844,842,841,839,838,837,835,834,833,831,830,829,827,
     *826,825,824,822,821,820,819,817,816,815,814,812,811,810,809,
     *808,806,805,804,803,802,801,800,798,797,796,795,794,793,792,
     *791,790,788,787,786,785/
      DATA (Md07(i),i=1973,2026) /854,852,851,850,848,847,846,844,
     *843,842,840,839,838,836,835,834,833,831,830,829,827,826,825,
     *824,823,821,820,819,818,817,815,814,813,812,811,810,808,807,
     *806,805,804,803,802,801,800,798,797,796,795,794,793,792,791,
     *790/
      DATA (Md07(i),i=2027,2073) /848,846,845,844,843,841,840,839,
     *837,836,835,834,832,831,830,829,828,826,825,824,823,822,821,
     *819,818,817,816,815,814,813,811,810,809,808,807,806,805,804,
     *803,802,801,800,799,801,796,795,794/

!.....Startadresse für Hoehenstufen 5 - 45 m:......................

      DATA (Add07( 1,i,1),i=5,45) /0,0,1,10,21,36,53,71,92,116,142,
     *168,198,231,266,302,339,379,422,469,527,613,699,784,868,951,
     *1033,1114,1194,1273,1351,1428,1503,1577,1649,1719,1786,1852,
     *1914,1973,2027/

!.....Unterer Durchmesser einer Hoehenstufe:.......................

      DATA (Add07( 1,i,2),i=5,45) /0,0,9,9,7,7,8,8,8,8,8,8,8,8,9,10,
     *11,12,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,29,31,34,
     *35,39,42,47,54/

!.....Oberer  Durchmesser einer Hoehenstufe:.......................

      DATA (Add07( 1,i,3),i=5,45) /0,0,17,19,21,23,25,28,31,33,33,
     *37,40,42,44,46,50,54,58,70,99,100,100,100,100,100,100,100,100,
     *100,100,100,100,100,100,100,100,100,100,100,100/

!     Tanne..................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=2074,2084) /634,511,470,430,400,400,400,400,
     *400,400,400/
      DATA (Md07(i),i=2085,2098) /919,794,707,633,586,544,520,492,
     *479,458,451,437,417,406/
      DATA (Md07(i),i=2099,2114) /934,840,771,737,702,681,657,630,
     *615,609,596,590,581,575,574,569/
      DATA (Md07(i),i=2115,2132) /937,865,826,796,765,746,725,712,
     *705,693,686,683,677,672,664,660,652,652/
      DATA (Md07(i),i=2133,2152) /944,881,846,819,791,788,767,764,
     *763,747,743,735,730,727,721,716,713,708,704,701/
      DATA (Md07(i),i=2153,2174) /942,893,845,835,822,814,805,808,
     *794,790,782,776,772,764,763,758,754,749,744,740,738,733/
      DATA (Md07(i),i=2175,2197) /935,888,858,833,829,825,830,829,
     *822,815,809,799,796,794,788,784,780,777,772,770,766,763,758/
      DATA (Md07(i),i=2198,2222) /938,894,867,852,848,840,842,846,
     *836,832,824,820,814,812,807,803,798,795,792,789,785,783,779,
     *776,772/
      DATA (Md07(i),i=2223,2249) /940,898,862,856,850,848,851,858,
     *851,845,839,833,830,825,821,816,815,811,807,804,801,796,794,
     *791,787,784,781/
      DATA (Md07(i),i=2250,2279) /938,901,875,860,852,850,856,865,
     *857,854,850,845,840,836,832,828,826,821,817,814,812,808,805,
     *802,799,796,793,791,787,785/
      DATA (Md07(i),i=2280,2312) /947,916,883,869,860,857,860,867,
     *863,861,855,852,848,844,840,836,832,830,826,823,819,817,814,
     *811,808,806,803,800,797,795,792,789,786/
      DATA (Md07(i),i=2313,2348) /923,900,882,871,866,867,871,869,
     *863,859,857,851,848,844,841,839,835,833,830,827,824,821,818,
     *815,813,810,808,805,803,800,797,795,793,790,788,786/
      DATA (Md07(i),i=2349,2386) /900,879,874,867,871,874,871,867,
     *862,858,856,851,848,845,843,841,836,834,832,829,827,824,821,
     *819,816,814,811,809,807,804,802,800,798,795,793,791,788,786/
      DATA (Md07(i),i=2387,2427) /900,884,872,869,872,873,869,867,
     *864,861,858,853,851,848,846,842,840,838,835,833,831,828,826,
     *824,821,819,817,814,813,810,808,806,804,801,799,797,795,793,
     *791,789,787/
      DATA (Md07(i),i=2428,2478) /884,877,869,870,871,870,867,865,
     *862,859,855,853,850,847,845,843,841,839,836,834,832,829,827,
     *825,823,821,819,817,815,813,810,809,807,805,803,801,799,797,
     *795,794,792,790,788,786,785,783,781,779,778,776,774/
      DATA (Md07(i),i=2479,2537) /878,872,873,873,870,867,866,862,
     *860,857,853,851,849,846,845,842,840,838,836,834,832,830,828,
     *826,824,822,820,818,816,815,813,811,809,808,806,804,802,801,
     *799,797,796,794,792,791,789,787,786,784,783,781,780,778,777,
     *775,774,772,771,769,768/
      DATA (Md07(i),i=2538,2613) /875,872,870,867,866,864,861,860,
     *856,853,851,850,847,845,843,841,839,837,836,833,832,830,828,
     *827,825,823,821,820,818,816,814,813,811,810,808,806,805,804,
     *802,801,799,798,796,795,793,792,790,789,787,786,785,783,782,
     *781,779,778,777,775,774,773,772,770,769,768,766,765,764,763,
     *761,760,759,758,756,755,754,753/
      DATA (Md07(i),i=2614,2720) /869,867,865,865,863,861,859,855,
     *853,851,849,847,845,843,842,840,838,837,835,833,831,830,828,
     *827,825,823,822,821,819,817,816,814,813,812,810,809,807,806,
     *805,803,802,801,799,798,797,796,794,793,792,790,789,788,787,
     *786,784,783,782,781,780,779,777,776,775,774,773,772,771,770,
     *768,767,766,765,764,763,762,761,760,759,758,757,755,754,753,
     *752,751,750,749,748,747,746,745,744,743,742,741,740,739,738,
     *736,735,734,733,732,731,730,729,728/
      DATA (Md07(i),i=2721,2826) /865,864,864,861,860,857,855,852,
     *850,849,847,845,844,842,840,839,837,835,834,832,831,829,828,
     *827,825,824,822,821,820,818,817,816,814,813,812,811,809,808,
     *807,806,804,803,802,801,800,799,798,796,795,794,793,792,791,
     *790,789,788,787,786,785,784,783,782,781,780,779,778,777,776,
     *775,774,773,772,771,770,769,768,767,766,765,765,764,763,762,
     *761,760,759,758,757,756,755,754,753,753,752,751,750,749,748,
     *747,746,745,744,743,742,742,741/
      DATA (Md07(i),i=2827,2931) /861,861,859,857,856,854,851,850,
     *848,846,845,843,841,840,838,837,836,834,833,831,830,829,827,
     *826,825,824,823,821,820,819,818,816,816,814,813,812,811,810,
     *809,808,807,806,805,804,803,802,801,800,799,798,797,796,795,
     *794,793,792,791,790,789,789,788,787,786,785,784,783,783,782,
     *781,780,779,778,778,777,776,775,774,774,773,772,771,770,770,
     *769,768,767,766,766,765,764,763,762,762,761,760,759,759,758,
     *757,756,755,755,754,753,752/
      DATA (Md07(i),i=2932,3034) /857,855,854,852,849,847,846,845,
     *843,842,840,839,838,837,835,834,833,831,831,829,828,827,826,
     *825,824,822,821,820,819,818,817,816,815,814,813,812,811,810,
     *809,808,808,807,806,805,804,803,802,802,801,800,799,798,797,
     *797,796,795,794,794,793,792,791,791,790,789,788,788,787,786,
     *786,785,784,783,783,782,781,781,780,779,779,778,777,777,776,
     *775,775,774,773,773,772,771,771,770,769,769,768,767,767,766,
     *765,765,764,763,763/
      DATA (Md07(i),i=3035,3136) /853,851,849,847,846,845,843,842,
     *840,839,838,837,836,835,833,832,832,830,829,828,827,826,825,
     *824,823,822,821,820,819,819,818,817,816,815,814,813,812,812,
     *811,810,809,809,808,807,806,806,805,804,803,803,802,801,801,
     *800,799,799,798,797,797,796,795,795,794,794,793,792,792,791,
     *791,790,789,789,788,788,787,787,786,785,785,784,784,783,783,
     *782,781,781,780,780,779,779,778,778,777,777,776,775,775,774,
     *774,770,773,772/
      DATA (Md07(i),i=3137,3237) /850,847,845,843,843,841,840,839,
     *838,837,836,835,834,833,832,831,830,829,828,827,826,825,824,
     *823,823,822,821,820,819,819,818,817,816,815,815,814,813,813,
     *812,811,811,810,809,809,808,807,807,806,806,805,804,804,803,
     *803,802,802,801,801,800,799,799,798,798,797,797,796,796,795,
     *795,795,794,794,793,793,792,792,791,791,790,790,790,789,789,
     *788,788,787,787,786,786,786,785,785,784,784,783,783,783,782,
     *782,781,781/
      DATA (Md07(i),i=3238,3336) /842,841,840,839,838,837,836,835,
     *834,833,833,832,831,830,829,828,827,827,826,825,824,823,823,
     *822,821,821,820,819,818,818,817,816,816,815,815,814,813,813,
     *812,812,811,811,810,810,809,809,808,808,807,807,806,806,805,
     *805,804,804,803,803,803,802,802,801,801,801,800,800,799,799,
     *799,798,798,798,797,797,796,796,796,795,795,795,794,794,794,
     *793,793,793,792,792,792,791,791,791,790,790,790,789,789,789,
     *788/
      DATA (Md07(i),i=3337,3433) /838,837,836,835,834,833,833,832,
     *831,830,829,829,828,827,826,826,825,824,824,823,822,822,821,
     *820,820,819,819,818,817,817,816,816,815,815,814,814,813,813,
     *813,812,812,811,811,810,810,810,809,809,808,808,808,807,807,
     *807,806,806,806,805,805,805,804,804,804,803,803,803,803,802,
     *802,802,802,801,801,801,800,800,800,800,799,799,799,799,798,
     *798,798,798,798,797,797,797,797,796,796,796,796,795,795/
      DATA (Md07(i),i=3434,3528) /833,833,832,831,830,830,829,829,
     *828,827,826,826,825,825,824,823,823,822,822,821,821,820,819,
     *819,819,818,818,817,817,816,816,815,815,815,814,814,813,813,
     *813,812,812,812,811,811,811,811,810,810,810,809,809,809,809,
     *808,808,808,808,807,807,807,807,807,806,806,806,806,806,805,
     *805,805,805,805,804,804,804,804,804,804,803,803,803,803,803,
     *803,803,802,802,802,802,802,802,802,801,801,801/
      DATA (Md07(i),i=3529,3621) /830,829,829,828,827,826,826,825,
     *825,824,824,823,823,822,822,821,821,820,820,819,819,819,818,
     *818,817,817,817,816,816,816,815,815,815,814,814,814,813,813,
     *813,813,812,812,812,812,812,811,811,811,811,811,810,810,810,
     *810,810,810,809,809,809,809,809,809,809,809,808,808,808,808,
     *808,808,808,808,808,808,807,807,807,807,807,807,807,807,807,
     *807,807,807,807,807,806,806,806,806,806/
      DATA (Md07(i),i=3622,3711) /826,825,824,824,823,823,823,822,
     *822,821,821,820,820,820,819,819,818,818,818,817,817,817,816,
     *816,816,816,815,815,815,815,814,814,814,814,814,813,813,813,
     *813,813,813,812,812,812,812,812,812,812,812,812,811,811,811,
     *811,811,811,811,811,811,811,811,811,811,811,811,811,811,811,
     *811,811,811,811,811,811,811,811,811,811,811,811,811,811,811,
     *811,811,811,811,811,811,811/
      DATA (Md07(i),i=3712,3800) /823,822,822,821,821,821,820,820,
     *819,819,819,818,818,818,817,817,817,817,816,816,816,816,815,
     *815,815,815,815,814,814,814,814,814,814,814,813,813,813,813,
     *813,813,813,813,813,813,813,813,813,813,812,813,812,812,812,
     *812,812,812,812,813,813,813,813,813,813,813,813,813,813,813,
     *813,813,813,813,813,813,813,813,813,813,814,814,814,814,814,
     *814,814,814,814,814,814/
      DATA (Md07(i),i=3801,3886) /819,819,818,818,818,818,817,817,
     *817,816,816,816,816,815,815,815,815,815,814,814,814,814,814,
     *814,814,814,814,813,813,813,813,813,813,813,813,813,813,813,
     *813,813,813,813,813,813,813,813,813,813,813,813,813,813,814,
     *814,814,814,814,814,814,814,814,814,814,814,815,815,815,815,
     *815,815,815,815,816,816,816,816,816,816,816,817,817,817,817,
     *817,817,817/
      DATA (Md07(i),i=3887,3970) /816,816,816,815,815,815,815,815,
     *814,814,814,814,814,814,813,813,813,813,813,813,813,813,813,
     *813,813,813,813,813,813,813,813,813,813,813,813,813,813,813,
     *813,813,813,813,813,813,813,814,814,814,814,814,814,814,814,
     *814,815,815,815,815,815,815,816,816,816,816,816,816,817,817,
     *817,817,817,817,818,818,818,818,818,819,819,819,819,819,820,
     *820/
      DATA (Md07(i),i=3971,4051) /813,813,813,813,813,812,812,812,
     *812,812,812,812,812,812,812,812,812,812,812,812,812,812,812,
     *812,812,812,812,812,812,812,812,812,812,812,812,813,813,813,
     *813,813,813,813,813,814,814,814,814,814,814,815,815,815,815,
     *815,816,816,816,816,816,817,817,817,817,818,818,818,818,819,
     *819,819,819,820,820,820,820,821,821,821,821,822,822/
      DATA (Md07(i),i=4052,4128) /810,810,810,810,810,810,810,810,
     *810,810,810,810,810,810,810,810,810,810,810,810,810,810,811,
     *811,811,811,811,811,811,811,812,812,812,812,812,813,813,813,
     *813,813,813,814,814,814,814,815,815,815,815,816,816,816,816,
     *817,817,817,817,818,818,818,818,819,819,819,820,820,820,820,
     *821,821,821,822,822,822,823,823,823/
      DATA (Md07(i),i=4129,4201) /808,808,808,808,808,808,808,808,
     *808,808,808,808,808,809,809,809,809,809,809,809,809,810,810,
     *810,810,810,811,811,811,811,811,812,812,812,812,813,813,813,
     *813,814,814,814,814,815,815,815,816,816,816,817,817,817,817,
     *818,818,818,819,819,819,820,820,820,821,821,821,822,822,822,
     *823,823,823,824,824/
      DATA (Md07(i),i=4202,4267) /806,806,806,806,806,807,807,807,
     *807,807,807,808,808,808,808,808,809,809,809,809,809,810,810,
     *810,810,811,811,811,812,812,812,813,813,813,813,814,814,814,
     *815,815,815,816,816,816,817,817,817,818,818,819,819,819,820,
     *820,820,821,821,822,822,822,823,823,824,824,824,825/
      DATA (Md07(i),i=4268,4325) /805,805,805,806,806,806,806,807,
     *807,807,807,808,808,808,809,809,809,809,810,810,810,811,811,
     *811,812,812,812,813,813,814,814,814,815,815,815,816,816,817,
     *817,817,818,818,819,819,819,820,820,821,821,821,822,822,823,
     *823,824,824,824,825/
      DATA (Md07(i),i=4326,4371) /806,807,807,807,807,808,808,809,
     *809,809,810,810,810,811,811,812,812,812,813,813,814,814,814,
     *815,815,816,816,816,817,817,818,818,819,819,820,820,820,821,
     *821,822,822,823,823,824,824,825/
      DATA (Md07(i),i=4372,4403) /810,810,810,811,811,812,812,813,
     *813,813,814,814,815,815,816,816,817,817,818,818,819,819,820,
     *820,821,821,821,822,822,823,823,824/
      DATA (Md07(i),i=4404,4420) /815,815,816,816,817,817,818,818,
     *819,819,820,820,821,821,822,823,823/

      DATA (Add07( 2,i,1),i=5,45) /0,0,2074,2085,2099,2115,2133,2153,
     *2175,2198,2223,2250,2280,2313,2349,2387,2428,2479,2538,2614,
     *2721,2827,2932,3035,3137,3238,3337,3434,3529,3622,3712,3801,
     *3887,3971,4052,4129,4202,4268,4326,4372,4404/

      DATA (Add07( 2,i,2),i=5,45) /0,0,9,8,8,8,8,8,8,8,8,8,8,9,10,
     *10,11,12,13,14,15,16,18,19,20,22,24,26,28,31,32,35,37,40,44,
     *48,55,63,75,89,104/

      DATA (Add07( 2,i,3),i=5,45) /0,0,19,21,23,25,27,29,30,32,34,
     *37,40,44,47,50,61,70,88,120,120,120,120,120,120,120,120,120,
     *120,120,120,120,120,120,120,120,120,120,120,120,120/

!     Douglasie..............................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=4421,4426) /451,487,480,446,401,400/
      DATA (Md07(i),i=4427,4435) /638,632,602,597,580,557,528,513,
     *494/
      DATA (Md07(i),i=4436,4448) /772,736,694,663,650,648,642,625,
     *605,596,581,574,563/
      DATA (Md07(i),i=4449,4464) /544,758,727,704,695,680,676,662,
     *660,651,639,632,621,614,604,599/
      DATA (Md07(i),i=4465,4481) /804,762,745,723,713,708,698,696,
     *688,675,667,664,656,645,638,633,631/
      DATA (Md07(i),i=4482,4510) /828,797,775,755,733,731,723,720,
     *711,698,691,685,683,677,668,662,661,662,664,653,651,644,641,
     *635,632,628,623,621,619/
      DATA (Md07(i),i=4511,4541) /848,815,798,776,756,755,747,743,
     *734,721,712,706,701,694,689,681,683,681,677,675,670,665,660,
     *655,654,651,648,642,641,636,628/
      DATA (Md07(i),i=4542,4574) /861,834,811,794,775,764,758,753,
     *744,737,728,720,714,706,704,698,697,694,692,688,682,680,676,
     *668,667,662,657,654,652,645,643,644,642/
      DATA (Md07(i),i=4575,4607) /881,846,827,809,791,779,772,767,
     *758,750,740,732,725,720,716,712,710,708,704,698,694,690,685,
     *680,675,671,669,666,662,659,653,655,651/
      DATA (Md07(i),i=4608,4642) /899,865,838,811,796,785,779,773,
     *764,757,751,742,738,731,726,724,722,719,713,709,705,700,695,
     *690,685,680,676,673,669,670,665,664,662,661,658/
      DATA (Md07(i),i=4643,4678) /874,848,825,807,792,784,778,770,
     *767,760,754,748,741,737,734,731,728,724,721,716,710,706,700,
     *696,691,687,686,685,681,680,679,678,673,675,674,673/
      DATA (Md07(i),i=4679,4716) /888,856,827,810,797,789,784,780,
     *775,768,761,755,750,748,743,741,737,734,729,725,719,714,709,
     *706,701,698,697,695,694,690,691,692,688,689,688,687,687,687/
      DATA (Md07(i),i=4717,4756) /901,863,837,816,801,795,793,787,
     *783,775,771,767,761,757,751,748,745,742,738,733,728,723,718,
     *715,711,708,707,705,703,700,698,698,697,698,696,697,688,696,
     *698,699/
      DATA (Md07(i),i=4757,4797) /871,844,825,810,803,801,795,790,
     *786,780,775,770,765,761,758,754,750,747,742,737,732,727,724,
     *720,717,716,713,711,708,706,705,705,706,705,704,704,706,705,
     *706,707,706/
      DATA (Md07(i),i=4798,4840) /879,851,834,818,811,808,802,797,
     *791,788,782,779,773,769,766,762,758,753,749,744,739,735,732,
     *729,726,723,721,718,716,714,713,711,710,709,708,708,709,707,
     *707,710,709,708,708/
      DATA (Md07(i),i=4841,4884) /859,840,826,818,814,808,803,799,
     *795,789,786,782,777,773,768,764,760,756,752,747,743,739,736,
     *733,731,728,726,724,721,720,718,716,716,716,714,713,714,713,
     *712,713,713,713,712,712/
      DATA (Md07(i),i=4885,4929) /847,833,825,821,814,807,804,801,
     *797,792,788,784,780,776,772,767,763,759,755,751,747,743,740,
     *738,735,733,731,728,726,725,723,721,720,719,718,716,717,716,
     *715,715,715,716,715,714,715/
      DATA (Md07(i),i=4930,4976) /852,840,831,826,819,815,811,807,
     *802,798,796,791,787,783,778,773,770,766,762,758,754,750,747,
     *744,741,739,737,734,732,730,728,727,725,723,723,721,722,720,
     *720,720,719,719,718,719,719,719,719/
      DATA (Md07(i),i=4977,5025) /836,831,826,821,816,813,809,806,
     *802,798,794,790,784,779,775,771,768,764,760,757,753,750,747,
     *745,743,740,738,736,734,732,730,729,727,726,727,724,723,723,
     *723,723,722,723,723,723,723,722,722,723,722/
      DATA (Md07(i),i=5026,5074) /832,827,821,817,813,811,807,803,
     *799,795,790,785,781,777,773,770,766,763,759,756,753,751,748,
     *746,744,742,740,738,735,734,732,731,730,729,728,728,728,726,
     *727,728,728,727,727,728,727,728,728,729,728/
      DATA (Md07(i),i=5075,5124) /832,827,822,819,816,812,808,804,
     *799,795,791,787,783,779,775,772,768,765,762,759,756,754,752,
     *750,747,745,743,740,739,737,736,735,734,734,733,733,734,734,
     *735,734,735,735,735,734,735,735,736,736,737,737/
      DATA (Md07(i),i=5125,5176) /827,824,820,816,813,808,804,800,
     *796,792,788,785,781,777,774,771,768,765,762,759,757,755,752,
     *750,747,745,743,742,741,740,739,739,738,738,739,739,739,740,
     *740,741,741,742,743,743,742,742,743,742,743,743,743,744/
      DATA (Md07(i),i=5177,5229) /828,824,821,817,813,808,804,800,
     *797,793,790,786,782,779,776,773,770,767,764,762,759,757,755,
     *753,750,749,747,746,745,744,744,743,743,744,744,744,744,745,
     *745,746,747,747,746,746,747,747,747,748,748,748,748,749,750/
      DATA (Md07(i),i=5230,5310) /832,829,825,821,817,813,809,805,
     *801,797,794,791,787,784,781,778,775,772,769,767,764,761,759,
     *757,755,754,752,751,750,749,749,748,748,748,748,748,748,749,
     *749,750,751,751,751,752,752,752,752,753,753,754,754,754,754,
     *755,755,756,757,757,758,758,759,759,760,761,761,762,762,763,
     *763,764,764,765,765,765,766,766,767,767,768,768,769/
      DATA (Md07(i),i=5311,5390) /834,829,825,821,817,813,809,806,
     *802,799,796,792,789,786,783,780,777,774,771,769,766,764,762,
     *760,759,757,756,755,754,753,753,753,752,752,752,752,753,753,
     *754,754,755,756,756,756,757,757,757,758,758,759,760,761,761,
     *762,762,763,764,764,765,764,765,765,766,766,767,767,767,767,
     *768,768,769,769,770,770,770,771,771,772,772,773/
      DATA (Md07(i),i=5391,5469) /834,829,825,821,818,814,811,807,
     *804,800,797,794,791,787,784,781,778,775,773,771,769,767,765,
     *764,762,761,759,758,758,757,757,756,756,756,756,757,757,758,
     *758,759,760,760,761,761,761,762,762,763,763,764,764,765,766,
     *766,767,767,768,768,768,769,769,769,770,770,771,771,772,772,
     *773,773,773,773,774,774,775,775,775,776,776/
      DATA (Md07(i),i=5470,5547) /833,829,826,822,818,815,811,808,
     *805,801,798,795,792,789,785,783,780,778,776,774,772,770,768,
     *767,765,764,763,762,761,761,760,760,760,760,760,761,761,762,
     *763,763,764,764,765,765,766,766,767,767,767,768,769,769,770,
     *770,771,771,772,772,773,773,774,774,775,775,776,776,777,777,
     *778,778,778,778,778,779,779,780,780,781/
      DATA (Md07(i),i=5548,5624) /833,830,826,823,820,816,812,809,
     *805,802,799,796,792,790,787,784,782,780,778,776,774,772,771,
     *769,768,767,766,765,765,764,764,764,764,764,765,765,766,766,
     *767,767,768,768,769,769,770,770,771,772,773,773,774,775,775,
     *776,776,777,777,777,778,779,779,780,780,780,781,782,782,783,
     *783,784,784,784,785,785,786,786,787/
      DATA (Md07(i),i=5625,5699) /830,827,824,820,817,813,810,806,
     *803,800,797,794,791,789,786,784,782,780,778,777,775,773,772,
     *771,770,769,768,768,767,767,767,768,768,768,769,770,770,771,
     *771,772,772,773,773,774,774,775,776,776,777,777,778,779,779,
     *780,780,781,781,782,782,783,783,784,784,784,785,785,786,786,
     *787,787,788,788,788,789,789/
      DATA (Md07(i),i=5700,5774) /835,831,828,824,821,818,814,811,
     *807,804,801,798,795,793,790,788,786,784,782,780,779,777,776,
     *774,773,773,772,771,771,771,771,771,771,772,772,772,773,774,
     *774,775,775,776,776,777,777,778,778,779,780,780,781,781,782,
     *782,783,783,784,784,785,785,786,786,787,787,788,788,788,789,
     *789,790,790,791,791,791,792/
      DATA (Md07(i),i=5775,5849) /841,837,833,829,826,822,819,815,
     *812,808,805,802,800,797,794,792,790,788,786,784,782,781,779,
     *778,777,776,775,775,774,774,774,774,774,775,775,775,776,776,
     *777,777,778,778,779,779,780,780,781,781,782,783,783,784,784,
     *785,785,786,786,787,787,788,788,788,789,789,790,790,791,791,
     *792,792,792,793,793,794,794/
      DATA (Md07(i),i=5850,5924) /846,843,838,834,831,827,823,820,
     *816,813,810,806,804,801,798,796,794,792,789,787,786,784,783,
     *781,782,779,779,778,778,777,777,777,777,778,778,778,779,779,
     *779,780,780,781,781,781,782,783,783,784,784,785,785,786,786,
     *787,787,788,788,789,789,790,790,791,791,792,792,792,793,793,
     *794,794,795,795,795,796,796/
      DATA (Md07(i),i=5925,5999) /851,846,843,839,836,832,828,824,
     *821,817,814,811,808,805,803,800,797,795,793,791,789,787,786,
     *784,783,782,782,781,781,780,780,780,780,780,781,781,781,782,
     *782,782,782,783,783,784,784,785,785,786,786,787,787,788,788,
     *789,789,790,790,791,791,792,792,793,793,794,794,794,795,795,
     *796,796,797,797,797,798,798/
      DATA (Md07(i),i=6000,6072) /848,843,840,836,833,829,826,822,
     *819,816,813,810,807,804,801,798,796,794,792,790,789,787,786,
     *785,785,784,784,783,783,783,783,783,783,783,784,784,784,784,
     *785,785,785,786,786,786,787,788,788,789,789,790,790,791,791,
     *792,792,793,793,794,794,795,795,795,796,796,797,797,798,798,
     *798,799,799,800,800/
      DATA (Md07(i),i=6073,6143) /845,841,838,834,831,827,824,820,
     *817,814,811,808,805,802,799,797,795,793,792,790,789,788,788,
     *787,787,786,786,786,786,786,786,786,786,786,786,786,787,787,
     *787,788,788,788,789,789,790,790,791,791,792,792,793,793,794,
     *794,795,795,796,796,797,797,798,798,798,799,799,800,800,800,
     *801,801,802/
      DATA (Md07(i),i=6144,6212) /845,841,837,833,829,826,822,819,
     *815,812,808,805,802,800,798,796,795,793,792,791,790,790,789,
     *789,789,789,789,788,788,788,788,788,788,789,789,789,790,790,
     *790,790,790,791,791,792,792,793,793,794,794,795,795,796,796,
     *797,797,798,798,799,799,800,800,800,801,801,802,802,802,803,
     *803/
      DATA (Md07(i),i=6213,6279) /843,838,835,831,828,824,820,816,
     *812,809,806,804,801,800,798,797,795,794,793,793,792,792,792,
     *791,791,791,791,791,791,791,791,791,791,791,792,792,792,793,
     *793,793,793,793,794,794,795,795,796,796,797,797,798,798,799,
     *799,800,800,800,801,801,802,802,803,803,803,804,804,805/
      DATA (Md07(i),i=6280,6344) /840,835,831,827,824,819,816,812,
     *811,807,805,803,801,800,798,797,796,795,795,795,794,794,794,
     *794,794,793,793,793,793,793,794,794,794,794,795,795,796,796,
     *797,797,797,798,797,798,798,798,798,799,799,800,800,800,801,
     *801,802,802,803,803,803,804,804,805,805,806,806/
      DATA (Md07(i),i=6345,6408) /838,834,829,826,822,819,815,814,
     *812,809,807,805,803,801,800,799,798,797,797,797,796,796,796,
     *796,796,796,796,796,796,796,796,796,797,797,798,798,798,798,
     *798,799,799,800,800,801,801,802,802,802,802,802,802,803,802,
     *803,803,804,804,805,805,805,806,806,807,807/

      DATA (Add07( 3,i,1),i=5,45) /0,0,4421,4427,4436,4449,4465,4482,
     *4511,4542,4575,4608,4643,4679,4717,4757,4798,4841,4885,4930,
     *4977,5026,5075,5125,5177,5230,5311,5391,5470,5548,5625,5700,
     *5775,5850,5925,6000,6073,6144,6213,6280,6345/

      DATA (Add07( 3,i,2),i=5,45) /0,0,10,9,7,7,8,8,8,8,8,8,9,9,9,
     *10,10,11,12,12,14,16,17,19,20,20,21,22,23,24,26,26,26,26,26,
     *28,30,32,34,36,37/

      DATA (Add07( 3,i,3),i=5,45) /0,0,15,17,19,22,24,36,38,40,40,
     *42,44,46,48,50,52,54,56,58,62,64,66,70,72,100,100,100,100,100,
     *100,100,100,100,100,100,100,100,100,100,100/

!     Kiefer..................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=6409,6420) /660,620,580,540,500,460,420,400,
     *400,400,400,400/
      DATA (Md07(i),i=6421,6432) /571,530,490,450,410,400,400,400,
     *400,400,400,400/
      DATA (Md07(i),i=6433,6445) /465,794,620,541,495,476,467,453,
     *444,432,429,420,430/
      DATA (Md07(i),i=6446,6459) /853,711,648,590,561,550,536,531,
     *524,532,535,540,539,549/
      DATA (Md07(i),i=6460,6478) /882,761,703,657,624,613,599,595,
     *597,595,604,606,617,622,627,634,640,647,654/
      DATA (Md07(i),i=6479,6505) /899,800,731,692,675,658,646,641,
     *641,646,651,658,664,671,672,679,685,691,697,703,708,714,716,
     *721,726,731,735/
      DATA (Md07(i),i=6506,6533) /910,826,765,725,704,688,677,673,
     *672,674,684,686,690,697,700,708,715,721,726,730,735,739,743,
     *746,752,755,760,762/
      DATA (Md07(i),i=6534,6568) /925,843,785,743,723,709,702,696,
     *694,700,701,707,713,717,722,728,731,737,742,746,751,757,759,
     *764,768,771,775,777,780,784,787,789,792,795,797/
      DATA (Md07(i),i=6569,6607) /857,794,764,737,724,711,705,710,
     *709,714,718,725,731,733,740,744,749,751,757,761,763,767,771,
     *776,779,783,786,789,791,795,798,800,803,805,807,810,812,814,
     *816/
      DATA (Md07(i),i=6608,6652) /863,807,780,755,736,729,720,717,
     *719,723,729,734,738,741,746,750,756,760,764,767,771,774,777,
     *781,784,787,790,793,795,797,801,804,806,808,811,813,814,817,
     *819,820,822,824,826,827,829/
      DATA (Md07(i),i=6653,6699) /819,785,769,752,738,724,726,728,
     *730,734,738,743,747,750,755,758,762,766,770,774,777,780,783,
     *786,789,792,794,797,799,802,804,806,808,810,812,815,817,819,
     *820,822,824,825,827,829,830,831,833/
      DATA (Md07(i),i=6700,6748) /829,789,769,749,737,728,730,731,
     *735,739,741,745,750,754,758,762,765,768,772,776,778,781,784,
     *787,789,792,795,798,800,802,804,806,808,810,812,813,815,817,
     *818,820,821,824,825,827,828,829,830,832,833/
      DATA (Md07(i),i=6749,6800) /837,797,773,756,735,733,733,734,
     *737,740,743,749,752,756,760,763,766,770,772,776,779,782,784,
     *786,789,791,794,796,798,800,802,804,806,808,810,811,813,815,
     *816,818,819,820,822,823,824,826,827,828,829,831,832,833/
      DATA (Md07(i),i=6801,6854) /797,773,757,744,735,736,736,739,
     *743,746,750,752,757,759,763,766,770,773,775,778,781,783,786,
     *788,790,792,795,796,798,801,802,804,806,807,809,810,812,813,
     *815,816,818,819,820,821,822,823,824,826,826,827,828,829,830,
     *831/
      DATA (Md07(i),i=6855,6910) /780,763,749,739,736,739,741,744,
     *746,751,755,757,761,764,767,769,772,775,778,780,783,784,787,
     *789,791,793,795,797,798,800,802,803,805,806,808,809,811,812,
     *813,814,815,817,818,819,820,821,822,823,824,825,825,826,827,
     *828,829,829/
      DATA (Md07(i),i=6911,6967) /761,749,740,740,739,743,746,749,
     *751,755,758,761,763,767,769,772,775,777,780,782,784,786,788,
     *790,792,793,795,797,799,800,801,803,804,806,807,808,809,811,
     *812,813,814,815,816,816,818,818,819,820,821,822,822,823,824,
     *825,825,826,827/
      DATA (Md07(i),i=6968,7028) /759,744,741,743,744,746,749,752,
     *755,759,761,764,767,770,772,774,776,779,781,783,785,787,789,
     *790,792,794,795,797,798,800,801,802,803,804,806,807,808,809,
     *810,811,812,813,814,815,815,816,817,818,818,819,820,820,821,
     *822,822,823,823,824,825,825,826/
      DATA (Md07(i),i=7029,7114) /745,742,744,744,748,750,753,756,
     *759,762,764,767,770,772,774,777,779,781,783,785,786,788,790,
     *791,793,794,796,797,798,799,800,802,803,804,805,806,807,808,
     *809,809,810,811,812,813,813,814,815,815,816,817,817,818,818,
     *819,820,820,821,821,822,822,822,823,823,824,824,825,825,825,
     *826,826,826,827,827,827,828,828,828,828,829,829,829,830,830,
     *830,830,831/
      DATA (Md07(i),i=7115,7199) /744,745,746,749,751,755,757,761,
     *763,765,768,770,772,775,777,779,781,782,784,786,788,789,791,
     *792,793,794,796,797,798,799,800,801,802,803,804,805,806,807,
     *808,808,809,810,811,811,812,812,813,814,814,815,815,816,816,
     *817,817,818,818,819,819,819,820,820,821,821,821,821,822,822,
     *822,823,823,823,824,824,824,824,825,825,825,825,825,826,826,
     *826,826/
      DATA (Md07(i),i=7200,7283) /748,748,751,754,756,759,762,764,
     *766,769,771,773,775,777,779,781,783,784,786,787,789,790,792,
     *793,794,795,796,797,799,799,800,801,802,803,804,805,806,806,
     *807,808,808,809,809,810,811,811,812,812,813,813,814,814,815,
     *815,815,816,816,816,817,817,817,818,818,818,819,819,819,819,
     *820,820,820,820,821,821,821,821,821,822,822,822,822,822,822,
     *823/
      DATA (Md07(i),i=7284,7366) /751,754,756,757,760,763,766,768,
     *770,772,774,776,778,780,782,783,785,786,788,789,790,792,793,
     *794,795,796,797,798,799,800,801,802,802,803,804,805,805,806,
     *806,807,808,808,809,809,810,810,811,811,812,812,812,813,813,
     *813,814,814,814,815,815,815,816,816,816,816,817,817,817,817,
     *817,818,818,818,818,818,818,819,819,819,819,819,819,819,819/
      DATA (Md07(i),i=7367,7448) /754,757,760,763,765,767,769,771,
     *773,775,778,779,781,783,784,785,787,788,790,791,792,793,794,
     *795,796,797,798,799,800,801,801,802,803,804,804,805,805,806,
     *806,807,807,808,808,809,809,810,810,811,811,811,811,812,812,
     *812,813,813,813,814,814,814,814,814,815,815,815,815,815,815,
     *816,816,816,816,816,816,816,816,817,817,817,817,817,817/
      DATA (Md07(i),i=7449,7529) /760,762,764,767,769,771,773,775,
     *777,779,781,782,784,785,787,788,789,791,792,793,794,795,796,
     *797,798,799,799,800,801,802,802,803,804,804,805,805,806,806,
     *807,807,808,808,808,809,809,810,810,810,810,811,811,811,812,
     *812,812,812,812,813,813,813,813,813,813,814,814,814,814,814,
     *814,814,814,814,815,815,815,815,815,815,815,815,815/
      DATA (Md07(i),i=7530,7609) /764,767,769,771,774,775,777,779,
     *781,782,784,786,787,788,789,791,792,793,794,795,796,797,798,
     *799,799,800,801,802,802,803,803,804,805,805,806,806,807,807,
     *807,808,808,808,809,809,809,810,810,810,810,811,811,811,811,
     *812,812,812,812,812,812,812,813,813,813,813,813,813,813,813,
     *813,813,813,813,813,814,814,814,814,814,814,814/
      DATA (Md07(i),i=7610,7688) /769,772,774,776,778,780,781,783,
     *784,786,787,789,790,791,793,793,795,796,797,797,798,799,800,
     *801,801,802,803,803,804,804,805,805,806,806,807,807,808,808,
     *808,809,809,809,809,810,810,810,810,811,811,811,811,811,811,
     *812,812,812,812,812,812,812,812,812,813,813,813,813,813,813,
     *813,813,813,813,813,813,813,813,813,813,813/
      DATA (Md07(i),i=7689,7766) /775,777,779,780,782,784,785,787,
     *788,790,791,792,793,794,795,796,798,798,799,800,801,801,802,
     *803,803,804,805,805,806,806,807,807,807,808,808,809,809,809,
     *809,810,810,810,810,811,811,811,811,811,811,812,812,812,812,
     *812,812,812,812,812,812,812,812,813,813,813,813,813,813,813,
     *813,813,813,812,812,812,812,812,812,812/
      DATA (Md07(i),i=7767,7843) /780,781,783,785,786,788,789,791,
     *792,793,794,795,796,798,799,799,800,801,802,803,803,804,804,
     *805,806,806,807,807,808,808,808,809,809,810,810,810,810,811,
     *811,811,811,811,812,812,812,812,812,812,812,813,813,813,813,
     *813,813,813,813,813,813,813,813,813,813,813,813,813,813,813,
     *813,813,813,813,812,812,812,812,812/
      DATA (Md07(i),i=7844,7919) /785,786,788,789,791,792,793,795,
     *796,797,798,799,800,801,802,803,803,804,805,805,806,807,807,
     *808,808,809,809,809,810,810,810,811,811,811,812,812,812,812,
     *813,813,813,813,813,813,813,813,814,814,814,814,814,814,814,
     *814,814,814,814,814,814,814,814,814,814,814,814,813,813,813,
     *813,813,813,813,813,813,813,813/
      DATA (Md07(i),i=7920,7994) /789,790,792,793,795,796,797,798,
     *799,801,801,803,803,804,805,806,806,807,808,808,809,809,810,
     *810,811,811,811,812,812,812,813,813,813,813,814,814,814,814,
     *814,814,815,815,815,815,815,815,815,815,815,815,815,815,815,
     *815,815,815,815,815,815,815,815,815,815,815,814,814,814,814,
     *814,814,814,814,813,813,813/
      DATA (Md07(i),i=7995,8068) /794,795,796,798,799,800,801,802,
     *803,804,805,806,807,808,808,809,810,810,811,811,812,812,813,
     *813,813,814,814,814,815,815,815,815,815,816,816,816,816,816,
     *816,816,817,817,817,817,817,817,817,817,817,817,817,817,817,
     *817,816,816,816,816,816,816,816,816,816,816,816,815,815,815,
     *815,815,815,815,814,814/
      DATA (Md07(i),i=8069,8141) /799,800,801,802,803,804,806,806,
     *807,808,809,810,810,811,812,812,813,813,814,814,815,815,815,
     *816,816,816,817,817,817,817,818,818,818,818,818,818,818,818,
     *819,819,819,819,819,819,819,819,819,819,819,819,818,818,818,
     *818,818,818,818,818,818,818,817,817,817,817,817,817,817,816,
     *816,816,816,816,816/
      DATA (Md07(i),i=8142,8213) /803,805,806,807,808,809,810,810,
     *811,812,813,813,814,815,815,816,816,817,817,818,818,818,819,
     *819,819,819,820,820,820,820,820,820,820,821,821,821,821,821,
     *821,821,821,821,821,821,821,821,821,821,821,820,820,820,820,
     *820,820,820,820,820,819,819,819,819,819,819,818,818,818,818,
     *818,818,817,817/
      DATA (Md07(i),i=8214,8278) /814,814,815,816,817,817,818,818,
     *819,819,820,820,820,821,821,821,822,822,822,822,822,823,823,
     *823,823,823,823,823,823,823,823,823,823,823,823,823,823,823,
     *823,823,823,823,823,823,822,822,822,822,822,822,822,821,821,
     *821,821,821,821,820,820,820,820,820,819,819,819/
      DATA (Md07(i),i=8279,8338) /820,821,822,822,822,823,823,823,
     *824,824,824,825,825,825,825,825,826,826,826,826,826,826,826,
     *826,826,826,826,826,826,826,826,826,826,826,825,825,825,825,
     *825,825,825,825,824,824,824,824,824,824,823,823,823,823,823,
     *822,822,822,822,821,821,821/
      DATA (Md07(i),i=8339,8390) /827,827,827,828,828,828,828,828,
     *829,829,829,829,829,829,829,829,829,829,829,829,829,829,828,
     *828,828,828,828,828,828,828,828,827,827,827,827,827,827,826,
     *826,826,826,826,825,825,825,825,824,824,824,824,823,823/
      DATA (Md07(i),i=8391,8436) /831,832,832,832,832,832,832,832,
     *832,832,832,832,832,832,832,832,831,831,831,831,831,831,831,
     *830,830,830,830,830,830,829,829,829,829,829,828,828,828,828,
     *827,827,827,827,826,826,826,826/
      DATA (Md07(i),i=8437,8470) /835,835,835,835,834,834,834,834,
     *834,834,834,833,833,833,833,833,832,832,832,832,832,831,831,
     *831,831,830,830,830,830,829,829,829,828,828/
      DATA (Md07(i),i=8471,8496) /837,837,837,836,836,836,836,836,
     *835,835,835,835,834,834,834,834,833,833,833,833,832,832,832,
     *831,831,831/
      DATA (Md07(i),i=8497,8512) /838,838,838,837,837,837,836,836,
     *836,836,835,835,835,834,834,834/

      DATA (Add07( 4,i,1),i=5,45) /0,0,6409,6421,6433,6446,6460,6479,
     *6506,6534,6569,6608,6653,6700,6749,6801,6855,6911,6968,7029,
     *7115,7200,7284,7367,7449,7530,7610,7689,7767,7844,7920,7995,
     *8069,8142,8214,8279,8339,8391,8437,8471,8497/

      DATA (Add07( 4,i,2),i=5,45) /0,0,7,7,7,8,8,8,8,8,9,9,10,10,10,
     *11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,36,
     *41,49,55,67,75,85/

      DATA (Add07( 4,i,3),i=5,45) /0,0,18,18,19,21,26,34,35,42,47,
     *53,56,58,61,64,67,69,74,100,100,100,100,100,100,100,100,100,
     *100,100,100,100,100,100,100,100,100,100,100,100,100/

!     Weymouthskiefer.........................(Aenderung 11.12.90 Kublin)

      DATA (Md07(i),i=8513,8521) /809,760,730,710,690,670,650,630,
     *610/
      DATA (Md07(i),i=8522,8528) /620,572,553,522,509,493,491/
      DATA (Md07(i),i=8529,8542) /465,794,730,695,701,704,694,673,
     *647,617,580,545,508,474/
      DATA (Md07(i),i=8543,8555) /881,818,775,764,748,744,731,708,
     *696,664,643,629,601/
      DATA (Md07(i),i=8556,8573) /943,868,810,800,787,768,757,743,
     *731,724,712,704,692,683,671,656,645,632/
      DATA (Md07(i),i=8574,8591) /990,921,854,826,812,797,788,772,
     *770,760,754,746,738,726,717,707,694,683/
      DATA (Md07(i),i=8592,8608) /941,884,841,832,817,808,792,796,
     *785,783,776,766,757,751,744,733,721/
      DATA (Md07(i),i=8609,8628) /878,864,838,828,812,807,802,798,
     *794,782,773,765,755,746,740,736,731,731,729,729/
      DATA (Md07(i),i=8629,8648) /916,888,870,849,838,825,818,808,
     *798,790,783,777,769,765,759,756,753,749,745,740/
      DATA (Md07(i),i=8649,8668) /948,915,893,870,853,843,825,818,
     *807,798,790,785,785,781,777,772,768,762,756,750/
      DATA (Md07(i),i=8669,8708) /986,956,918,898,878,861,847,829,
     *817,810,806,804,802,796,791,786,779,773,768,763,760,757,755,
     *754,752,751,749,749,749,750,752,753,753,752,751,749,747,750,
     *752,753/
      DATA (Md07(i),i=8709,8747) /981,939,909,887,870,855,841,834,
     *827,825,821,817,811,803,795,789,782,774,769,764,762,760,757,
     *754,754,753,753,755,751,751,756,754,752,755,751,753,754,754,
     *754/
      DATA (Md07(i),i=8748,8782) /883,871,861,851,845,841,834,826,
     *817,808,799,793,787,780,775,771,768,765,762,759,758,759,759,
     *758,757,755,758,761,757,758,758,758,757,756,759/
      DATA (Md07(i),i=8783,8827) /905,895,885,875,864,855,843,834,
     *824,816,807,801,794,788,783,777,773,769,767,764,764,764,762,
     *766,763,759,761,761,761,761,760,762,760,758,759,760,760,759,
     *759,761,762,763,764,764,764/
      DATA (Md07(i),i=8828,8872) /931,917,904,890,879,867,854,840,
     *828,818,812,805,799,794,789,784,780,777,773,771,770,767,770,
     *766,767,768,767,767,765,767,765,766,763,763,763,762,761,763,
     *764,764,765,765,765,766,765/
      DATA (Md07(i),i=8873,8912) /873,860,849,838,829,820,813,807,
     *802,796,792,788,784,781,779,775,777,778,772,772,771,769,771,
     *769,770,770,770,769,768,766,765,765,766,766,766,767,766,768,
     *766,766/
      DATA (Md07(i),i=8913,8957) /886,873,861,849,840,831,824,817,
     *810,804,800,796,793,789,785,780,780,779,778,776,778,775,776,
     *776,775,774,773,771,772,770,770,770,769,768,769,770,768,768,
     *768,768,769,768,769,769,769/
      DATA (Md07(i),i=8958,9002) /898,884,872,860,849,841,832,826,
     *818,813,807,802,799,796,791,790,788,786,784,785,781,781,780,
     *779,778,776,777,774,774,773,773,774,772,773,773,773,772,772,
     *771,771,772,772,772,771,771/
      DATA (Md07(i),i=9003,9047) /911,896,882,870,860,851,842,835,
     *827,821,816,811,810,803,801,799,796,793,789,789,788,787,785,
     *783,783,780,780,779,778,777,778,778,778,778,777,776,775,775,
     *775,775,775,774,775,774,774/
      DATA (Md07(i),i=9048,9097) /921,905,892,880,869,860,852,846,
     *838,832,827,820,818,815,811,808,803,799,798,797,795,792,789,
     *786,786,785,784,782,783,783,783,782,781,780,781,779,779,779,
     *778,779,778,778,778,778,777,778,779,779,779,779/
      DATA (Md07(i),i=9098,9147) /930,913,899,886,877,868,861,853,
     *846,839,837,828,825,821,816,812,807,805,803,801,798,798,794,
     *793,792,790,788,788,788,787,786,785,785,785,785,784,784,783,
     *783,783,783,782,782,781,781,781,781,782,781,781/
      DATA (Md07(i),i=9148,9187) /846,842,837,831,826,820,814,811,
     *808,808,804,803,799,797,798,795,795,794,793,791,790,790,789,
     *789,788,788,787,787,786,786,786,785,785,785,785,784,785,785,
     *785,785/
      DATA (Md07(i),i=9188,9227) /855,849,843,837,835,828,821,818,
     *814,813,811,809,807,804,803,800,799,798,798,798,795,794,793,
     *792,792,792,792,791,790,790,790,790,789,789,789,789,789,788,
     *788,788/
      DATA (Md07(i),i=9228,9267) /868,861,854,847,843,836,828,824,
     *822,820,818,815,812,811,809,808,806,804,803,802,801,800,798,
     *798,797,796,795,795,795,795,794,793,794,794,793,793,793,792,
     *792,792/
      DATA (Md07(i),i=9268,9292) /813,811,810,808,809,807,805,804,
     *803,802,801,801,800,799,800,799,798,798,798,798,798,798,797,
     *797,797/
      DATA (Md07(i),i=9293,9317) /819,817,816,814,814,811,810,809,
     *808,807,807,806,805,805,805,804,804,803,803,803,802,802,802,
     *802,802/
      DATA (Md07(i),i=9318,9342) /826,824,822,821,820,817,816,814,
     *814,813,812,812,811,811,810,809,809,809,808,808,807,807,806,
     *806,806/
      DATA (Md07(i),i=9343,9357) /818,817,817,816,816,815,815,814,
     *813,813,812,812,812,811,811/
      DATA (Md07(i),i=9358,9372) /823,822,822,822,821,820,820,819,
     *819,818,817,817,817,817,816/
      DATA (Md07(i),i=9373,9387) /830,828,828,827,827,827,826,825,
     *825,824,823,823,822,822,821/

      DATA (Add07( 5,i,1),i=5,45) /0,0,8513,8522,8529,8543,8556,8574,
     *8592,8609,8629,8649,8669,8709,8748,8783,8828,8873,8913,8958,
     *9003,9048,9098,9148,9188,9228,9268,9293,9318,9343,9358,9373,
     *0,0,0,0,0,0,0,0,0/

      DATA (Add07( 5,i,2),i=5,45) /0,0,7,9,7,8,8,8,9,11,11,11,11,12,
     *16,16,16,21,21,21,21,21,21,31,31,31,46,46,46,56,56,56,0,0,0,
     *0,0,0,0,0,0/

      DATA (Add07( 5,i,3),i=5,45) /0,0,15,15,20,20,25,25,25,30,30,
     *30,50,50,50,60,60,60,65,65,65,70,70,70,70,70,70,70,70,70,70,
     *70,0,0,0,0,0,0,0,0,0/

!     L?rche..................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=9388,9400) /758,710,670,630,590,550,510,490,
     *460,430,400,400,400/
      DATA (Md07(i),i=9401,9414) /634,549,475,410,400,400,400,400,
     *400,400,400,400,400,400/
      DATA (Md07(i),i=9415,9432) /522,831,739,678,609,544,487,449,
     *401,400,400,400,400,400,400,400,400,400/
      DATA (Md07(i),i=9433,9452) /910,860,766,695,641,586,554,523,
     *489,456,423,400,400,400,400,400,400,400,400,400/
      DATA (Md07(i),i=9453,9472) /910,827,782,736,689,654,623,591,
     *560,529,506,482,465,447,434,420,405,400,400,400/
      DATA (Md07(i),i=9473,9494) /976,940,893,844,804,759,727,695,
     *661,631,608,585,561,542,523,508,492,481,468,454,441,427/
      DATA (Md07(i),i=9495,9517) /967,925,894,850,810,780,747,722,
     *690,665,640,620,600,584,567,550,537,522,508,496,483,473,463/
      DATA (Md07(i),i=9518,9540) /937,887,858,825,793,767,739,713,
     *687,665,648,630,611,596,581,566,553,539,528,517,505,493,483/
      DATA (Md07(i),i=9541,9564) /892,859,828,802,778,755,731,707,
     *689,669,649,633,619,607,595,580,568,555,544,532,523,513,501,
     *491/
      DATA (Md07(i),i=9565,9590) /921,885,860,838,814,788,766,745,
     *724,706,688,673,657,641,632,616,602,590,579,568,557,546,537,
     *527,518,510/
      DATA (Md07(i),i=9591,9617) /911,889,865,840,818,795,776,756,
     *740,723,706,691,676,663,649,636,624,612,602,591,580,571,561,
     *552,545,537,529/
      DATA (Md07(i),i=9618,9646) /911,886,862,842,821,800,782,766,
     *751,735,721,706,692,679,666,653,642,632,621,611,601,592,583,
     *574,565,557,549,541,533/
      DATA (Md07(i),i=9647,9676) /903,880,860,842,822,806,791,776,
     *761,746,732,718,704,692,680,668,658,648,638,629,620,611,602,
     *593,585,577,569,562,556,550/
      DATA (Md07(i),i=9677,9709) /911,895,878,861,843,827,812,797,
     *783,768,755,742,728,716,704,693,682,672,662,653,645,636,628,
     *619,611,603,595,590,579,569,564,560,555/
      DATA (Md07(i),i=9710,9744) /919,905,892,876,860,845,830,816,
     *802,789,776,763,751,740,728,717,706,696,686,677,668,660,651,
     *643,635,625,620,614,608,602,596,590,584,578,576/
      DATA (Md07(i),i=9745,9779) /903,889,874,860,846,832,819,807,
     *795,783,772,760,750,739,728,719,709,700,691,682,674,668,662,
     *655,648,641,634,627,619,612,610,603,600,592,589/
      DATA (Md07(i),i=9780,9816) /914,900,886,873,859,846,835,823,
     *812,801,791,780,770,759,749,740,730,721,712,703,694,686,678,
     *670,667,659,656,648,645,637,633,625,621,613,609,601,597/
      DATA (Md07(i),i=9817,9853) /897,885,871,859,848,838,827,817,
     *807,797,787,777,768,759,750,740,728,719,710,702,698,689,685,
     *681,676,667,663,658,653,645,640,635,630,625,620,615,610/
      DATA (Md07(i),i=9854,9894) /907,895,882,872,861,851,841,831,
     *822,813,803,794,785,777,768,760,750,740,730,720,715,710,705,
     *700,694,689,683,677,672,666,660,655,652,646,640,635,632,626,
     *623,620,616/
      DATA (Md07(i),i=9895,9937) /915,904,892,882,873,864,854,844,
     *835,826,818,809,801,792,787,775,769,758,752,742,735,729,723,
     *717,710,704,701,698,691,685,679,672,669,663,659,655,652,645,
     *642,638,634,630,628/
      DATA (Md07(i),i=9938,9980) /900,891,882,874,866,857,848,839,
     *831,823,815,808,801,794,786,779,772,764,757,750,743,735,728,
     *721,717,713,709,702,698,691,687,683,678,674,670,665,661,656,
     *654,650,647,643,640/
      DATA (Md07(i),i=9981,10025) /909,900,891,883,876,868,860,853,
     *845,838,830,822,814,806,797,789,785,777,773,765,760,752,747,
     *740,735,730,725,720,715,711,706,701,696,691,688,683,680,675,
     *672,669,666,663,660,657,655/
      DATA (Md07(i),i=10026,10072) /917,908,900,892,885,878,871,864,
     *857,848,839,830,825,816,811,802,797,792,787,778,773,767,762,
     *756,751,745,740,734,731,726,723,717,714,708,705,701,698,694,
     *691,687,684,680,678,674,672,668,666/
      DATA (Md07(i),i=10073,10121) /909,902,894,887,881,871,866,856,
     *851,841,836,830,824,814,808,802,796,790,784,778,775,769,765,
     *759,755,749,746,742,738,734,730,726,722,718,714,710,707,703,
     *699,695,693,690,688,684,681,678,675,673,670/
      DATA (Md07(i),i=10122,10172) /918,911,903,895,886,876,870,864,
     *858,851,845,838,835,828,822,815,808,801,797,791,787,783,778,
     *774,770,765,761,756,752,747,745,740,738,733,730,726,723,718,
     *714,711,708,705,702,698,695,691,688,685,682,679,676/
      DATA (Md07(i),i=10173,10223) /910,903,896,889,882,875,868,861,
     *857,850,845,838,834,826,822,814,809,805,800,795,790,785,782,
     *778,775,770,766,762,758,753,750,747,744,740,737,732,729,725,
     *722,718,715,712,708,705,702,698,696,694,691,689,688/
      DATA (Md07(i),i=10224,10274) /901,893,890,882,877,869,865,860,
     *855,850,844,839,834,828,823,817,812,807,801,798,794,789,785,
     *782,778,774,771,767,763,760,757,753,750,746,742,738,735,731,
     *728,725,722,719,716,713,711,708,706,703,701,700,699/
      DATA (Md07(i),i=10275,10327) /909,901,896,891,886,881,875,869,
     *863,858,854,848,845,839,833,827,823,817,813,809,805,801,797,
     *793,789,785,782,778,775,771,768,764,761,757,755,750,748,745,
     *742,739,736,733,730,727,724,722,720,717,715,713,711,708,706/
      DATA (Md07(i),i=10328,10384) /913,908,902,896,894,888,884,878,
     *874,868,863,859,855,848,844,838,833,829,824,820,815,811,806,
     *802,799,796,793,789,786,783,780,777,774,770,767,764,761,757,
     *754,751,749,746,743,740,738,736,733,731,729,726,723,721,719,
     *716,714,712,710/
      DATA (Md07(i),i=10385,10443) /916,910,908,905,901,897,893,888,
     *884,879,874,869,864,859,854,849,844,839,834,829,826,821,818,
     *813,809,806,802,799,797,793,791,787,785,781,779,775,773,769,
     *767,763,761,758,755,753,750,748,746,743,741,738,736,733,731,
     *728,726,724,723,721,719/
      DATA (Md07(i),i=10444,10502) /916,912,908,903,898,893,890,885,
     *882,876,873,867,864,858,855,849,845,840,836,831,827,823,820,
     *816,814,810,807,803,801,798,795,792,789,787,784,781,778,775,
     *772,770,768,765,763,760,758,755,753,750,748,745,743,741,739,
     *737,735,733,731,729,727/
      DATA (Md07(i),i=10503,10561) /924,919,914,909,903,898,894,888,
     *885,881,877,873,869,865,861,856,854,849,847,842,839,835,832,
     *828,825,820,817,814,811,808,805,802,799,796,794,791,789,786,
     *783,780,778,776,773,771,769,766,764,762,759,757,755,752,750,
     *748,746,744,742,740,738/
      DATA (Md07(i),i=10562,10620) /920,914,910,906,902,898,894,889,
     *885,880,877,873,870,865,862,857,854,850,847,844,840,837,833,
     *830,826,823,821,817,815,812,810,807,805,802,799,797,794,792,
     *789,787,785,782,780,778,776,773,771,769,767,764,762,760,758,
     *755,754,752,750,748,746/
      DATA (Md07(i),i=10621,10685) /925,921,917,912,910,905,902,897,
     *894,889,885,880,876,873,869,865,862,858,856,852,849,845,843,
     *839,836,833,830,828,825,822,819,816,814,811,809,807,805,802,
     *800,797,795,793,791,789,787,785,782,780,778,776,774,772,770,
     *767,766,764,762,760,758,757,755,754,752,750,749/
      DATA (Md07(i),i=10686,10748) /925,920,917,913,910,904,900,896,
     *892,888,884,880,876,872,869,866,863,860,857,854,851,848,845,
     *842,839,836,833,830,828,826,823,821,819,816,814,812,809,807,
     *804,802,800,798,796,794,792,790,788,786,784,782,780,778,777,
     *775,773,771,770,768,767,766,764,762,761/
      DATA (Md07(i),i=10749,10811) /933,929,925,921,917,912,908,904,
     *899,895,892,887,884,881,878,874,871,867,865,862,859,856,853,
     *850,847,845,842,840,837,835,832,829,827,825,823,820,818,816,
     *814,811,810,808,806,804,802,800,798,796,794,792,791,789,787,
     *785,784,782,781,779,778,776,775,773,772/
      DATA (Md07(i),i=10812,10872) /931,926,923,919,914,909,906,902,
     *899,895,891,888,885,882,879,875,872,870,867,864,861,858,856,
     *853,851,848,845,842,840,838,836,833,831,829,827,825,822,820,
     *819,816,815,812,811,809,807,805,804,802,800,798,797,795,794,
     *792,791,789,788,786,785,783,784/
      DATA (Md07(i),i=10873,10933) /939,933,930,925,921,917,913,909,
     *906,902,899,895,892,889,886,883,880,877,875,872,870,866,864,
     *861,858,855,853,850,848,846,844,841,840,837,835,833,831,829,
     *827,825,823,821,819,817,815,814,812,811,809,808,806,805,803,
     *802,800,799,798,797,795,794,792/
      DATA (Md07(i),i=10934,10988) /919,916,912,909,906,903,900,897,
     *895,891,888,885,883,879,876,874,871,868,866,863,860,858,856,
     *854,852,850,848,846,843,841,839,837,835,833,831,829,827,826,
     *824,822,821,819,818,816,815,813,812,810,809,808,806,805,804,
     *802,801/

      DATA (Add07( 6,i,1),i=5,45) /0,0,9388,9401,9415,9433,9453,9473,
     *9495,9518,9541,9565,9591,9618,9647,9677,9710,9745,9780,9817,
     *9854,9895,9938,9981,10026,10073,10122,10173,10224,10275,10328,
     *10385,10444,10503,10562,10621,10686,10749,10812,10873,10934/

      DATA (Add07( 6,i,2),i=5,45) /0,0,7,9,7,8,9,8,9,11,13,13,14,15,
     *16,16,16,18,18,20,20,20,22,22,22,24,24,26,28,28,28,28,30,30,
     *32,32,34,34,36,36,42/

      DATA (Add07( 6,i,3),i=5,45) /0,0,19,22,24,27,28,29,31,33,36,
     *38,40,43,45,48,50,52,54,56,60,62,64,66,68,72,74,76,78,80,84,
     *86,88,88,90,96,96,96,96,96,96/

!     Jap.L?rche..............................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=10989,10995) /854,810,770,730,690,650,630/
      DATA (Md07(i),i=10996,10999) /450,445,410,400/
      DATA (Md07(i),i=11000,11007) /667,619,558,523,487,465,435,404/
      DATA (Md07(i),i=11008,11019) /607,812,743,698,655,622,586,569,
     *551,528,503,477/
      DATA (Md07(i),i=11020,11034) /594,818,786,747,712,687,659,641,
     *623,602,579,556,563,556,546/
      DATA (Md07(i),i=11035,11049) /479,855,830,784,752,732,707,690,
     *673,652,631,608,612,610,604/
      DATA (Md07(i),i=11050,11063) /876,851,810,780,763,742,728,709,
     *688,668,652,652,647,644/
      DATA (Md07(i),i=11064,11079) /892,863,829,815,795,777,762,744,
     *722,700,684,682,679,673,665,654/
      DATA (Md07(i),i=11080,11097) /902,876,847,838,820,800,787,770,
     *749,732,713,709,703,700,690,681,683,684/
      DATA (Md07(i),i=11098,11119) /896,867,860,839,822,806,790,775,
     *757,737,734,727,721,712,702,704,703,702,699,694,697,697/
      DATA (Md07(i),i=11120,11141) /901,885,875,857,843,829,810,795,
     *777,757,755,749,741,734,723,724,721,720,715,710,712,713/
      DATA (Md07(i),i=11142,11163) /912,898,888,875,860,845,832,816,
     *798,777,772,768,761,752,742,741,738,735,730,724,726,726/
      DATA (Md07(i),i=11164,11185) /913,908,906,891,880,863,849,832,
     *817,796,793,786,777,769,758,757,754,751,745,740,740,739/
      DATA (Md07(i),i=11186,11220) /921,922,919,910,898,883,867,849,
     *833,815,809,801,793,784,773,772,770,765,759,753,752,750,747,
     *743,740,738,736,734,732,703,727,725,723,720,717/
      DATA (Md07(i),i=11221,11252) /926,916,898,884,865,849,833,825,
     *818,808,800,788,786,783,777,772,764,762,760,757,754,749,748,
     *745,743,740,737,735,734,731,729,726/
      DATA (Md07(i),i=11253,11284) /945,932,915,899,882,864,848,842,
     *832,823,813,802,798,794,789,782,775,773,770,767,763,758,757,
     *754,751,748,745,744,742,740,737,735/
      DATA (Md07(i),i=11285,11316) /967,949,930,911,895,879,863,855,
     *846,836,826,815,811,805,800,793,786,783,780,776,772,766,765,
     *762,759,756,752,752,750,748,745,742/
      DATA (Md07(i),i=11317,11357) /965,946,927,909,892,874,867,858,
     *848,837,827,822,816,810,803,795,793,789,784,780,775,773,770,
     *767,764,760,759,757,755,753,750,751,750,749,748,747,747,747,
     *746,746,745/
      DATA (Md07(i),i=11358,11398) /982,960,938,919,903,887,880,871,
     *860,850,839,833,827,821,813,806,802,797,793,788,782,780,777,
     *774,771,767,766,764,762,760,757,757,757,756,755,754,755,754,
     *754,753,752/
      DATA (Md07(i),i=11399,11441) /970,951,932,914,897,891,882,872,
     *861,849,843,837,830,823,815,810,806,800,795,789,787,784,780,
     *777,773,772,770,768,766,764,764,764,763,762,761,761,761,761,
     *760,759,760,760,760/
      DATA (Md07(i),i=11442,11480) /907,901,892,882,871,860,853,846,
     *840,831,824,819,814,809,803,797,794,791,787,783,780,778,777,
     *775,773,770,770,770,770,769,768,768,768,768,767,766,767,767,
     *767/
      DATA (Md07(i),i=11481,11519) /918,910,901,891,880,869,862,855,
     *848,840,832,827,822,816,810,804,801,798,794,790,786,785,783,
     *781,779,776,776,776,776,775,774,775,774,774,773,773,774,774,
     *774/
      DATA (Md07(i),i=11520,11553) /879,872,864,856,848,840,835,829,
     *824,818,812,808,805,801,796,792,791,789,787,785,782,782,782,
     *782,781,780,781,781,780,780,779,780,780,780/
      DATA (Md07(i),i=11554,11587) /888,880,872,864,856,847,842,837,
     *831,825,819,816,811,807,803,798,797,795,793,791,788,788,788,
     *788,787,786,787,787,786,786,785,786,786,787/
      DATA (Md07(i),i=11588,11621) /898,889,880,872,864,855,850,844,
     *839,833,827,823,819,814,809,805,803,801,799,797,794,794,794,
     *794,793,792,792,792,792,791,791,791,792,792/
      DATA (Md07(i),i=11622,11655) /906,897,888,880,871,862,857,852,
     *846,840,834,830,826,821,816,811,809,807,805,803,800,800,800,
     *799,798,797,798,797,797,797,796,797,797,798/
      DATA (Md07(i),i=11656,11679) /842,838,833,828,823,818,816,814,
     *812,809,807,806,806,805,804,803,803,803,803,802,802,802,803,
     *803/
      DATA (Md07(i),i=11680,11703) /850,845,840,835,830,825,823,821,
     *818,816,813,812,812,811,810,809,809,809,808,808,807,808,808,
     *808/
      DATA (Md07(i),i=11704,11714) /817,816,815,815,814,814,813,812,
     *813,813,813/
      DATA (Md07(i),i=11715,11725) /823,822,820,820,820,819,818,817,
     *818,818,818/

      DATA (Add07( 7,i,1),i=5,45) /0,0,10989,10996,11000,11008,11020,
     *11035,11050,11064,11080,11098,11120,11142,11164,11186,11221,
     *11253,11285,11317,11358,11399,11442,11481,11520,11554,11588,
     *11622,11656,11680,11704,11715,0,0,0,0,0,0,0,0,0/

      DATA (Add07( 7,i,2),i=5,45) /0,0,7,10,9,7,7,7,8,8,8,9,9,9,9,
     *9,12,12,12,13,13,14,18,18,23,23,23,23,33,33,46,46,0,0,0,0,0,
     *0,0,0,0/

      DATA (Add07( 7,i,3),i=5,45) /0,0,13,13,16,18,21,21,21,23,25,
     *30,30,30,30,43,43,43,43,53,53,56,56,56,56,56,56,56,56,56,56,
     *56,0,0,0,0,0,0,0,0,0/

!     Buche...................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=11726,11732) /455,460,430,400,400,400,400/
      DATA (Md07(i),i=11733,11742) /669,596,588,569,536,526,552,511,
     *485,472/
      DATA (Md07(i),i=11743,11754) /726,681,675,637,641,632,651,622,
     *602,604,583,579/
      DATA (Md07(i),i=11755,11769) /786,767,759,743,725,716,701,709,
     *704,680,674,660,651,644,640/
      DATA (Md07(i),i=11770,11787) /776,790,785,766,754,748,750,742,
     *738,716,709,703,692,691,684,678,673,669/
      DATA (Md07(i),i=11788,11808) /515,830,794,785,772,771,771,766,
     *760,749,739,732,721,716,713,710,703,701,698,696,691/
      DATA (Md07(i),i=11809,11835) /498,837,807,798,790,784,786,790,
     *777,764,755,747,745,739,733,728,727,722,721,719,715,713,711,
     *708,706,705,703/
      DATA (Md07(i),i=11836,11866) /836,819,808,803,797,790,784,787,
     *783,772,763,759,756,749,745,745,741,737,736,734,732,730,727,
     *726,725,723,721,718,717,716,714/
      DATA (Md07(i),i=11867,11899) /855,834,829,807,806,806,808,802,
     *791,785,775,770,769,764,761,756,753,753,749,747,745,745,742,
     *740,740,737,736,734,733,732,730,729,727/
      DATA (Md07(i),i=11900,11934) /857,839,832,816,812,813,815,809,
     *802,791,789,783,779,776,772,768,766,765,762,760,758,757,755,
     *753,752,751,750,748,747,745,744,743,741,740,738/
      DATA (Md07(i),i=11935,11975) /959,856,843,838,822,821,819,819,
     *818,807,801,798,793,788,783,781,780,776,774,773,771,769,767,
     *767,764,764,762,761,759,759,757,756,755,754,753,751,750,749,
     *747,746,745/
      DATA (Md07(i),i=11976,12021) /867,853,840,828,818,823,823,821,
     *815,808,803,799,796,792,790,788,785,784,782,780,779,777,776,
     *775,774,772,771,770,769,768,767,766,764,763,762,761,760,759,
     *758,757,756,755,753,752,751,750/
      DATA (Md07(i),i=12022,12069) /857,843,831,822,829,826,829,823,
     *818,813,807,805,800,799,795,794,792,789,789,787,786,785,784,
     *782,781,781,779,778,777,776,776,774,774,773,772,771,770,769,
     *768,767,766,765,764,763,762,761,759,759/
      DATA (Md07(i),i=12070,12165) /858,849,836,828,831,833,832,828,
     *823,817,813,810,806,804,802,800,799,798,796,795,794,793,792,
     *791,790,789,788,787,786,785,785,784,783,782,781,780,779,778,
     *778,777,776,775,774,773,772,771,770,769,768,767,766,765,764,
     *763,762,761,760,759,758,757,756,754,753,752,751,750,748,747,
     *746,745,743,742,741,740,738,737,736,734,733,732,730,729,727,
     *726,725,723,722,720,719,717,716,714,713,711,710,708/
      DATA (Md07(i),i=12166,12260) /863,844,831,834,837,836,832,828,
     *824,819,817,814,811,809,807,806,805,803,803,802,801,799,799,
     *797,796,796,795,794,794,793,792,792,791,790,789,788,788,787,
     *786,785,785,784,783,782,781,780,779,779,778,777,776,775,774,
     *773,772,771,770,769,768,767,766,765,764,763,762,761,760,759,
     *758,756,755,754,753,752,751,749,748,747,746,744,743,742,741,
     *739,738,737,735,734,733,731,730,728,727,726,724/
      DATA (Md07(i),i=12261,12354) /846,839,836,839,842,837,832,830,
     *826,823,819,818,815,814,812,811,809,809,808,807,806,805,805,
     *804,803,802,802,801,801,800,799,798,798,797,796,796,795,795,
     *794,793,792,792,791,790,790,789,788,787,786,786,785,784,783,
     *782,781,781,780,779,778,777,776,775,774,773,772,771,770,769,
     *768,767,766,765,764,763,762,761,760,759,758,756,755,754,753,
     *752,751,750,748,747,746,745,743,742,741,740/
      DATA (Md07(i),i=12355,12448) /850,842,846,845,847,840,837,832,
     *830,827,825,823,821,819,818,817,816,815,814,813,813,812,811,
     *810,810,809,809,808,808,807,806,806,805,805,804,804,803,802,
     *802,801,801,800,800,799,798,798,797,796,796,795,794,793,793,
     *792,791,790,790,789,788,787,786,786,785,784,783,782,781,780,
     *780,779,778,777,776,775,774,773,772,771,770,769,768,767,766,
     *765,764,763,762,760,759,758,757,756,755,754/
      DATA (Md07(i),i=12449,12541) /844,845,846,850,844,840,837,835,
     *832,829,827,826,824,823,822,821,820,820,819,818,818,817,817,
     *816,816,815,815,814,814,813,813,813,812,811,811,810,810,810,
     *809,808,808,807,807,806,806,805,804,804,803,803,802,801,801,
     *800,800,799,798,797,797,796,795,795,794,793,792,792,791,790,
     *789,788,788,787,786,785,784,783,782,782,781,780,779,778,777,
     *776,775,774,773,772,771,770,769,768,767/
      DATA (Md07(i),i=12542,12633) /851,853,853,847,845,841,840,838,
     *834,833,831,830,829,828,827,826,825,825,824,824,823,823,823,
     *822,822,821,821,821,820,820,819,819,818,818,818,817,817,816,
     *816,816,815,815,814,814,813,813,812,812,811,811,810,810,809,
     *809,808,807,807,806,806,805,804,804,803,802,802,801,800,800,
     *799,798,797,797,796,795,794,794,793,792,791,790,790,789,788,
     *787,786,785,784,784,783,782,781,780/
      DATA (Md07(i),i=12634,12725) /852,853,856,852,849,847,844,841,
     *839,837,835,834,833,833,832,831,831,830,830,829,829,828,828,
     *828,827,827,827,826,826,826,825,825,825,824,824,824,823,823,
     *823,822,822,822,821,821,820,820,820,819,819,818,818,817,817,
     *816,816,815,815,814,814,813,813,812,812,811,811,810,809,809,
     *808,808,807,806,806,805,804,804,803,802,801,801,800,799,799,
     *798,797,796,796,795,794,793,792,792/
      DATA (Md07(i),i=12726,12816) /860,859,855,851,850,847,845,843,
     *842,840,839,838,837,836,836,835,835,835,834,834,834,833,833,
     *832,832,832,832,832,831,831,831,831,830,830,830,830,829,829,
     *829,829,828,828,828,827,827,827,826,826,826,825,825,825,824,
     *824,823,823,822,822,821,821,820,820,820,819,819,818,817,817,
     *816,816,815,815,814,814,813,812,812,811,811,810,809,809,808,
     *807,807,806,805,805,804,803,802/
      DATA (Md07(i),i=12817,12906) /862,859,855,853,851,849,847,846,
     *844,843,842,842,841,840,840,839,839,839,838,839,838,838,838,
     *838,837,837,837,837,837,836,836,836,836,836,836,835,835,835,
     *835,835,834,834,834,833,833,833,833,832,832,832,831,831,831,
     *830,830,830,829,829,829,828,828,827,827,827,826,826,825,825,
     *824,824,823,823,822,822,821,821,820,820,819,819,818,817,817,
     *816,816,815,815,814,813,813/
      DATA (Md07(i),i=12907,12995) /861,858,857,854,852,851,850,848,
     *847,846,845,845,844,844,844,844,843,843,843,843,843,842,842,
     *842,842,842,842,842,842,842,842,841,841,841,841,841,841,841,
     *840,840,840,840,840,840,839,839,839,839,838,838,838,838,837,
     *837,837,837,836,836,836,835,835,835,834,834,833,833,833,832,
     *832,832,831,831,830,830,829,829,828,828,827,827,827,826,825,
     *825,824,824,823,823,822/
      DATA (Md07(i),i=12996,13084) /864,861,859,857,856,854,853,851,
     *850,850,850,849,848,848,848,848,848,847,847,847,847,847,847,
     *847,847,847,847,847,847,847,847,847,847,847,846,846,846,846,
     *846,846,846,846,846,845,845,845,845,845,845,845,844,844,844,
     *844,844,843,843,843,843,842,842,842,841,841,841,841,840,840,
     *840,839,839,839,838,838,837,837,837,836,836,835,835,835,834,
     *834,833,833,832,832,831/
      DATA (Md07(i),i=13085,13172) /864,863,861,859,858,856,856,855,
     *854,853,853,853,852,852,852,852,852,852,851,852,852,851,852,
     *852,852,852,852,852,852,852,852,852,852,852,852,852,852,852,
     *851,851,851,851,851,851,851,851,851,851,851,850,850,850,850,
     *850,850,850,849,849,849,849,849,848,848,848,848,844,847,847,
     *847,846,846,846,846,845,845,845,844,844,844,843,843,842,842,
     *842,841,841,841,840/
      DATA (Md07(i),i=13173,13259) /865,864,863,861,860,859,858,857,
     *857,857,856,856,856,856,856,856,856,856,856,856,856,856,856,
     *856,856,856,856,856,856,857,857,857,857,857,857,857,857,857,
     *857,857,857,857,857,857,857,857,856,856,856,856,856,856,856,
     *856,856,856,856,855,855,855,855,855,855,854,854,854,854,854,
     *853,853,853,853,852,852,852,852,851,851,851,851,850,850,850,
     *849,849,849,848/
      DATA (Md07(i),i=13260,13345) /866,865,864,863,862,861,861,860,
     *860,860,860,859,859,860,860,860,860,860,860,860,860,860,860,
     *860,861,861,861,861,861,861,861,861,861,862,862,862,862,862,
     *862,862,862,862,862,862,862,862,862,862,862,862,862,862,862,
     *862,862,862,862,861,861,861,861,861,861,861,861,861,860,860,
     *860,860,860,859,859,859,859,859,858,858,858,858,858,857,857,
     *857,857,856/
      DATA (Md07(i),i=13346,13430) /868,867,867,866,865,864,864,864,
     *863,863,863,863,863,863,863,863,864,864,864,864,864,864,865,
     *865,865,865,865,865,866,866,866,866,866,866,866,866,867,867,
     *867,867,867,867,867,867,867,867,867,867,867,867,867,867,867,
     *867,867,867,867,867,867,867,867,867,867,867,867,867,867,866,
     *866,866,866,866,866,866,866,865,865,865,865,865,865,864,864,
     *864,864/
      DATA (Md07(i),i=13431,13513) /870,869,868,867,867,866,867,866,
     *867,867,867,867,867,867,867,868,868,868,868,868,869,869,869,
     *869,869,870,870,870,870,870,871,871,871,871,871,871,872,872,
     *872,872,872,872,872,872,872,873,873,873,873,873,873,873,873,
     *873,873,873,873,873,873,873,873,873,873,873,873,873,873,873,
     *873,872,872,872,872,872,872,872,872,872,871,871,871,871,871/
      DATA (Md07(i),i=13514,13593) /870,870,870,870,870,870,870,870,
     *870,871,871,871,871,871,872,872,872,872,873,873,873,873,874,
     *874,874,874,875,875,875,875,875,876,876,876,876,876,872,877,
     *877,877,877,877,877,878,878,878,878,878,878,878,878,878,878,
     *878,879,879,879,879,879,879,879,879,879,879,879,879,879,879,
     *879,878,878,878,878,878,878,878,878,878,878,878/
      DATA (Md07(i),i=13594,13670) /873,873,873,873,873,874,874,874,
     *874,875,875,875,875,876,876,876,877,877,877,878,878,878,878,
     *879,879,879,879,880,880,880,881,881,881,881,881,881,882,882,
     *882,882,882,883,883,883,883,883,883,883,884,884,884,884,884,
     *884,884,884,884,884,884,884,884,884,884,884,885,885,885,885,
     *885,884,884,884,884,884,884,884,884/
      DATA (Md07(i),i=13671,13744) /876,877,877,877,877,878,878,878,
     *879,879,879,880,880,880,881,881,881,882,882,882,883,883,883,
     *884,884,884,884,885,885,885,885,886,886,886,886,887,887,887,
     *887,887,888,888,888,888,888,888,889,889,889,889,889,889,889,
     *890,890,890,890,890,890,890,890,890,890,890,890,890,890,890,
     *890,890,890,890,890,890/
      DATA (Md07(i),i=13745,13808) /883,883,884,884,884,885,885,885,
     *886,886,886,887,887,887,888,888,888,889,889,889,890,890,890,
     *890,891,891,891,892,892,892,892,892,893,893,893,893,894,894,
     *894,894,894,894,895,895,895,895,895,895,895,895,896,896,896,
     *896,896,896,896,896,896,896,896,896,896,896/
      DATA (Md07(i),i=13809,13859) /891,892,892,892,893,893,893,894,
     *894,894,895,895,895,896,896,896,896,897,897,897,897,898,898,
     *898,898,899,899,899,899,899,900,900,900,900,900,900,901,901,
     *901,901,899,901,901,902,902,902,902,902,902,902,902/
      DATA (Md07(i),i=13860,13895) /900,901,901,901,902,902,902,902,
     *903,903,903,903,904,904,904,904,905,905,905,905,905,906,906,
     *906,906,906,907,907,907,907,907,907,907,908,908,908/
      DATA (Md07(i),i=13896,13916) /909,909,910,910,910,910,910,911,
     *911,911,911,912,912,912,912,912,912,913,913,913,913/

      DATA (Add07( 8,i,1),i=5,45) /0,0,11726,11733,11743,11755,11770,
     *11788,11809,11836,11867,11900,11935,11976,12022,12070,12166,
     *12261,12355,12449,12542,12634,12726,12817,12907,12996,13085,
     *13173,13260,13346,13431,13514,13594,13671,13745,13809,13860,
     *13896,0,0,0/

      DATA (Add07( 8,i,2),i=5,45) /0,0,10,9,9,7,7,7,7,8,8,8,7,8,9,
     *9,10,11,11,12,13,13,14,15,16,16,17,18,19,20,22,25,28,31,41,54,
     *69,84,0,0,0/

      DATA (Add07( 8,i,3),i=5,45) /0,0,16,18,20,21,24,27,33,38,40,
     *42,47,53,56,104,104,104,104,104,104,104,104,104,104,104,104,
     *104,104,104,104,104,104,104,104,104,104,104,0,0,0/

!     Eiche.............................................................

      DATA (Md07(i),i=13917,13917) /437/
      DATA (Md07(i),i=13918,13926) /416,522,604,632,663,673,667,665,
     *653/
      DATA (Md07(i),i=13927,13945) /557,501,502,507,523,557,614,669,
     *703,718,720,712,718,714,721,711,704,712,713/
      DATA (Md07(i),i=13946,13970) /811,701,644,612,618,620,633,677,
     *717,738,746,744,743,743,744,744,745,746,752,752,748,748,752,
     *755,755/
      DATA (Md07(i),i=13971,14000) /868,771,724,691,687,684,691,714,
     *744,759,763,767,770,766,768,770,771,776,776,776,776,772,777,
     *778,776,776,774,771,771,770/
      DATA (Md07(i),i=14001,14033) /887,827,762,732,729,731,734,740,
     *764,774,775,781,786,785,784,786,787,792,792,791,787,786,789,
     *791,789,786,785,785,785,786,789,793,795/
      DATA (Md07(i),i=14034,14073) /916,846,802,766,750,751,765,767,
     *777,784,789,792,793,794,794,794,799,800,800,802,798,796,798,
     *800,798,795,794,794,794,796,797,799,802,802,801,799,797,794,
     *793,791/
      DATA (Md07(i),i=14074,14129) /928,876,832,798,772,767,777,778,
     *788,791,794,799,802,801,802,802,805,808,808,807,805,803,804,
     *805,803,801,799,799,799,801,802,805,807,807,806,804,803,801,
     *800,798,798,796,794,793,792,790,788,787,785,784,782,781,779,
     *778,776,774/
      DATA (Md07(i),i=14130,14189) /953,895,857,815,792,777,787,794,
     *794,801,801,804,806,806,808,809,811,812,812,811,811,811,809,
     *809,807,805,805,804,804,806,807,810,812,811,810,808,808,806,
     *805,804,802,800,798,797,797,796,795,793,792,790,788,786,784,
     *783,782,781,780,778,775,775/
      DATA (Md07(i),i=14190,14254) /965,912,876,827,806,787,795,805,
     *805,805,808,812,811,810,813,814,816,817,817,814,814,813,813,
     *812,811,809,809,807,809,811,812,814,815,815,815,813,813,811,
     *810,809,806,804,802,802,802,801,801,799,799,797,794,793,790,
     *789,788,787,785,784,783,781,780,779,777,776,775/
      DATA (Md07(i),i=14255,14323) /925,884,840,818,799,800,810,808,
     *808,810,816,816,816,819,818,820,821,819,819,818,817,817,815,
     *814,812,813,811,814,815,815,818,819,819,819,818,817,816,814,
     *813,811,809,807,806,807,807,806,805,804,803,800,799,796,795,
     *794,793,791,789,787,785,784,783,781,780,779,778,777,775,774,
     *773/
      DATA (Md07(i),i=14324,14396) /898,848,823,809,808,813,811,814,
     *814,818,820,821,822,822,824,825,823,822,822,820,821,819,818,
     *816,816,815,817,819,819,821,822,822,821,820,820,819,818,817,
     *815,813,811,810,811,811,810,809,808,807,806,804,802,800,799,
     *798,797,795,793,791,791,789,788,787,785,784,783,782,780,779,
     *779,780,780,781,782/
      DATA (Md07(i),i=14397,14471) /856,831,819,818,815,813,819,818,
     *821,823,822,825,825,827,827,827,824,824,824,823,823,821,819,
     *819,818,821,822,822,824,824,824,824,823,823,823,821,820,818,
     *816,815,814,815,815,814,813,812,812,810,809,807,806,805,803,
     *801,800,798,797,796,795,794,793,792,790,789,788,787,786,785,
     *785,798,785,786,786,786,785/
      DATA (Md07(i),i=14472,14549) /838,827,824,819,819,821,821,823,
     *826,826,829,828,829,830,828,827,827,826,826,825,824,823,822,
     *821,823,825,825,826,826,825,825,825,826,825,824,823,821,820,
     *819,819,819,818,818,817,816,816,814,813,811,810,808,807,806,
     *805,803,802,802,801,800,798,797,796,795,794,793,792,791,791,
     *790,789,788,787,786,785,784,783,782,781/
      DATA (Md07(i),i=14550,14629) /830,821,821,824,825,826,830,830,
     *831,830,830,831,830,829,829,828,829,829,828,827,827,826,825,
     *827,829,829,828,826,825,826,827,826,826,825,824,823,822,822,
     *821,820,821,820,820,819,818,816,815,814,812,811,810,809,808,
     *807,806,806,805,804,802,801,801,800,799,798,797,796,795,794,
     *793,792,791,791,790,789,788,787,786,787,787,787/
      DATA (Md07(i),i=14630,14711) /826,823,825,828,830,832,832,833,
     *833,833,834,833,832,832,832,832,832,831,831,831,830,828,830,
     *831,831,830,829,827,828,828,828,828,827,826,825,824,824,822,
     *823,823,823,823,821,820,819,818,816,815,814,813,813,812,812,
     *811,810,810,809,808,807,806,805,804,804,803,802,801,800,799,
     *799,798,797,796,796,795,794,794,794,795,796,797,798,799/
      DATA (Md07(i),i=14712,14795) /829,830,833,835,835,836,835,836,
     *837,836,835,834,834,834,835,835,834,834,833,832,833,834,834,
     *833,832,830,830,830,830,830,829,828,828,827,827,826,825,826,
     *826,826,824,824,822,821,820,820,819,818,818,817,817,816,816,
     *815,814,813,813,812,811,811,810,809,808,807,817,806,806,805,
     *804,804,803,803,802,801,802,802,803,805,806,807,808,809,809,
     *809/
      DATA (Md07(i),i=14796,14878) /833,837,837,837,837,838,837,839,
     *839,839,839,837,837,838,838,837,837,837,835,836,837,837,836,
     *834,833,833,833,833,832,831,831,831,830,830,830,829,829,829,
     *828,827,826,825,824,824,824,823,823,823,822,822,821,821,820,
     *820,819,807,818,817,817,816,815,815,814,813,813,812,812,811,
     *811,810,810,809,809,809,810,810,811,812,813,814,814,815,815/
      DATA (Md07(i),i=14879,14960) /840,840,838,838,838,838,839,840,
     *842,841,840,841,840,841,840,840,840,839,839,840,840,839,838,
     *837,837,837,836,835,835,835,835,834,834,834,833,833,833,833,
     *832,832,830,830,830,829,828,828,828,827,827,827,826,826,826,
     *826,825,825,824,823,822,822,821,820,820,820,819,819,819,818,
     *818,817,817,817,818,818,818,818,818,819,820,820,821,821/
      DATA (Md07(i),i=14961,15041) /841,840,839,838,840,840,841,842,
     *842,842,842,842,843,843,843,843,842,842,842,843,842,841,841,
     *840,840,840,839,839,839,838,838,838,838,838,837,837,837,837,
     *837,836,835,834,834,834,833,833,833,833,832,832,832,832,831,
     *831,831,830,830,829,828,827,827,827,826,826,826,825,825,825,
     *824,824,824,824,825,825,825,824,824,824,825,825,825/
      DATA (Md07(i),i=15042,15121) /840,840,839,840,841,841,841,842,
     *843,844,843,845,845,845,846,845,845,845,846,846,844,844,844,
     *844,844,843,843,844,843,842,842,842,842,842,842,842,842,841,
     *840,840,839,839,838,838,838,838,838,838,838,837,837,837,837,
     *837,836,836,835,834,834,833,832,832,832,831,831,831,831,830,
     *831,831,831,830,830,830,830,829,829,829,829,830/
      DATA (Md07(i),i=15122,15200) /842,840,841,840,841,840,842,844,
     *845,845,847,847,847,849,849,848,848,849,849,847,847,847,847,
     *847,847,847,847,847,846,846,846,846,846,846,846,846,845,845,
     *844,844,844,844,844,843,843,843,843,843,843,843,842,842,842,
     *841,841,840,840,839,838,838,838,838,838,838,837,837,837,837,
     *837,837,836,836,835,835,835,834,834,834,833/
      DATA (Md07(i),i=15201,15278) /840,842,840,840,840,843,846,847,
     *847,848,848,849,850,851,851,851,852,852,851,850,850,850,850,
     *850,850,850,850,850,850,850,849,849,849,849,849,849,849,849,
     *849,849,849,849,848,848,847,847,847,847,847,847,847,846,846,
     *846,845,845,844,843,843,843,843,843,843,842,842,841,841,841,
     *841,841,841,841,841,841,840,840,840,839/
      DATA (Md07(i),i=15279,15355) /841,840,839,840,844,847,848,850,
     *850,850,851,852,853,853,853,854,855,854,853,853,853,853,853,
     *853,853,853,853,852,852,853,853,853,853,853,853,853,852,852,
     *852,852,863,851,850,850,851,851,850,851,851,851,851,850,850,
     *849,849,849,849,849,848,848,848,848,848,847,847,847,847,847,
     *847,847,846,846,846,846,845,845,845/
      DATA (Md07(i),i=15356,15431) /840,838,840,844,848,850,851,852,
     *852,853,854,855,855,855,857,857,857,856,856,856,856,856,856,
     *857,857,857,856,856,856,856,857,857,857,856,856,856,856,856,
     *856,856,855,855,855,855,855,855,855,855,855,855,855,854,854,
     *854,853,853,854,853,853,854,853,853,853,852,852,851,851,852,
     *852,852,852,851,851,851,851,851/
      DATA (Md07(i),i=15432,15506) /838,840,845,849,851,853,852,854,
     *855,856,856,856,858,859,859,859,858,859,859,859,859,860,860,
     *859,859,859,860,860,860,860,860,860,860,860,860,860,860,859,
     *859,859,859,859,859,859,859,858,858,858,859,859,859,859,858,
     *857,857,858,858,858,858,858,858,857,857,856,857,857,857,857,
     *857,857,856,856,856,856,856/
      DATA (Md07(i),i=15507,15580) /841,845,850,853,855,855,856,857,
     *857,858,858,859,860,861,861,861,862,862,862,863,863,863,862,
     *862,862,862,862,863,863,863,864,864,864,864,864,863,863,863,
     *862,862,862,861,861,861,862,862,862,862,863,862,862,862,861,
     *862,862,862,862,862,862,862,862,862,861,861,861,861,861,862,
     *862,862,861,861,862,861/
      DATA (Md07(i),i=15581,15653) /846,851,854,856,856,858,858,858,
     *859,859,860,861,862,862,862,863,863,863,864,865,864,865,865,
     *864,865,865,865,866,866,866,866,867,867,867,866,866,866,866,
     *866,866,866,865,865,865,866,866,867,867,867,867,866,866,866,
     *866,866,866,866,867,867,866,866,866,866,866,866,866,866,866,
     *866,866,866,866,866/
      DATA (Md07(i),i=15654,15725) /852,856,857,857,859,859,859,860,
     *860,861,862,863,863,863,863,864,864,865,865,865,866,867,866,
     *867,867,868,868,868,869,869,869,870,870,870,870,870,870,870,
     *870,870,870,869,869,870,870,871,871,871,871,871,871,870,871,
     *871,871,870,870,870,870,870,870,870,870,870,870,870,870,870,
     *870,870,870,870/
      DATA (Md07(i),i=15726,15796) /857,858,858,860,859,860,861,861,
     *862,863,864,864,864,864,865,866,867,867,868,868,868,869,869,
     *870,870,870,871,871,872,872,872,872,872,873,873,873,874,874,
     *874,874,873,874,874,874,875,875,875,876,875,875,875,875,875,
     *875,875,874,874,874,874,874,874,874,874,874,874,874,873,873,
     *873,873,873/
      DATA (Md07(i),i=15797,15866) /859,859,861,860,862,862,862,863,
     *864,865,865,865,866,867,867,868,869,870,871,871,871,871,872,
     *872,873,873,874,874,875,875,875,876,876,876,876,876,876,876,
     *876,876,877,877,878,878,878,879,879,878,878,878,878,878,878,
     *878,878,878,878,878,878,877,877,877,877,877,877,877,877,876,
     *876,876/
      DATA (Md07(i),i=15867,15935) /859,860,861,863,863,864,864,865,
     *866,866,866,867,867,868,868,870,871,871,871,872,873,873,873,
     *874,875,875,876,877,877,878,878,878,878,879,878,878,878,879,
     *879,880,880,880,880,881,881,881,881,880,880,885,880,880,881,
     *881,881,881,880,880,880,880,880,880,880,880,880,880,879,879,
     *878/
      DATA (Md07(i),i=15936,16003) /861,862,864,865,865,866,867,867,
     *866,867,867,868,869,871,872,872,873,873,874,875,875,876,876,
     *877,877,878,878,879,879,880,880,881,881,881,881,881,882,882,
     *882,882,882,882,882,882,882,882,882,881,881,881,881,880,880,
     *879,879,879,879,879,879,879,879,879,878,878,878,878,878,877/
      DATA (Md07(i),i=16004,16070) /864,865,866,867,867,868,868,869,
     *870,870,871,872,873,874,874,875,875,876,877,877,878,878,879,
     *880,880,881,881,882,882,883,883,884,884,884,885,885,886,886,
     *886,886,886,886,886,886,886,886,886,886,885,885,885,884,884,
     *884,884,884,883,883,883,883,883,883,883,883,883,883,882/

      DATA (Add07( 9,i,1),i=5,45) /0,0,13917,13918,13927,13946,13971,
     *14001,14034,14074,14130,14190,14255,14324,14397,14472,14550,
     *14630,14712,14796,14879,14961,15042,15122,15201,15279,15356,
     *15432,15507,15581,15654,15726,15797,15867,15936,16004,0,0,0,
     *0,0/

      DATA (Add07( 9,i,2),i=5,45) /0,0,17,14,9,8,8,8,8,8,8,8,9,10,
     *11,12,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
     *33,34,0,0,0,0,0/

      DATA (Add07( 9,i,3),i=5,45) /0,0,17,22,27,32,37,40,47,63,67,
     *72,77,82,85,89,93,96,100,100,100,100,100,100,100,100,100,100,
     *100,100,100,100,100,100,100,100,0,0,0,0,0/

!     Roteiche...............................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=16071,16074) /659,620,580,540/
      DATA (Md07(i),i=16075,16078) /527,490,460,430/
      DATA (Md07(i),i=16079,16081) /616,526,474/
      DATA (Md07(i),i=16082,16087) /813,676,608,566,555,530/
      DATA (Md07(i),i=16088,16094) /874,737,679,642,612,598,585/
      DATA (Md07(i),i=16095,16101) /892,780,723,680,656,643,633/
      DATA (Md07(i),i=16102,16112) /917,824,758,716,697,675,675,680,
     *677,690,700/
      DATA (Md07(i),i=16113,16128) /932,841,785,745,722,709,701,704,
     *708,716,727,727,719,711,701,705/
      DATA (Md07(i),i=16129,16145) /935,869,802,768,746,735,729,729,
     *731,737,745,743,737,732,730,724,722/
      DATA (Md07(i),i=16146,16162) /890,831,798,773,757,750,750,750,
     *758,759,756,752,746,745,741,739,737/
      DATA (Md07(i),i=16163,16180) /905,852,811,789,780,773,771,770,
     *771,771,770,765,760,757,754,753,753,752/
      DATA (Md07(i),i=16181,16200) /919,870,830,810,794,787,785,782,
     *783,781,779,776,770,768,766,765,763,761,762,759/
      DATA (Md07(i),i=16201,16226) /927,883,849,827,811,803,800,797,
     *792,790,789,783,780,778,776,774,773,771,771,768,767,766,764,
     *764,761,760/
      DATA (Md07(i),i=16227,16251) /897,863,841,827,817,814,806,800,
     *797,796,791,788,786,783,782,781,780,779,777,775,774,774,771,
     *770,769/
      DATA (Md07(i),i=16252,16280) /910,874,854,839,830,821,813,811,
     *804,804,799,794,794,791,789,788,787,786,784,783,782,781,780,
     *778,776,776,774,772,771/
      DATA (Md07(i),i=16281,16312) /833,825,820,815,811,809,806,802,
     *800,797,795,794,794,792,791,790,790,789,787,785,784,784,783,
     *781,779,778,777,775,773,772,771,769/
      DATA (Md07(i),i=16313,16345) /837,832,826,820,816,814,810,807,
     *805,804,803,800,799,799,798,797,795,794,792,791,790,790,789,
     *787,786,785,784,782,780,779,780,778,776/
      DATA (Md07(i),i=16346,16375) /823,822,820,815,812,812,809,809,
     *807,805,804,802,801,799,800,798,796,796,796,795,793,792,791,
     *790,788,787,786,785,784,782/
      DATA (Md07(i),i=16376,16411) /825,824,821,817,817,815,814,812,
     *811,810,808,807,805,805,803,801,801,799,800,799,798,797,796,
     *794,793,792,792,790,789,789,789,787,786,785,784,784/
      DATA (Md07(i),i=16412,16443) /822,820,818,817,816,815,813,813,
     *810,810,808,807,807,805,805,803,803,802,801,800,799,798,797,
     *796,796,795,794,793,792,791,790,789/
      DATA (Md07(i),i=16444,16483) /826,824,823,821,821,818,817,816,
     *815,814,813,812,811,810,809,808,808,807,806,805,804,803,803,
     *802,801,800,800,799,798,797,797,796,795,794,793,793,792,791,
     *790,790/
      DATA (Md07(i),i=16484,16527) /830,828,827,825,825,823,822,821,
     *820,819,816,815,815,814,813,813,811,810,810,809,808,807,806,
     *805,804,804,803,802,801,800,799,798,798,797,796,795,794,793,
     *791,791,790,788,787,787/
      DATA (Md07(i),i=16528,16568) /829,830,826,825,824,823,822,821,
     *819,819,818,818,815,815,815,814,812,811,811,811,809,808,808,
     *806,805,805,803,802,802,800,800,799,797,797,796,794,794,793,
     *791,791,790/
      DATA (Md07(i),i=16569,16611) /835,830,829,828,826,825,824,824,
     *822,822,821,820,819,819,818,816,816,815,814,813,813,811,811,
     *810,809,808,808,807,806,805,805,803,803,802,800,800,800,798,
     *798,797,795,795,794/
      DATA (Md07(i),i=16612,16648) /828,827,826,825,825,823,823,821,
     *821,819,819,819,817,817,816,815,814,813,812,811,811,809,809,
     *807,807,805,805,803,803,801,801,799,799,799,797,797,796/
      DATA (Md07(i),i=16649,16681) /827,827,825,825,824,823,822,821,
     *820,820,819,818,817,816,815,815,813,813,811,811,810,809,808,
     *808,806,806,804,804,803,802,801,801,799/
      DATA (Md07(i),i=16682,16704) /822,821,820,820,818,817,817,815,
     *815,814,813,812,812,811,810,809,809,807,807,806,805,804,804/

      DATA (Add07(10,i,1),i=5,45) /0,0,16071,16075,16079,16082,16088,
     *16095,16102,16113,16129,16146,16163,16181,16201,16227,16252,
     *16281,16313,16346,16376,16412,16444,16484,16528,16569,16612,
     *16649,16682,0,0,0,0,0,0,0,0,0,0,0,0/

      DATA (Add07(10,i,2),i=5,45) /0,0,7,7,9,8,8,8,8,8,8,9,9,9,9,10,
     *10,14,14,17,18,22,22,22,25,26,32,36,46,0,0,0,0,0,0,0,0,0,0,0,
     *0/

      DATA (Add07(10,i,3),i=5,45) /0,0,10,10,11,13,14,14,18,23,24,
     *25,26,28,34,34,38,45,46,46,53,53,61,65,65,68,68,68,68,0,0,0,
     *0,0,0,0,0,0,0,0,0/

!     Esche.............................................................

      DATA (Md07(i),i=16705,16706) /560,599/
      DATA (Md07(i),i=16707,16711) /613,632,649,661,634/
      DATA (Md07(i),i=16712,16718) /912,714,693,677,676,676,669/
      DATA (Md07(i),i=16719,16726) /801,778,737,725,706,698,692,681/
      DATA (Md07(i),i=16727,16735) /498,794,758,733,721,712,710,701,
     *688/
      DATA (Md07(i),i=16736,16744) /818,777,750,742,733,723,715,703,
     *694/
      DATA (Md07(i),i=16745,16755) /840,796,777,753,741,732,726,721,
     *718,715,709/
      DATA (Md07(i),i=16756,16768) /857,817,786,767,754,747,741,737,
     *733,728,725,718,711/
      DATA (Md07(i),i=16769,16782) /824,796,778,768,759,753,749,744,
     *740,736,731,727,720,686/
      DATA (Md07(i),i=16783,16798) /827,802,787,776,769,764,759,754,
     *750,744,740,734,729,726,723,718/
      DATA (Md07(i),i=16799,16815) /808,794,784,779,772,768,763,758,
     *753,750,746,742,740,740,740,736,730/
      DATA (Md07(i),i=16816,16832) /816,802,792,785,780,775,770,765,
     *762,758,756,754,752,752,752,753,754/
      DATA (Md07(i),i=16833,16850) /812,803,796,791,785,780,774,771,
     *768,765,763,763,763,763,764,766,768,771/
      DATA (Md07(i),i=16851,16870) /821,812,804,799,794,788,784,780,
     *777,775,774,773,774,774,775,775,777,780,782,785/
      DATA (Md07(i),i=16871,16891) /826,813,805,799,795,791,789,787,
     *784,781,782,782,782,783,784,786,787,789,791,794,796/
      DATA (Md07(i),i=16892,16913) /821,812,805,802,799,796,794,792,
     *790,790,791,791,791,792,793,794,796,797,799,801,803,804/
      DATA (Md07(i),i=16914,16936) /818,811,807,804,803,802,800,799,
     *798,797,798,798,798,800,800,802,803,805,806,807,809,810,811/
      DATA (Md07(i),i=16937,16960) /817,812,810,807,807,806,806,805,
     *804,804,804,804,805,806,807,808,809,810,811,812,813,814,815,
     *816/
      DATA (Md07(i),i=16961,16985) /817,814,813,811,810,811,810,810,
     *809,809,809,810,810,811,812,813,814,815,815,816,817,817,818,
     *818,819/
      DATA (Md07(i),i=16986,17011) /818,816,816,816,815,815,814,814,
     *814,813,813,814,814,815,816,816,817,818,818,819,819,819,820,
     *820,820,820/
      DATA (Md07(i),i=17012,17037) /819,819,819,819,818,818,818,817,
     *817,817,817,817,818,819,819,819,820,820,820,821,821,821,821,
     *821,821,821/
      DATA (Md07(i),i=17038,17064) /823,822,822,821,821,821,821,821,
     *820,820,820,821,821,821,821,821,822,822,822,822,822,822,822,
     *822,821,821,820/
      DATA (Md07(i),i=17065,17091) /824,823,824,824,824,823,823,823,
     *823,823,823,823,823,823,823,823,823,823,823,823,822,822,821,
     *821,820,820,819/
      DATA (Md07(i),i=17092,17119) /825,826,826,826,826,826,826,825,
     *825,825,825,824,824,824,824,824,824,824,823,823,822,822,821,
     *820,820,819,818,817/
      DATA (Md07(i),i=17120,17148) /827,827,828,828,828,827,827,826,
     *826,826,826,825,825,825,825,825,824,824,823,822,822,821,820,
     *820,819,818,817,816,815/
      DATA (Md07(i),i=17149,17178) /829,829,829,829,829,829,828,828,
     *827,827,827,827,826,826,825,825,824,823,823,822,822,820,820,
     *819,818,817,816,815,814,813/
      DATA (Md07(i),i=17179,17209) /830,830,831,830,830,830,829,828,
     *828,828,827,827,826,826,825,824,823,823,822,821,820,819,818,
     *817,816,815,814,813,812,811,810/
      DATA (Md07(i),i=17210,17240) /831,832,832,831,831,830,830,829,
     *829,828,827,827,826,825,824,823,823,822,821,820,819,818,817,
     *816,815,814,813,812,811,809,808/
      DATA (Md07(i),i=17241,17272) /833,833,832,832,831,831,830,829,
     *829,828,827,826,825,824,824,822,822,821,820,819,818,816,815,
     *814,813,812,811,810,809,808,807,806/
      DATA (Md07(i),i=17273,17304) /833,833,833,832,831,830,829,828,
     *828,827,826,825,824,822,821,821,819,818,817,816,815,814,813,
     *812,811,810,808,807,806,805,804,803/
      DATA (Md07(i),i=17305,17337) /834,834,833,832,831,830,829,828,
     *827,826,825,824,823,821,820,819,818,817,816,815,814,812,811,
     *810,809,808,807,805,804,803,802,800,799/
      DATA (Md07(i),i=17338,17369) /833,832,831,830,829,828,827,826,
     *825,824,823,822,820,819,818,817,816,814,813,812,811,810,808,
     *807,806,805,804,802,801,800,798,797/

      DATA (Add07(11,i,1),i=5,45) /0,0,0,0,16705,16707,16712,16719,
     *16727,16736,16745,16756,16769,16783,16799,16816,16833,16851,
     *16871,16892,16914,16937,16961,16986,17012,17038,17065,17092,
     *17120,17149,17179,17210,17241,17273,17305,17338,0,0,0,0,0/

      DATA (Add07(11,i,2),i=5,45) /0,0,0,0,10,9,7,7,7,8,8,8,9,9,10,
     *10,11,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,29,30,
     *32,0,0,0,0,0/

      DATA (Add07(11,i,3),i=5,45) /0,0,0,0,11,13,13,14,15,16,18,20,
     *22,24,26,26,28,30,32,34,36,38,40,42,43,45,47,49,51,53,55,56,
     *58,60,62,63,0,0,0,0,0/

!     Erle....................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=17370,17372) /981,940,900/
      DATA (Md07(i),i=17373,17376) /596,784,692,624/
      DATA (Md07(i),i=17377,17382) /408,841,759,699,668,626/
      DATA (Md07(i),i=17383,17389) /871,807,748,707,690,696,690/
      DATA (Md07(i),i=17390,17397) /908,831,785,745,734,737,744,750/
      DATA (Md07(i),i=17398,17406) /919,854,811,778,763,767,772,777,
     *778/
      DATA (Md07(i),i=17407,17417) /937,871,832,803,787,790,792,794,
     *796,794,783/
      DATA (Md07(i),i=17418,17431) /935,882,850,816,806,806,808,808,
     *808,805,797,790,786,786/
      DATA (Md07(i),i=17432,17450) /941,890,860,833,820,820,821,820,
     *817,813,809,804,798,793,788,782,777,771,766/
      DATA (Md07(i),i=17451,17475) /947,909,864,852,842,831,830,833,
     *829,825,819,815,808,802,797,795,794,792,788,787,783,779,776,
     *771,767/
      DATA (Md07(i),i=17476,17499) /869,856,840,843,843,839,833,829,
     *825,817,813,810,805,804,803,799,794,791,789,788,788,788,788,
     *788/
      DATA (Md07(i),i=17500,17524) /856,855,853,853,848,841,833,828,
     *821,817,814,813,810,807,804,800,796,794,793,792,790,787,785,
     *782,778/
      DATA (Md07(i),i=17525,17547) /863,861,854,847,833,828,823,818,
     *817,814,811,808,805,803,799,797,795,794,792,789,786,783,779/
      DATA (Md07(i),i=17548,17570) /860,860,852,834,828,823,820,817,
     *814,810,809,806,804,801,799,796,795,793,790,787,784,781,779/
      DATA (Md07(i),i=17571,17593) /861,857,849,834,828,823,819,816,
     *813,810,809,806,804,801,800,797,795,794,791,788,785,783,782/
      DATA (Md07(i),i=17594,17617) /854,847,834,828,823,819,816,812,
     *810,808,806,804,802,800,798,796,794,792,789,787,785,784,782,
     *780/
      DATA (Md07(i),i=17618,17640) /835,829,824,819,815,812,810,807,
     *806,804,803,801,798,796,794,792,790,789,787,786,784,782,781/
      DATA (Md07(i),i=17641,17664) /835,829,824,819,816,813,811,808,
     *806,805,803,800,799,797,794,792,791,790,788,787,785,783,782,
     *781/
      DATA (Md07(i),i=17665,17689) /836,830,826,821,817,814,812,809,
     *807,805,803,801,799,797,795,793,793,791,790,788,787,785,784,
     *782,781/
      DATA (Md07(i),i=17690,17714) /831,827,823,819,816,813,810,808,
     *806,804,802,800,798,797,795,794,793,792,790,789,787,786,784,
     *783,781/
      DATA (Md07(i),i=17715,17740) /829,826,821,817,814,812,809,807,
     *805,804,801,800,798,797,796,795,794,792,791,790,788,787,785,
     *784,782,780/
      DATA (Md07(i),i=17741,17765) /828,823,819,816,813,810,809,808,
     *806,804,802,800,800,799,798,796,795,794,793,791,790,788,787,
     *785,784/
      DATA (Md07(i),i=17766,17789) /826,822,819,816,812,811,810,809,
     *807,804,803,802,801,800,799,798,797,796,795,793,792,791,789,
     *788/

      DATA (Add07(12,i,1),i=5,45) /0,0,0,17370,17373,17377,17383,17390,
     *17398,17407,17418,17432,17451,17476,17500,17525,17548,17571,
     *17594,17618,17641,17665,17690,17715,17741,17766,0,0,0,0,0,0,
     *0,0,0,0,0,0,0,0,0/

      DATA (Add07(12,i,2),i=5,45) /0,0,0,7,7,7,8,8,8,8,8,8,8,11,12,
     *14,15,15,16,18,18,18,19,20,21,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     *0/

      DATA (Add07(12,i,3),i=5,45) /0,0,0,9,10,12,14,15,16,18,21,26,
     *32,34,36,36,37,37,39,40,41,42,43,45,45,45,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,0,0/

!     Birke...................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=17790,17791) /431,400/
      DATA (Md07(i),i=17792,17795) /402,400,400,400/
      DATA (Md07(i),i=17796,17796) /494/
      DATA (Md07(i),i=17797,17801) /953,678,699,502,522/
      DATA (Md07(i),i=17802,17807) /914,666,677,676,661,643/
      DATA (Md07(i),i=17808,17830) /891,849,816,786,759,843,907,847,
     *791,739,693,649,662,619,579,583,583,579,574,565,556,544,532/
      DATA (Md07(i),i=17831,17852) /823,785,755,828,786,835,785,735,
     *691,649,662,669,627,627,624,618,609,599,588,576,563,570/
      DATA (Md07(i),i=17853,17879) /933,761,732,787,748,784,738,694,
     *710,718,719,715,671,664,655,644,632,619,605,591,596,599,582,
     *582,579,576,571/
      DATA (Md07(i),i=17880,17906) /906,850,798,840,791,815,765,718,
     *725,727,722,715,705,693,679,665,650,635,619,622,623,622,603,
     *600,595,589,583/
      DATA (Md07(i),i=17907,17942) /823,775,805,826,839,787,738,737,
     *733,725,715,702,689,674,659,665,647,649,630,629,625,621,615,
     *608,601,604,595,595,594,593,590,586,582,584,578,579/
      DATA (Md07(i),i=17943,17978) /890,829,781,795,802,754,753,748,
     *740,728,745,728,711,693,675,677,658,657,653,648,642,635,628,
     *631,621,621,620,618,615,611,614,609,610,603,602,601/
      DATA (Md07(i),i=17979,18014) /864,807,818,769,774,772,767,757,
     *744,731,742,724,706,709,689,688,684,679,673,666,657,648,650,
     *640,639,637,634,630,633,628,629,622,621,619,617,614/
      DATA (Md07(i),i=18015,18050) /914,850,850,799,795,788,778,765,
     *750,733,741,722,703,703,701,697,691,685,676,668,670,660,659,
     *657,654,650,646,648,642,642,641,639,637,633,630,625/
      DATA (Md07(i),i=18051,18085) /888,881,822,813,802,787,771,782,
     *761,763,741,720,717,712,705,698,689,680,670,671,670,667,664,
     *660,662,656,657,656,655,652,649,645,646,641,641/
      DATA (Md07(i),i=18086,18119) /906,843,829,812,795,777,784,761,
     *761,738,734,729,721,713,704,694,695,683,682,679,675,678,672,
     *673,673,671,669,666,662,658,658,653,652,650/
      DATA (Md07(i),i=18120,18152) /860,842,791,773,782,761,762,759,
     *754,747,739,730,720,709,709,697,695,691,687,690,683,684,683,
     *681,678,675,670,671,671,670,668,666,663/
      DATA (Md07(i),i=18153,18184) /821,802,781,786,786,783,776,768,
     *759,748,737,726,714,712,709,705,700,694,695,695,694,691,688,
     *684,685,680,679,678,676,673,670,670/
      DATA (Md07(i),i=18185,18215) /811,789,790,788,782,775,765,755,
     *743,732,731,729,725,720,715,708,709,708,706,703,699,695,690,
     *690,689,687,685,681,682,682,681/
      DATA (Md07(i),i=18216,18245) /795,794,789,782,773,762,751,752,
     *738,736,733,728,722,715,715,714,712,709,705,701,701,701,699,
     *697,694,691,691,690,689,687/
      DATA (Md07(i),i=18246,18274) /818,809,799,787,774,761,759,745,
     *741,736,730,723,715,715,713,710,712,708,708,707,705,703,700,
     *701,700,699,698,695,696/
      DATA (Md07(i),i=18275,18302) /809,797,784,771,757,755,751,745,
     *739,732,732,723,721,718,714,715,715,714,712,710,711,707,706,
     *705,703,704,701,701/
      DATA (Md07(i),i=18303,18328) /783,768,754,751,746,740,733,734,
     *733,730,727,723,724,724,722,720,717,718,718,717,716,713,711,
     *711,707,706/
      DATA (Md07(i),i=18329,18351) /758,751,744,744,743,741,737,733,
     *733,733,731,729,726,726,726,725,723,721,718,718,717,715,713/
      DATA (Md07(i),i=18352,18372) /748,747,745,741,737,738,737,736,
     *733,730,731,730,729,727,725,725,725,724,722,720,718/

      DATA (Add07(13,i,1),i=5,45) /0,0,17790,17792,17796,17797,17802,
     *17808,17831,17853,17880,17907,17943,17979,18015,18051,18086,
     *18120,18153,18185,18216,18246,18275,18303,18329,18352,0,0,0,
     *0,0,0,0,0,0,0,0,0,0,0,0/

      DATA (Add07(13,i,2),i=5,45) /0,0,7,7,11,8,8,8,9,9,9,10,10,10,
     *10,11,12,13,14,15,16,17,18,20,23,25,0,0,0,0,0,0,0,0,0,0,0,0,
     *0,0,0/

      DATA (Add07(13,i,3),i=5,45) /0,0,8,10,11,12,13,30,30,35,35,45,
     *45,45,45,45,45,45,45,45,45,45,45,45,45,45,0,0,0,0,0,0,0,0,0,
     *0,0,0,0,0,0/

!     Pappel..................................(Aenderung: 11.12.90 Kublin)

      DATA (Md07(i),i=18373,18381) /911,870,830,790,750,710,670,630,
     *590/
      DATA (Md07(i),i=18382,18389) /456,420,408,400,400,400,400,400/
      DATA (Md07(i),i=18390,18400) /616,591,549,519,477,443,416,400,
     *400,400,400/
      DATA (Md07(i),i=18401,18414) /813,676,652,604,586,557,520,520,
     *490,479,464,448,436,427/
      DATA (Md07(i),i=18415,18429) /825,737,679,656,626,610,585,573,
     *559,542,529,521,509,500,495/
      DATA (Md07(i),i=18430,18446) /853,762,709,680,667,643,624,612,
     *598,591,579,571,559,555,543,538,530/
      DATA (Md07(i),i=18447,18465) /872,777,733,706,688,666,649,649,
     *635,621,616,606,599,589,585,578,570,565,559/
      DATA (Md07(i),i=18466,18485) /872,789,752,717,713,693,679,677,
     *658,649,643,634,626,620,614,607,603,597,592,585/
      DATA (Md07(i),i=18486,18507) /882,800,762,734,722,714,696,693,
     *681,672,665,656,651,644,637,632,627,622,618,612,607,603/
      DATA (Md07(i),i=18508,18531) /896,807,769,743,738,725,714,712,
     *701,691,683,677,671,663,659,652,649,643,640,635,630,626,622,
     *618/
      DATA (Md07(i),i=18532,18557) /892,815,778,755,744,740,730,727,
     *716,707,699,692,689,683,677,672,667,662,658,654,650,646,642,
     *637,633,630/
      DATA (Md07(i),i=18558,18585) /897,820,786,758,757,752,743,740,
     *729,724,715,708,703,697,693,687,683,679,675,670,667,663,659,
     *656,652,648,645,641/
      DATA (Md07(i),i=18586,18616) /900,824,786,770,767,762,754,752,
     *745,735,730,721,716,711,706,702,699,694,691,686,683,679,675,
     *671,668,664,661,659,655,652,648/
      DATA (Md07(i),i=18617,18650) /899,829,796,778,776,772,764,765,
     *754,748,739,734,730,724,720,714,711,707,704,700,697,694,689,
     *686,684,680,677,674,670,668,665,661,658,655/
      DATA (Md07(i),i=18651,18686) /905,829,802,783,784,780,773,773,
     *763,757,750,746,740,735,731,728,723,719,716,713,709,706,703,
     *700,697,694,691,688,685,682,680,676,674,671,668,665/
      DATA (Md07(i),i=18687,18724) /913,832,805,791,790,787,783,780,
     *773,767,760,755,751,745,742,738,734,731,728,724,722,718,715,
     *713,709,707,704,701,699,695,693,690,687,685,683,680,677,675/
      DATA (Md07(i),i=18725,18763) /839,804,794,797,793,790,790,782,
     *776,771,766,760,756,752,749,745,741,738,735,733,729,726,724,
     *722,719,716,714,710,708,706,703,701,698,696,693,691,688,686,
     *684/
      DATA (Md07(i),i=18764,18804) /813,800,803,802,795,798,788,784,
     *779,773,769,766,761,758,755,752,748,745,743,740,738,734,732,
     *730,727,724,722,720,718,715,713,710,708,706,704,701,699,697,
     *695,693,690/
      DATA (Md07(i),i=18805,18848) /809,810,808,804,802,796,791,786,
     *782,779,773,770,767,763,760,758,755,752,750,747,745,743,740,
     *738,736,733,731,729,726,724,722,720,718,716,713,711,709,707,
     *705,703,701,699,697,695/
      DATA (Md07(i),i=18849,18938) /814,816,815,812,810,803,798,794,
     *790,786,781,779,775,772,769,767,765,762,759,757,755,752,750,
     *748,746,743,741,739,737,735,733,731,729,727,725,723,721,719,
     *717,715,713,711,709,708,706,704,702,700,698,696,694,692,691,
     *689,687,685,683,682,680,678,676,674,672,671,669,667,665,664,
     *662,660,658,656,655,653,651,649,648,646,644,642,641,639,637,
     *635,634,632,630,628,627,625/
      DATA (Md07(i),i=18939,19027) /823,823,815,816,809,806,802,797,
     *794,790,786,784,781,777,775,773,770,768,766,763,761,760,757,
     *755,753,751,749,747,745,744,742,739,738,736,734,732,730,729,
     *727,725,723,721,719,718,716,714,712,711,709,707,705,704,702,
     *700,699,697,695,694,692,690,688,687,685,683,682,680,679,677,
     *675,674,672,670,669,667,665,664,662,660,659,657,655,654,652,
     *651,649,647,646,644,642/
      DATA (Md07(i),i=19028,19115) /825,822,822,817,812,807,805,800,
     *797,793,790,788,786,783,781,778,776,774,772,770,768,766,764,
     *763,761,759,757,755,753,751,750,748,746,744,743,741,739,737,
     *736,734,733,731,729,728,726,724,723,721,719,718,716,714,713,
     *711,710,708,707,705,704,702,700,699,697,696,694,693,691,689,
     *688,686,685,683,682,680,679,677,676,674,673,671,670,668,667,
     *665,664,662,661,659/
      DATA (Md07(i),i=19116,19201) /828,822,818,813,810,806,804,801,
     *797,795,792,791,788,786,784,782,781,778,776,775,773,771,769,
     *768,766,764,763,761,759,757,756,754,753,751,750,748,746,745,
     *743,742,740,739,737,735,734,733,731,729,728,726,725,723,722,
     *721,719,718,716,715,713,712,710,709,707,706,705,703,702,700,
     *699,697,696,695,693,692,690,689,687,686,685,683,682,680,679,
     *678,676,675/
      DATA (Md07(i),i=19202,19286) /827,823,819,816,813,810,808,805,
     *802,800,797,796,794,792,790,788,786,785,783,781,780,778,776,
     *775,773,771,770,768,767,765,764,762,761,759,758,756,755,753,
     *752,751,749,748,746,745,743,742,741,739,738,736,735,734,732,
     *731,730,728,727,725,724,723,721,720,719,717,716,715,713,712,
     *711,709,708,707,705,704,703,702,700,699,698,696,695,694,692,
     *691,690/
      DATA (Md07(i),i=19287,19369) /826,822,820,817,814,812,809,807,
     *804,802,801,799,797,795,794,792,791,789,787,786,784,783,781,
     *780,778,777,776,774,773,771,770,769,767,766,765,763,762,760,
     *759,758,757,755,754,752,751,750,749,747,746,745,744,742,741,
     *740,738,737,736,735,733,732,731,730,728,727,726,725,724,722,
     *721,720,719,717,716,715,714,713,711,710,709,708,707,705,704/
      DATA (Md07(i),i=19370,19450) /826,823,820,817,815,814,812,810,
     *808,806,804,803,801,800,798,797,795,794,792,791,790,788,787,
     *785,784,783,781,780,779,778,776,775,774,773,771,770,769,767,
     *766,765,764,763,761,760,759,758,757,755,754,753,752,751,750,
     *748,747,746,745,744,743,741,740,739,738,737,736,735,733,732,
     *731,730,729,728,727,726,724,723,722,721,720,719,718/
      DATA (Md07(i),i=19451,19529) /826,824,822,820,818,816,815,813,
     *811,809,808,807,805,804,802,801,800,799,797,796,795,794,792,
     *791,790,789,787,786,785,784,783,781,780,779,778,777,776,774,
     *773,772,771,770,769,768,767,766,764,763,762,761,760,759,758,
     *757,756,755,754,753,752,751,750,749,747,746,745,744,743,742,
     *741,740,739,738,737,736,735,734,733,732,731/
      DATA (Md07(i),i=19530,19607) /830,828,826,824,822,821,819,818,
     *816,815,814,812,811,809,808,807,806,805,804,802,801,800,799,
     *798,797,796,794,793,792,791,790,789,788,787,786,785,784,783,
     *782,781,779,778,777,776,775,774,773,772,771,770,769,768,767,
     *766,765,764,763,762,762,761,760,759,758,757,756,755,754,753,
     *752,751,750,749,748,747,746,745,744,743/
      DATA (Md07(i),i=19608,19683) /832,830,828,827,826,824,823,821,
     *820,819,818,817,815,814,813,812,811,810,809,808,806,805,804,
     *803,802,801,800,799,798,797,796,795,794,793,792,791,790,789,
     *789,788,787,786,785,784,783,782,781,780,779,778,777,776,776,
     *775,774,773,772,771,770,769,768,768,767,766,765,764,763,762,
     *761,761,760,759,758,757,756,755/
      DATA (Md07(i),i=19684,19757) /834,833,832,830,829,828,826,825,
     *824,823,822,821,820,819,818,817,816,815,814,813,812,811,810,
     *809,808,807,806,805,804,803,803,802,801,800,799,798,797,796,
     *796,795,794,793,792,791,790,789,789,788,787,786,785,784,784,
     *783,782,781,780,780,779,778,777,776,776,775,774,773,772,772,
     *771,770,769,768,768,767/
      DATA (Md07(i),i=19758,19830) /839,837,836,835,834,833,832,831,
     *829,829,828,826,825,825,824,823,822,821,820,819,818,817,816,
     *816,815,814,813,812,811,811,810,809,808,807,806,806,805,804,
     *803,802,802,801,800,799,799,798,797,796,796,795,794,793,792,
     *792,791,790,789,789,788,787,787,786,785,784,784,783,782,781,
     *781,780,779,779,778/
      DATA (Md07(i),i=19831,19901) /842,841,840,839,838,837,836,835,
     *834,833,832,831,830,829,829,828,827,826,825,824,824,823,822,
     *821,821,820,819,818,818,817,816,815,815,814,813,812,812,811,
     *810,809,809,808,807,807,806,805,804,804,803,802,802,801,800,
     *800,799,798,798,797,796,796,795,794,794,793,792,792,791,790,
     *790,789,788/
      DATA (Md07(i),i=19902,19970) /846,845,844,843,842,841,840,839,
     *838,838,837,836,835,834,834,833,832,831,831,830,829,829,828,
     *827,826,826,825,824,824,823,822,822,821,820,820,819,818,818,
     *817,816,816,815,814,814,813,813,812,811,811,810,809,809,808,
     *808,807,806,806,805,805,804,803,803,802,802,801,800,800,799,
     *799/
      DATA (Md07(i),i=19971,20037) /849,848,848,847,846,845,845,844,
     *843,842,842,841,840,840,839,838,837,837,836,836,835,834,834,
     *833,832,832,831,830,830,829,829,828,827,827,826,826,825,824,
     *824,823,823,822,822,821,820,820,819,819,818,818,817,817,816,
     *815,815,814,814,813,813,812,812,811,811,810,809,809,808/
      DATA (Md07(i),i=20038,20102) /853,853,852,851,850,850,849,848,
     *848,847,847,846,845,845,844,844,843,842,842,841,841,840,839,
     *839,838,838,837,837,836,835,835,834,834,833,833,832,832,831,
     *831,830,830,829,829,828,828,827,827,826,826,825,825,824,824,
     *823,823,822,822,821,821,820,820,819,819,818,818/
      DATA (Md07(i),i=20103,20165) /858,857,856,856,855,854,854,853,
     *853,852,852,851,851,850,849,849,848,848,847,847,846,846,845,
     *845,844,844,843,843,842,842,841,841,840,840,839,839,838,838,
     *838,837,837,836,836,835,835,834,834,834,833,833,832,832,831,
     *831,831,830,830,829,829,828,828,828,827/
      DATA (Md07(i),i=20166,20226) /862,862,861,860,860,859,859,858,
     *858,857,857,856,856,855,855,854,854,853,853,852,852,852,851,
     *851,850,850,849,849,849,848,848,847,847,846,846,846,845,845,
     *844,844,844,843,843,842,842,842,841,841,840,840,840,839,839,
     *839,838,838,838,837,837,836,836/
      DATA (Md07(i),i=20227,20284) /866,866,865,865,864,864,863,863,
     *862,862,861,861,861,860,860,859,859,859,858,858,857,857,857,
     *856,856,856,855,855,854,854,854,853,853,853,852,852,852,851,
     *851,851,850,850,850,849,849,849,848,848,848,847,847,847,846,
     *846,846,845,845,845/
      DATA (Md07(i),i=20285,20333) /868,868,867,867,866,866,866,865,
     *865,865,864,864,864,863,863,863,862,862,862,861,861,861,860,
     *860,860,860,859,859,859,858,858,858,857,857,857,857,856,856,
     *856,855,855,855,855,854,854,854,854,853,853/

      DATA (Add07(14,i,1),i=5,45) /0,0,18373,18382,18390,18401,18415,
     *18430,18447,18466,18486,18508,18532,18558,18586,18617,18651,
     *18687,18725,18764,18805,18849,18939,19028,19116,19202,19287,
     *19370,19451,19530,19608,19684,19758,19831,19902,19971,20038,
     *20103,20166,20227,20285/

      DATA (Add07(14,i,2),i=5,45) /0,0,7,10,9,8,8,8,8,8,8,8,8,8,8,
     *8,8,8,9,10,11,11,12,13,15,16,18,20,22,23,25,27,28,30,32,34,36,
     *38,40,43,52/

      DATA (Add07(14,i,3),i=5,45) /0,0,15,17,19,21,22,24,26,27,29,
     *31,33,35,38,41,43,45,47,50,54,100,100,100,100,100,100,100,100,
     *100,100,100,100,100,100,100,100,100,100,100,100/

      END
