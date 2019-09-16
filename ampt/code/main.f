c.....driver program for A Multi-Phase Transport model
      PROGRAM AMPT
c
      PARAMETER     (MAXX   =   20,  MAXZ  =    24) !xj dencon
      double precision xmp, xmu, alpha, rscut2, cutof2, dshadow
      double precision smearp,smearh,dpcoal,drcoal,ecritl
      CHARACTER FRAME*8, PROJ*8, TARG*8
      character*25 amptvn
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON /HPARNT/HIPR1(100), IHPR2(50), HINT1(100), IHNT2(50)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /AROUT/ IOUT
      COMMON /AREVT/ IAEVT, IARUN, MISS
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common/anim/nevent,isoft,isflag,izpc
c     parton coalescence radii in case of string melting:
      common /coal/dpcoal,drcoal,ecritl
      common/snn/efrm,npart1,npart2,epsiPz,epsiPt,PZPROJ,PZTARG
c     initialization value for parton cascade:
      common /para2/ xmp, xmu, alpha, rscut2, cutof2
      common /para7/ ioscar,nsmbbbar,nsmmeson
      common /para8/ idpert,npertd,idxsec
      common /rndm3/ iseedp
c     initialization value for hadron cascade:
      COMMON /RUN/ NUM
	COMMON /EVENT/ NEVNT !xj dencon
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON /INPUT2/ ILAB, MANYB, NTMAX, ICOLL, INSYS, IPOT, MODE, 
     &   IMOMEN, NFREQ, ICFLOW, ICRHO, ICOU, KPOTEN, KMUL
      common/oscar1/iap,izp,iat,izt
      common/oscar2/FRAME,amptvn
      common/resdcy/NSAV,iksdcy
clin-6/2009:
c      common/phidcy/iphidcy
      common/phidcy/iphidcy,pttrig,ntrig,maxmiss
      common/embed/iembed,nsembd,pxqembd,pyqembd,xembd,yembd,
     1     psembd,tmaxembd,phidecomp
clin-7/2009:
      common/cmsflag/dshadow,ishadow

	COMMON  /POT/UB(2001,2001),UBBAR(2001,2001),SCA(2001,2001) !xj
     &,PIM(2001,2001),PI0(2001,2001),PIP(2001,2001),DRHO,RHOMAX

	COMMON/DenConxj/ DenConxz(30,-MAXX:MAXX,-MAXX:MAXX)
     &,DenaConxz(30,-MAXX:MAXX,-MAXX:MAXX)
     &,DenConxy(30,-MAXX:MAXX,-MAXX:MAXX)
     &,DenaConxy(30,-MAXX:MAXX,-MAXX:MAXX) !xj
      common/grid1/nx, ny, nz
	COMMON/DenConnjl/ DenConnjlxz(10,100,40)
     &,DenaConnjlxz(10,100,40),DenConnjlxy(10,100,100)
     &,DenaConnjlxy(10,100,100) !xj
	DOUBLE PRECISION DenConnjlxz,DenaConnjlxz,DenConnjlxy,DenaConnjlxy

	COMMON/HIS/HISDP(20),HISDR(20)
     &,HISMDP(20),HISMDR(20),HISBDP(20),HISBDR(20)
      COMMON/DHIS/DPHIS,DRHIS !xj: histogram of dp and dr for coalescence

clin-10/2011 ctest on ppbar: read in 3 factors:
c     1st factor is being multiplied to cross sections of
c     lightB+lightBbar -> M+M & lightB+strangeB -> M+K, 
c     which decrease the baryon# and the anti-baryon#.
c     2nd factor is being multiplied to cross sections of
c     lightB+lightBbar <- M+M & lightB+strangeB <- M+K, 
c     which increase the baryon# and the anti-baryon#.
c     3rd factor is being multiplied to cross sections of
c     lightB+M <-> strangeB+K, lightB+K <-> strangeB+M, 
c     strangeB+K <-> moreStrangeB+M, lightB+phi <-> strangeB+K,
c     and lightB+lightB -> lightB+strangeB+K, 
c     which do not change the baryon# and the anti-baryon# 
c     but trade light-(anti)baryon# with strange-(anti)baryon#.
      common/ppbar/xbbann,xbbcre,xbstrg

      EXTERNAL HIDATA, PYDATA, LUDATA, ARDATA, PPBDAT, zpcbdt
      SAVE   
c****************
      OPEN (24, FILE = '../input/input.ampt'
     &, STATUS = 'UNKNOWN')
      OPEN (12, FILE = 'ana/version', STATUS = 'UNKNOWN')
      READ (24, *) EFRM
c     format-read characters (for ALPHA compilers):
      READ (24, 111) FRAME
      READ (24, 111) PROJ
      READ (24, 111) TARG
      READ (24, *) IAP
      READ (24, *) IZP
      READ (24, *) IAT
      READ (24, *) IZT
      READ (24, *) NEVNT
      READ (24, *) BMIN
      READ (24, *) BMAX
c     flag to select default AMPT or string melting:
      READ (24, *) isoft
c     read initialization value for hadron cascade:
      READ (24, *) NTMAX
      READ (24, *) DT
c     parj(41) and (42) are a and b parameters in Lund string fragmentation:
      READ (24, *) PARJ(41)
      READ (24, *) PARJ(42)
c     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and 
c     increase the net-baryon stopping in rapidity (value HIJING is 1):
      READ (24, *) ipop
      if(ipop.eq.1) IHPR2(11)=3
c     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
      READ (24, *) PARJ(5)
c     shadowing flag in HIJING:
      READ (24, *) IHPR2(6)
c     quenching flag in HIJING:
      READ (24, *) IHPR2(4)
c     quenching rate when quenching flag is on (=1.0 GeV/fm):
      READ (24, *) HIPR1(14)
c     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV. 
      READ (24, *) HIPR1(8)
c     read initialization value for parton cascade:
      READ (24, *) xmu
      READ (24, *) izpc
      READ (24, *) alpha
c     quark coalescence radii in momentum and space for string melting:
      READ (24, *) dpcoal
      READ (24, *) drcoal
c     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
      READ (24, *) ihjsed
c     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
      READ (24, *) nseed
      READ (24, *) iseedp
      READ (24, *) iksdcy
      READ (24, *) iphidcy
c     flag for OSCAR output for final partons and hadrons:
      READ (24, *) ioscar
clin-5/2008     flag for perturbative treatment of deuterons:
      READ (24, *) idpert
      READ (24, *) npertd
      READ (24, *) idxsec
clin-6/2009 To select events that have at least 1 high-Pt minijet parton:
      READ (24, *) pttrig
      READ (24, *) maxmiss
      READ (24, *) IHPR2(2)
      READ (24, *) IHPR2(5)
clin-6/2009 To embed a back-to-back q/qbar pair into each event:
      READ (24, *) iembed
      READ (24, *) pxqembd, pyqembd
      READ (24, *) xembd, yembd
      READ (24, *) nsembd,psembd,tmaxembd
clin-7/2009 Allow modification of nuclear shadowing:
      READ (24, *) ishadow
      READ (24, *) dshadow

clin-10/2011 ctest on ppbar:
      READ (24, *) xbbann,xbbcre,xbstrg
c
      CLOSE (24)
 111  format(a8)
clin-6/2009 ctest off turn on jet triggering:
c      IHPR2(3)=1
c     Trigger Pt of high-pt jets in HIJING:
c      HIPR1(10)=7.
c
      if(isoft.eq.1) then
         amptvn = '1.25t7 (Default)'
      elseif(isoft.eq.4) then
         amptvn = '2.25t7 (StringMelting)'
      else
         amptvn = 'Test-Only'
      endif
      WRITE(6,50) amptvn
      WRITE(12,50) amptvn
 50   FORMAT(' '/
     &11X,'##################################################'/1X,
     &10X,'#      AMPT (A Multi-Phase Transport) model      #'/1X,
     &10X,'#          Version ',a25,                  '     #'/1X,
     &10X,'#                09/08/2011                      #'/1X,
     &10X,'##################################################'/1X,
     &10X,' ')
c     when ihjsed=11: use environment variable at run time for HIJING nseed:
      if(ihjsed.eq.11) then
         PRINT *,
     1 '# Read in NSEED in HIJING at run time (e.g. 20030819):'
      endif
      READ (*, *) nseedr
      if(ihjsed.eq.11) then
         nseed=nseedr
      endif
c     an odd number is needed for the random number generator:
      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
      if(ihjsed.eq.11) then      
         PRINT *, '#   read in: ', nseed
         WRITE(12,*) '# Read in NSEED in HIJING at run time:',nseed
      endif
      CLOSE(12)
c     9/26/03 random number generator for f77 compiler:
      CALL SRAND(NSEED)
c
c.....turn on warning messages in nohup.out when an event is repeated:
      IHPR2(10) = 1
c     string formation time:
      ARPAR1(1) = 0.7
c     smearp is the smearing halfwidth on parton z0, 
c     set to 0 for now to avoid overflow in eta.
c     smearh is the smearing halfwidth on string production point z0.
      smearp=0d0
      IAmax=max(iap,iat)
      smearh=1.2d0*IAmax**0.3333d0/(dble(EFRM)/2/0.938d0)
      nevent=NEVNT
c
c     AMPT momentum and space info at freezeout:
      OPEN (16, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      OPEN (14, FILE = 'ana/zpc.dat', STATUS = 'UNKNOWN')
cxj2011begin
	OPEN(10,FILE='ana/initial_ampt.dat')
	OPEN(97,FILE='den.txt')
	OPEN(98,FILE='den_NJL.txt')
        OPEN(11,FILE='Ach_d.dat')

	OPEN(1,FILE='../input/BARYON.TXT')
	OPEN(2,FILE='../input/ANTIBARYON.TXT')
	OPEN(3,FILE='../input/SCALAR.TXT')
	OPEN(4,FILE='../input/PIONM.TXT')
	OPEN(5,FILE='../input/PION0.TXT')
	OPEN(7,FILE='../input/PIONP.TXT')
	DO I = 1,2001
	DO J = 1,2001
	READ(1,*) UB(I,J)
	READ(2,*) UBBAR(I,J)
	READ(3,*) SCA(I,J)
	READ(4,*) PIM(I,J)
	READ(5,*) PI0(I,J)
	READ(7,*) PIP(I,J)
        ENDDO
	ENDDO
	CLOSE(1)
	CLOSE(2)
	CLOSE(3)
	CLOSE(4)
	CLOSE(5)
	CLOSE(7)
	!test_begin
c	OPEN(1,FILE='BARYON_DISPLAY.TXT')
c	OPEN(2,FILE='ANTIBARYON_DISPLAY.TXT')
c	OPEN(3,FILE='SCALAR_DISPLAY.TXT')
c	OPEN(4,FILE='PIONM_DISPLAY.TXT')
c	OPEN(5,FILE='PION0_DISPLAY.TXT')
c	OPEN(6,FILE='PIONP_DISPLAY.TXT')
	RHO0 = 0.16
	DRHO = 0.001*RHO0
	RHOMAX = 2.*RHO0
	!xj coalescence histogram
	DPHIS = 0.1
	DRHIS = 0.5
	DO I = 1,20
	HISMDP(I)=0
	HISMDR(I)=0
	HISBDP(I)=0
	HISBDR(I)=0
	HISDP(I)=0
	HISDR(I)=0
	ENDDO
c	DO RHO1 = 0.1*RHO0,RHOMAX,0.1*RHO0
c		DO RHO2 = 0.1*RHO0,RHOMAX,0.1*RHO0
c	I = RHO1/DRHO + 1
c	J = RHO2/DRHO + 1
c	WRITE(1,'(F10.3,$)') UB(I,J)
c	WRITE(2,'(F10.3,$)') UBBAR(I,J)
c	WRITE(3,'(F10.3,$)') SCA(I,J)
c	WRITE(4,'(F10.3,$)') PIM(I,J)
c	WRITE(5,'(F10.3,$)') PI0(I,J)
c	WRITE(6,'(F10.3,$)') PIP(I,J)
c		ENDDO
c	WRITE(1,*)
c	WRITE(2,*)
c	WRITE(3,*)
c	WRITE(4,*)
c	WRITE(5,*)
c	WRITE(6,*)
c	ENDDO
c	CLOSE(1)
c	CLOSE(2)
c	CLOSE(3)
c	CLOSE(4)
c	CLOSE(5)
c	CLOSE(6)
c	stop
	!test_end
cxj2011end
ctest off for resonance (phi, K*) studies:
c      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN')
c      OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN')
      CALL HIJSET(EFRM, FRAME, PROJ, TARG, IAP, IZP, IAT, IZT)
      CALL ARTSET
	NEVNT = NEVNT/NUM !xj
      CALL INIZPC
clin-5/2009 ctest off:
c      call flowp(0)
c      call flowh0(NEVNT,0)
c      call iniflw(NEVNT,0)
c      call frztm(NEVNT,0)
c
       DO 2000 J = 1, NEVNT
          IAEVT = J
c          DO 1000 K = 1, NUM
c             IARUN = K
c             IF (IAEVT .EQ. NEVNT .AND. IARUN .EQ. NUM) THEN
c                IOUT = 1
c             END IF
             PRINT *, ' EVENT ', J !, ', RUN ', K
c             imiss=0
             write(11,*) 'EVENT',J
 100         CALL HIJING(EFRM,FRAME, BMIN, BMAX)
c             IAINT2(1) = NATT !xj use MULTI1         
clin-6/2009 ctest off
c           if(J.eq.-2) then 
c              write(98,*) HIPR1
c              write(98,*) ' '
c              write(98,*) IHPR2
c              write(98,*) ' '
c              write(98,*) (HINT1(i),i=1,20)
c              write(98,*) ' '
c              write(98,*) (HINT1(i),i=21,40)
c              write(98,*) ' '
c              write(98,*) (HINT1(i),i=41,60)
c              write(98,*) ' '
c              write(98,*) (HINT1(i),i=61,80)
c              write(98,*) ' '
c              write(98,*) (HINT1(i),i=81,100)
c              write(98,*) ' '
c              write(98,*) IHNT2
c           endif

c     evaluate Npart (from primary NN collisions) for both proj and targ:
c             call getnp
c     switch for final parton fragmentation:
c             IF (IHPR2(20) .EQ. 0) GOTO 2000
c     In the unlikely case of no interaction (even after loop of 20 in HIJING),
c     still repeat the event to get an interaction 
c     (this may have an additional "trigger" effect):
c             if(NATT.eq.0) then
c                imiss=imiss+1
c                if(imiss.le.20) then
c                   write(6,*) 'repeated event: natt=0,j,imiss=',j,imiss
c                   goto 100
c                else
c                   write(6,*) 'missed event: natt=0,j=',j
c                   goto 2000
c                endif
c             endif
c.....ART initialization and run
	   DO 1000 K = 1,NUM
             CALL ARINI(K)
             CALL ARINI2(K)
1000     CONTINUE
c
C          CALL ARTAN1
C          CALL HJANA3
          CALL ARTMN
C          CALL HJANA4
C          CALL ARTAN2
 2000  CONTINUE
c
C       CALL ARTOUT(NEVNT)
clin-5/2009 ctest off:
c       call flowh0(NEVNT,2)
c       call flowp(2)
c       call iniflw(NEVNT,2)
c       call frztm(NEVNT,2)
c
cxj2011begin
	CLOSE(10)
	CLOSE(97)
	CLOSE(98)
	CLOSE(70)
	CLOSE(71)
        CLOSE(11)

cxj2011end

      do ntest = 1,30
	IF(ntest.eq.1) then
      OPEN (1, FILE = 'ana/RHXZt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (1, FILE = 'ana/RHXZt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (1, FILE = 'ana/RHXZt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (1, FILE = 'ana/RHXZt4.dat', STATUS = 'UNKNOWN')
      ELSE IF(ntest.eq.5) then
	OPEN (1, FILE = 'ana/RHXZt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then      
	OPEN (1, FILE = 'ana/RHXZt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then      
	OPEN (1, FILE = 'ana/RHXZt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (1, FILE = 'ana/RHXZt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (1, FILE = 'ana/RHXZt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (1, FILE = 'ana/RHXZt10.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.11) then	   
      OPEN (1, FILE = 'ana/RHXZt11.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.12) then	   
      OPEN (1, FILE = 'ana/RHXZt12.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.13) then	  
      OPEN (1, FILE = 'ana/RHXZt13.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.14) then	   
      OPEN (1, FILE = 'ana/RHXZt14.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.15) then	   
      OPEN (1, FILE = 'ana/RHXZt15.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.16) then	   
      OPEN (1, FILE = 'ana/RHXZt16.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.17) then	   
      OPEN (1, FILE = 'ana/RHXZt17.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.18) then	  
      OPEN (1, FILE = 'ana/RHXZt18.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.19) then	   
      OPEN (1, FILE = 'ana/RHXZt19.dat',STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.20) then	    
      OPEN (1, FILE = 'ana/RHXZt20.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.21) then	   
      OPEN (1, FILE = 'ana/RHXZt21.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.22) then	  
      OPEN (1, FILE = 'ana/RHXZt22.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.23) then	   
      OPEN (1, FILE = 'ana/RHXZt23.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.24) then	   
      OPEN (1, FILE = 'ana/RHXZt24.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.25) then	  
      OPEN (1, FILE = 'ana/RHXZt25.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.26) then	   
      OPEN (1, FILE = 'ana/RHXZt26.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.27) then	   
      OPEN (1, FILE = 'ana/RHXZt27.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.28) then	  
      OPEN (1, FILE = 'ana/RHXZt28.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.29) then	   
      OPEN (1, FILE = 'ana/RHXZt29.dat',STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.30) then	    
      OPEN (1, FILE = 'ana/RHXZt30.dat',STATUS = 'UNKNOWN') 
	ENDIF
	
	IF(ntest.eq.1) then 
      OPEN (2, FILE = 'ana/RHAXZt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (2, FILE = 'ana/RHAXZt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (2, FILE = 'ana/RHAXZt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (2, FILE = 'ana/RHAXZt4.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.5) then
      OPEN (2, FILE = 'ana/RHAXZt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then
      OPEN (2, FILE = 'ana/RHAXZt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then
      OPEN (2, FILE = 'ana/RHAXZt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (2, FILE = 'ana/RHAXZt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (2, FILE = 'ana/RHAXZt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (2, FILE = 'ana/RHAXZt10.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.11) then	  
      OPEN (2, FILE = 'ana/RHAXZt11.dat',STATUS = 'UNKNOWN')   
	ELSE IF(ntest.eq.12) then	 
      OPEN (2, FILE = 'ana/RHAXZt12.dat',STATUS = 'UNKNOWN')   
	ELSE IF(ntest.eq.13) then	 
      OPEN (2, FILE = 'ana/RHAXZt13.dat',STATUS = 'UNKNOWN')   
	ELSE IF(ntest.eq.14) then	 
      OPEN (2, FILE = 'ana/RHAXZt14.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.15) then	  
      OPEN (2, FILE = 'ana/RHAXZt15.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.16) then	  
      OPEN (2, FILE = 'ana/RHAXZt16.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.17) then	   
      OPEN (2, FILE = 'ana/RHAXZt17.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.18) then	   
      OPEN (2, FILE = 'ana/RHAXZt18.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.19) then	   
      OPEN (2, FILE = 'ana/RHAXZt19.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.20) then	   
      OPEN (2, FILE = 'ana/RHAXZt20.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.21) then	   
      OPEN (2, FILE = 'ana/RHAXZt21.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.22) then	  
      OPEN (2, FILE = 'ana/RHAXZt22.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.23) then	   
      OPEN (2, FILE = 'ana/RHAXZt23.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.24) then	   
      OPEN (2, FILE = 'ana/RHAXZt24.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.25) then	   
      OPEN (2, FILE = 'ana/RHAXZt25.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.26) then	  
      OPEN (2, FILE = 'ana/RHAXZt26.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.27) then	  
      OPEN (2, FILE = 'ana/RHAXZt27.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.28) then	   
      OPEN (2, FILE = 'ana/RHAXZt28.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.29) then	   
      OPEN (2, FILE = 'ana/RHAXZt29.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.30) then	  
      OPEN (2, FILE = 'ana/RHAXZt30.dat',STATUS = 'UNKNOWN')  
	ENDIF
        do ixCLW=-20,20 !-20 ~ 20 fm
          do izCLW=-20,20 !-20 ~ 20 fm
	  RhoCont=DenConxz(ntest,ixCLW,izCLW)/0.16 !rho/rho0
	  RhoaCont=DenaConxz(ntest,ixCLW,izCLW)/0.16 !rhoa/rho0
	WRITE(1,'(G18.8,$)') RHOCont !xj
	WRITE(2,'(G18.8,$)') RHOaCont !xj
	    end do
	WRITE(1,*) !xj
	WRITE(2,*) !xj
	  end do
	CLOSE(1)
	CLOSE(2)
	enddo

      do ntest = 1,30
	IF(ntest.eq.1) then
      OPEN (1, FILE = 'ana/RHXYt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (1, FILE = 'ana/RHXYt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (1, FILE = 'ana/RHXYt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (1, FILE = 'ana/RHXYt4.dat', STATUS = 'UNKNOWN')
      ELSE IF(ntest.eq.5) then
	OPEN (1, FILE = 'ana/RHXYt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then      
	OPEN (1, FILE = 'ana/RHXYt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then      
	OPEN (1, FILE = 'ana/RHXYt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (1, FILE = 'ana/RHXYt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (1, FILE = 'ana/RHXYt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (1, FILE = 'ana/RHXYt10.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.11) then	   
      OPEN (1, FILE = 'ana/RHXYt11.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.12) then	   
      OPEN (1, FILE = 'ana/RHXYt12.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.13) then	  
      OPEN (1, FILE = 'ana/RHXYt13.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.14) then	   
      OPEN (1, FILE = 'ana/RHXYt14.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.15) then	   
      OPEN (1, FILE = 'ana/RHXYt15.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.16) then	   
      OPEN (1, FILE = 'ana/RHXYt16.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.17) then	   
      OPEN (1, FILE = 'ana/RHXYt17.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.18) then	  
      OPEN (1, FILE = 'ana/RHXYt18.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.19) then	   
      OPEN (1, FILE = 'ana/RHXYt19.dat',STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.20) then	    
      OPEN (1, FILE = 'ana/RHXYt20.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.21) then	   
      OPEN (1, FILE = 'ana/RHXYt21.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.22) then	  
      OPEN (1, FILE = 'ana/RHXYt22.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.23) then	   
      OPEN (1, FILE = 'ana/RHXYt23.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.24) then	   
      OPEN (1, FILE = 'ana/RHXYt24.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.25) then	  
      OPEN (1, FILE = 'ana/RHXYt25.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.26) then	   
      OPEN (1, FILE = 'ana/RHXYt26.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.27) then	   
      OPEN (1, FILE = 'ana/RHXYt27.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.28) then	  
      OPEN (1, FILE = 'ana/RHXYt28.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.29) then	   
      OPEN (1, FILE = 'ana/RHXYt29.dat',STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.30) then	    
      OPEN (1, FILE = 'ana/RHXYt30.dat',STATUS = 'UNKNOWN') 
	ENDIF
	
	IF(ntest.eq.1) then 
      OPEN (2, FILE = 'ana/RHAXYt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (2, FILE = 'ana/RHAXYt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (2, FILE = 'ana/RHAXYt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (2, FILE = 'ana/RHAXYt4.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.5) then
      OPEN (2, FILE = 'ana/RHAXYt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then
      OPEN (2, FILE = 'ana/RHAXYt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then
      OPEN (2, FILE = 'ana/RHAXYt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (2, FILE = 'ana/RHAXYt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (2, FILE = 'ana/RHAXYt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (2, FILE = 'ana/RHAXYt10.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.11) then	  
      OPEN (2, FILE = 'ana/RHAXYt11.dat',STATUS = 'UNKNOWN')   
	ELSE IF(ntest.eq.12) then	 
      OPEN (2, FILE = 'ana/RHAXYt12.dat',STATUS = 'UNKNOWN')   
	ELSE IF(ntest.eq.13) then	 
      OPEN (2, FILE = 'ana/RHAXYt13.dat',STATUS = 'UNKNOWN')   
	ELSE IF(ntest.eq.14) then	 
      OPEN (2, FILE = 'ana/RHAXYt14.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.15) then	  
      OPEN (2, FILE = 'ana/RHAXYt15.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.16) then	  
      OPEN (2, FILE = 'ana/RHAXYt16.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.17) then	   
      OPEN (2, FILE = 'ana/RHAXYt17.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.18) then	   
      OPEN (2, FILE = 'ana/RHAXYt18.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.19) then	   
      OPEN (2, FILE = 'ana/RHAXYt19.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.20) then	   
      OPEN (2, FILE = 'ana/RHAXYt20.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.21) then	   
      OPEN (2, FILE = 'ana/RHAXYt21.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.22) then	  
      OPEN (2, FILE = 'ana/RHAXYt22.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.23) then	   
      OPEN (2, FILE = 'ana/RHAXYt23.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.24) then	   
      OPEN (2, FILE = 'ana/RHAXYt24.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.25) then	   
      OPEN (2, FILE = 'ana/RHAXYt25.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.26) then	  
      OPEN (2, FILE = 'ana/RHAXYt26.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.27) then	  
      OPEN (2, FILE = 'ana/RHAXYt27.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.28) then	   
      OPEN (2, FILE = 'ana/RHAXYt28.dat',STATUS = 'UNKNOWN') 
	ELSE IF(ntest.eq.29) then	   
      OPEN (2, FILE = 'ana/RHAXYt29.dat',STATUS = 'UNKNOWN')  
	ELSE IF(ntest.eq.30) then	  
      OPEN (2, FILE = 'ana/RHAXYt30.dat',STATUS = 'UNKNOWN')  
	ENDIF
        do ixCLW=-20,20 !-20 ~ 20 fm
          do iyCLW=-20,20 !-20 ~ 20 fm
	  RhoCont=DenConxy(ntest,ixCLW,iyCLW)/0.16 !rho/rho0
	  RhoaCont=DenaConxy(ntest,ixCLW,iyCLW)/0.16 !rhoa/rho0
	WRITE(1,'(G18.8,$)') RHOCont !xj
	WRITE(2,'(G18.8,$)') RHOaCont !xj
	    end do
	WRITE(1,*) !xj
	WRITE(2,*) !xj
	  end do
	CLOSE(1)
	CLOSE(2)
	enddo

      do ntest = 1,10
	IF(ntest.eq.1) then
      OPEN (1, FILE = 'ana/RHXYnjlt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (1, FILE = 'ana/RHXYnjlt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (1, FILE = 'ana/RHXYnjlt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (1, FILE = 'ana/RHXYnjlt4.dat', STATUS = 'UNKNOWN')
      ELSE IF(ntest.eq.5) then
	OPEN (1, FILE = 'ana/RHXYnjlt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then      
	OPEN (1, FILE = 'ana/RHXYnjlt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then      
	OPEN (1, FILE = 'ana/RHXYnjlt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (1, FILE = 'ana/RHXYnjlt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (1, FILE = 'ana/RHXYnjlt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (1, FILE = 'ana/RHXYnjlt10.dat',STATUS = 'UNKNOWN') 	   
	ENDIF
	
	IF(ntest.eq.1) then 
      OPEN (2, FILE = 'ana/RHAXYnjlt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (2, FILE = 'ana/RHAXYnjlt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (2, FILE = 'ana/RHAXYnjlt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (2, FILE = 'ana/RHAXYnjlt4.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.5) then
      OPEN (2, FILE = 'ana/RHAXYnjlt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then
      OPEN (2, FILE = 'ana/RHAXYnjlt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then
      OPEN (2, FILE = 'ana/RHAXYnjlt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (2, FILE = 'ana/RHAXYnjlt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (2, FILE = 'ana/RHAXYnjlt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (2, FILE = 'ana/RHAXYnjlt10.dat',STATUS = 'UNKNOWN')  
	ENDIF
        do iyCLW=1,ny !-12.5 ~ 12.5 fm
          do ixCLW=1,nx !-11 ~ 11 fm
	  RhoCont=DenConnjlxy(ntest,ixCLW,iyCLW)
	  RhoaCont=DenaConnjlxy(ntest,ixCLW,iyCLW)
	WRITE(1,'(G18.8,$)') RHOCont !xj
	WRITE(2,'(G18.8,$)') RHOaCont !xj
	    end do
	WRITE(1,*) !xj
	WRITE(2,*) !xj
	  end do
	CLOSE(1)
	CLOSE(2)
	enddo

      do ntest = 1,10
	IF(ntest.eq.1) then
      OPEN (1, FILE = 'ana/RHXZnjlt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (1, FILE = 'ana/RHXZnjlt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (1, FILE = 'ana/RHXZnjlt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (1, FILE = 'ana/RHXZnjlt4.dat', STATUS = 'UNKNOWN')
      ELSE IF(ntest.eq.5) then
	OPEN (1, FILE = 'ana/RHXZnjlt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then      
	OPEN (1, FILE = 'ana/RHXZnjlt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then      
	OPEN (1, FILE = 'ana/RHXZnjlt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (1, FILE = 'ana/RHXZnjlt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (1, FILE = 'ana/RHXZnjlt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (1, FILE = 'ana/RHXZnjlt10.dat',STATUS = 'UNKNOWN') 	   
	ENDIF
	
	IF(ntest.eq.1) then 
      OPEN (2, FILE = 'ana/RHAXZnjlt1.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.2) then
      OPEN (2, FILE = 'ana/RHAXZnjlt2.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.3) then
      OPEN (2, FILE = 'ana/RHAXZnjlt3.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.4) then
      OPEN (2, FILE = 'ana/RHAXZnjlt4.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.5) then
      OPEN (2, FILE = 'ana/RHAXZnjlt5.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.6) then
      OPEN (2, FILE = 'ana/RHAXZnjlt6.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.7) then
      OPEN (2, FILE = 'ana/RHAXZnjlt7.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.8) then
      OPEN (2, FILE = 'ana/RHAXZnjlt8.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.9) then
      OPEN (2, FILE = 'ana/RHAXZnjlt9.dat', STATUS = 'UNKNOWN')
	ELSE IF(ntest.eq.10) then
      OPEN (2, FILE = 'ana/RHAXZnjlt10.dat',STATUS = 'UNKNOWN')  
	ENDIF
        do ixCLW=1,nx !-11 ~ 11 fm
          do izCLW=1,nz !-5 ~ 5 fm
	  RhoCont=DenConnjlxz(ntest,ixCLW,izCLW)
	  RhoaCont=DenaConnjlxz(ntest,ixCLW,izCLW)
	WRITE(1,'(G18.8,$)') RHOCont !xj
	WRITE(2,'(G18.8,$)') RHOaCont !xj
	    end do
	WRITE(1,*) !xj
	WRITE(2,*) !xj
	  end do
	CLOSE(1)
	CLOSE(2)
	enddo

	OPEN(1,FILE='his_dp.txt')
	OPEN(2,FILE='his_dr.txt')

	DO I = 1,20
	WRITE(1,*) (I-0.5)*DPHIS,HISDP(I)/NUM/NEVNT
     &,HISMDP(I)/NUM/NEVNT,HISBDP(I)/NUM/NEVNT
	WRITE(2,*) (I-0.5)*DRHIS,HISDR(I)/NUM/NEVNT
     &,HISMDR(I)/NUM/NEVNT,HISBDR(I)/NUM/NEVNT
	ENDDO

	CLOSE(1)
	CLOSE(2)

       STOP
       END
c     FYI: taken file unit numbers are 12-14, 16-88, 91-93; 
c     so free file unit numbers are 11,15,89,97-99.
