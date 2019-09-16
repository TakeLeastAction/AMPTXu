	subroutine getnjl(efrm,ixj)
      PARAMETER (MAXPTN=400001,MAXR=100)
	implicit real*8 (a-h,o-z)
	common/particle/ par(9,20000,100), q(9,20000,100), q0(3,20000,100)
     &, ihd(20000,100)
      common/totoal/ itotal(1000), ntotal(1000) 
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
      COMMON /ilistnjl/ LSTRGnjl(MAXPTN,MAXR), LPARTnjl(MAXPTN,MAXR)
      COMMON /PARA1/ MUL
	COMMON /NJLMUL/ NJLMUL(MAXR)
      real  efrm
! number of partitions (max=10)
	npart = 1 
C      do npara=1, nevent !<----------------------------------------------------- new event
	npara = ixj
c       iparthd(npara)=0 !<======this is the array recording the # of hadronized partons in each event.
       !<=========Initially, all the partons aren't hadronized, so set all the components to zero.
       k=0
       do npartition=1, npart
c          read(17,*) iev, i2, ipart(npara), bim, i5, i6, i7, i8 
c		  read(17,*) nsg(npara) 

c          do i=1, ipart(npara)
           do i=1,mul
c             read(17,*) id,lstring,lpart,xx,yy,zz,qx,qy,qz,xmass,time 
             k=k+1 !k=i
                 par(1,k,npara)=GX0(I)
                 par(2,k,npara)=GY0(I)
!xj17	         par(3,k,npara)=GZ0(I)
                 dxj = 14.*0.938/EFRM !xj17
                 par(3,k,npara) = (RANART(NSEED)-0.5)*2.*dxj !xj17
                 par(4,k,npara)=FT0(I)
                 par(5,k,npara)=PX0(I)
                 par(6,k,npara)=PY0(I)
                 par(7,k,npara)=PZ0(I)
                 par(8,k,npara)=XMASS0(I)
                 par(9,k,npara)=float(ITYP0(I))

          if(par(9,k,npara).gt.-2.5.and.par(9,k,npara).lt.-1.5)then
           xnum_u_bar=xnum_u_bar+1.
          elseif(par(9,k,npara).gt.-1.5.and.par(9,k,npara).lt.-0.5) then
           x num_d_bar=xnum_d_bar+1.
          elseif(par(9,k,npara).gt.0.5.and.par(9,k,npara).lt.1.5) then
           xnum_d=xnum_d+1.
          elseif(par(9,k,npara).gt.1.5.and.par(9,k,npara).lt.2.5) then
           xnum_u=xnum_u+1.
          endif
c			 ls(k,npara) = LSTRG0(I)
c			 lp(k,npara) = LPART0(I) 
	!propagate to formation time
	E = SQRT(par(5,k,npara)**2+par(6,k,npara)**2+par(7,k,npara)**2
     &+par(8,k,npara)**2)
	par(1,k,npara) = par(1,k,npara) + par(5,k,npara)/E* par(4,k,npara)
	par(2,k,npara) = par(2,k,npara) + par(6,k,npara)/E* par(4,k,npara)
	par(3,k,npara) = par(3,k,npara) + par(7,k,npara)/E* par(4,k,npara)

	LSTRGnjl(i,ixj) = LSTRG0(i)
	LPARTnjl(i,ixj) = LPART0(i)

	      enddo ! for particles in one partition
        delta_parton=((xnum_d-xnum_d_bar)-(xnum_u-xnum_u_bar))/
     &               ((xnum_d-xnum_d_bar)+(xnum_u-xnum_u_bar))
        write(11,163) delta_parton,xnum_d,xnum_d_bar,xnum_u,xnum_u_bar
163     format(1x,5E15.6)
       enddo ! for particles in one event
!    itotal(npara)=ntotal
       ntotal(npara)=k
c	enddo ! for events
	NJLMUL(npara) = mul
	RETURN
	END
		
	subroutine NJL
      PARAMETER (MAXPTN=400001,MAXR=100)
	implicit real*8 (a-h,o-z)
	common/particle/ par(9,20000,100), q(9,20000,100), q0(3,20000,100)
     &, ihd(20000,100)
      common/scalar/ Vu(100,100,40), Vd(100,100,40), Vs(100,100,40)! scalar potential
	common/dennjl/ denqnjl(100,100,40),denqbarnjl(100,100,40) !quark density
      common/vector/ Ut(100,100,40), Ux(100,100,40), Uy(100,100,40)
     &, Uz(100,100,40) ! vector poential
      common/isovector/ Utiv(100,100,40), Uxiv(100,100,40)
     &, Uyiv(100,100,40), Uziv(100,100,40) ! vector-isovector poential
      common/grid1/nx, ny, nz
      common/grid2/dx, dy, dz

      common/totoal/ itotal(1000), ntotal(1000) 

      dimension pxs(100,100,40), pys(100,100,40), pzs(100,100,40)
      dimension vxs(100,100,40), vys(100,100,40), vzs(100,100,40)
      dimension pxvu(100,100,40), pyvu(100,100,40), pzvu(100,100,40)
     &, chau(100,100,40)
      dimension pxvd(100,100,40), pyvd(100,100,40), pzvd(100,100,40)
     &, chad(100,100,40)
      dimension pxv(100,100,40), pyv(100,100,40), pzv(100,100,40)
     &, chab(100,100,40)
      dimension Edd(100,100,40), Euu(100,100,40), Ess(100,100,40)
      dimension dd0(100,100,40), uu0(100,100,40), ss0(100,100,40)

c	dimension ipart(100),nsg(100),ls(20000,100),lp(20000,100)
c     &,iparthd(100) 

	dimension imp(20000,100) !keep initial order

      COMMON/RNDF77/NSEED
	COMMON /NJLMUL/ NJLMUL(MAXR)
      COMMON /precnjl/GXnjl(MAXPTN,MAXR),GYnjl(MAXPTN,MAXR)
     &,GZnjl(MAXPTN,MAXR),FTnjl(MAXPTN,MAXR)
     &,PXnjl(MAXPTN,MAXR), PYnjl(MAXPTN,MAXR), PZnjl(MAXPTN,MAXR)
     &,Enjl(MAXPTN,MAXR),XMASSnjl(MAXPTN,MAXR), ITYPnjl(MAXPTN,MAXR)
      COMMON /ilistnjl/ LSTRGnjl(MAXPTN,MAXR), LPARTnjl(MAXPTN,MAXR)
      COMMON  /RUN/     NUM

	COMMON /EVENT/ NEVNT !xj dencon
      COMMON /AREVT/ IAEVT, IARUN, MISS !xj dencon
	COMMON/DenConnjl/ DenConnjlxz(10,100,40)
     &,DenaConnjlxz(10,100,40),DenConnjlxy(10,100,100)
     &,DenaConnjlxy(10,100,100) !xj
		
c	open (unit=99, file='den_NJL.txt',status='unknown')
c      open (unit=14, file='zpc_NJL.dat', status='unknown')
	! vector coupling strength (gv/Gs)
! 1./3. for Fierz transformation
! 2.2/3. for Weise's paper
	ratio = 2.2/3. 

! number of events (max=100)
	nevent = NUM

! number of partitions (max=10)
	npart = 1 

! parton cross section (mb)
	cs = 1. 
        sigma=cs/10./npart ! cross section is 10 mb = 1 fm^2

! seed # for random numbers
c	nseed = irun*10+1

	iglobal = 1 !1: global hadronization, 0: local hadronization 
	thd = 1.6 !2. !time after which hadronization is allowed
	tend = 2. !8. !end time

c	do ixj = 1,1 
c	write(99,*) ixj
c	write(*,*) ixj 

	pi=3.14159265
	nc=3

	nt=60 ! iteration number for quark mass

        tau=0.05 ! initial time

        nx=44 
	ny=50 
	nz=40 
        dx=5.d-1 !xjdencon: x=-11 ~ 11
	dy=5.d-1 !xjdencon: y=-12.5 ~ 12.5
	dz=2.5d-1 !xjdencon: z=-5 ~ 5
	dV=dx*dy*dz*5.07**3.

        nxh=int(nx/2.+1.d-8)
	nyh=int(ny/2.+1.d-8)
	nzh=int(nz/2.+1.d-8)

! case Weise
       xm0=3.6d-3
       xms0=87.d-3
       xlam=0.75 ! unit is GeV
       gs=3.64/xlam**2  
       gg=1.82/xlam**2. ! in Weise's original paper, gg=3.6
       gk=8.9/xlam**5.
       gv=ratio*gs
       gis = 1.d-5 !gg !scalar-isovector 
       giv = 2.*gs !vector-isovector 

! this is for vacuum
	xmu=xm0
	xmd=xm0
	xms=xms0
      xlam2=xlam**2.
      do i=1, 100
       alpha=log(xlam/xmu +sqrt((xlam/xmu)**2.+1.))
       xmu2=xmu**2.
	   uu=-xmu*((2.*nc)/(2.*pi)**3.)*pi*xmu2*(sinh(2.*alpha)-2.*alpha)
       alpha=log(xlam/xmd +sqrt((xlam/xmd)**2.+1.))
       xmd2=xmd**2.
	   dd=-xmd*((2.*nc)/(2.*pi)**3.)*pi*xmd2*(sinh(2.*alpha)-2.*alpha)

       alpha=log(xlam/xms +sqrt((xlam/xms)**2.+1.))
       xms2=xms**2.
	   ss=-xms*((2.*nc)/(2.*pi)**3.)*pi*xms2*(sinh(2.*alpha)-2.*alpha)

	   xmu=xm0 -4.*gg*uu +2.*gk*dd*ss - 2.*gis*(uu-dd)
	   xmd=xm0 -4.*gg*dd +2.*gk*uu*ss + 2.*gis*(uu-dd)
	   xms=xms0-4.*gg*ss +2.*gk*uu*dd 
      enddo
	uuv=uu ! uu condensate in vacuum
      ddv=dd
      ssv=ss ! ss condensate in vacuum
    !<========== record the constitute quark mass in vacuum
      xmuv=xmu
	xmdv=xmd
      xmsv=xms
    !<==========
c      alpha=log(xlam/xm0 +sqrt((xlam/xm0)**2.+1.))
c      ubottom=((2.*nc)/(2.*pi)**3.)*pi*xm0*xm0*(sinh(2.*alpha)-2.*alpha) ! for u, d

c      alpha=log(xlam/xms0 +sqrt((xlam/xms0)**2.+1.))
c      sbottom=((2.*nc)/(2.*pi)**3.)*pi*xms0*xms0
c     &*(sinh(2.*alpha)-2.*alpha) ! for s


! this is for vacuum energy density
      deg=2.*nc
      xmu2=xmu**2.
	xmd2=xmd**2.
      xms2=xms**2.
	Ed0=sqrt(xmu2 +xlam2)
      Eu0=sqrt(xmd2 +xlam2)
	Es0=sqrt(xms2 +xlam2)
	veu1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmu2)*Eu0 
     &-xmu2*xmu2*log((Eu0+xlam)/xmu))
	ved1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmd2)*Ed0 
     &-xmd2*xmd2*log((Ed0+xlam)/xmd))
	ves1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xms2)*Es0 
     &-xms2*xms2*log((Es0+xlam)/xms))

      ve2=2.*gg*(uuv*uuv +ddv*ddv +ssv*ssv) -4.*gk*uuv*ddv*ssv 
     &+gis*(uuv-ddv)**2

      vedn=(veu1 +ved1 +ves1) +ve2 

c      call SRAND(nseed) 
c      do npara=1, nevent !<----------------------------------------------------- new event
c       iparthd(npara)=0 !<======this is the array recording the # of hadronized partons in each event.
       !<=========Initially, all the partons aren't hadronized, so set all the components to zero.
c       k=0
c       do npartition=1, npart
c          read(17,*) iev, i2, ipart(npara), bim, i5, i6, i7, i8 
c		  read(17,*) nsg(npara) 

c          do i=1, ipart(npara)
c		 do i=1,mul
c             read(17,*) id,lstring,lpart,xx,yy,zz,qx,qy,qz,xmass,time 
c             k=k+1 !k=i
c	         par(1,k,npara)=GX0(I)
c	         par(2,k,npara)=GY0(I)
c	         par(3,k,npara)=GZ0(I)
c               par(4,k,npara)=FT0(I)
c	         par(5,k,npara)=PX0(I)
c	         par(6,k,npara)=PY0(I)
c	         par(7,k,npara)=PZ0(I)
c               par(8,k,npara)=XMASS0(I)
c               par(9,k,npara)=float(ITYP0(I))
c			 ls(k,npara) = LSTRG0(I)
c			 lp(k,npara) = LPART0(I) 
	!propagate to formation time
c	E = SQRT(par(5,k,npara)**2+par(6,k,npara)**2+par(7,k,npara)**2
c     &+par(8,k,npara)**2)
c	par(1,k,npara) = par(1,k,npara) + par(5,k,npara)/E* par(4,k,npara)
c	par(2,k,npara) = par(2,k,npara) + par(6,k,npara)/E* par(4,k,npara)
c	par(3,k,npara) = par(3,k,npara) + par(7,k,npara)/E* par(4,k,npara)
c	      enddo ! for particles in one partition
c       enddo ! for particles in one event
!    itotal(npara)=ntotal
c       ntotal(npara)=k
c	enddo ! for events

!goto 99
!************************************************************************ particle generation

      t=tau
      do i=1, nx
	do j=1, ny
	do k=1, nz
       Vu(i,j,k)=xmu
       Vd(i,j,k)=xmd
       Vs(i,j,k)=xms
	   Ut(i,j,k)=0.d0 !<--- phi
	   Ux(i,j,k)=0.d0 !<--- Ax
	   Uy(i,j,k)=0.d0 !<--- Ay
	   Uz(i,j,k)=0.d0 !<--- Az
	enddo
	enddo
	enddo

      do i=1, nevent
	   itotal(i)=0
	enddo

      itotal2=0
!   dt=0.02*exp(t/1.9)
      dt=2.d-2 

	!xjbegin
	if(IAEVT.eq.1) then
	      DO I=1,10
	   DO II=1,100
	     DO III=1,40
         DenConnjlxz(I,II,III)=0.
         DenaConnjlxz(I,II,III)=0.
	     END DO
	   END DO
		   DO II=1,100
	     DO III=1,100
         DenConnjlxy(I,II,III)=0.
         DenaConnjlxy(I,II,III)=0.
	     END DO
	   END DO
	END DO
	endif
	!xjend
	ntest = 0
 30   continue !<----------------------------------------------------- new time step

      do npara=1, nevent
       mp=itotal(npara)
       np=0
	   do i=1, ntotal(npara)
          if ((par(4,i,npara) .gt. t-dt) 
     &.and. (par(4,i,npara) .le. t)) then
             mp=mp+1
			 np=np+1
       E = SQRT(par(5,i,npara)**2+par(6,i,npara)**2+par(7,i,npara)**2
     &+par(8,i,npara)**2)
	q(1,mp,npara) = par(1,i,npara) + par(5,i,npara)/E*(t-par(4,i,npara))        !propogate to now time t
	q(2,mp,npara) = par(2,i,npara) + par(6,i,npara)/E*(t-par(4,i,npara))
	q(3,mp,npara) = par(3,i,npara) + par(7,i,npara)/E*(t-par(4,i,npara))
               q(4,mp,npara)=par(4,i,npara)
		     q(5,mp,npara)=par(5,i,npara)
		     q(6,mp,npara)=par(6,i,npara)
		     q(7,mp,npara)=par(7,i,npara)
               q(8,mp,npara)=par(8,i,npara)
		     q(9,mp,npara)=par(9,i,npara)
		     ihd(mp,npara)=0  !<======It's an array denoting whether a parton is hadronized or not. ihd=1 means hadronized.
			 imp(mp,npara)=i !keep initial order
	      endif
	   enddo
       itotal(npara)=itotal(npara)+np
       itotal2=itotal2+np
      enddo
	dt=0.02*exp(t/1.9)
      if (itotal2 .lt. 1) then
	   t=t+dt
         ntest = ntest + 1
	   goto 30
	endif

! goto 77
!********************************************************************** potential
	ntest = ntest + 1

	do i=1, nx
	   do j=1, ny
	      do k=1, nz
             pxs(i,j,k)=0.d0
             pys(i,j,k)=0.d0
             pzs(i,j,k)=0.d0

             pxv(i,j,k)=0.d0
             pyv(i,j,k)=0.d0
             pzv(i,j,k)=0.d0

             pxvu(i,j,k)=0.d0 !han
             pyvu(i,j,k)=0.d0 !han
             pzvu(i,j,k)=0.d0 !han
	       pxvd(i,j,k)=0.d0 !han
             pyvd(i,j,k)=0.d0 !han
             pzvd(i,j,k)=0.d0 !han

		   chau(i,j,k)=0.d0 !han
		   chad(i,j,k)=0.d0 !han

		   chab(i,j,k)=0.d0
	denqnjl(i,j,k) = 0.
	denqbarnjl(i,j,k) = 0.
          enddo
       enddo
	enddo

      do npara=1, nevent
	   do i=1, itotal(npara)

        if(ihd(i,npara).eq.1) goto 2015 !xjden

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 2015
          if ((my .le. 1) .or. (my .ge. ny)) goto 2015
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 2015

	if(q(9,i,npara).gt.0) then
      denqnjl(mx,my,mz) = denqnjl(mx,my,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/3.  	 
      denqnjl(mx+1,my,mz) = denqnjl(mx+1,my,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqnjl(mx-1,my,mz) = denqnjl(mx-1,my,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqnjl(mx,my+1,mz) = denqnjl(mx,my+1,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqnjl(mx,my-1,mz) = denqnjl(mx,my-1,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqnjl(mx,my,mz+1) = denqnjl(mx,my,mz+1) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqnjl(mx,my,mz-1) = denqnjl(mx,my,mz-1) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
	endif

	if(q(9,i,npara).lt.0) then
      denqbarnjl(mx,my,mz) = denqbarnjl(mx,my,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/3.   
      denqbarnjl(mx+1,my,mz) = denqbarnjl(mx+1,my,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqbarnjl(mx-1,my,mz) = denqbarnjl(mx-1,my,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqbarnjl(mx,my+1,mz) = denqbarnjl(mx,my+1,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqbarnjl(mx,my-1,mz) = denqbarnjl(mx,my-1,mz) 
     &+ 1./(dx*dy*dz)/npart/nevent/9. 
      denqbarnjl(mx,my,mz+1) = denqbarnjl(mx,my,mz+1) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.  
      denqbarnjl(mx,my,mz-1) = denqbarnjl(mx,my,mz-1) 
     &+ 1./(dx*dy*dz)/npart/nevent/9.
      endif
2015   continue

       enddo
	enddo

	!xjbegin
	If(ntest.eq.1.or.ntest.eq.10.or.ntest.eq.19.or.ntest.eq.28
     &.or.ntest.eq.37.or.ntest.eq.46.or.ntest.eq.55.or.ntest.eq.64
     &.or.ntest.eq.73.or.ntest.eq.82) then !xj
	idencont = (ntest-1)/9 + 1
	write(*,*) 'dencon:',idencont

        do ix0=1,nx
        do iy0=1,ny
        do iz0=1,nz
        if(denqnjl(ix0,iy0,iz0).eq.0..and.denqbarnjl(ix0,iy0,iz0).eq.0.)
     &go to 8888

	if(iy0.ge.ny/2-1.and.iy0.le.ny/2+1) then
	DenConnjlxz(idencont,ix0,iz0)=DenConnjlxz(idencont,ix0,iz0)
     &+denqnjl(ix0,iy0,iz0)/3./NEVNT
	DenaConnjlxz(idencont,ix0,iz0)=DenaConnjlxz(idencont,ix0,iz0)
     &+denqbarnjl(ix0,iy0,iz0)/3./NEVNT
	end if

	if(iz0.ge.nz/2-2.and.iz0.le.nz/2+2) then
	DenConnjlxy(idencont,ix0,iy0)=DenConnjlxy(idencont,ix0,iy0)
     &+denqnjl(ix0,iy0,iz0)/5./NEVNT
	DenaConnjlxy(idencont,ix0,iy0)=DenaConnjlxy(idencont,ix0,iy0)
     &+denqbarnjl(ix0,iy0,iz0)/5./NEVNT
	end if

8888	CONTINUE
	enddo
	enddo
	enddo

	end if
	!xjend

! step 0 : sum momentum in each fluid cell
      do npara=1, nevent
	   do i=1, itotal(npara)

        if(ihd(i,npara).eq.1) goto 32 !xjfeng

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 32
          if ((my .le. 1) .or. (my .ge. ny)) goto 32
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 32

		     pxs(mx,my,mz)= pxs(mx,my,mz) +q(5,i,npara)/3.
		     pxs(mx+1,my,mz)= pxs(mx+1,my,mz) +q(5,i,npara)/9.
		     pxs(mx,my+1,mz)= pxs(mx,my+1,mz) +q(5,i,npara)/9.
		     pxs(mx,my,mz+1)= pxs(mx,my,mz+1) +q(5,i,npara)/9.
		     pxs(mx-1,my,mz)= pxs(mx-1,my,mz) +q(5,i,npara)/9.
		     pxs(mx,my-1,mz)= pxs(mx,my-1,mz) +q(5,i,npara)/9.
		     pxs(mx,my,mz-1)= pxs(mx,my,mz-1) +q(5,i,npara)/9.

             pys(mx,my,mz)= pys(mx,my,mz) +q(6,i,npara)/3.
             pys(mx+1,my,mz)= pys(mx+1,my,mz) +q(6,i,npara)/9.
             pys(mx,my+1,mz)= pys(mx,my+1,mz) +q(6,i,npara)/9.
             pys(mx,my,mz+1)= pys(mx,my,mz+1) +q(6,i,npara)/9.
             pys(mx-1,my,mz)= pys(mx-1,my,mz) +q(6,i,npara)/9.
             pys(mx,my-1,mz)= pys(mx,my-1,mz) +q(6,i,npara)/9.
             pys(mx,my,mz-1)= pys(mx,my,mz-1) +q(6,i,npara)/9.

             pzs(mx,my,mz)= pzs(mx,my,mz) +q(7,i,npara)/3.
             pzs(mx+1,my,mz)= pzs(mx+1,my,mz) +q(7,i,npara)/9.
             pzs(mx,my+1,mz)= pzs(mx,my+1,mz) +q(7,i,npara)/9.
             pzs(mx,my,mz+1)= pzs(mx,my,mz+1) +q(7,i,npara)/9.
             pzs(mx-1,my,mz)= pzs(mx-1,my,mz) +q(7,i,npara)/9.
             pzs(mx,my-1,mz)= pzs(mx,my-1,mz) +q(7,i,npara)/9.
             pzs(mx,my,mz-1)= pzs(mx,my,mz-1) +q(7,i,npara)/9.

 32   continue

       enddo
	enddo



!****************************************
! this is for scalar mean-fields

      it=0
 31   continue


! step 1 : initializing
	do i=1, nx
	   do j=1, ny
	      do k=1, nz

			 Edd(i,j,k)=0.d0
			 Euu(i,j,k)=0.d0
			 Ess(i,j,k)=0.d0

             dd0(i,j,k)=0.d0
             uu0(i,j,k)=0.d0
             ss0(i,j,k)=0.d0

          enddo
       enddo
	enddo



! step 2 : sum energy in each fluid cell
      do npara=1, nevent
	   do i=1, itotal(npara)

          if(ihd(i,npara).eq.1) goto 33 !<===== if a parton is hadronized, it won't contribute to the mean field.
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 33
          if ((my .le. 1) .or. (my .ge. ny)) goto 33
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 33

          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.

          if (abs(q(9,i,npara)) .lt. 1.5) then		       ! for d quarks
             Edd(mx,my,mz)= Edd(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)/3.
          Edd(mx+1,my,mz)= Edd(mx+1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my+1,mz)= Edd(mx,my+1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my,mz+1)= Edd(mx,my,mz+1) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx-1,my,mz)= Edd(mx-1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my-1,mz)= Edd(mx,my-1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my,mz-1)= Edd(mx,my,mz-1) +sqrt(q(8,i,npara)**2.+q2)/9.
          else if (abs(q(9,i,npara)) .lt. 2.5) then		       ! for u quarks
             Euu(mx,my,mz)= Euu(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)/3.
          Euu(mx+1,my,mz)= Euu(mx+1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my+1,mz)= Euu(mx,my+1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my,mz+1)= Euu(mx,my,mz+1) +sqrt(q(8,i,npara)**2.+q2)/9.           
		Euu(mx-1,my,mz)= Euu(mx-1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my-1,mz)= Euu(mx,my-1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my,mz-1)= Euu(mx,my,mz-1) +sqrt(q(8,i,npara)**2.+q2)/9.
		  else                                             ! for s quarks
             Ess(mx,my,mz)= Ess(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)/3.
          Ess(mx+1,my,mz)= Ess(mx+1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Ess(mx,my+1,mz)= Ess(mx,my+1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Ess(mx,my,mz+1)= Ess(mx,my,mz+1) +sqrt(q(8,i,npara)**2.+q2)/9.
          Ess(mx-1,my,mz)= Ess(mx-1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Ess(mx,my-1,mz)= Ess(mx,my-1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Ess(mx,my,mz-1)= Ess(mx,my,mz-1) +sqrt(q(8,i,npara)**2.+q2)/9.
          endif

 33   continue

       enddo
	enddo



! step 3 : find fluid velocity
	do i=1, nx
	   do j=1, ny
	      do k=1, nz

             if (Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k) .gt. 1.d-9) then
              vxs(i,j,k)= pxs(i,j,k)/(Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k))
              vys(i,j,k)= pys(i,j,k)/(Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k))
              vzs(i,j,k)= pzs(i,j,k)/(Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k))
             endif

          enddo
       enddo
	enddo



! step 4 : move to CM frame
      do npara=1, nevent
	   do i=1, itotal(npara)

        if(ihd(i,npara).eq.1) goto 34 !<=====if a parton is hadronized, it won't contribute to the mean field.

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 34
          if ((my .le. 1) .or. (my .ge. ny)) goto 34
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 34

          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.
          E1=sqrt(q(8,i,npara)**2.+q2)

             bx=vxs(mx,my,mz)
		     by=vys(mx,my,mz)
		     bz=vzs(mx,my,mz)

          call xLorentz(E1,q(5,i,npara),q(6,i,npara),q(7,i,npara)
     &, bx,by,bz, E2,px,py,pz)

          q2= px*px +py*py +pz*pz
		  if (q2 .gt. xlam2) goto 34

          if (abs(q(9,i,npara)) .lt. 1.5) then
             dd0(mx,my,mz)= dd0(mx,my,mz) +1./E1/dV/3.
             dd0(mx+1,my,mz)= dd0(mx+1,my,mz) +1./E1/dV/9.
             dd0(mx,my+1,mz)= dd0(mx,my+1,mz) +1./E1/dV/9.
             dd0(mx,my,mz+1)= dd0(mx,my,mz+1) +1./E1/dV/9.
             dd0(mx-1,my,mz)= dd0(mx-1,my,mz) +1./E1/dV/9.
             dd0(mx,my-1,mz)= dd0(mx,my-1,mz) +1./E1/dV/9.
             dd0(mx,my,mz-1)= dd0(mx,my,mz-1) +1./E1/dV/9.
          else if (abs(q(9,i,npara)) .lt. 2.5) then
             uu0(mx,my,mz)= uu0(mx,my,mz) +1./E1/dV/3.
             uu0(mx+1,my,mz)= uu0(mx+1,my,mz) +1./E1/dV/9.
             uu0(mx,my+1,mz)= uu0(mx,my+1,mz) +1./E1/dV/9.
             uu0(mx,my,mz+1)= uu0(mx,my,mz+1) +1./E1/dV/9.
             uu0(mx-1,my,mz)= uu0(mx-1,my,mz) +1./E1/dV/9.
             uu0(mx,my-1,mz)= uu0(mx,my-1,mz) +1./E1/dV/9.
             uu0(mx,my,mz-1)= uu0(mx,my,mz-1) +1./E1/dV/9.
          else
             ss0(mx,my,mz)= ss0(mx,my,mz) +1./E1/dV/3.
             ss0(mx+1,my,mz)= ss0(mx+1,my,mz) +1./E1/dV/9.
             ss0(mx,my+1,mz)= ss0(mx,my+1,mz) +1./E1/dV/9.
             ss0(mx,my,mz+1)= ss0(mx,my,mz+1) +1./E1/dV/9.
             ss0(mx-1,my,mz)= ss0(mx-1,my,mz) +1./E1/dV/9.
             ss0(mx,my-1,mz)= ss0(mx,my-1,mz) +1./E1/dV/9.
             ss0(mx,my,mz-1)= ss0(mx,my,mz-1) +1./E1/dV/9.
		  endif

 34   continue

       enddo
	enddo


	xmuucen = 0.d0 
	xmudcen = 0.d0 
	xmuscen = 0.d0 
      totalS=0.
	centerS=0.
! step 5 : determine condensate and scalar mean-fields
	do i=1, nx
	   do j=1, ny
	      do k=1, nz

             if (Edd(i,j,k) .gt. 1.d-9) then
			  alpha=log(xlam/Vd(i,j,k) +sqrt((xlam/Vd(i,j,k))**2.+1.))
                dmass2=Vd(i,j,k)**2.
			    constI=((2.*nc)/(2.*pi)**3.)*pi*dmass2*(sinh(2.*alpha)
     &-2.*alpha)
		        dd= -Vd(i,j,k) *(constI -dd0(i,j,k)/npart/nevent)
                if (dd .gt. 0.d0) dd=0.d0
			 else
			    dd=ddv
			 endif

             if (Euu(i,j,k) .gt. 1.d-9) then
			  alpha=log(xlam/Vu(i,j,k) +sqrt((xlam/Vu(i,j,k))**2.+1.))
                umass2=Vu(i,j,k)**2.
			    constI=((2.*nc)/(2.*pi)**3.)*pi*umass2*(sinh(2.*alpha)
     &-2.*alpha)
		        uu= -Vu(i,j,k) *(constI -uu0(i,j,k)/npart/nevent)
                if (uu .gt. 0.d0) uu=0.d0
			 else
			    uu=uuv
			 endif

             if (Ess(i,j,k) .gt. 1.d-9) then
			  alpha=log(xlam/Vs(i,j,k) +sqrt((xlam/Vs(i,j,k))**2.+1.))
                smass2=Vs(i,j,k)**2.
			    constI=((2.*nc)/(2.*pi)**3.)*pi*smass2*(sinh(2.*alpha)
     &-2.*alpha)
		        ss= -Vs(i,j,k) *(constI -ss0(i,j,k)/npart/nevent)
                if (ss .gt. 0.d0) ss=0.d0
			 else
			    ss=ssv
			 endif

	         Vd(i,j,k)= xm0 -4.*gg*dd +2.*gk*uu*ss - 2.*gis*(uu-dd)
	         Vu(i,j,k)= xm0 -4.*gg*uu +2.*gk*dd*ss + 2.*gis*(uu-dd)
	         Vs(i,j,k)= xms0-4.*gg*ss +2.*gk*uu*dd

         if ((i .ge. nxh-1) .and. (i .le. nxh+2)) then
	      if ((j .ge. nyh-1) .and. (j .le. nyh+2)) then
	      if ((k .ge. nzh) .and. (k .le. nzh+1)) then
			 xmudcen = xmudcen + Vd(i,j,k) !d quark mass in central cells
			 xmuucen = xmuucen + Vu(i,j,k) !u quark mass in central cells
			 xmuscen = xmuscen + Vs(i,j,k) !s quark mass in central cells
	      endif
	      endif
	      endif

! this is to calculate energy density from scalar mean fields

             deg=2.*nc
			 xmd2=Vd(i,j,k)**2.
			 xmu2=Vu(i,j,k)**2.
			 xms2=Vs(i,j,k)**2.
	         Ed0=sqrt(xmd2 +xlam2)
	         Eu0=sqrt(xmu2 +xlam2)
             Es0=sqrt(xms2 +xlam2)
	         ved1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmd2)*Ed0 
     &-xmd2*xmd2*log((Ed0+xlam)/sqrt(xmd2)))
	         veu1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmu2)*Eu0 
     &-xmu2*xmu2*log((Eu0+xlam)/sqrt(xmu2)))
	         ves1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xms2)*Es0 
     &-xms2*xms2*log((Es0+xlam)/sqrt(xms2)))

             ve2=2.*gg*(uu*uu +dd*dd +ss*ss) -4.*gk*uu*dd*ss 
     &+gis*(uuv-ddv)**2

             pot=(ved1 +veu1 +ves1) +ve2 -vedn
             totalS=totalS +pot

             if ((i .ge. nxh-3) .and. (i .le. nxh+4)) then
	         if ((j .ge. nyh-3) .and. (j .le. nyh+4)) then
	         if ((k .ge. nzh-1) .and. (k .le. nzh+2)) then
                centerS=centerS +pot
	         endif
	         endif
	         endif

          enddo
       enddo
	enddo

	xmudcen = xmudcen/32. !d quark mass in central cells
	xmuucen = xmuucen/32. !u quark mass in central cells
	xmuscen = xmuscen/32. !s quark mass in central cells

! step 6 : update particle mass
      summas=0.d0
	sumnum=0.d0
      do npara=1, nevent
	   do i=1, itotal(npara)

          if(ihd(i,npara).eq.1) goto 35 !<=====if a parton is hadronized, its mass won't be affected by mean field.

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 35
          if ((my .le. 1) .or. (my .ge. ny)) goto 35 
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 35 

          if (abs(q(9,i,npara)) .lt. 1.5) then
		     q(8,i,npara)= Vd(mx,my,mz)
			 summas=summas+Vd(mx,my,mz)
			 sumnum=sumnum+1.d0
          else if (abs(q(9,i,npara)) .lt. 2.5) then
		     q(8,i,npara)= Vu(mx,my,mz)
			 summas=summas+Vu(mx,my,mz)
			 sumnum=sumnum+1.d0
		  else
		     q(8,i,npara)= Vs(mx,my,mz)
		  endif

 35   continue

       enddo
	enddo


      it=it+1
	if (it .le. nt) goto 31

! step 7 : update total momentum & energy in each cell AFTER interation
      do npara=1, nevent
	   do i=1, itotal(npara)

          if(ihd(i,npara).eq.1) goto 39 !<=====if a parton is hadronized, it won't contribute to the mean field.
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 39
          if ((my .le. 1) .or. (my .ge. ny)) goto 39
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 39

          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.

          if (abs(q(9,i,npara)) .lt. 1.5) then		       ! for u, d quarks
             Edd(mx,my,mz)= Edd(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)/3.
          Edd(mx+1,my,mz)= Edd(mx+1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my+1,mz)= Edd(mx,my+1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my,mz+1)= Edd(mx,my,mz+1) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx-1,my,mz)= Edd(mx-1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my-1,mz)= Edd(mx,my-1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Edd(mx,my,mz-1)= Edd(mx,my,mz-1) +sqrt(q(8,i,npara)**2.+q2)/9.
          else if (abs(q(9,i,npara)) .lt. 2.5) then		   ! for u, d quarks
             Euu(mx,my,mz)= Euu(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)/3.
          Euu(mx+1,my,mz)= Euu(mx+1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my+1,mz)= Euu(mx,my+1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my,mz+1)= Euu(mx,my,mz+1) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx-1,my,mz)= Euu(mx-1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my-1,mz)= Euu(mx,my-1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
          Euu(mx,my,mz-1)= Euu(mx,my,mz-1) +sqrt(q(8,i,npara)**2.+q2)/9.
		  else                                             ! for s quarks
             Ess(mx,my,mz)= Ess(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)/3.
         Ess(mx+1,my,mz)= Ess(mx+1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
         Ess(mx,my+1,mz)= Ess(mx,my+1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
         Ess(mx,my,mz+1)= Ess(mx,my,mz+1) +sqrt(q(8,i,npara)**2.+q2)/9.
         Ess(mx-1,my,mz)= Ess(mx-1,my,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
         Ess(mx,my-1,mz)= Ess(mx,my-1,mz) +sqrt(q(8,i,npara)**2.+q2)/9.
         Ess(mx,my,mz-1)= Ess(mx,my,mz-1) +sqrt(q(8,i,npara)**2.+q2)/9.
          endif

 39   continue

       enddo
	enddo

!****************************************
! this is for vector mean-fields

! step 1 : find to charge flow velocity
!	do i=1, nx
!	   do j=1, ny
!	      do k=1, nz

!             vxv(i,j,k)= pxv(i,j,k)/(Eud(i,j,k)+Ess(i,j,k))
!             vyv(i,j,k)= pyv(i,j,k)/(Eud(i,j,k)+Ess(i,j,k))
!             vzv(i,j,k)= pzv(i,j,k)/(Eud(i,j,k)+Ess(i,j,k))

!          enddo
!       enddo
!	enddo



! step 2 : move to no current frame
      do npara=1, nevent
	   do i=1, itotal(npara)

          if(ihd(i,npara).eq.1) goto 36 !<=====if a parton is hadronized, it won't contribute to the mean field.

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .le. 1) .or. (mx .ge. nx)) goto 36
          if ((my .le. 1) .or. (my .ge. ny)) goto 36 
          if ((mz .le. 1) .or. (mz .ge. nz)) goto 36  

          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.
          E1=sqrt(q(8,i,npara)**2.+q2)

          bx=vxs(mx,my,mz) 
		  by=vys(mx,my,mz) 
		  bz=vzs(mx,my,mz) 

          call xLorentz(E1,q(5,i,npara),q(6,i,npara),q(7,i,npara)
     &, bx,by,bz, E2,px,py,pz)

          q2= px*px +py*py +pz*pz
		  if (q2 .gt. xlam2) goto 36

!<======count the charge current density in the fireball frame.
		  if(abs(q(9,i,npara)).lt.1.5) then

      chad(mx,my,mz)= chad(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))/3.
      chad(mx+1,my,mz)= chad(mx+1,my,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chad(mx,my+1,mz)= chad(mx,my+1,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chad(mx,my,mz+1)= chad(mx,my,mz+1) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chad(mx-1,my,mz)= chad(mx-1,my,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chad(mx,my-1,mz)= chad(mx,my-1,mz)
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chad(mx,my,mz-1)= chad(mx,my,mz-1) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.

          pxvd(mx,my,mz)= pxvd(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pxvd(mx+1,my,mz)= pxvd(mx+1,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvd(mx,my+1,mz)= pxvd(mx,my+1,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvd(mx,my,mz+1)= pxvd(mx,my,mz+1) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvd(mx-1,my,mz)= pxvd(mx-1,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvd(mx,my-1,mz)= pxvd(mx,my-1,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvd(mx,my,mz-1)= pxvd(mx,my,mz-1) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

          pyvd(mx,my,mz)= pyvd(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pyvd(mx+1,my,mz)= pyvd(mx+1,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvd(mx,my+1,mz)= pyvd(mx,my+1,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvd(mx,my,mz+1)= pyvd(mx,my,mz+1) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvd(mx-1,my,mz)= pyvd(mx-1,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvd(mx,my-1,mz)= pyvd(mx,my-1,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvd(mx,my,mz-1)= pyvd(mx,my,mz-1) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

          pzvd(mx,my,mz)= pzvd(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pzvd(mx+1,my,mz)= pzvd(mx+1,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvd(mx,my+1,mz)= pzvd(mx,my+1,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvd(mx,my,mz+1)= pzvd(mx,my,mz+1) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvd(mx-1,my,mz)= pzvd(mx-1,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvd(mx,my-1,mz)= pzvd(mx,my-1,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvd(mx,my,mz-1)= pzvd(mx,my,mz-1) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

		  else if(abs(q(9,i,npara)).lt.2.5) then

      chau(mx,my,mz)= chau(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))/3.
      chau(mx+1,my,mz)= chau(mx+1,my,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chau(mx,my+1,mz)= chau(mx,my+1,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chau(mx,my,mz+1)= chau(mx,my,mz+1) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chau(mx-1,my,mz)= chau(mx-1,my,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chau(mx,my-1,mz)= chau(mx,my-1,mz)
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chau(mx,my,mz-1)= chau(mx,my,mz-1) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.

          pxvu(mx,my,mz)= pxvu(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pxvu(mx+1,my,mz)= pxvu(mx+1,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvu(mx,my+1,mz)= pxvu(mx,my+1,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvu(mx,my,mz+1)= pxvu(mx,my,mz+1) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvu(mx-1,my,mz)= pxvu(mx-1,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvu(mx,my-1,mz)= pxvu(mx,my-1,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxvu(mx,my,mz-1)= pxvu(mx,my,mz-1) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

          pyvu(mx,my,mz)= pyvu(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pyvu(mx+1,my,mz)= pyvu(mx+1,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvu(mx,my+1,mz)= pyvu(mx,my+1,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvu(mx,my,mz+1)= pyvu(mx,my,mz+1) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvu(mx-1,my,mz)= pyvu(mx-1,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvu(mx,my-1,mz)= pyvu(mx,my-1,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyvu(mx,my,mz-1)= pyvu(mx,my,mz-1) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

          pzvu(mx,my,mz)= pzvu(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pzvu(mx+1,my,mz)= pzvu(mx+1,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvu(mx,my+1,mz)= pzvu(mx,my+1,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvu(mx,my,mz+1)= pzvu(mx,my,mz+1) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvu(mx-1,my,mz)= pzvu(mx-1,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvu(mx,my-1,mz)= pzvu(mx,my-1,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzvu(mx,my,mz-1)= pzvu(mx,my,mz-1) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

		  endif

      chab(mx,my,mz)= chab(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))/3.
      chab(mx+1,my,mz)= chab(mx+1,my,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chab(mx,my+1,mz)= chab(mx,my+1,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chab(mx,my,mz+1)= chab(mx,my,mz+1) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chab(mx-1,my,mz)= chab(mx-1,my,mz) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chab(mx,my-1,mz)= chab(mx,my-1,mz)
     &+q(9,i,npara)/abs(q(9,i,npara))/9.
      chab(mx,my,mz-1)= chab(mx,my,mz-1) 
     &+q(9,i,npara)/abs(q(9,i,npara))/9.

          pxv(mx,my,mz)= pxv(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pxv(mx+1,my,mz)= pxv(mx+1,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxv(mx,my+1,mz)= pxv(mx,my+1,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxv(mx,my,mz+1)= pxv(mx,my,mz+1) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxv(mx-1,my,mz)= pxv(mx-1,my,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxv(mx,my-1,mz)= pxv(mx,my-1,mz) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pxv(mx,my,mz-1)= pxv(mx,my,mz-1) +q(5,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

          pyv(mx,my,mz)= pyv(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pyv(mx+1,my,mz)= pyv(mx+1,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyv(mx,my+1,mz)= pyv(mx,my+1,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyv(mx,my,mz+1)= pyv(mx,my,mz+1) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyv(mx-1,my,mz)= pyv(mx-1,my,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyv(mx,my-1,mz)= pyv(mx,my-1,mz) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pyv(mx,my,mz-1)= pyv(mx,my,mz-1) +q(6,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

          pzv(mx,my,mz)= pzv(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/3.
        pzv(mx+1,my,mz)= pzv(mx+1,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzv(mx,my+1,mz)= pzv(mx,my+1,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzv(mx,my,mz+1)= pzv(mx,my,mz+1) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzv(mx-1,my,mz)= pzv(mx-1,my,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzv(mx,my-1,mz)= pzv(mx,my-1,mz) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.
        pzv(mx,my,mz-1)= pzv(mx,my,mz-1) +q(7,i,npara)/E1*q(9,i,npara)
     &/abs(q(9,i,npara))/9.

 36   continue

       enddo
	enddo

        totalV=0.
        centerV=0.
        cent_Ut=0.
        cent_Utiv=0.
! step 3 : go back to Lab. frame and determine vector-mean fields
        do i=1, nx
          do j=1, ny
            do k=1, nz

             comp0=chab(i,j,k)/dV/npart/nevent
             comp1= pxv(i,j,k)/dV/npart/nevent
             comp2= pyv(i,j,k)/dV/npart/nevent
             comp3= pzv(i,j,k)/dV/npart/nevent

             Ut(i,j,k)= gv*comp0
             Ux(i,j,k)= gv*comp1
             Uy(i,j,k)= gv*comp2
             Uz(i,j,k)= gv*comp3

             comp0iv=(chau(i,j,k)-chad(i,j,k))/dV/npart/nevent
             comp1iv=(pxvu(i,j,k)-pxvd(i,j,k))/dV/npart/nevent
             comp2iv=(pyvu(i,j,k)-pyvd(i,j,k))/dV/npart/nevent
             comp3iv=(pzvu(i,j,k)-pzvd(i,j,k))/dV/npart/nevent

             Utiv(i,j,k)= giv*comp0iv
             Uxiv(i,j,k)= giv*comp1iv
             Uyiv(i,j,k)= giv*comp2iv
             Uziv(i,j,k)= giv*comp3iv

! this is to calculate energy density from vector mean fields

             totalV=totalV +Ut(i,j,k)**2./(2.*gv) +Utiv(i,j,k)**2./giv

             if ((i .ge. nxh-3) .and. (i .le. nxh+4)) then
             if ((j .ge. nyh-3) .and. (j .le. nyh+4)) then
             if ((k .ge. nzh-1) .and. (k .le. nzh+2)) then
             centerV=centerV +Ut(i,j,k)**2./(2.*gv) +Utiv(i,j,k)**2./giv
             cent_Ut=cent_Ut+Ut(i,j,k)
             cent_Utiv=cent_Utiv+Utiv(i,j,k)
             endif
             endif
             endif

            enddo
         enddo
      enddo

!********************************************************************** cascade
 77   continue

	do npara=1, nevent
       do i=1, itotal(npara)-1
	           if(ihd(i,npara).eq.1) goto 70 
	      E1=sqrt(q(8,i,npara)**2. +q(5,i,npara)**2. +q(6,i,npara)**2. 
     &+q(7,i,npara)**2.)

	      do iq=i+1, itotal(npara)
		          if(ihd(iq,npara).eq.1) goto 60 
             dis=sqrt((q(1,i,npara)-q(1,iq,npara))**2.+(q(2,i,npara)
     &-q(2,iq,npara))**2.+(q(3,i,npara)-q(3,iq,npara))**2.)
             if (dis .ge. 3.*sqrt(sigma/pi)) goto 60
             if (dis .lt. 1.d-9) goto 60 !<--- to avoid the successive collision of two particles
	         E2=sqrt(q(8,iq,npara)**2. +q(5,iq,npara)**2. 
     &+q(6,iq,npara)**2. +q(7,iq,npara)**2.)


             E=E1+E2
	         qx=q(5,i,npara)+q(5,iq,npara)
	         qy=q(6,i,npara)+q(6,iq,npara)
	         qz=q(7,i,npara)+q(7,iq,npara)

             bx=qx/E
	         by=qy/E
	         bz=qz/E

             beta2=bx*bx +by*by +bz*bz
	         gamma=1./sqrt(1.-beta2)

             call xLorentz(E1,q(5,i,npara), q(6,i,npara), q(7,i,npara)
     &,  bx,by,bz, Ei, px, py, pz)
	         call xLorentz(t, q(1,i,npara), q(2,i,npara), q(3,i,npara)
     &,  bx,by,bz, ti, xi, yi, zi )
             call xLorentz(t, q(1,iq,npara),q(2,iq,npara),q(3,iq,npara)
     &, bx,by,bz, tiq,xiq,yiq,ziq)


             s=E*E -qx*qx -qy*qy -qz*qz
	         pr2=(s-(q(8,i,npara)+q(8,iq,npara))**2.)*(s-(q(8,i,npara)
     &-q(8,iq,npara))**2.)/(4.*s)
             pr2=max(pr2,1.d-9)
	         pr=sqrt(pr2)

             vxi=px/Ei
	         vyi=py/Ei
	         vzi=pz/Ei

	         Eiq=sqrt(q(8,iq,npara)**2.+px*px +py*py +pz*pz)
	         vxiq=-px/Eiq
	         vyiq=-py/Eiq
	         vziq=-pz/Eiq

	         xi0=xi-vxi*ti
	         yi0=yi-vyi*ti
	         zi0=zi-vzi*ti

	         xiq0=xiq-vxiq*tiq
	         yiq0=yiq-vyiq*tiq
	         ziq0=ziq-vziq*tiq

	         xnum=(vxi-vxiq)*(xi0-xiq0) +(vyi-vyiq)*(yi0-yiq0) 
     &+(vzi-vziq)*(zi0-ziq0)
	         deno=(vxi-vxiq)**2. +(vyi-vyiq)**2. +(vzi-vziq)**2.
             tc=-xnum/deno

             d2=(xi0-xiq0)**2.+(yi0-yiq0)**2.+(zi0-ziq0)**2. 
     &-xnum**2./deno

	         xic=xi0+vxi*tc
	         yic=yi0+vyi*tc
	         zic=zi0+vzi*tc

	         xiqc=xiq0+vxiq*tc
	         yiqc=yiq0+vyiq*tc
	         ziqc=ziq0+vziq*tc

             xc=(xic+xiqc)/2.
             yc=(yic+yiqc)/2.
             zc=(zic+ziqc)/2.

	         call xLorentz(tc,xc,yc,zc, -bx,-by,-bz, tc0,xc0,yc0,zc0)

             if ((d2 .lt. sigma/pi) .and. (tc0 .lt. t+dt) 
     &.and. (tc0 .gt. t)) then

c	            t1=rand()
c		        t2=rand()
				t1 = RANART(NSEED)
				t2 = RANART(NSEED)
		        cs=(t1-0.5)*2.
		        sn=sqrt(1.-cs*cs)
		        ph=t2*(2.*pi)

                Ecm1=sqrt(q(8,i,npara)**2. +pr2)
                Ecm2=sqrt(q(8,iq,npara)**2. +pr2)
                px=pr*sn*cos(ph)
		        py=pr*sn*sin(ph)
		        pz=pr*cs

		     call xLorentz(Ecm1,px,py,pz, -bx,-by,-bz, E1,qx1,qy1,qz1)

                q(5,i,npara)=qx1
 	            q(6,i,npara)=qy1
		        q(7,i,npara)=qz1
		        q(1,i,npara)=xc0 
		        q(2,i,npara)=yc0 
		        q(3,i,npara)=zc0 

		  call xLorentz(Ecm2,-px,-py,-pz, -bx,-by,-bz, E2,qx2,qy2,qz2)

		        q(5,iq,npara)=qx2
		        q(6,iq,npara)=qy2
		        q(7,iq,npara)=qz2
		        q(1,iq,npara)=xc0 
		        q(2,iq,npara)=yc0 
		        q(3,iq,npara)=zc0 

		        ncoll=ncoll+1
		        goto 70

             endif
 60   continue
	      enddo
 70   continue
       enddo
      enddo

!********************************************************************** next time step
!77 continue

      denq = 0.d0            !quark density in central cells
      denqbar = 0.d0           !antiquark density in central cells
      den_d=0.d0
      den_u=0.d0
      den_d_bar=0.d0
      den_u_bar=0.d0
      totalE=0.d0
      energyden=0.d0
      do npara=1, nevent
          do i=1, itotal(npara)

          if(ihd(i,npara).eq.1) goto 79 !<=====if a parton is hadronized, it won't be affected by mean field and no propagation.

          px2=q(5,i,npara)
	      py2=q(6,i,npara)
	      pz2=q(7,i,npara)

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

!	goto 78

          if ((mx .le. 1) .or. (my .le. 1) .or. (mz .le. 1)) then
!		     print *, "lower"
!print *, i, q(1,i,npara), q(2,i,npara), q(3,i,npara)
	         goto 78
		  endif

	      if ((mx .ge. nx) .or. (my .ge. ny) .or. (mz .ge. nz)) then
!		     print *, "upper"
!print *, i, q(1,i,npara), q(2,i,npara), q(3,i,npara)
			 goto 78
		  endif

	      Eqst=sqrt(q(5,i,npara)**2. +q(6,i,npara)**2. 
     &+q(7,i,npara)**2. +q(8,i,npara)**2.)

! momentum update

             bx=vxs(mx,my,mz)
             by=vys(mx,my,mz)
             bz=vzs(mx,my,mz)


              call xLorentz(Eqst,px2,py2,pz2, bx,by,bz, s0,s1,s2,s3)

		  q2=s1**2.+s2**2.+s3**2.
		  if (q2 .gt. xlam**2.) goto 78 

!	goto 78

          if (abs(q(9,i,npara)) .lt. 1.5) then
	         q(5,i,npara)=q(5,i,npara)-dt*(Vd(mx+1,my,mz)
     &-Vd(mx-1,my,mz))/(2.*dx)*q(8,i,npara)/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Vd(mx,my+1,mz)
     &-Vd(mx,my-1,mz))/(2.*dy)*q(8,i,npara)/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Vd(mx,my,mz+1)
     &-Vd(mx,my,mz-1))/(2.*dz)*q(8,i,npara)/Eqst
          else if (abs(q(9,i,npara)) .lt. 2.5) then
	         q(5,i,npara)=q(5,i,npara)-dt*(Vu(mx+1,my,mz)
     &-Vu(mx-1,my,mz))/(2.*dx)*q(8,i,npara)/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Vu(mx,my+1,mz)
     &-Vu(mx,my-1,mz))/(2.*dy)*q(8,i,npara)/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Vu(mx,my,mz+1)
     &-Vu(mx,my,mz-1))/(2.*dz)*q(8,i,npara)/Eqst
          else
	         q(5,i,npara)=q(5,i,npara)-dt*(Vs(mx+1,my,mz)
     &-Vs(mx-1,my,mz))/(2.*dx)*q(8,i,npara)/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Vs(mx,my+1,mz)
     &-Vs(mx,my-1,mz))/(2.*dy)*q(8,i,npara)/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Vs(mx,my,mz+1)
     &-Vs(mx,my,mz-1))/(2.*dz)*q(8,i,npara)/Eqst
		  endif

          if (q(9,i,npara) .gt. 0.d0) then

	         q(5,i,npara)=q(5,i,npara)-dt*(Ut(mx+1,my,mz)
     &-Ut(mx-1,my,mz))/(2.*dx)
	         q(6,i,npara)=q(6,i,npara)-dt*(Ut(mx,my+1,mz)
     &-Ut(mx,my-1,mz))/(2.*dy)
	         q(7,i,npara)=q(7,i,npara)-dt*(Ut(mx,my,mz+1)
     &-Ut(mx,my,mz-1))/(2.*dz)

             q(5,i,npara)=q(5,i,npara)-dt*(Ux(mx,my+1,mz)
     &-Ux(mx,my-1,mz))/(2.*dy) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)+dt*(Uy(mx+1,my,mz)
     &-Uy(mx-1,my,mz))/(2.*dx) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)+dt*(Uz(mx+1,my,mz)
     &-Uz(mx-1,my,mz))/(2.*dx) *pz2/Eqst
             q(5,i,npara)=q(5,i,npara)-dt*(Ux(mx,my,mz+1)
     &-Ux(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(6,i,npara)=q(6,i,npara)+dt*(Ux(mx,my+1,mz)
     &-Ux(mx,my-1,mz))/(2.*dy) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Uy(mx+1,my,mz)
     &-Uy(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)+dt*(Uz(mx,my+1,mz)
     &-Uz(mx,my-1,mz))/(2.*dy) *pz2/Eqst
             q(6,i,npara)=q(6,i,npara)-dt*(Uy(mx,my,mz+1)
     &-Uy(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(7,i,npara)=q(7,i,npara)+dt*(Ux(mx,my,mz+1)
     &-Ux(mx,my,mz-1))/(2.*dz) *px2/Eqst
             q(7,i,npara)=q(7,i,npara)-dt*(Uz(mx+1,my,mz)
     &-Uz(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(7,i,npara)=q(7,i,npara)+dt*(Uy(mx,my,mz+1)
     &-Uy(mx,my,mz-1))/(2.*dz) *py2/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Uz(mx,my+1,mz)
     &-Uz(mx,my-1,mz))/(2.*dy) *py2/Eqst

			 if(abs(q(9,i,npara)) .lt. 1.5) then

	         q(5,i,npara)=q(5,i,npara)+dt*(Utiv(mx+1,my,mz)
     &-Utiv(mx-1,my,mz))/(2.*dx)
	         q(6,i,npara)=q(6,i,npara)+dt*(Utiv(mx,my+1,mz)
     &-Utiv(mx,my-1,mz))/(2.*dy)
	         q(7,i,npara)=q(7,i,npara)+dt*(Utiv(mx,my,mz+1)
     &-Utiv(mx,my,mz-1))/(2.*dz)

             q(5,i,npara)=q(5,i,npara)+dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)-dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)-dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *pz2/Eqst
             q(5,i,npara)=q(5,i,npara)+dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(6,i,npara)=q(6,i,npara)-dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)+dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *pz2/Eqst
             q(6,i,npara)=q(6,i,npara)+dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(7,i,npara)=q(7,i,npara)-dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *px2/Eqst
             q(7,i,npara)=q(7,i,npara)+dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *py2/Eqst
	         q(7,i,npara)=q(7,i,npara)+dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *py2/Eqst

             else if (abs(q(9,i,npara)) .lt. 2.5) then

	         q(5,i,npara)=q(5,i,npara)-dt*(Utiv(mx+1,my,mz)
     &-Utiv(mx-1,my,mz))/(2.*dx)
	         q(6,i,npara)=q(6,i,npara)-dt*(Utiv(mx,my+1,mz)
     &-Utiv(mx,my-1,mz))/(2.*dy)
	         q(7,i,npara)=q(7,i,npara)-dt*(Utiv(mx,my,mz+1)
     &-Utiv(mx,my,mz-1))/(2.*dz)

             q(5,i,npara)=q(5,i,npara)-dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)+dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)+dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *pz2/Eqst
             q(5,i,npara)=q(5,i,npara)-dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(6,i,npara)=q(6,i,npara)+dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)+dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *pz2/Eqst
             q(6,i,npara)=q(6,i,npara)-dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(7,i,npara)=q(7,i,npara)+dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *px2/Eqst
             q(7,i,npara)=q(7,i,npara)-dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(7,i,npara)=q(7,i,npara)+dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *py2/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *py2/Eqst

			 endif

          else

	         q(5,i,npara)=q(5,i,npara)+dt*(Ut(mx+1,my,mz)
     &-Ut(mx-1,my,mz))/(2.*dx)
	         q(6,i,npara)=q(6,i,npara)+dt*(Ut(mx,my+1,mz)
     &-Ut(mx,my-1,mz))/(2.*dy)
	         q(7,i,npara)=q(7,i,npara)+dt*(Ut(mx,my,mz+1)
     &-Ut(mx,my,mz-1))/(2.*dz)

             q(5,i,npara)=q(5,i,npara)+dt*(Ux(mx,my+1,mz)
     &-Ux(mx,my-1,mz))/(2.*dy) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)-dt*(Uy(mx+1,my,mz)
     &-Uy(mx-1,my,mz))/(2.*dx) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)-dt*(Uz(mx+1,my,mz)
     &-Uz(mx-1,my,mz))/(2.*dx) *pz2/Eqst
             q(5,i,npara)=q(5,i,npara)+dt*(Ux(mx,my,mz+1)
     &-Ux(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(6,i,npara)=q(6,i,npara)-dt*(Ux(mx,my+1,mz)
     &-Ux(mx,my-1,mz))/(2.*dy) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)+dt*(Uy(mx+1,my,mz)
     &-Uy(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Uz(mx,my+1,mz)
     &-Uz(mx,my-1,mz))/(2.*dy) *pz2/Eqst
             q(6,i,npara)=q(6,i,npara)+dt*(Uy(mx,my,mz+1)
     &-Uy(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(7,i,npara)=q(7,i,npara)-dt*(Ux(mx,my,mz+1)
     &-Ux(mx,my,mz-1))/(2.*dz) *px2/Eqst
             q(7,i,npara)=q(7,i,npara)+dt*(Uz(mx+1,my,mz)
     &-Uz(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Uy(mx,my,mz+1)
     &-Uy(mx,my,mz-1))/(2.*dz) *py2/Eqst
	         q(7,i,npara)=q(7,i,npara)+dt*(Uz(mx,my+1,mz)
     &-Uz(mx,my-1,mz))/(2.*dy) *py2/Eqst

			 if(abs(q(9,i,npara)) .lt. 1.5) then
             
	         q(5,i,npara)=q(5,i,npara)-dt*(Utiv(mx+1,my,mz)
     &-Utiv(mx-1,my,mz))/(2.*dx)
	         q(6,i,npara)=q(6,i,npara)-dt*(Utiv(mx,my+1,mz)
     &-Utiv(mx,my-1,mz))/(2.*dy)
	         q(7,i,npara)=q(7,i,npara)-dt*(Utiv(mx,my,mz+1)
     &-Utiv(mx,my,mz-1))/(2.*dz)

             q(5,i,npara)=q(5,i,npara)-dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)+dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)+dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *pz2/Eqst
             q(5,i,npara)=q(5,i,npara)-dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(6,i,npara)=q(6,i,npara)+dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)+dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *pz2/Eqst
             q(6,i,npara)=q(6,i,npara)-dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(7,i,npara)=q(7,i,npara)+dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *px2/Eqst
             q(7,i,npara)=q(7,i,npara)-dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(7,i,npara)=q(7,i,npara)+dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *py2/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *py2/Eqst			 

			 else if (abs(q(9,i,npara)) .lt. 2.5) then
			 
	         q(5,i,npara)=q(5,i,npara)+dt*(Utiv(mx+1,my,mz)
     &-Utiv(mx-1,my,mz))/(2.*dx)
	         q(6,i,npara)=q(6,i,npara)+dt*(Utiv(mx,my+1,mz)
     &-Utiv(mx,my-1,mz))/(2.*dy)
	         q(7,i,npara)=q(7,i,npara)+dt*(Utiv(mx,my,mz+1)
     &-Utiv(mx,my,mz-1))/(2.*dz)

             q(5,i,npara)=q(5,i,npara)+dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)-dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *py2/Eqst
	         q(5,i,npara)=q(5,i,npara)-dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *pz2/Eqst
             q(5,i,npara)=q(5,i,npara)+dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(6,i,npara)=q(6,i,npara)-dt*(Uxiv(mx,my+1,mz)
     &-Uxiv(mx,my-1,mz))/(2.*dy) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)+dt*(Uyiv(mx+1,my,mz)
     &-Uyiv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *pz2/Eqst
             q(6,i,npara)=q(6,i,npara)+dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *pz2/Eqst

	         q(7,i,npara)=q(7,i,npara)-dt*(Uxiv(mx,my,mz+1)
     &-Uxiv(mx,my,mz-1))/(2.*dz) *px2/Eqst
             q(7,i,npara)=q(7,i,npara)+dt*(Uziv(mx+1,my,mz)
     &-Uziv(mx-1,my,mz))/(2.*dx) *px2/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(Uyiv(mx,my,mz+1)
     &-Uyiv(mx,my,mz-1))/(2.*dz) *py2/Eqst
	         q(7,i,npara)=q(7,i,npara)+dt*(Uziv(mx,my+1,mz)
     &-Uziv(mx,my-1,mz))/(2.*dy) *py2/Eqst	
			 			 
			 endif

             endif

              Eqst=sqrt(q(5,i,npara)**2. +q(6,i,npara)**2. 
     &+q(7,i,npara)**2. +q(8,i,npara)**2.)
          totalE=totalE +Eqst

              if ((mx .ge. nxh-3) .and. (mx .le. nxh+4)) then
              if ((my .ge. nyh-3) .and. (my .le. nyh+4)) then
              if ((mz .ge. nzh-1) .and. (mz .le. nzh+2)) then
                energyden=energyden+Eqst
                if(q(9,i,npara).gt.0) then
                  denq = denq + 1. 
                  if(abs(q(9,i,npara)).lt.1.5) then
                  den_d=den_d+1.
                  elseif(abs(q(9,i,npara)).lt.2.5) then
                  den_u=den_u+1.
                  endif
                elseif(q(9,i,npara).lt.0) then
                  denqbar = denqbar + 1.
                  if(abs(q(9,i,npara)).lt.1.5) then
                  den_d_bar=den_d_bar+1.
                  elseif(abs(q(9,i,npara)).lt.2.5) then
                  den_u_bar=den_u_bar+1.
                  endif
                endif
 
              endif
              endif
              endif

! position update

 78    continue
             Eqst=sqrt(q(5,i,npara)**2. +q(6,i,npara)**2.
     &+q(7,i,npara)**2. +q(8,i,npara)**2.)
	      vx=q(5,i,npara)/Eqst
	      vy=q(6,i,npara)/Eqst
	      vz=q(7,i,npara)/Eqst
	      q(1,i,npara)=q(1,i,npara) +dt*vx
          q(2,i,npara)=q(2,i,npara) +dt*vy
          q(3,i,npara)=q(3,i,npara) +dt*vz


! for p, not for p*

	      q0(1,i,npara)= Ux(mx,my,mz)
          q0(2,i,npara)= Uy(mx,my,mz)
          q0(3,i,npara)= Uz(mx,my,mz)

!The output is p*

!<============Determine the hadronization mass by flavor============
          if (abs(q(9,i,npara)) .lt. 1.5) then ! for d quarks
            xmhd = 0.5*xmdv
          else if (abs(q(9,i,npara)).lt.2.5) then !for u quarks
            xmhd=0.5*xmuv
          else !for s quarks
            xmhd=0.5*xmsv
          end if
!<==================================================================

!<=============Once a parton is hadronized, relabel it and record its freeze out time=====

        if(iglobal.eq.0) then 
          if((q(8,i,npara).gt.xmhd).and.(ihd(i,npara).eq.0)
     &.and.t.gt.thd) then
c            iparthd(npara)=iparthd(npara)+1 !<====# of hadronized partons +1
            ihd(i,npara)=1
            q(4,i,npara)=t
          end if
        endif
 79   continue
!<==============================================================
	   enddo
      enddo

      t=t+dt

      write(*,*) "t=", t, dt, ntest

      centerEd=(energyden/(8.*32.*dx*dy*dz)/npart/nevent 
     &+(centerS+centerV)/(8.*32.)*5.07**3.) ! unit is GeV/fm^3
      totalE=(totalE/npart/nevent +(totalS+totalV)*dV) ! unit is GeV
      centdenq = denq/(8.*32.*dx*dy*dz)/npart/nevent 
      centdenqbar = denqbar/(8.*32.*dx*dy*dz)/npart/nevent 
      centdenu = (den_u-den_u_bar)/((8.*32.)*dx*dy*dz)/npart/nevent
      centdend = (den_d-den_d_bar)/((8.*32.)*dx*dy*dz)/npart/nevent
      xmucen = xmudcen/2. + xmuucen/2.
      cent_Ut=cent_Ut/(8.*32.)
      cent_Utiv=cent_Utiv/(8.*32)
      write(98,2014) t,centdenq,centdenqbar,centerEd,cent_Ut,cent_Utiv
2014  format(1X,6F10.4)
      if(t.lt.tend) then
        if(iglobal.eq.1) then
           if((t .lt. thd) .or. (xmucen .lt. 0.5*xmu)) goto 30
        else if(iglobal.eq.0) then
           goto 30
        endif
      endif

 15   format(1x,3(i7,1x))
 16   format(1x,3(i5,1x),8(f9.6,1x))
2013  FORMAT(3I6,8(1X,F10.3))
      do npara=1, nevent
c       write(14,15) npara,ipart(npara),nsg(npara)
c	   do i=1, ipart(npara)

	do i=1,itotal(npara)
	      Eqs=sqrt(q(5,i,npara)**2. +q(6,i,npara)**2. 
     &+q(7,i,npara)**2. +q(8,i,npara)**2.)	 
	k = imp(i,npara) !keep the initial order
	GXnjl(k,npara) = q(1,i,npara)
	GYnjl(k,npara) = q(2,i,npara)
	GZnjl(k,npara) = q(3,i,npara)
	if(ihd(i,npara).eq.1) then
	FTnjl(k,npara) = q(4,i,npara)
	else
	FTnjl(k,npara) = t-dt
	endif
	PXnjl(k,npara) = q(5,i,npara)
	PYnjl(k,npara) = q(6,i,npara)
	PZnjl(k,npara) = q(7,i,npara)
	XMASSnjl(k,npara) = q(8,i,npara)
	Enjl(k,npara) = Eqs
	ITYPnjl(k,npara) = int(q(9,i,npara))
	enddo

	do i=1,njlmul(npara)
	if(par(4,i,npara).gt.t-dt) then
	      Eqs=sqrt(par(5,i,npara)**2. +par(6,i,npara)**2. 
     &+par(7,i,npara)**2. +par(8,i,npara)**2.)  
	GXnjl(i,npara) = par(1,i,npara)
	GYnjl(i,npara) = par(2,i,npara)
	GZnjl(i,npara) = par(3,i,npara)
	FTnjl(i,npara) = par(4,i,npara)
	PXnjl(i,npara) = par(5,i,npara)
	PYnjl(i,npara) = par(6,i,npara)
	PZnjl(i,npara) = par(7,i,npara)
	XMASSnjl(i,npara) = par(8,i,npara)
	Enjl(i,npara) = Eqs
	ITYPnjl(i,npara) = int(par(9,i,npara))
	endif
	enddo
		
	enddo

c	enddo 
c	close(99) 
c      end program
	end

	  subroutine xLorentz(t1,x1,y1,z1, vx,vy,vz, t2,x2,y2,z2)
	  implicit real*8 (a-h,o-z)

      beta2=vx*vx +vy*vy +vz*vz
      if (beta2.gt.1.d-9) then  
      gam=1./sqrt(1.-beta2)

      t2=gam*(t1 -vx*x1 -vy*y1 -vz*z1)
      x2=x1 -vx*gam*t1 +(gam-1.)*vx*vx/beta2*x1 +(gam-1.)*vx*vy/beta2*y1 
     &+(gam-1.)*vx*vz/beta2*z1
      y2=y1 -vy*gam*t1 +(gam-1.)*vy*vx/beta2*x1 +(gam-1.)*vy*vy/beta2*y1 
     &+(gam-1.)*vy*vz/beta2*z1
      z2=z1 -vz*gam*t1 +(gam-1.)*vz*vx/beta2*x1 +(gam-1.)*vz*vy/beta2*y1 
     &+(gam-1.)*vz*vz/beta2*z1
      else
      t2=t1
      x2=x1
      y2=y1
      z2=z1
      end if

      return
      end

      subroutine potential(i,npara,i1,i2,i3,h)
	  implicit real*8 (a-h,o-z)
	common/particle/ par(9,20000,100), q(9,20000,100), q0(3,20000,100)
     &, ihd(20000,100)
      common/scalar/ Vu(100,100,40), Vd(100,100,40), Vs(100,100,40)! scalar potential
      common/vector/ Ut(100,100,40), Ux(100,100,40), Uy(100,100,40)
     &, Uz(100,100,40) ! vector poential
      common/grid1/nx, ny, nz
      common/grid2/dx, dy, dz

      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	  my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	  mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

      if (q(9,i,npara).ge.0.) then
          sign=1.
      else
          sign=-1.
      endif

!      xmass=q(8,i,npara)
      if (abs(q(9,i,npara)).lt.1.5) then
          xmass=Vd(i1,i2,i3)
      else if (abs(q(9,i,npara)).lt.2.5) then
          xmass=Vu(i1,i2,i3)
      else
          xmass=Vs(i1,i2,i3)
      endif

      px=q(5,i,npara)+sign*(Ux(mx,my,mz)-Ux(i1,i2,i3)) !q(5-7,i,npara) is p*, not p !!!
      py=q(6,i,npara)+sign*(Uy(mx,my,mz)-Uy(i1,i2,i3))
      pz=q(7,i,npara)+sign*(Uz(mx,my,mz)-Uz(i1,i2,i3))
      pot=sign*Ut(i1,i2,i3)

      h=sqrt(xmass**2. +px*px +py*py +pz*pz) +pot

      return
      end

c      SUBROUTINE SRAND(ISEED)
!C
!C  This subroutine sets the integer seed to be used with the
!C  companion RAND function to the value of ISEED.  A flag is
!C  set to indicate that the sequence of pseudo-random numbers
!C  for the specified seed should start from the beginning.
!C
c      COMMON /SEED/JSEED,IFRST
!C
c      JSEED = ISEED
c      IFRST = 0
!C
c      RETURN
c      END

c      REAL FUNCTION RAND()
!C
!C  This function returns a pseudo-random number for each invocation.
!C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
!C  standard number generator whose Pascal code appears in the article:
!C
!C     Park, Steven K. and Miller, Keith W., "Random Number Generators:
!C     Good Ones are Hard to Find", Communications of the ACM,
!C     October, 1988.
!C
c      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773, 
c     & MOMDMP=2836)
!C
c      COMMON  /SEED/JSEED,IFRST
c      INTEGER HVLUE, LVLUE, TESTV, NEXTN
c      SAVE    NEXTN
!C
c      IF (IFRST .EQ. 0) THEN
c        NEXTN = JSEED
c        IFRST = 1
c      ENDIF
!C
c      HVLUE = NEXTN / MOBYMP
c      LVLUE = MOD(NEXTN, MOBYMP)
c      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
c      IF (TESTV .GT. 0) THEN
c        NEXTN = TESTV
c      ELSE
c        NEXTN = TESTV + MODLUS
c      ENDIF
c      RAND = REAL(NEXTN)/REAL(MODLUS)
!C
c      RETURN
c      END
c      BLOCKDATA RANDBD
c      COMMON /SEED/JSEED,IFRST
!C
c      DATA JSEED,IFRST/123456789,0/
!C
c      END
