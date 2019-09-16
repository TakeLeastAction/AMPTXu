    PROGRAM YieldCollect
	implicit real*8 (a-h,o-z) 
	integer,parameter :: Njobmax= 10
        character(50) :: datatype
        character(3) :: tit          ! title   
        character(50):: DirP
	dimension xpv(Njobmax),xdv(Njobmax),xppv(Njobmax)
        dimension a(3),b(3),c(3)
        datatype='/yield.dat'
        DirP = 'data0/ampt-'    
       
        a(1)=1.0
        a(2)=2.0
        a(3)=3.0
        b(1)=1.0
        b(2)=2.0
        b(3)=3.0

        c = 1/2*b

         write(*,*)'a=',a
          write(*,*)'c=',c

       xp = 0.
	   xn = 0.
	   xd = 0.
           xd2 = 0.
           xt2 = 0.
	   xt=0.
	   xpp=0.
	   xpm=0.
	   xkp=0.
	   xkm=0.
	open (unit=20, file='yield.dat', status='unknown')
	do njob=1,njobmax
           write(*,*)"job=", njob
	   write(tit,"(I3.3)")njob   
           open (116, file=trim(DirP)//tit//trim(datatype), status='unknown')
           read(116,*)xpionp,xpionm,xkaonp,xkaonm,xnp,xnn,xnd,xnt
       write(20,101)xpionp,xpionm,xkaonp,xkaonm,xnp,xnn,xnd,xnt
       xp = xp+xnp
	   xn = xn +xnn
	   xd = xd+xnd
	   xt=xt+xnt
           xd2 = xd2+xnd**2.
           xt2 = xt2+xnt**2.
	   xpp = xpp+xpionp
	   xpm = xpm+xpionm
	   xkp = xkp+xkaonp
	   xkm = xkm+xkaonm
	   
	   xpv(njob) = xnp
	   xppv(njob) = xpionp
	   xdv(njob) = xnd
	   
	 enddo
	 xp = xp/Njobmax
	 xn = xn/njobmax
	 xd = xd/njobmax
	 xt = xt/njobmax
         xd2 = xd2/njobmax
         xt2 = xt2/njobmax
	 xpp = xpp/njobmax
	 xpm = xpm/njobmax
	 xkp = xkp/njobmax
	 xkm = xkm/njobmax
	 erd = sqrt(xd2-xd**2.)
         ert = sqrt(xt2-xt**2.)
	 write(*,*)"number in central rapidity:",xpp,xpm,xkp,xkm,xp,xn,xd,xt
	 write(*,*)"charged particle mult. = ",xpp+xpm+xkp+xkm
	 write(*,*)"he3/p ratio:",xt/xp,ert/xp
	 write(*,*)"d/p ratio:",xd/xp,erd/xp
	 write(*,*)"p/pion ratio:",xp/xpp
         write(*,*)"dv = ",xdv	
101 format(8(1x,F16.10))
    Stop
    end
