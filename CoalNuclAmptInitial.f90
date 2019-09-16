    PROGRAM CoalNucl
	implicit real*8 (a-h,o-z) 
	integer,parameter :: Njobmax= 1
	integer,parameter :: NumPar= 1000, Ntest=1000
	dimension xpro(8,NumPar,Ntest),xneu(8,NumPar,Ntest)
	dimension ipro(Ntest),ineu(Ntest),xndv(Ntest)
	dimension xpv(Njobmax),xdv(Njobmax),xppv(Njobmax),xtv(Njobmax)
	
       xp = 0.
	   xn = 0.
	   xd = 0.
	   xhe3=0.
	   xd2 = 0.
	   xt=0.
	   xt2 = 0.
	   xpp=0.
	   xpm=0.
	   xkp=0.
	   xkm=0.
	   nevent=80
	open (unit=20, file='yield.dat', status='unknown')
	!open (unit=21, file='yieldv.dat', status='unknown')
	do njob=1,njobmax
	   CALL DataRead(nevent,njob,xpro,xneu,ipro,ineu,xnp,xnn,xpionp,xpionm,xkaonp,xkaonm)
	   write(*,*)"dear master: ",ineu(1:nevent)
	   CALL COALdeu(nevent,njob,xpro,xneu,ipro,ineu,xnd) 
	   write(*,*)"d/p ratio:",xnd/xnp
	   rms = 1.76
	   nevent=40
	   CALL COALhe3(nevent,njob,xpro,xneu,ipro,ineu,rms,xhe3)
       write(20,101)xpionp,xpionm,xkaonp,xkaonm,xnp,xnn,xnd,xhe3
	   !do k=1,nevent
	   !write(21,102)xndv(1:20)
	   !enddo
       xp = xp+xnp
	   xn = xn +xnn
	   xd = xd+xnd
	   xt = xt+xhe3
	   xpp = xpp+xpionp
	   xpm = xpm+xpionm
	   xkp = xkp+xkaonp
	   xkm = xkm+xkaonm
	   xpv(njob) = xnp
	   xppv(njob) = xpionp
	   xdv(njob) = xnd
	   xtv(njob) = xhe3
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
	 write(*,*)"number in central rapidity:",xpp,xpm,xkp,xkm,xp,xn,xd,xt
	 write(*,*)"charged particle mult. = ",xpp+xpm+xkp+xkm
	 write(*,*)"proton: = ",xp," neuteron: =",xn
	 write(*,*)"p/pion ratio:",xp/xpp
	 write(*,*)"d/p ratio:",xd/xp	 
	 write(*,*)"he-3/p ratio:",xt/xp
	 write(*,*)"he-3 p/d^2 ratio:",xt*xp/xd**2.*xn/xp
	 write(20,103)xd/xp,xt/xp,xt*xp/xd**2.*xn/xp
102 format(20(1x,F16.10))	 
101 format(8(1x,F16.10))
103 format(3(1x,F16.10))
       STOP
       END
	   
	   
	  SUBROUTINE DataRead(Nevent,njob,xpro,xneu,ipro,ineu,avp,avn,avpionp,avpionm,avkaonp,avkaonm)
	  implicit real*8 (a-h,o-z)
	  integer,parameter :: NumPar= 1000, Ntest=1000
	  dimension xpro(8,NumPar,Ntest),xneu(8,NumPar,Ntest)
	  dimension ipro(Ntest),ineu(Ntest)
	  character(50) :: datatype
      character(3) :: tit          ! title   
      character(50):: DirP
	  !integer :: npara,nevent,i,kp,kn,id,num
	  !integer :: i1,i2,i3,i4,i5,i6,i7,i8,k
	  !real*8 :: xx,yy,zz,time,px,py,pz,xmass,f1,f2
	  
      datatype='ana/initial_ampt.dat' 
      DirP = 'ampt-'	  
	  write(tit,"(I3.3)")njob 
	  !OPEN (116, FILE = 'ana/ampt.dat', STATUS = 'UNKNOWN')
      open (116, file=trim(datatype), status='unknown')
	  open (117, file="neutron.dat", status='unknown')
	  open (118, file="proton.dat", status='unknown')
	  open (119, file="inucleon.dat", status='unknown')
	  open (120, file="ipion.dat", status='unknown')
	  open (121, file="pion.dat", status='unknown')	 
	  avp = 0.
	  avn = 0.
	  avpionp = 0.
	  avpionm = 0.
	  avkaonp = 0.
	  avkaonm = 0.
	  do npara=1, Nevent 
       !print *, "event = ", npara   
          read(116,*) i1,i2,num,f1,i3,i4,i5,i6,i7,i8,f2
		  
		  	 kp=0
		     kn=0
			 kpcut=0
			 kncut=0
			 kpi = 0
          do k=1, num
             read(116,*) id,px,py,pz,xmass,xx,yy,zz,time
			 
			 P = (px**2.+py**2.+pz**2.)**0.5
             E = (P**2.+xmass**2.)**0.5
				 
			 Pt = (px**2.+py**2.)**0.5
	         PL = pz
	         yita = 1./2.*log((E+PL)/(E-PL))
			 if((id.eq.211).and.(abs(yita)<0.5))then
			    avpionp = avpionp+1
			 endif
			 
			 if((id.eq.-211).and.(abs(yita)<0.5))then
			    avpionm = avpionm+1
			 endif
			 if((id.eq.321).and.(abs(yita)<0.5))then
			    avkaonp = avkaonp+1
			 endif	

			 if((id.eq.-321).and.(abs(yita)<0.5))then
			    avkaonm = avkaonm+1
			 endif				 
			 
             if((id.eq.2212).and.(time>0.5).and.(time<100).and.(abs(pz)<100.)) then
			 kp=kp+1
			 xpro(1,kp,npara)=xx
	         xpro(2,kp,npara)=yy
	         xpro(3,kp,npara)=zz
             xpro(4,kp,npara)=time
	         xpro(5,kp,npara)=px
	         xpro(6,kp,npara)=py
	         xpro(7,kp,npara)=pz
             xpro(8,kp,npara)=xmass
			 write(118,101) id,px,py,pz,xmass,xx,yy,zz,time
			 if(abs(yita)<0.5)then
                kpcut=kpcut+1			 
			 endif
			 
			 !write(*,*) id,px,py,pz,xmass,xx,yy,zz,time
			 end if
			 
			 if((id.eq.211).and.(time>0.5).and.(time<100).and.(abs(pz)<100.)) then
			 kpi=kpi+1
			 write(121,101) id,px,py,pz,xmass,xx,yy,zz,time
			 end if
			 
			 if((id.eq.2112).and.(time>0.5).and.(time<100).and.(abs(pz)<100.)) then
			 kn=kn+1
			 xneu(1,kn,npara)=xx
	         xneu(2,kn,npara)=yy
	         xneu(3,kn,npara)=zz
             xneu(4,kn,npara)=time
	         xneu(5,kn,npara)=px
	         xneu(6,kn,npara)=py
	         xneu(7,kn,npara)=pz
             xneu(8,kn,npara)=xmass
			 write(117,101) id,px,py,pz,xmass,xx,yy,zz,time
			 if(abs(yita)<0.5)then
                kncut=kncut+1			 
			 endif
			 end if
	      enddo ! for particles 
		  
		  ipro(npara) = kp
		  ineu(npara) = kn
		  avp = avp+kpcut
		  avn = avn+kncut
		  !write(*,*)kp,kn	
		  write(119,*) kn,kp
		  write(120,*) kpi
		  
		
		enddo ! for events
       
	   avp = avp/float(Nevent)
	   avn = avn/float(Nevent)
	   avpionp = avpionp/float(Nevent)
	   avpionm = avpionm/float(Nevent)
	   avkaonp = avkaonp/float(Nevent)
	   avkaonm = avkaonm/float(Nevent)
	   !write(*,*)"deuteron:",xndv
	   write(*,*)"njob = ", njob
	   write(*,*)"proton number:",avp," neutron number: ", avn
	   
	  
	  
101 format((1x,I6),8(1x,F10.5))	  
	  return
	  end
	  
	   
	  SUBROUTINE COALdeu(Nevent,njob,xpro,xneu,ipro,ineu,avnd)
	  implicit real*8 (a-h,o-z)
	  integer,parameter :: NumPar= 1000, Ntest=1000
	  dimension xpro(8,NumPar,Ntest),xneu(8,NumPar,Ntest)
	  dimension ipro(Ntest),ineu(Ntest)
	  character(50) :: datatype
      character(3) :: tit          ! title   
      character(50):: DirP
	  !integer :: npara,nevent,i,kp,kn,id,num
	  !integer :: i1,i2,i3,i4,i5,i6,i7,i8,k
	  !real*8 :: xx,yy,zz,time,px,py,pz,xmass,f1,f2
	  
      
	  avnd = 0.
       nparak = 0
       do npara=1,Nevent
	   write(*,*)"dear master: deuteron coalescence: ",npara
	   do npara1=1,Nevent
	      nparak = nparak+1
              !write(*,*)"dear master: testing deu",npara,npara1
	      cut_dr = 5.
		  cut_dp = 5.
		  hbc = 0.19733
		  hbc2=hbc**2.
		  sf = 3.2
		  xnd=0.
		  gc = 3./4.
		  nnp = ipro(npara)
        do k1=1,nnp
			xx1 = xpro(1,k1,npara)
			yy1 = xpro(2,k1,npara)
			zz1 = xpro(3,k1,npara)
			t1 = xpro(4,k1,npara)
			px1 = xpro(5,k1,npara)
			py1=xpro(6,k1,npara)
			pz1 = xpro(7,k1,npara)
			xm1 = xpro(8,k1,npara)
			xE1 = sqrt(xm1**2.+px1**2.+py1**2.+pz1**2.)
		    !write(*,*) "pro: ",k1,t1,xx1,xm1
			
			nnn = ineu(npara1)
            !write(*,*) "dear master testing",npara,npara1,nnp,nnn,ineu(npara1),ineu(34)
             
		   do k2=1,nnn
		    xx2 = xneu(1,k2,npara1)
			yy2 = xneu(2,k2,npara1)
			zz2 = xneu(3,k2,npara1)
			t2 = xneu(4,k2,npara1)
			px2 = xneu(5,k2,npara1)
			py2=xneu(6,k2,npara1)
			pz2 = xneu(7,k2,npara1)
			xm2 = xneu(8,k2,npara1)
		   	xE2 = sqrt(xm2**2.+px2**2.+py2**2.+pz2**2.)

			
			!px=px1+px2
			 !py=py1+py2
			 pz=pz1+pz2
			 xmass=xm1+xm2
			 Ed = xE1+xE2
			 !Pt = (px**2.+py**2.)**0.5
	         PL = pz
	         yitad = 1./2.*log((Ed+PL)/(Ed-PL))
			 if(abs(yitad)>0.5) goto 10
			
            betax = px/Ed
			betay = py/Ed
			betaz = pz/Ed
			
			xx1new = xx1
			yy1new = yy1
			zz1new = zz1
            t1new = t1
			xx2new = xx2
			yy2new = yy2
			zz2new = zz2
			t2new = t2
			
			px1new = px1
			py1new = py1
			pz1new = pz1
            xE1new = sqrt(xm1**2.+px1new**2.+py1new**2.+pz1new**2.)
			px2new = px2
			py2new = py2
			pz2new = pz2
			xE2new = sqrt(xm2**2.+px2new**2.+py2new**2.+pz2new**2.)
			
			!call lorentzboost(betax,betay,betaz,t1,xx1,yy1,zz1,t1new,xx1new,yy1new,zz1new)
			!call lorentzboost(betax,betay,betaz,xE1,px1,py1,pz1,xE1new,px1new,py1new,pz1new)
			!call lorentzboost(betax,betay,betaz,t2,xx2,yy2,zz2,t2new,xx2new,yy2new,zz2new)
			!call lorentzboost(betax,betay,betaz,xE2,px2,py2,pz2,xE2new,px2new,py2new,pz2new)
			
			
			
			
			 
			!write(*,*)npara,k1,k2, t1,t2,tmax
			!write(*,*) xE1,xE2
			
			tmax=max(t1new,t2new)
			  xx2new = xx2new+(tmax-t2new)*px2new/xE2new
			  yy2new = yy2new+(tmax-t2new)*py2new/xE2new
			  zz2new = zz2new+(tmax-t2new)*pz2new/xE2new
			  
			 
			  xx1new = xx1new+(tmax-t1new)*px1new/xE1new
			  yy1new = yy1new+(tmax-t1new)*py1new/xE1new
			  zz1new = zz1new+(tmax-t1new)*pz1new/xE1new
			  
			  
			  
			
			
		    dx1 = (xx1new-xx2new) !/2.**0.5
			dx2 = (yy1new-yy2new) !/2.**0.5
			dx3 = (zz1new-zz2new) !/2.**0.5
			dp1 = (px1new-px2new)/2. !**0.5
			dp2 = (py1new-py2new)/2. !**0.5
			dp3 = (pz1new-pz2new)/2. !**0.5
			disx = (dx1**2.+dx2**2.+dx3**2.)**0.5
 			if(disx >cut_dr*sf) cycle
			disp = (dp1**2.+dp2**2.+dp3**2.)**0.5
			if(disp>cut_dp*hbc/sf) cycle
	 
     !call wigner_d_CLW(dx,dp,w) 
			w = gc*8.*exp(-disx**2./sf**2.-disp**2.*sf**2./hbc2)
     
			xnd = xnd+w

10      continue			
		   enddo
		   
		enddo
		
		!xndv(nparak) = xnd
		avnd = avnd+xnd
	   
	   enddo
	   write(*,*)"Events:",npara,nparak,"deuteron number:",avnd/float(nparak)
	   enddo

	   avnd = avnd/float(Nevent)**2.
	   !write(*,*)"deuteron:",xndv
	   write(*,*)"njob = ", njob
	   !write(*,*)"dv: ",xndv(1:nevent)
	   write(*,*) "deuteron number:",avnd
	   
	  
	  
101 format((1x,I6),8(1x,F10.5))	  
	  return
	  end	  
	  
	  SUBROUTINE COALhe3(Nevent,njob,xpro,xneu,ipro,ineu,rms,avt)
	  implicit real*8 (a-h,o-z)
	  integer,parameter :: NumPar= 1000, Ntest=1000
	  dimension xpro(8,NumPar,Ntest),xneu(8,NumPar,Ntest)
	  dimension ipro(Ntest),ineu(Ntest)
	  character(50) :: datatype
      character(3) :: tit          ! title   
      character(50):: DirP
      datatype='ana/ampt.dat' 
      DirP = 'ampt-'	  
	  write(tit,"(I3.3)")njob 
      !open (117, file=trim(datatype), status='unknown')
	  avnt = 0.
	  avp = 0.
	  avn = 0.
	  sig1 = rms
	  sig2 = rms
	  hbc = 0.1973
	  hbc2=hbc**2.
	  cut_dr = 5.
	  cut_dp = 5.
	  gc = 1./4.
	  nparak=0
       do npara=1,Nevent-2
	   do npara1=npara+1,Nevent-1
	   do npara2=npara1+1,Nevent
		  xnt=0.
		  nparak=nparak+1
		if(ipro(npara)<2) cycle
        do k1=1,ipro(npara)
			xx1 = xpro(1,k1,npara)
			yy1 = xpro(2,k1,npara)
			zz1 = xpro(3,k1,npara)
			t1 = xpro(4,k1,npara)
			px1 = xpro(5,k1,npara)
			py1=xpro(6,k1,npara)
			pz1 = xpro(7,k1,npara)
			xm1 = xpro(8,k1,npara)
			xE1 = sqrt(xm1**2.+px1**2.+py1**2.+pz1**2.)
			do k2=1,ipro(npara1)
		    xx2 = xpro(1,k2,npara1)
			yy2 = xpro(2,k2,npara1)
			zz2 = xpro(3,k2,npara1)
			t2 = xpro(4,k2,npara1)
			px2 = xpro(5,k2,npara1)
			py2=xpro(6,k2,npara1)
			pz2 = xpro(7,k2,npara1)
			xm2 = xpro(8,k2,npara1)
			xE2 = sqrt(xm2**2.+px2**2.+py2**2.+pz2**2.)
			
		   do k3=1,ineu(npara2)
		    xx3 = xneu(1,k3,npara2)
			yy3 = xneu(2,k3,npara2)
			zz3 = xneu(3,k3,npara2)
			t3 = xneu(4,k3,npara2)
			px3 = xneu(5,k3,npara2)
			py3=xneu(6,k3,npara2)
			pz3 = xneu(7,k3,npara2)
			xm3 = xneu(8,k3,npara2)
			xE3 = sqrt(xm3**2.+px3**2.+py3**2.+pz3**2.)
			
			px=px1+px2+px3
			 py=py1+py2+py3
			 pz=pz1+pz2+pz3
			 xmass=xm1+xm2+xm3
			 !P = (px**2.+py**2.+pz**2.)**0.5
             !E = (P**2.+xmass**2.)**0.5
			  Et = xE1+xE2+xE3	 
			 Pt = (px**2.+py**2.)**0.5
	         PL = pz
	         yitat = 1./2.*log((Et+PL)/(Et-PL))
            if(abs(yitat)>0.5) cycle
			
			
			tmax = max(t3,t1,t2)
			xx1new = xx1+(tmax-t1)*px1/xE1
			yy1new = yy1+(tmax-t1)*py1/xE1
			zz1new = zz1+(tmax-t1)*pz1/xE1
			xx2new = xx2+(tmax-t2)*px2/xE2
			yy2new = yy2+(tmax-t2)*py2/xE2
			zz2new = zz2+(tmax-t2)*pz2/xE2
			xx3new = xx3+(tmax-t3)*px3/xE3
			yy3new = yy3+(tmax-t3)*py3/xE3
			zz3new = zz3+(tmax-t3)*pz3/xE3
			
		    dx1 = (xx1new-xx2new)/2.**0.5
			dx2 = (yy1new-yy2new)/2.**0.5
			dx3 = (zz1new-zz2new)/2.**0.5
			dp1 = (px1-px2)/2.**0.5
			dp2 = (py1-py2)/2.**0.5
			dp3 = (pz1-pz2)/2.**0.5
			
			disx = (dx1**2.+dx2**2.+dx3**2.)**0.5
			disp = (dp1**2.+dp2**2.+dp3**2.)**0.5
			
			dxx1 = (0.5*xx1new+0.5*xx2new-xx3new)*(2./3.)**0.5
			dxx2 = (0.5*yy1new+0.5*yy2new-yy3new)*(2./3.)**0.5
			dxx3 = (0.5*zz1new+0.5*zz2new-zz3new)*(2./3.)**0.5
	 
			dpp1 = (px1+px2-2.*px3)/6.**0.5
			dpp2 = (py1+py2-2.*py3)/6.**0.5
			dpp3 = (pz1+pz2-2.*pz3)/6.**0.5

			disx1 = (dxx1**2.+dxx2**2.+dxx3**2.)**0.5
			!if(disx1 >cut_dr*sig2) cycle
			disp1 = (dpp1**2.+dpp2**2.+dpp3**2.)**0.5
			!if(disp1>cut_dp*hbc/sig2) cycle
	 
			w = gc*8.**2.*exp(-disx**2./sig1**2.-disp**2.*sig1**2./hbc2-disx1**2./sig2**2.-disp1**2.*sig2**2./hbc2)
	 
			
			xnt = xnt+w
		   
		   enddo
		enddo
		enddo
		
		!xntv(nparak) = xnt
		avnt = avnt+xnt
		
	   enddo
	   enddo
	   write(*,*)"Events:",npara,nparak,"helium-3 number:",avnt/float(nparak)
	   enddo

	   avt = avnt/float(nparak)
	   
	   !write(*,*)"deuteron:",xndv
	   write(*,*)"njob = ", njob
	   write(*,*) "helium-3 number:",avt
	  
	  
	  
	  return
	  end
	  
	  
	subroutine lorentzboost(betax,betay,betaz,t,x,y,z,tp,xp,yp,zp)
!       -------------------------
!       loretnz transformation
!       -------------------------
        implicit double precision (a-h,o-z) 
     	implicit integer*4 (i-n)
	tp=t
	xp=x
	yp=y
	zp=z
	beta2 = betax**2+betay**2+betaz**2
	if (beta2.gt.1.0e-5) then 
	gamma = 1./sqrt(1.-beta2)
	xlam00= gamma
	xlam01 = -gamma*betax
	xlam02 = -gamma*betay
	xlam03 = -gamma*betaz
	xlam10 = xlam01
	xlam11 = 1.+(gamma-1.)*betax**2/beta2
	xlam12 = (gamma-1.)*betax*betay/beta2
	xlam13 = (gamma-1.)*betax*betaz/beta2
	xlam20 = xlam02
	xlam21 = xlam12
	xlam22 = 1.+(gamma-1.)*betay**2/beta2
	xlam23 = (gamma-1.)*betay*betaz/beta2
	xlam30 = xlam03
	xlam31 = xlam13
	xlam32 = xlam23
	xlam33 = 1.+(gamma-1.)*betaz**2/beta2
	tp = t*xlam00+x*xlam01+y*xlam02+z*xlam03
	xp = t*xlam10+x*xlam11+y*xlam12+z*xlam13
	yp = t*xlam20+x*xlam21+y*xlam22+z*xlam23
	zp = t*xlam30+x*xlam31+y*xlam32+z*xlam33
	
	else
	endif   
	return
	end
