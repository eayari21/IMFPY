c File name: isd-1992.f (modified from isd-new.f according to the Date/Year)
c Input  file: istart.in
c Output file: Coord.dat (t,x,y,z,vx,vy,vz,p,q)
c current sheet is inculded (see Miyake, 2005, Eq.2)
c HCS by Miyake, tilt(time) Hoeksema 
c tilt varies as observed (Hoeksema) 
c Tilt(time) in subroutine fields starts from 1964
c User input: the Date of observation (e.g. 2005)

c modeling ISD 
c magnetic filed model: see ISD project
c forces: solar gravity, radiation pressure, Lorentz force, Poynting-Robertson force
c dust charge: +4V (this can be changed) /Kimura & Mann, The Astrophys.
c Journal, 499, 5454-462, 1998, Fig.2
c
c (created/modified from nanodust-circ.f - ISD project)
c
c Written by Antal Juhasz (22 March 2013)
c
c Last modified from isd.f (LASP, 04 September, 2013)
c added: solar wind velocity radial dependence w(r)
c added: Q(Date) instead of using a constant potential (Mihaly did it)

c=================================================================c
c  CGS units are used
c 
c The calculations are in Helicentric Equatorial Coordinate System )HCI):
c However, the dust initial coordinates and the results are in 
c Heliocentric Ecliptic Inertial (HAE) coord. system

c                       z ^        .Y
c                         i      .
c               north     i    .
c                         i  .
c                         i.
c                         SUN------------> x 
c
c
c  parameter 'n' : # of particles
c=================================================================c
      program risd_2022

c Last modified by Antal Juhasz (02 May, 2017)
c N_start particles released initially (input parameter) on a CIRCLE
c with keplerian start (in HEI coord. system)
c B_0 is set to 3.5nT !

c only ONE particle size is used (INPUT parameter)

c currents sheet position testing --> OK (see also : tilt.f )
c New CS model (Pei, 2012) (cs-test.f, cs-test.pro) 
c this is good for all tilt angles, while Miyake model is good
c for low tilt<30 deg
c calculation of dust longitude (fid) corrected! 

      parameter (n=250000)

      parameter (ns=1)

      implicit real*8 (a-h,o-z)

	  dimension date(41),potenc(41)

      common /data/  xx(n),yy(n),zz(n),vxx(n),vyy(n),vzz(n),
     *               qq(n),gr0(ns),rg(n),tlife(n),lost(n),pol(n)
      common /tempo/ x,y,z,vx,vy,vz,q
      common /const/ eq,em,pm,bk,rc,comm,gc,c,pi,sqpi,elcc,eqpbk,
     *               omega,dau,prcc1,prcc2,eqgsp4,oxcc1,oxcc2,pj2,pln
      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
      common /servi/ igone,icrash,r1,r2,r3,nochar,mflag
      common /currs/ ecc,hic,pec,secc
      common /plasma/ dion,tion,dec,tec
      common /rzeroparams/ fimin,fimax,eta,maxfun
      common /field/ ex,ey,ez,bx,by,bz,vplasx,vplasy,vplasz
      common /flags/ ifgr,ifl,ifrp
      common /dusttype/ type,f1
      common /start/ deltat,nstart,modstep,kstart
      common /solarwind/ densw,vsw,bsw,fi0
      common /eforce/ fex,fey,fez,fbx,fby,fbz
      common /sun/ rs,sm,oneau,frp,fsgrz,tilt0
      common /startgrid/ ystmax,zstmax
      common /averb/ p,br,bfi,tilt,p0,t0,thetad,theta
      common /force/ fx,fy,fz,fsg


c T_0 determines the Tilt_0 (see cs-2.ps)
c T_0=0 --> tilt=0
c T_0=5.5 --> tilt=90
       print*,'Enter observational date [year] (e.g. 2005) '
       read(*,*) tobs
c Observation date (near Sun) (2022)
c       tobs=2022.
	   print*,''
       print*,'T_obs=',tobs
       print*,'' 
       if(tobs.ge.2011.) tobs=tobs-22.

c output file for coordinates and velocities
       open(unit=55, file='Coord.dat', status='unknown')
       print*,''
       print*,'Output file: Coord.dat (t,x,y,z,vx,vy,vz)'
c       print*,''

c output file for E and B filed components
c       open(unit=65, file='Fields.dat', status='unknown')
c       print*,''
c       print*,'Output file: Fields.dat (t,ex,ey,ez,bx,by,bz)'
c       print*,''

c output file for F_E and F_B force components
c       open(unit=75, file='Forces.dat', status='unknown')
c       print*,''
c       print*,'Output file: Forces.dat (t,fex,fey,fez,fbx,fby,fbz)'
c       print*,''

c output file for C.S. and dust latitudes an "p" B polarity
c       open(unit=85, file='CSlat.dat', status='unknown')
c       print*,''
c       print*,'Output file: CSlat.dat (t,cslat, dustlat, p)'
c       print*,''


c START:
      call startup(tplan)

c read charge.dat (Date,Potential)
c potential is time (date dependent due to I_ph(t) )
	  open(unit=95, file='charge.dat', status='unknown')
	  do 1963 iq=1,41
	   read(95,*) date(iq),potenc(iq)
1963  continue	    
      close(95)
       
	  tdatemin=1970.
	  deltatd=1.
	  
c Parameters:
      lmax=0
      idum=-1

c set init=1 to restart
      init=0
c influx: the number of new particle input/deltat
      influx=0
c kstart is a flag for printing info on the screen only once
      kstart=1

c initial number of lost particles
      ngone=0
      ncrash=0

c Initial conditions for the cloud at t=0
      time=0.

c solar rotation frequency [1/sec] (T=25.88 days, sideral rot.)
      omegasun=2.86533e-6
c B_r=B_fi field at 1AU (3.5 nT ) 1 Tesla = 10^4 Gauss --> 3.5nT=3.5e-5 G
c B_0 is the magnitude of the radial component at 1AU
      b0=3.5e-5
c the CS position and dynamics in case of ISD is not sensitive to the 
c actual value of fi0 therefore fi0 can be arbitrary
      fi0deg=0.
      fi0=fi0deg*pi/180.

c initial ISD velocity [km/s] (CHANGED: earlier it was 26 km/s) 
      visd=27.
c initial distance [AU]
      rstau=80.
c Sun's tilt
      alf=7.25*pi/180.
      xstau=rstau*cos(alf)
c toff [years] time is necessary fo the dust to travel rst distance
c      rstkm=rstau*oneau*1.e-5 
c      toff=rstkm/visd/86400./365.
      xstkm=rstau*oneau*1.e-5 
      toff=xstkm/visd/86400./365.
c calculate the initial time for the Tilt
      t0=tobs-1964.-toff
c      print*,'toff=',toff
      print*,'Particles start in ',t0+1964.
c      print*,'T_0=',t0

c number of particles in X,Z
      nxx=20
      nzz=20
c (x,z) lattice for start [AU]
c      xstmin=-15.+xh0
c      xstmax=15.+xh0
c      xstmin=-25.
      xstmin=-20.
      xstmax=-5.
c      zstmin=-15.+zh0
c      zstmax=15.+zh0
      zstmin=0.
c      zstmax=20.
      zstmax=10.
      
c total number of particles at start
      nstart=nxx*nzz

      npar=nstart

      print*,''
      print*,' Number of particles at start= ',nstart
      print*,''
      print*,' Initial distance of ISD particles:',rstau,' AU'

c      beta=0.1
      print*,''
      print*,'   beta=',beta
      print*,''

c Sun's tilt and (x_HAE - x_HCI) angle
      alf=7.25*pi/180.
c        alf=0.
      if(abs(alf).lt.0.01) print*,'The tilt of Sun is set to Zero !'   
      slamb=76.*pi/180.

c ISD initial coordinates in HAE
c the ISD arriving direction in HAE (longitude=259 deg., latitude=8 deg.)
c the ISD arriving direction in HAE (longitude=259 deg., latitude=5 deg.)
      dustlong=259.*pi/180.
      dustlat=5.*pi/180.
c ISD initial direction [AU]
      rst=rstau*oneau          
      xh0=rstau*cos(dustlat)*cos(dustlong)
      yh0=rstau*cos(dustlat)*sin(dustlong)
      zh0=rstau*sin(dustlat)
      print*,'x0,y0,z0=',xh0,yh0,zh0
c ISD initial velocity [cm/s] in HAE   
      v0=visd*1.e5
      vxh=v0*cos(dustlat)*cos(dustlong-pi)
      vyh=v0*cos(dustlat)*sin(dustlong-pi)
      vzh=-v0*sin(dustlat)
c velocity components in HCI (=HSE)
c HAE=HEI ---> Heliocentric Solar Equatorial (=HGI) transformation
c two rotations 
      vxhsep=vxh*cos(slamb)+vyh*sin(slamb)
      vyhsep=-vxh*sin(slamb)+vyh*cos(slamb)
      vxhse=vxhsep
      vyhse=vzh*sin(alf)+vyhsep*cos(alf)
      vzhse=vzh*cos(alf)-vyhsep*sin(alf)

c initial coordinates, velocities, charges and particle sizes
      i=1
      isize=1
c start from a slab at rst
      do 20 ix=1,nxx
       do 25 iz=1,nzz

        rg(i)=gr0(isize)

        xh=(xstmin+ix*(xstmax-xstmin)/nxx)*oneau
        yh=(yh0-ran2(idum))*oneau
        zh=(zstmin+iz*(zstmax-zstmin)/nzz)*oneau

c HAE=HEI ---> Heliocentric Solar Equatorial (=HGI) transformation
c two rotations 
	xhsep=xh*cos(slamb)+yh*sin(slamb)
	yhsep=-xh*sin(slamb)+yh*cos(slamb)
        xhse=xhsep
        yhse=zh*sin(alf)+yhsep*cos(alf)
        zhse=zh*cos(alf)-yhsep*sin(alf)

	xx(i)=xhse
        yy(i)=yhse
        zz(i)=zhse

        vxx(i)=vxhse
        vyy(i)=vyhse
        vzz(i)=vzhse

c set dut potential to +7 V  --> Q=Pot*rg/300.
	    qq(i)=7.*rg(i)/300.

        i=i+1

c Q/m ratio expressed in (e/m_p)
c	qperm=qq(i)/gm
c        epermp=eq/pm
c        rqm=qperm/epermp
c        rqe=qq(i)/eq
c        rmm=gm/pm
c        if(i.eq.2) print*,'size=',i,rg(i),qq(i),gm,qperm

c        if(i.eq.2) then
c         print*,''
c         print*,'   rg= ',rg(i),gm,qq(i)
c         print*,'' 
c	 print*,'   Q/m= ',rqm,' [e/m_p]',eq
c	 print*,''
c	 print*,'   Q/e= ',rqe,' [# electron]'      
c	 print*,''
c	 print*,'   m_dust/m_proton= ',int(rmm)      
c	 print*,''
c        endif

25    continue 	
20    continue 	



      l=0
      nstep=0
      timest=time
      ncount=0

      par=1.

c INTEGRATION begins
 100  continue
c T_plan is in years
c      tend=time+deltat*(2.-ran2(idum))
      tend=time+deltat
      ncount=ncount+1

c actual date (year)
      tdate=t0+1964.+tend/86400./365.

        if(tend.gt.(365.*86400.*tplan)) then
           print*,'       '
           print*,' Time > T_plan'
           go to 1000
        endif

        times=time
        j=0
        npar=npar+influx

       do 200 i=1,npar

          igone=0
          icrash=0
c     no integration is required for lost particles
            if (l.ne.0) then
               do k=1,l
                if(i.eq.lost(k)) go to 200
               enddo
            endif

         x=xx(i)
         y=yy(i)
         z=zz(i)
         vx=vxx(i)
         vy=vyy(i)
         vz=vzz(i)
         q=qq(i)
c       print*,q/gr*300.,tcq*qq(i)/gr*300.
         gr=rg(i)
         time=times

         gr=gr0(i)

c set up size dependent parameters
c         call dust

c the integration step
         call merson(time,tend)

c print info on the screen only once
	 if(i*ncount.eq.1) then

c calculate the gyrofrequency
          bsun=sqrt(bx**2+by**2+bz**2)
          omg=q*bsun/gm/c
          tg=(2.*pi/omg)/86400.
c Larmor radius
          vpar=(vx*bx+vy*by+vz*bz)/bsun
          vabs=sqrt(vx**2.+vy**2.+vz**2.)
          vperp=sqrt(vabs**2.-vpar**2.) 
          rlarcm=vperp/omg
	  rlarau=rlarcm/oneau
          rlarrs=rlarcm/rs
c particle energy per unit mass
          v=sqrt(vx**2.+vy**2.+vz**2.)
          r=sqrt(x**2.+y**2.+z**2.)
          sintheta=z/r
          ekin=v**2./2.
          epot=-gc*sm*(1.-beta)/r
          eem=-q*b0*omegasun*oneau**2./gm*sintheta

c          print*,'' 
c          print*,'   Dust gyrofrequency=',omg,' [1/sec]  '
c          print*,'   Dust gyroperiod   =',tg,' [day]'
          print*,''
c          print*,'   Dust orbital frequency=',omk,' [1/sec]  '
c          print*,'        orbital period   =',torb,' [day]'
          print*,''
c	  print*,'   ---->  Omega_g/Omega_K=',omg/omk
          print*,''
c          print*,'   Dust Larmor radius=',rlarau,' [AU]'
c          print*,'                      ',rlarrs,' [R_s]' 
          print*,''
          print*,'CS tilt=',asin(abs(sin(tilt)))*180./pi
          print*,'---------------------------------------'
      endif

	 if(igone.eq.0.and.icrash.eq.0) then
           j=i
c           j=j+1
           xx(j)=x
           yy(j)=y
           zz(j)=z
           vxx(j)=vx
           vyy(j)=vy
           vzz(j)=vz
c NEW -----------------------------------           
c the dust potential varies with date and radial distance
c Wang and Richardson, JGR 106, A12, p29401-29407, 2001, Fig.4
c           rxy=sqrt(x**2.+y**2.)
c Mihaly calculates the Q(date,r) matrix with resolution of (delta_r,delta_td)
c I need rxymin, deltar, tdatemin, deltatd
c           ir=(rxy-rxymin)/deltar
            itd=(tdate-tdatemin)/deltatd             
            qdate=potenc(itd)*gr/300. 
            qq(j)=qdate
c           qq(j)=7.*tcq*gr/300.                
c ---------------------------------------
c          qq(j)=q
           
c         print*,'pot=',qq(j)/gr*300.
           rg(j)=gr
           pol(j)=p

	 else
           l=l+1
           lost(l)=i
           tlife(l)=time/86400./365.
c lost particles
          if(igone.eq.1) then
             ngone=ngone+igone
          endif
          if(icrash.eq.1) then
             ncrash=ncrash+icrash
          endif

         endif
 200   continue

       nstep=nstep+1
c       npar=j

       if ((ngone+ncrash).eq.nstart) go to 600

c Print TIME on the screen yearly:
cc deltat=3600 s (1 hour), 1 year=8760 hours
c       modstep=8760
c deltat=86400/2 s (=0.5 day), 1 year=365 days
       modstep=365

       if(mod(nstep,modstep).eq.0) then
        print*,''
        write(*,15) time/86400./365.
c        print*,'vsw=',vsw*1.e-5
       endif
 15    format('Time [years]:',f5.1)

c deltat=0.5 days
c Write into files according to modstep2*delatat=20*0.5 days=10 days 
       modstep2=20
       kiir=modstep2

       if(mod(nstep,kiir).eq.0) then

        
        do 500 k=1,npar

c HCI (=HGI,HSE)---> HAE (=HEI) transformation
        yheip=zz(k)*sin(-alf)+yy(k)*cos(-alf)
        zheip=zz(k)*cos(-alf)-yy(k)*sin(-alf)
	    xhei=xx(k)*cos(-slamb)+yheip*sin(-slamb)
	    yhei=-xx(k)*sin(-slamb)+yheip*cos(-slamb)
        zhei=zheip
  
        vyheip=vzz(k)*sin(-alf)+vyy(k)*cos(-alf)
        vzhei=vzz(k)*cos(-alf)-vyy(k)*sin(-alf)
	    vxhei=vxx(k)*cos(-slamb)+vyheip*sin(-slamb)
	    vyhei=-vxx(k)*sin(-slamb)+vyheip*cos(-slamb)


         xau=xhei/oneau
         yau=yhei/oneau
         zau=zhei/oneau

c in HCI
c         xau=xx(k)/oneau
c         yau=yy(k)/oneau
c         zau=zz(k)/oneau

c Write HAE coordinates and velocities into file: "Coord-Year.dat"
c only for the |y|<3 AU region
c        if(abs(yau).lt.2..and.abs(xau).lt.2..and.abs(zau).lt.2.) then   
          write(55,32) time/86400.,xau,yau,zau,
     *     vxhei*1.e-5,vyhei*1.e-5,vzhei*1.e-5,pol(k),qq(k)/gr*300.
c        endif
c     *     vxhei*1.e-5,vyhei*1.e-5,vzhei*1.e-5,pol(k),qq(k)/gr*300.,k
c     *                thetacs(k)*180./pi,thetadu(k)*180./pi,pol(k)
c dust es C.S. latitudes in HCI
          cslat=90.-theta*180./pi
          dustlat=90.-thetad*180./pi
c         if(nstart.eq.1)  write(85,33) time/86400.,cslat,dustlat,p
c write out Fields and Forces if only Nstart=1
c         if(nstart.eq.1) write(65,31) time/86400.,ex,ey,ez,bx,by,bz
c         if(nstart.eq.1) write(75,32) time/86400.,fex,fey,fez,
c     *                   fbx,fby,fbz,fsg

 500    continue
       endif

c 31    format(7e16.6)
c 32    format(9e16.6,i8)
 32    format(9e16.6)
c 33    format(3e16.6,1f16.1)


 400   continue

       go to 100


c Check for remaining particle number

 1000  continue


       print*,''
 600   print 27,ngone
 27    format('NGONE=',i6)
       print*,''
       print 28,ncrash
 28    format('NCRASH=',i6)


1999  continue 

      print*,' '
      print*,'          Program completed successfully, Bye!'

      close(55)
c      close(65)
c      close(75)
c      close(85)

      STOP

      END



c================================STARTUP================================

      subroutine startup(tplan)

      parameter (n=250000)

      parameter (ns=1)


      implicit real*8 (a-h,o-z)

      dimension a0(ns)
      common /data/  xx(n),yy(n),zz(n),vxx(n),vyy(n),vzz(n),
     *               qq(n),gr0(ns),rg(n),tlife(n),lost(n),pol(n)
      common /const/ eq,em,pm,bk,rc,comm,gc,c,pi,sqpi,elcc,eqpbk,
     *               omega,dau,prcc1,prcc2,eqgsp4,oxcc1,oxcc2,pj2,pln
      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
      common /merso/ acc,h,hmin,jtest,nmer
      common /servi/ igone,icrash,r1,r2,r3,nochar,mflag
      common /flags/ ifgr,ifl,ifrp
      common /dusttype/ type,f1
      common /start/ deltat,nstart,modstep,kstart
      common /rzeroparams/ fimin,fimax,eta,maxfun
      common /sun/ rs,sm,oneau,frp,fsgrz,tilt0

c      data a0/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/

c only ONE particle size is used here
c       data a0/0.005/
c       print*,'Input rg [micron]'
c       read(*,*) rgmicron
c       rgmicron=0.25
       rgmicron=0.3
       a0(1)=rgmicron
       print*,'rg [micron]= ',a0(1)

c to neglect charging set 'nochar=1' (otherwise set nochar=0)
c MODSTEP: output frequency
c Input for start:
      open(unit=35,file='istart.in',status='unknown')
c
c      read(35,*) tpl,ifrp,ifl,nochar
      read(35,*) tpl,ifgr,ifrp,ifl,nochar
      close(35)

c Constants:
      eq=4.803e-10
      em=9.108e-28
      pm=1.660e-24
      gc=6.668e-8
      c=2.998e10
      bk=1.380e-16
      eqpbk=eq/bk
      prcc1=sqrt(2.*bk/pm)
      prcc2=sqrt(2.*eq/pm)
c Solar radius [cm], mass [g] and 1 AU [cm]
      rs=6.96e10
      sm=1.98892e33
      oneau=149598.e8
c for heavy ions (H2O) m=18m_p
      oxcc1=prcc1/sqrt(18.)
      oxcc2=prcc2/sqrt(18.)
      pi=acos(-1.)
      sqpi=sqrt(pi)

c papameters for rzero (equilibrium potential search)
      fimin=-100./300.
      fimax=100./300.
      eta=1.e-4
      maxfun=1000

c integration timestep 0.5 day (1 day = 86400 sec)
      deltat=86400./2.

c Type of grain  1:conducting         2:dielectric
      type=2.

c Grain's size: micron-->cm conversion
      do 10 i=1,ns
 10      gr0(i)=a0(i)*1.e-4
      gr=gr0(1)

c Parameters for secondary current
      delsec=1.5
      eamxse0=400.
      delmax=2.*3.7*delsec
c      delmax=3.7*delsec
      emaxse=emaxse0*1.6e-12
      emaxp4k=emaxse/4./bk
      tsetk=2.5*1.6e-12

c particles will be followed for 'tplan'  years
      tplan=tpl

c set parameters for merson
       nmer=7
       acc=1.e-5
       h=1.
       hmin=1.e-10
       jtest=0

c Photoemission parameters:
      hi=1.0
      if(type.eq.2) hi=0.1
      tp=2.0*1.6e-12
      eqpktp=eq/tp
      f1=2.5e10*hi/dau**2

c set up size dependent quantities:
      call dust

c Print initial info on the scrren:
      if(kstart.eq.1) then
        print*,'           '
        print 20, nstart
 20     format('      N_start=',i5)
        print*,''
        print 21, tpl
 21     format('       T_plan=',f5.1,' [days]')
        print*,''
c        print 22,grmin0
c 22     format('       rg_min=',f5.3,' [micron]')
        print*,''
c        print 23,r0
c 23     format('       R_init=',f4.2,' [R_S]')
        print*,''
        print 24,vej*1.e-5
 24     format('         V_ej=',f4.2,' [km/s]')
        print*,''
        print 25,rkiir
 25     format('        Output=',f4.1,' hourly')
        print*,'           '
c        print*,'           J2=',ij2
c        print*,'     ion drag=',iond
c        print*,' neutral drag=',neud
c        print*,'Lorentz force=',ifl
        print*,'       nochar=',nochar
        print*,'Rad. pressure=',ifrp
c        print*,'   Sputtering=',itsp
        print*,'           '
      endif
c print the above info on the screen only once
      kstart=kstart+1

      return
      end

c================================DUST================================

      subroutine dust

c Set up size dependent parameters
c physical properties of the grain :
c
c    type=1      conducting magnetite
c    type=2      dielectric olivine
c                differences:  photoelectron efficiency
c                              light scattering efficiency
c                              density
c                              secondary yield

      implicit real*8 (a-h,o-z)

      common /const/ eq,em,pm,bk,rc,comm,gc,c,pi,sqpi,elcc,eqpbk,
     *               omega,dau,prcc1,prcc2,eqgsp4,oxcc1,oxcc2,pj2,pln
      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
c      common /poten/ pot,qp,w
      common /sece/ delsec
      common /dusttype/ type,f1

      dimension r(21),ros(21),q1(21),q2(21)

      data r/0.,1.e-5,1.44e-5,1.77e-5,2.04e-5,2.98e-5,4.51e-5,6.63e-5,
     1  1.02e-4,2.36e-4,5.57e-4,1.3e-3,2.94e-3,6.51e-3,1.42e-2,3.08e-2,
     2  6.66e-2,1.44e-1,3.1e-1,6.68e-1,1.0/

      data ros/3.,2.9,2.85,2.82,2.8,2.72,2.59,2.45,2.26,1.8,1.38,1.09,
     1  .94,.86,.83,.81,.8,.8,.8,.8,.8/

      data q1/0.00,1.54,1.88,2.03,2.07,1.97,1.73,1.52,1.38,1.18,1.06,
     1     10*1.0/

      data q2/0.00,0.31,0.59,0.76,0.88,1.05,1.18,1.10,0.96,0.74,0.63,
     1        0.65,0.74,0.87,7*1.0/


      do 1 i=2,18
      if(gr.lt.r(i)) go to 2
 1     continue
      rr=1.0
      go to 3
 2     rr=(gr-r(i-1))/(r(i)-r(i-1))
 3     continue
c     ro=ros(i-1)+(ros(i)-ros(i-1))*rr

c Density of the nanodust grain [g/cm^3]:
      ro=2.5

c size dependent quantities
      gs=4*pi*gr**2
      eqgsp4=eq*gs/4.
      elcc=eq*gr**2*sqrt(8*bk*pi/em)
      gm=4.*pi/3.*gr**3*ro
      gcplmm=gc*comm*gm

c beta (=F_rp/F_grav) calculation (using Landgraf's curve)
      call betacalc
c      beta=1.
C Q_pr determination:
      if(i.ne.2) then
        qpr=q1(i-1)+(q1(i)-q1(i-1))*rr
        if(type.eq.2) qpr=q2(i-1)+(q2(i)-q2(i-1))*rr
      else
        dr=r(2)-gr
        qpr=(q1(2)-0.3)*exp(-dr*1.e5)+0.3
        if(type.eq.2) qpr=q2(2)*exp(-dr*5.5e5)
      endif

c Constant for the light pressure force (at 1 AU):
      dau=1.
      clpf=pi*gr**2/c*1.36e6*qpr/dau**2

c Photoelectron emission:
      cphe=gs/4.*eq*f1


      return
      end

c==============================MERSON================================

      subroutine merson (x,xend)
      implicit real*8 (a-h,o-z)
      common /merso/  acc,h,hmin,jtest,n
      common /local/  yz(7) , a(7) , b(7)
      common /solut/  w(7)
      common /prime/  f(7)
      common /tempo/  y(7)
      common /servi/ igone,icrash,r1,r2,r3,nochar,mflag
      logical ok
c
c     rzero is a number with a magnitude of order equal to the noise
c     level of the machine i.e. in the range of the rounding off errors.
c
      data rzero / 1.e-13 /
c
      ok=.true.
c
c     store internally parameters in list
c
      nn=n
      do 1 k=1,nn
    1 w(k)=y(k)
      z=x
      zend=xend
      bcc=acc
      zmin=hmin
      itest=jtest
      s   =h
      iswh=0
c
    2 hsv=s
      cof=zend-z
      if (abs(s).lt.abs(cof)) go to 8
      s=cof
      if (abs(cof/hsv).lt.rzero) go to 50
      iswh=1
c
c     if iswh=1 then  s is equal to maximum possible steplength
c     within the remaining part of the domain of integration.
c
    8 do 10 k=1,nn
   10 yz(k)=w(k)
   12 ht=.3333333333333*s
c
      call diff(z)
c
      z=z+ht
      do 20 k=1,nn
      a(k)=ht*f(k)
   20 w(k)=a(k)+yz(k)
c
      call diff(z)
c
      do 22 k=1,nn
      a(k)=.5*a(k)
   22 w(k)=.5*ht*f(k)+a(k)+yz(k)
c
      call diff(z)
c
      z=z+.5*ht
      do 24 k=1,nn
      b(k)=4.5*ht*f(k)
   24 w(k)=.25*b(k)+.75*a(k)+yz(k)
c
      call diff(z)
c
      z=z+.5*s
      do 26 k=1,nn
      a(k)=2.*ht*f(k)+a(k)
   26 w(k)=3.*a(k)-b(k)+yz(k)
c
      call diff(z)
c
      do 28 k=1,nn
      b(k)=-.5*ht*f(k)-b(k)+2.*a(k)
      w(k)=w(k)-b(k)
      a(k)=abs(5.*bcc*w(k))
      b(k)=abs(b(k))
      if (abs(w(k)).le.rzero) go to 28
      if(b(k).gt.a(k)) go to 60
   28 continue
c
c     required accuracy obtained for all computed function values.
c
      if (iswh.eq.1) go to 50
c
c     test if steplength doubling is possible.
c
   40 do 42 k=1,nn
      if(b(k).gt. .03125*a(k)) go to 2
   42 continue
      s=s+s
      go to 2
c
c     calculation finished.replace input function values with the func-
c     tion values computed for the output point xend.replace input step-
c     length h with new computed steplength.
c
   50 h=hsv
      x=z
      do 52 k=1,nn
   52 y(k)=w(k)
c
      return
c
c     required accuracy not obtained.
c
c
   60 cof=.5*s
      if (abs(cof).ge.zmin) go to 80
      if (itest.eq.0) go to 84
c
c     jtest=1,continue with constant steplength equal hmin.
c
      s=zmin
      if (hsv.lt.0.) s=-s
      if (iswh.eq.1) go to 50
      go to 2
c
c     do calculations related to halving of steplength.
c
   80 do 82 k=1,nn
   82 w(k)=yz(k)
      z=z-s
      s=cof

      go to 2
c
c     jtest=0 and abs(s).lt.hmin.print error message,set ok=.false. and
c     return to calling program.
c
   84 print 88 , itest , s , zmin , z
      ok=.false.
      go to 50
c
   86 print 90 , n
      stop
c
   88 format(//,5x,31h*** subroutine merson error ***,2x,8hjtest = ,i2,
     12x,4hh = ,e12.5,2x,7hhmin = ,e12.5,2x,4hx = ,e12.5,//)
   90 format(//,5x,31h*** subroutine merson error ***,4h  n=,i4,1x,55hgr
     1eater than the maximum number of equations permitted.,//)
c
      end

c===============================DIFF===================================

      subroutine diff(t)

c      parameter (nx=160, ny=160, nz=160)

      implicit real*8 (a-h,o-z)
      common /solut/  x,y,z,vx,vy,vz,q
      common /prime/  xp,yp,zp,vxp,vyp,vzp,qp
      common /force/ fx,fy,fz,fsg
      common /grain/  cphe,clpf,eqpktp,gr,gs,gm,
     *                delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
      common /servi/ igone,icrash,r1,r2,r3,nochar,mflag
      common /rzeroparams/ fimin,fimax,eta,maxfun
      common /plasma/ dion,tion,dec,tec
      common /potenc/ pot
      common /averb/ p,br,bfi,tilt,p0,t0,thetad,theta
      common /field/ ex,ey,ez,bx,by,bz,vplasx,vplasy,vplasz

c      external qzero

      if(nochar.eq.0) call fields(t)


      call forces

      if(icrash.eq.1.or.igone.eq.1) go to 100
      xp=vx
      yp=vy
      zp=vz
      vxp=fx/gm
      vyp=fy/gm
      vzp=fz/gm
      if(nochar.eq.1) then
       qp=0.
c calculates with constant potential !
      else
        qp=0.
c           print*,tcq*qq(i) 
c           qq(j)=q*tcq             
c        call rzero(fimin,fimax,pot,fv,eta,maxfun,qzero)
c	q=pot*gr
c       else
c        call currents(qp)
c       endif
      endif
      return
100   xp=0.
      yp=0.
      zp=0.
      vxp=0.
      vzp=0.
      qp=0.

      return
      end


c=============================FORCES================================

      subroutine forces

      implicit real*8 (a-h,o-z)

      common /const/ eq,em,pm,bk,rc,comm,gc,c,pi,sqpi,elcc,eqpbk,
     *               omega,dau,prcc1,prcc2,eqgsp4,oxcc1,oxcc2,pj2,pln
      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
      common /solut/ x,y,z,vx,vy,vz,q
      common /field/ ex,ey,ez,bx,by,bz,vplasx,vplasy,vplasz
      common /servi/ igone,icrash,r1,r2,r3,nochar,mflag
      common /force/ fx,fy,fz,fsg
      common /eforce/ fex,fey,fez,fbx,fby,fbz
c      common /plasma/ dion,tion,dec,tec
      common /flags/ ifgr,ifl,ifrp
      common /startgrid/ ystmax,zstmax
      common /sun/ rs,sm,oneau,frp,fsgrz,tilt0

c cevt: eV---> K conversion factor
c      data cevt/11605./

      r2=x**2.+y**2.+z**2.
      rxy=sqrt(x**2.+y**2.)
      r=sqrt(r2)
      r3=r1**3.

      icrash=0
c crash=1 if close to Sun
      if(r.gt.2.5*rs) go to 100
       icrash=1
      return
 100   continue

      igone=0.
c igone=1 if lost
c cm--> AU conversion first
      xau=x/oneau
      axau=abs(x)/oneau
      ayau=abs(y)/oneau
      azau=abs(z)/oneau
      rau=r/oneau
      rxyau=rxy/oneau
	
c      print*,'x,y,z=',axau,ayau,azau

c      if((xau.gt.-10.).and.
c     *  (ayau.le.ystmax).and.(azau.le.zmaxst)) go to 200
c      igone=1
c      print*,'GONE'
c      return

      if(yau.gt.10.) igone=1
c      print*,igone
       
 200  continue

c-------------------------------------------------------------
c Sun's gravity and light pressure force:
      fsgr=gc*sm*gm/r2
      fsgrx=-fsgr*(1.-ifrp*beta)*x/r
      fsgry=-fsgr*(1.-ifrp*beta)*y/r
      fsgrz=-fsgr*(1.-ifrp*beta)*z/r
      fsg=sqrt(fsgrx**2.+fsgry**2.+fsgrz**2.)

      fx=ifgr*fsgrx
      fy=ifgr*fsgry
      fz=ifgr*fsgrz

c Poynting-Robertson force (see Krivov, Icarus 134,1998, page 314, Eq.4)
      vr=(x*vx+y*vy+z*vz)/r
      fprx=-fsgr*beta*(vr*x/r+vx)/c
      fpry=-fsgr*beta*(vr*y/r+vy)/c
      fprz=-fsgr*beta*(vr*z/r+vz)/c

      fx=fx+ifrp*fprx
      fy=fy+ifrp*fpry
      fz=fz+ifrp*fprz
	

c------------------------------------
c Lorentz Force   q(e + v/c x b)
      if(nochar.eq.1) return

       fex=q*ex
       fey=q*ey
       fez=q*ez	
       fbx= q*(vy*bz-vz*by)/c
       fby= q*(vz*bx-vx*bz)/c
       fbz= q*(vx*by-vy*bx)/c

       fxl=fex+fbx
       fyl=fey+fby
       fzl=fez+fbz

c for testing
c       fxl=fex
c       fyl=fey
c       fzl=fez 
 	
       fl=sqrt(fxl**2.+fyl**2.+fzl**2.)

      fx=fx+ifl*fxl
      fy=fy+ifl*fyl
      fz=fz+ifl*fzl

c      print*,'r,Frgm,Frp,FL=',r,fsgr,frpx,fl

      return
      end



c===========================QZERO====================================

      function qzero(pot,k)

      implicit real*8 (a-h,o-z)

      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
      common /solut/ x,y,z,vx,vy,vz,q

      q=gr*pot
      call currents(qp)
      qzero=qp

      return
      end

c============================CURRENTS===============================

      subroutine currents(qp)

      implicit real*8 (a-h,o-z)
      common /const/ eq,em,pm,bk,rc,comm,gc,c,pi,sqpi,elcc,eqpbk,
     *               omega,dau,prcc1,prcc2,eqgsp4,oxcc1,oxcc2,pj2,pln
      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta
      common /currs/ ecc,hic,pec,secc
      common /field/ ex,ey,ez,bx,by,bz,vplasx,vplasy,vplasz
      common /plasma/ dion,tion,dec,tec
      common /solut/ x,y,z,vx,vy,vz,q
      common /servi/ igone,icrash,r1,r2,r3,nochar,mflag
      common /potenc/ pot

c      if(nochar.eq.0) call plasma_mhd
c      if(nochar.eq.0) call fields

       pot=q/gr

c (vplasx,vplasy,vplasz) comes from subroutine fields
       v=sqrt((vx-vplasx)**2+(vy-vplasy)**2+(vz-vplasz)**2)

c -------------------------- negative potential
      if(pot.le.0.) then

c cold electron current
	hec=eqpbk*pot/tec
        ec0c=dec*elcc*sqrt(tec)
        ecc=ec0c*exxp(hec)

c photoelectron current
        pec=cphe

c heavy ion current
        alfa=oxcc1*sqrt(tion)
        if(v.lt.alfa*1.e-3) v=alfa*1.e-3
        vpalfa=v/alfa
        alfapv=alfa/v
        u=oxcc2*sqrt(abs(pot))
        hic=dion*eqgsp4*v*((1.+0.5*alfapv**2+(u/v)**2)*erf1(vpalfa)+
     *      alfapv/sqpi*exp(-vpalfa**2))

c secondary electrons induced by cold electrons
        b1c =emaxp4k/tec
        secc=delmax*ec0c*exxp(hec)*f5(b1c)

      else
c --------------------------- positive potential

c cold electron current
        hec=eqpbk*pot/tec
        ec0c=dec*elcc*sqrt(tec)
        ecc=ec0c*(1+hec)

c photoelectron current
        pec=cphe*exxp(-eqpktp*pot)

c modified heavy ion currents
        alfa=oxcc1*sqrt(tion)
        if(v.lt.alfa*1.e-3) v=alfa*1.e-3
        vpalfa=v/alfa
        alfapv=alfa/v
        u=oxcc2*sqrt(abs(pot))
        hic=dion*eqgsp4*0.5*v*((1.+0.5*alfapv**2-(u/v)**2)*
     *         (erf1((v+u)/alfa)+erf1((v-u)/alfa))+
     *   alfapv/sqpi*((u/v+1)*exp(-((v-u)/alfa)**2)-
     *                (u/v-1)*exp(-((v+u)/alfa)**2)))


c secondary electrons induced by cold electrons
        arg =  -eq*pot/tsetk
        hs    =   (1.0-arg)*exxp(arg)
        b1c =emaxp4k/tec
        b10c  = sqrt(hec/b1c)
        secc=delmax*ec0c*exxp(hec)*fb5(b10c,b1c)*hs

      endif

      qp=hic-ecc+pec+secc

      return
      end


c=============================F5(x)=================================

      function f5(xx)
      implicit real*8 (a-h,o-z)
c to calculate the f5 function as defined in :
c Meyer-Vernet: Astron. Astrophys. 105,98-106,1982
c
      data pi/3.1415926/
      x=xx
      if(x.le.1.e-10) x=1.e-10

      f1=pi/x
      f2=0.25/x
      f3=0.5/sqrt(x)
      f0=0.5*x**2*sqrt(f1)*exxp(f2)*erfc1(f3)
      f1=  (x-f0/x)*.5
      f2=   (f0-f1)*.5/x
      f3= (2*f1-f2)*.5/x
      f4= (3*f2-f3)*.5/x
      f5= (4*f3-f4)*.5/x
      if(f5.le.0.) f5=0.0
      return
      end


c=============================F5B(x)================================

      function fb5(b,x)
      implicit real*8 (a-h,o-z)
c to calculate the fb5 function as defined in :
c n.meyer-vernet: astron.astrophys.105,98-106,1982
      data pi/3.14159654/,sqpi/1.772453851/

      if(x.le.1.e-10) x=1.e-10
      f1=sqrt(x)
      f2=0.5/f1
      f3=b*f1+f2
      f4=0.25/x
      fb5=erfc1(f3)
      f0=fb5*0.5*x*x*sqpi/f1*exxp(f4)
      fexxp=-x*b*b-b
      fexxp=exxp(fexxp)
      f1=          (x*fexxp-f0/x)*.5
      f2=         (f0-f1+x**2*b*fexxp)*.5/x
      f3=     (2*f1-f2+x**2*b*b*fexxp)*.5/x
      f4=   (3*f2-f3+x**2*b*b*b*fexxp)*.5/x
      fb5=(4*f3-f4+x**2*b*b*b*b*fexxp)*.5/x
      if(fb5.le.0.) fb5=0.0
      return
      end

c=============================ERFC1(x)================================

      function erfc1(x)
      implicit real*8 (a-h,o-z)

      z=abs(x)
      t=1./(1.+0.5*z)
      arg=-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     *    t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     *    t*(1.48851587+t*(-.82215223+t*.17087277))))))))
      erfc1=t*exxp(arg)
      if (x.lt.0.) erfc1=2.-erfc1
      return
      end

c==============================ERF1(x)================================

      function erf1(x)
      implicit real*8 (a-h,o-z)

      z=abs(x)
      t=1./(1.+0.5*z)
      arg=-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     *    t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     *    t*(1.48851587+t*(-.82215223+t*.17087277))))))))
      erfc1=t*exxp(arg)
      if (x.lt.0.) erfc1=2.-erfc1
      erf1=1-erfc1
      return
      end

c=============================EXXP(x)===============================

      function exxp(x)
      implicit real*8 (a-h,o-z)
      arg=x
      xx=abs(x)
      if(xx.gt.70.) arg=x/xx*70.
      exxp=exp(arg)
      return
      end

c==============================RAN2=================================



      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
 11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software K#+%.


c==============================RZERO=================================

      SUBROUTINE RZERO(A,B,X,R,ETA,MAXFUN,FCN)
      implicit real*8 (a-h,o-z)
      external fcn
C
      EPSI=ETA
      IF(EPSI.LE.1.E-6 ) EPSI=1.E-6
      FLOW=1.E30
      E=1.
      MC=0
c      XA=AMIN1(A,B)
      XA=DMIN1(A,B)
c      XB=AMAX1(A,B)
      XB=DMAX1(A,B)
      I=1
      FA=FCN(XA,I)
      MC=MC+1
      I=2
      FB=FCN(XB,I)
      IF(FA*FB.GT.0.) GO TO 16
      MC=MC+1
C
    4 X=0.5*(XA+XB)
      R=X-XA
      EE=ABS(X)+E
      IF(R.LE.EE*EPSI) GO TO 18
      F1=FA
      X1=XA
      F2=FB
      X2=XB
    1 CONTINUE
      MC=MC+1
      IF(MC.GT.MAXFUN) GO TO 17
      FX=FCN(X,I)
C
      IF(FX*FA.GT.0) GO TO 2
      FB=FX
      XB=X
      GO TO 3
    2 XA=X
      FA=FX
    3 CONTINUE
C
C     PARABOLA ITERATION
C
      F3=FX
      X3=X
      IF(ABS(F1-F2).GE.FLOW*ABS(X1-X2)) GO TO 4
      U1=(F1-F2)/(X1-X2)
      IF(ABS(F2-FX).GE.FLOW*ABS(X2-X)) GO TO 4
      U2=(F2-FX)/(X2-X)
      CA=U1-U2
      CB=(X1+X2)*U2-(X2+X)*U1
      CC=(X1-X)*F1-X1*(CA*X1+CB)
      IF(ABS(CB).GE.FLOW*ABS(CA)) GO TO 8
      U3=0.5*CB/CA
      IF(ABS(CC).GE.FLOW*ABS(CA)) GO TO 4
      U4=U3**2-CC/CA
      IF(U4.LT.0.) GO TO4
      U5=SQRT(U4)
      IF(X.GE.-U3) GO TO 10
      X=-U3-U5
      GO TO9
   10 X=-U3+U5
      GO TO 9
    8 IF(ABS(CC).GE.FLOW*ABS(CB)) GO TO 4
      X=-CC/CB
    9 CONTINUE
      IF(X.LT.XA) GO TO 4
      IF(X.GT.XB) GO TO 4
C
C     TEST FOR OUTPUT
C
      R=ABS(X-X3)
      R1=ABS(X-X2)
      IF(R.GT.R1) R=R1
      EE=ABS(X)+E
      IF(R/EE.GT.EPSI) GO TO 5
      MC=MC+1
      IF(MC.GT.MAXFUN) GO TO 17
      FX=FCN(X,I)
      IF(FX.EQ.0.) GO TO 18
      IF(FX*FA.LT.0.) GO TO 7
      XX=X+EPSI*EE
      IF(XX.GE.XB) GO TO 18
      MC=MC+1
      IF(MC.GT.MAXFUN) GO TO 17
      FF=FCN(XX,I)
      FA=FF
      XA=XX
      GO TO 6
    7 XX=X-EPSI*EE
      IF(XX.LE.XA) GO TO 18
      MC=MC+1
      FF=FCN(XX,I)
      FB=FF
      XB=XX
    6 IF(FX*FF.GT.0.) GO TO 14
   18 CONTINUE
      R=EPSI*EE
      I=3
      FF=FCN(X,I)
      RETURN
   14 F1=F3
      X1=X3
      F2=FX
      X2=X
      X=XX
      FX=FF
      GO TO 3
C
    5 CONTINUE
      F1=F2
      X1=X2
      F2=F3
      X2=X3
      GO TO 1
C
   16 PRINT 301
  301 FORMAT(5X,51HRZERO    FCN(A,I)  AND FCN(B,I)  HAVE THE SAME SIGN )
      R=-2.*(XB-XA)
      RETURN
C
   17 PRINT 300,MC
  300 FORMAT(10X,7HRZERO  ,I5,22H CALLS OF THE FUNCTION/10X,19HLINE LIMI
     1T EXCEEDED///)
      R=-0.5*ABS(XB-XA)
      RETURN
      END


c------------------------------------------------------------------------
      subroutine fields(t)


      implicit real*8 (a-h,o-z)

      common /solut/ x,y,z,vx,vy,vz,q
      common /sun/ rs,sm,oneau,frp,fsgrz,tilt0
      common /averb/ p,br,bfi,tilt,p0,t0,thetad,theta
      common /const/ eq,em,pm,bk,rc,comm,gc,c,pi,sqpi,elcc,eqpbk,
     *               omega,dau,prcc1,prcc2,eqgsp4,oxcc1,oxcc2,pj2,pln
      common /field/ ex,ey,ez,bx,by,bz,vplasx,vplasy,vplasz
      common /solarwind/ densw,vsw,bsw,fi0

 
c solar wind plasma velocity [cm/s]
      r=sqrt(x**2.+y**2.+z**2.)
      rxy=sqrt(x**2.+y**2.)
      
      sinfi=y/rxy
      cosfi=x/rxy
      sinth=rxy/r
      costh=z/r

c NEW --------------------------------
c  radial dependent solar wind velocity [cm/s]
c       ccc=1.1
c      print*,rxy/oneau
c      vr=(410.-r/oneau)*1.e5
      vr=400.e5
c we replace the constant velocity with a radial dependent one
c Mihaly will give me the F(rxy) function
c      vr=F(rxy)       
c ------------------------------------      
      vth=0.
      vfi=0.
c polar ---> cartesian coord transformation
      vswx=vr*sinth*cosfi+vth*costh*cosfi-vfi*sinfi
      vswy=vr*sinth*sinfi+vth*costh*sinfi+vfi*cosfi
      vswz=vr*costh-vth*sinth

      vplasx=vswx
      vplasy=vswy
      vplasz=vswz
      vsw=sqrt(vswx**2+vswy**2+vswz**2)
      
c solar rotation (sideral time=25.38 days) frequency [1/sec]
      omegasun=2.86533e-6
c B_r=B_fi field at 1AU (3.5 nT ) 1 Tesla = 10^4 Gauss --> 3.5nT=3e-5 G
c B_0 is the magnitude of the radial component at 1AU
      b0=3.5e-5
c radius of the Sun [km]
c      rs=6.96e5
c radius of the source surface rss=2.5*Rs (WSO B field model)
       rss=2.5*rs
 
c solar cycle (22 years)
      tsc=22.*365.*86400.
      omegasc=2.*pi/tsc
c half solar cycle
      t11=11.*365.*86400.
    
c T_0 determines the tilt at t=0
      t0sec=t0*365.*86400.

      r0=oneau

c solrar magnetic field [in Gauss]
c conversions: 1 Tesla=10^4 Gauss, 1 gamma=1 nT= 10^-5 Gauss
c for now we use Parker field (see Burger, 2005)

      cthetad=z/r
      sthetad=rxy/r
      thetad=acos(cthetad)
      fid=atan2(y,x)
      if(fid.lt.0.) fid=fid+2.*pi

c incorrect:  fid=asin(y/rxy)

c the CS position and dynamics is sensitive to the actual value of fi0
c therefore fi0 should be determined carefully!
c for every start date there is a corresponding fi0 value
c that is fi0=fi0(date) function !
c      fi0=0.

c the solar wind travels about tv=4 days to reach Earth (r=1AU)
       tv=(r-rss)/vsw

c tilt time variation 4. version
       rm1=pi/8.
       rm2=pi/14.
       b2=pi/2.-4.*rm2
       b3=pi-11.*rm1
       b4=3.*pi/2.-15.*rm2
       b5=2.*pi-22.*rm1
       b6=5.*pi/2.-26*rm2
       b7=3.*pi-33.*rm1
       b8=7.*pi/2.-37.*rm2
       b9=4.*pi-44.*rm1
       b10=9.*pi/2.-48.*rm2

c       tyear=t/365./86400.
       
       tt=t0+(t-tv)/365./86400.
       tyear=tt
         
c tyear=0  --> T_date=1964
c tyear=48 --> T_date=2012
       if(tyear.ge.0.and.tyear.lt.4.) then
        tilt=rm1*tt
        p0=-1.
       endif
       if(tyear.ge.4.and.tyear.lt.11.) then
        tilt=rm2*tt+b2
        p0=1.
       endif
       if(tyear.ge.11.and.tyear.lt.15.) then
        tilt=rm1*tt+b3
        p0=1.
       endif
       if(tyear.ge.15.and.tyear.lt.22.) then
        tilt=rm2*tt+b4   
        p0=-1.
       endif
       if(tyear.ge.22.and.tyear.lt.26.) then
        tilt=rm1*tt+b5   
        p0=-1.
       endif
       if(tyear.ge.26.and.tyear.lt.33.) then
        tilt=rm2*tt+b6   
        p0=1.
       endif
       if(tyear.ge.33.and.tyear.lt.37.) then
        tilt=rm1*tt+b7   
        p0=1.
       endif
       if(tyear.ge.37.and.tyear.lt.44.) then
        tilt=rm2*tt+b8   
        p0=-1.
       endif
       if(tyear.ge.44..and.tyear.lt.48.) then
        tilt=rm1*tt+b9
        p0=-1.
       endif
       if(tyear.ge.48..and.tyear.lt.55.) then
        tilt=rm2*tt+b10
        p0=1.
       endif

c the Current Sheet position (theta) from Miyake (2005) 
c theta is measured from the z axis !
c this CS model is good for small tilts (tilt<30 deg)
       sintilt=abs(sin(tilt))
       fi00=fid-fi0-omegasun*t+tv*omegasun
       sinfi00=sin(fi00)
c       theta=pi/2.-asin(sintilt*sinfi00)
c HCS position from Pei [2012] gives the same result as Miyake
c for small tilts but this is the correct formula for all tilt angles! 
       tgtilt=abs(tan(tilt))
       theta=pi/2.-atan(tgtilt*sinfi00)

c if particle above the C.S. then p=p0 below it p=-p0
       if(thetad.lt.theta) p=p0
       if(thetad.gt.theta) p=-p0

c       print*,tyear,dasin(szog1)*180./pi,dasin(szog2)*180./pi
c       print*,1964.+tyear,theta*180./pi,asin(szog1)*180./pi

       br=p*b0*(r0/r)**2.
c B_fi from Burger, Adv. Space Res. 35. p637, Eq(1), 2005
      if(r.gt.rss) then
       bfi=-p*b0*(r0/r)**2.*omegasun*(r-rss)*sthetad/vsw
      else
       bfi=0.
      endif
      bth=0.

c polar ---> cartesian coord transformation
      bx=br*sinth*cosfi+bth*costh*cosfi-bfi*sinfi
      by=br*sinth*sinfi+bth*costh*sinfi+bfi*cosfi
      bz=br*costh-bth*sinth


c electric field (convective, E=-v x B/c)
      ex=-(vplasy*bz-vplasz*by)/c
      ey=-(vplasz*bx-vplasx*bz)/c
      ez=-(vplasx*by-vplasy*bx)/c

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine betacalc


      implicit real*8 (a-h,o-z)

      common /grain/ cphe,clpf,eqpktp,gr,gs,gm,
     *               delmax,emaxse,emaxp4k,tsetk,gcplmm,beta

      parameter (n=100)

      dimension tm(n),pbeta(n)


c open input data file (beta-2.dat)
	open(35,file='beta-2.dat',status='old')

c read data file
        do 50 i=1,n
 50	   read(35,*,end=55) tm(i),pbeta(i)
 55	   continue
        close(35)

	if(gm.lt.1.02e-17) then
	 beta=0.1
        else
c find the appropriate mass bin
        do 100 k=1,n-1
 	   if(gm.gt.tm(k).and.gm.lt.tm(k+1)) then
            kfound=k
            go to 200
           endif
 100	continue
   
        
 200	rr=(gm-tm(kfound))/(tm(kfound+1)-tm(kfound))

        beta=pbeta(kfound)+rr*(pbeta(kfound+1)-pbeta(kfound)) 
c 	print*,'size, mass, beta=',gr,gm,beta
        endif

	return
	END

c============================= END ====================================



