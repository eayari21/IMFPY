c
c     MAIN CODE
c
c-------------------------------------------------------------------------
c This code requires the iput file file 'in.dat ' and generates the output file 'Out.dat'
c
c     Variables:
c     Ud is the current position, stored instantaneously as x,y,z 
c     Vd is the current velocity 
c     t  is the current time in the simulation
c     dm is dust mass 
c     d is solar distance in [AU]
c     c is speed of light
c     Sd is the solar constant
c
c     B0, far-field Magnetic field oriented along x 
c     v0, plasma speed along z (a positive value)
c     E0 = 1/c B0 v0 far field Electric field oriented along y   
c     Ex, Ey, Ez are the electric field in x, y and z
c     z0, distance between heliopause nucleus and and apex
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
c   
      common /variables/ c,  B0, Z0, v0
      common /pos/Ux,Uy,Uz
      common /vel/Vx,Vy,Vz
      common /C1C2/ C1,C2
      common /solution_status/ Istatus, JS
      common /Solution_data/ t(300000),ud(300000,3),vd(300000,3)
c      
      character*20 material
      dimension yd(10)
      
c
      open(unit=5,file='in.dat', status='unknown')
      open(unit=6,file='out.dat', status='unknown')
c
c initialise solution step to 0 and solution status to 0
      Js=0      
      Istatus=0
c
c     read name of material type, density and diameter
      read(5,*) material,rho, a
      read(5,*) c,d,Sd,B0,z0,v0
c read initial time, initial position and initial velocity     
      read(5,*) t(1),(Ud(1,i),i=1,3),(Vd(1,j),j=1,3)
c
      Ux=ud(1,1)
      Uy=ud(1,2)
      Uz=ud(1,3)   
      Vx=Vd(1,1)
      Vy=Vd(1,2)
      Vz=Vd(1,3)  
c Output results one axis variable at a time 
      write(6,*)'------I N I T I A L   C O N D I T I O N S-------'
      write(6,10)
      write(6,20)t(1),Vx,Ux,Vy,Uy,Vz,Uz
10    format(5x,'time',9x,'Vx',9x,'Ux',9x,'Vy',9x,'Uy',9x,'Vz',9x,'Uz')
20    format(7(1pe12.5))

c     
c read simulation time step length and the number of simulation steps
      read(5,*) simulation_time_step
      write(6,30)simulation_time_step
      read(5,*) N_simulation_steps
30    format('Simulation time step:',1pe12.5)
      write(6,40)N_simulation_steps
40    format('Number of simulation steps:',I8)
c
c compute mass     
      dm=rho * 3.14159 * a**3 / 6.d0
      write(6,60) dm
60    format('Particle Mass:',1pe12.5)
c compute scattering coefficient, Qpr      
      call Qprf(Qpr,a,material)
c compute Chi
      call Chif(Chi,a,material)
c compute coefficients C1 and C2 for equation of motion      
      C1=3.14159*a**2/c/d**2*Qpr*Sd/dm
      C2=C1*Chi
c Output Qpr and Chi values
      write(6,70) Qpr,Chi
70    format('Qpr:', 1pe12.5,'  Chi: ', 1pe12.5) 
      write(6,80 )C1,C2
80    format(' C1: ',1pe12.5, '  C2: ', 1pe12.5)
c
c ### Begin simulation Loop ###
      do Js=1,N_simulation_steps 
       t(Js+1)=t(Js) + simulation_time_step
       do i=1,3
        Yd(I)  =Vd(Js,i)
        YD(I+3)=UD(Js,i)
       end do
       call get_A(Yd)
       do i=1,3
        Ud(Js+1,i)=0.5*Yd(i)*t(js+1)**2+Vd(1,i)*t(js+1)+Ud(1,i)
        Vd(Js+1,i)=Yd(i)*t(js+1)+vd(1,i)       
       end do  
      end do
c### End simulation loop.
c
c  Output results and exit  
      call exit_simulation
c
      stop
      end
c 
c-----------------------------------------------------------------------------------            
      subroutine get_A(Y)
c-----------------------------------------------------------------------------------            
      implicit real*8 (a-h,o-z)   
      common /Efield/  Ex, Ey, Ez, Emod
      common /C1C2/ C1,C2
c
      dimension Y(10)
c compute Ex Ey and Ez      
      call Ef(Y(4),y(5),y(6))
c      
      Y(1)= C2*Ex 
      Y(2)= C2*Ey
      Y(3)= C2*Ez -C1
c
c      
      return
      end
c
      

               
c***************************************************************************
      Subroutine Ef(ux,uy,uz)
c***************************************************************************
c
c Generate the electric field at a point (ux,uy,uz)
c
c INPUT: ux,uy,uz and c: speed of light, B0: the far-field magnetic strength B0, and 
c Z0, the distance between nucleus and ionopause apex
c OUTPUT: Ex, Ey and Ez, and the modulus of E
c
      implicit real*8 (A-h,o-z)
c
      common /variables/ c,  B0, Z0, v0
      common /Efield/     Ex,Ey,Ez,Emod
      common /solution_status/ Istatus, JS
c
      r=dsqrt(ux**2.+uy**2.)
      r1=dsqrt(r**2.+uz**2.)
      S=(r**2.+2.*z0**2*(uz/r1-1.))
c  check if the particle is inside the IONOPAUSE; Exit when it is the case.
      if(S.le.0.0) then 
       Istatus=1
       call exit_simulation
      else
       S=dsqrt(S)
       S2=-S/r**2+1./S-z0**2*uz/r1**(3.)/S
       Coef=B0*v0/c/r
c
       Ex=coef*S2*ux*uy
       Ey=coef*(S2*uy**2.+S)
       Ez=coef*z0**2*r**2/r1**(3.) /S *uy
c
       Emod= (Ex**2+ Ey**2+ Ez**2)**(0.5)
c       Emod2=Coef*sqrt((r**4*(1/S-z0**2*z/r1**3/S)**2-S**2)*(y/r)**2
c     *      +S**2+ (z0**2*r**2/r1**3/S)**2)
c
      endif
c     
      Return
      end
c------------------------------------------------------------------------
      subroutine exit_simulation
      implicit real*8 (a-h,o-z)
c
      common /solution_status/ Istatus, JS
       if(istatus.eq.1) then
        write(6,*)'Premature end -Particule inside Ionopause'
        call output_data   
        STOP
       else if(istatus.eq.0)then  
        write(6,*)  '     '
        write (6,*) '******Simulation ended normally*****'
        write(6,*) 
        call output_data
       end if
       return
       end
c-------------------------------------------------------------------------
      subroutine Output_data
      implicit real*8 (a-h,o-z)
c
      common /solution_status/ Istatus, JS
      common /Solution_data/ t(300000),ud(300000,3),vd(300000,3)
c
      write(6,10) 
10    format(5x,'time',9x,'Vx',9x,'Ux',9x,'Vy',9x,'Uy',9x,'Vz',9x,'Uz')
      write(6,20)(t(i),(vd(i,j),ud(i,j),j=1,3),i=1,Js)
20    format(7(1pe12.5,2x))
      return 
      end

c***************************************************************************
      subroutine Qprf(Qpr,a,material)
c***************************************************************************      
c function returns Qpr, given a dust sie of diameter a [m]
c Only two materials are currently supported Olivine and Magnetite
      implicit real*8 (a-h,o-z)
      real*8 a, Qpr
      character*20 material
c      
c ----Magnetite fits expect a in meter
      if(material(1:9).eq.'Magnetite')then 
       if(a.le. 0.1D-6) then 
        Qpr= 17.89 * (1.e-6*a)
       else 
        Qpr=-0.8 * (1.e+6*a) +2.
       end if 
c Olivine fits expect a in meter

      else if (material(1:7).eq.'olivine')then
       if (a.le. 0.1e-6)then 
         Qpr= 3.0275 * (1.e+6*a)
       else
         Qpr=-0.4 * (1.e+6*a) + 1.4
       end if
      else
       write(*,*)'no material assigned'
       STOP
      endif
c
      RETURN
      END

      
c***************************************************************************
      subroutine Chif(Chi,a,material)
c***************************************************************************      
c  function returns Chi, given a dust particle of diam a
      implicit real*8 (a-h,o-z)
      real*8 a, Chi
      character*20 material
c      
c ----Magnetite fits expect a in meter

      if(material(1:9).eq.'Magnetite')then 
         Chi= 0.0292*(1.e+6*a)**(-1.238)
c ----Olivine fits expect a in meter
      else if(material(1:7).eq.'olivine')then
         Chi=0.0357*(1.e+6*a)**(-1.65)
      else
       write(*,*)'no material assigned'
       STOP
      end if
c
      RETURN
      END
