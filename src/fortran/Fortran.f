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
c
C     Yd(1), Yd(2), Yd(3) INITIAL VALUE OF THE DEPENDENT VARIABLES; THEY are              
C              Overridden WITH THE NEW VALUES in RUNG4 
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      external get_dydx
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
c read number of substeps within a simulation time step        
      read(5,*)nstep
      time_step=simulation_time_step/nstep
      write(6,50)nstep
50    Format('Number of iterations in a step:',I10)
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
c assign the number of dependent varaibles to 6 (3 velocities and 3 positions)
      neq=6
c
c
c
c ### Begin simulation Loop ###
      do Js=1,N_simulation_steps                 
       do i=1,3
        Yd(I)  =Vd(Js,i)
        YD(I+3)=UD(Js,i)
       end do
       t_init=t(Js)
       call runge4(neq,t_init,Yd,nstep,time_step,get_dydx) 
       do i=1,3
        Vd(Js+1,i)=Yd(i)
        Ud(Js+1,i)=Yd(i+3)       
       end do  
c 
       Vx=Yd(1)
       Vy=Yd(2)
       Vz=Yd(3)
       Ux=Yd(4)
       Uy=Yd(5)
       Uz=Yd(6)
c advance time to the next time station      
       t(Js+1)=t(Js) + simulation_time_step
      end do
c### End simulation loop.
c
c
c  Output results and exit  
      call exit_simulation
c
      stop
      end
c 
c-----------------------------------------------------------------------------------            
      subroutine get_dydx(t,Y,DYDX)
c-----------------------------------------------------------------------------------            
      implicit real*8 (a-h,o-z)   
      common /Efield/  Ex, Ey, Ez, Emod
      common /C1C2/ C1,C2
c
      dimension Y(10),DYDX(10)
c      
      call Ef(Y)
c      
      DYDX(1)= C2*Ex
      DYDX(2)= C2*Ey
      DYDX(3)= C2*Ez -C1
c
      DYDX(4)= Y(1)
      DYDX(5)= Y(2)
      DYDX(6)= Y(3)
c      
      return
      end
c
c-----------------------------------------------------------------------------------                                                                 
      SUBROUTINE runge4(N,X,Y,K,H,get_dydx)         
c-----------------------------------------------------------------------------------      
C-----THIS SUBROUTINE ADVANCES THE SOLUTION TO THE SYSTEM OF N ODE'S  
C-----DY(i)/DX = f(X,Y(1),Y(2),,.....Y(N)), i=1,...,N       
C-----USING THE FOURTH ORDER RUNGE-KUTTA METHOD.                        
C-----N = Number of differential equations in the system. N<=10
C-----X = INITIAL VALUE OF THE INDEPENDENT VARIABLE; IT IS INCREASED          
C-----    BY THE SUBROUTINE.                                                  
C-----Y = INITIAL VALUE OF THE DEPENDENT VARIABLE(S); THE ROUTINE             
C-----    OVERWRITES THIS WITH THE NEW VALUE(S).    
C-----H = time step
C-----K = number of steps to advance the equations
C-----get_dydx:  subroutine  get_dydx(X,Y,DYDX)

      IMPLICIT REAL*8(A-H,O-Z)                     
      EXTERNAL get_dydx      
      DIMENSION Y(10),YS(10),YSS(10),YSSS(10),T1(10),T2(10),                  
     1          T3(10),T4(10)                                                 
C-----THE MAIN LOOP                      
C      print*, X, Y(1), Y(2), Y(3), Y(4), Y(5)

      DO 1 I=1,K                                                              
C-----TEMPORARY ARRAYS ARE NEEDED FOR THE FUNCTIONS TO SAVE THEM              
C-----FOR THE FINAL CORRECTOR STEP.                                           
C-----FIRST (HALF STEP) PREDICTOR                                             
      call get_dydx(X,Y,T1)
      DO 2 J=1,N                                                              
      YS(J) = Y(J) + .5 * H * T1(J)
    2 CONTINUE                                                                
C-----SECOND STEP (HALF STEP CORRECTOR)                                       
      X = X + .5 * H                                                          
      call get_dydx(X,YS,T2)
      DO 3 J=1,N                                                              
      YSS(J) = Y(J) + .5 * H * T2(J)                                          
    3 CONTINUE                                                                
C-----THIRD STEP (FULL STEP MIDPOINT PREDICTOR)                               
      call get_dydx(X,YSS,T3)
      DO 4 J=1,N                                                              
      YSSS(J) = Y(J) + H * T3(J)                                              
    4 CONTINUE                                                                
C-----FINAL STEP (SIMPSON'S RULE CORRECTOR)                                   
      X = X + .5 * H                                                          
      call get_dydx(X,YSSS,T4)
      DO 5 J=1,N                                                              
      Y(J) = Y(J) + (H/6.)*(T1(J) + 2.*(T2(J)+T3(J)) + T4(J))                 
    5 CONTINUE                                                                
    1 CONTINUE                                                                
      RETURN                                                                  
      END         

      
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
               
c***************************************************************************
      Subroutine Ef(Y)
c***************************************************************************
c
c Generate the electric field at a point (x,y,z)
c INPUT: x,y,z and c: speed of light, B0 far field magnetic strength B0 and 
c distance between nucleus and ionopause apex
c OUTPUT Ex, Ey and Ez, and the modulus of E
c
      implicit real*8 (A-h,o-z)
c
      common /variables/ c,  B0, Z0, v0
      common /Efield/     Ex,Ey,Ez,Emod
      common /solution_status/ Istatus, JS
c
      dimension Y(10)
 
      ux=Y(4)
      uy=Y(5)
      uz=Y(6)
c
      r=dsqrt(ux**2.+uy**2.)
      r1=dsqrt(r**2.+uz**2.)
      S=(r**2.+2.*z0**2*(uz/r1-1.))
c  check if the particle inside the IONOPAUSE; Exit when it is the case.
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
       Ez=coef*z0**2*r**2/r1**(3.) /S * uy
c
       Emod= (Ex**2+ Ey**2+ Ez**2)**(0.5)
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
