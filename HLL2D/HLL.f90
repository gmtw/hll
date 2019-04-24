!=======================================================================
!   This program solves the Euler equations (2D) with the Lax Method
!=======================================================================
!un cambio diferente en el programa
!otro cambio
!mas cambios
!=======================================================================
!   This module contains global variables
      module globals
      implicit none

!------------------------------------------------------------------------------
!   nx = number of points used to discretize X
!   ny = number of points used to discretize Y
!   neq = number of equations
!   xmax = extent of X
!   ymax = extent of Y
!   gamma = 4/3 (1.3333) or 7/5 (1.4) or 5/3 (1.6666)
!
!   tmax = max time for the simulation
!   dtprint = outputs every dtprint
!
!   rhoin, pin, vxin, vyin = density, pressure, velocity at x<rad, y<rad
!   rhoout, pout, vxout, vyout = density, pressure, velocity at x>rad, y>rad
!
!   Co = courant number
!
!   **we can modify these values**
!------------------------------------------------------------------------------
      integer, parameter :: nx=200, ny=200, neq=4
      real, parameter :: xmax=3.0, dx=xmax/float(nx)
      real, parameter :: ymax=3.0, dy=ymax/float(ny)
      real, parameter :: gamma=4./3.

      real, parameter :: tmax= 7.     ! maximum integration time
      real, parameter :: dtprint=tmax/20.       ! interval between outputs

      real, parameter :: rhoin = 12.0
      real, parameter :: rhoout = 1.0
      real, parameter :: pin = 10.
      real, parameter :: pout = 1.0
      real, parameter :: vxin = 0.0
      real, parameter :: vyin = 0.0
      real, parameter :: vxout = 0.0
      real, parameter :: vyout = 0.0

      real, parameter :: xc = 3.
      real, parameter :: yc = 3.

      real, parameter :: Co = 0.7

!--------------------------------------
! Boundary Conditions
! 0.0 = open boundary conditions
! 10.0 = Jet injection
!--------------------------------------
      real, parameter :: bound=10.0
 
!--------------------------------------
! Election of fluxes
! Lax-Friederich = 1
! Harten-Lax-van-Leer = 2
!--------------------------------------
     integer, parameter :: choose_f = 2
!-----------------------------------------------------------------------------
!Mode relativistic
!Relativistic = 1
!Newton = 0

      integer, parameter :: relativity = 1
!-----------------------------------------------------------------------------



!------------------------------------------------------------------------------
!   This is a vector that contains u(x,y) (do not touch)
!------------------------------------------------------------------------------
      real :: u(neq,0:nx+1,0:ny+1)

!-----------------------------------------------------------------------------
!   This varial help us to mode relativistic!

      real :: qu(neq)
      real :: qp(neq)
      real :: qpp(neq,0:nx+1,0:ny+1)
!-----------------------------------------------------------------------------

        end module globals
!------------------------------------------------------------------------------
! End module globals
!==============================================================================

!--------

!==============================================================================
!   This module is the main program
!------------------------------------------------------------------------------
      program euler2D_lax
      use globals
      implicit none
      real      :: time, dt         ! t, $\Delta t$
      real      :: tprint           ! time of next output
      integer   :: itprint          ! number of current output

!------------------------------------------------------------------------------
! This subroutine generates the initial conditions main loop,
! iterate until maximum time is reached
!------------------------------------------------------------------------------
      call initconds(time, tprint, itprint)     ! conditions at t=0
      do while (time.lt.tmax)                   ! stops when t = tmax
        if(time.ge.tprint) then                 ! prints data in terminal
          write (*,*) itprint+1,time,tmax,dt
          call output (itprint)                 ! print data output
          tprint=tprint+dtprint
          itprint=itprint+1
        end if

!------------------------------------------------------------------------------
! Calculates and uses the CFL dt criterium
!------------------------------------------------------------------------------
        call courant(dt)
!------------------------------------------------------------------------------
! Calculates U(time) and F => U(time+dt)
!------------------------------------------------------------------------------
        call ulax(dt,time)
        time=time+dt

      end do

      end program euler2D_lax
!------------------------------------------------------------------------------
! end of main program module
!==============================================================================

!--------

!==============================================================================
! In this module we set the initial condition
!------------------------------------------------------------------------------
      subroutine initconds(time,tprint,itprint)
      use globals
      implicit none
      real, intent(out) :: time, tprint
      integer, intent (out) :: itprint
      integer ::i,j
      real :: x,y, rad, lorin,lorout, hin,hout, lorleft,lorright,hleft,hright

!------------------------------------------------------------------------------
! For the 2D circular blast:
! u(1,i,j) = rho(i,j)
! u(2,i,j) = vx(i,j)
! u(3,i,j) = vy(i,j)
! u(4,i,j) = etot(i,j) = eint + ekin = P/(gamma-1)
!------------------------------------------------------------------------------
    if(relativity==0)then
        
        do i=0,nx+1
          do j=0,ny+1
            x=float(i)*dx          ! obtain the position $x_i$
            y=float(j)*dy          ! obtain the position $y_j$
            rad=sqrt((x-xc)**2+(y-yc)**2)
            if (rad < 0. ) then
              u(1,i,j)=rhoin
              u(2,i,j)=rhoin*vxin
              u(3,i,j)=rhoin*vyin
              u(4,i,j)=pin/(gamma-1.) + 0.5*u(2,i,j)*u(2,i,j)/u(1,i,j) + 0.5*u(3,i,j)*u(3,i,j)/u(1,i,j)          
            else
              u(1,i,j)=rhoout
              u(2,i,j)=rhoout*vxout
              u(3,i,j)=rhoout*vyout
              u(4,i,j)=pout/(gamma-1.) + 0.5*u(2,i,j)*u(2,i,j)/u(1,i,j) + 0.5*u(3,i,j)*u(3,i,j)/u(1,i,j)
            end if

            
          end do
        end do

      elseif(relativity==1)then

         do i=0,nx+1
          do j=0,ny+1
            x=float(i)*dx          ! obtain the position $x_i$
            y=float(j)*dy          ! obtain the position $y_j$
           rad=abs(sqrt((x-xc)**2+(y-yc)**2))

            if (rad < 0.) then
             
              lorin=1/sqrt(1-(vxin**2+vyin**2))
              hin=1.+gamma/(gamma-1.)*pin/rhoin
             
              u(1,i,j)=rhoin*lorin
              u(2,i,j)=rhoin*vxin*lorin**2*hin
              u(3,i,j)=rhoin*vyin*lorin**2*hin
              u(4,i,j)=rhoin*lorin**2*hin-pin

            else
              lorout=1./sqrt(1.-(vxout**2+vyout**2))
              hout=1.+gamma/(gamma-1.)*pout/rhoout

              u(1,i,j)=rhoout*lorout
              u(2,i,j)=rhoout*vxout*lorout**2*hout
              u(3,i,j)=rhoout*vyout*lorout**2*hout
              u(4,i,j)=rhoout*lorout**2*hout-pout

           

             end if
            end do
           end do


    endif

!------------------------------------------------------------------------------
! end of the 2D circular blast initial condition
! reset the counters and time to 0
!------------------------------------------------------------------------------
      time=0
      tprint=0
      itprint=0

      return
      end subroutine initconds
!------------------------------------------------------------------------------
! end of the init condition module
!==============================================================================

!--------

!==============================================================================
! output files module
!------------------------------------------------------------------------------
      subroutine output(itprint)
      use globals
      implicit none
      integer, intent(in) :: itprint
      character (len=20) filerho, fileE, fileVel
      integer :: i,j
      real :: rho,vx,vy,P

!------------------------------------------------------------------------------
! the name of the output files are:
!------------------------------------------------------------------------------
      write(filerho,'(a,i2.2,a)') '2D_HLL-',itprint,'.dat'
      open(unit=10,file=filerho,status='unknown')

!------------------------------------------------------------------------------
! writes the output datafiles (ascii)
! writes x, y, rho, vx, vy and P
!------------------------------------------------------------------------------
      if (relativity==0)then

          do j=1,ny
            do i=1,nx
              rho=u(1,i,j)
              vx=u(2,i,j)/rho
              vy=u(3,i,j)/rho
              P=(u(4,i,j)-0.5*rho*(vx**2)*(vy**2))*(gamma-1.)
              write(10,'(7es12.5)') float(i)*dx,float(j)*dy,rho, vx, vy, P
            end do
            write(10,*)
          end do
        close(10)


      elseif(relativity==1)then

        do j=1,ny
          do i=1,nx
           ! rho=u(1,i,j)
           ! vx=u(2,i,j)/rho
           ! vy=u(3,i,j)/rho
           ! P=(u(4,i,j)-0.5*rho*(vx**2)*(vy**2))*(gamma-1.)
            qu(1) = u(1,i,j)
            qu(2) = u(2,i,j)
            qu(3) = u(3,i,j)
            qu(4) = u(4,i,j)

            call uprim(qu,qp)

            qpp(:,i,j)=qp

            rho=qpp(1,i,j)
            vx=qpp(2,i,j)
            vy=qpp(3,i,j)
            P=qpp(4,i,j)
            
            write(10,'(7es12.5)') float(i)*dx,float(j)*dy,rho, vx, vy, P
          end do
          write(10,*)
        end do
      close(10)



  

      endif


      return
      end subroutine output
!------------------------------------------------------------------------------
! end of the output files module
!==============================================================================

!--------

!==============================================================================
! CFL criterium module
!------------------------------------------------------------------------------
      subroutine courant(dt)
      use globals
      implicit none
      real, intent(out) ::dt
      real :: rho, vx, vy, P, cs
      integer :: i,j

!------------------------------------------------------------------------------
! Calculate the CFL criterium
!------------------------------------------------------------------------------
      dt=1E30

      if(relativity==0)then
          do i=0,nx+1
            do j=0,ny+1
              rho=u(1,i,j)
              vx=u(2,i,j)/rho
              vy=u(3,i,j)/rho
              P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)
              cs=sqrt(gamma*P/rho)
              dt=min( dt,Co*dx/(abs(vx)+cs) )
              dt=min( dt,Co*dy/(abs(vy)+cs) )

            end do
          end do

      elseif(relativity==1)then
          do i=0,nx+1
            do j=0,ny+1
             ! rho=u(1,i,j)
             ! vx=u(2,i,j)/rho
             ! vy=u(3,i,j)/rho
             ! P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)
              qu(1) = u(1,i,j)
              qu(2) = u(2,i,j)
              qu(3) = u(3,i,j)
              qu(4) = u(4,i,j)

              call uprim(qu,qp)

              qpp(:,i,j)=qp
              
              rho=qpp(1,i,j)
              vx=qpp(2,i,j)
              vy=qpp(3,i,j)
              P=qpp(4,i,j)

              cs=sqrt(gamma*P/rho)
              dt=min( dt,Co*dx/(abs(vx)+cs) )
              dt=min( dt,Co*dy/(abs(vy)+cs) )
              
              
        end do
      end do


        endif

      return
      end subroutine courant
!------------------------------------------------------------------------------
! end of the CFL criterium module
!==============================================================================

!--------

!==============================================================================
!Lax integration method module to calculate the conservative variables
!------------------------------------------------------------------------------
      subroutine ulax(dt,time)
      use globals
      implicit none
      real, intent(in) :: dt, time
      real :: up(neq,0:nx+1,0:ny+1), f(neq,0:nx+1,0:ny+1), g(neq,0:nx+1,0:ny+1)
      real :: fhll(neq,0:nx+1,0:ny+1), ghll(neq,0:nx+1,0:ny+1)
      real :: dtx, dty
      integer :: i,j

!------------------------------------------------------------------------------
! obtain the fluxes
!------------------------------------------------------------------------------
 

      

!------------------------------------------------------------------------------
!   Here is the Lax method, notice that the values at the extremes can not be
!   calculated, we need to enter then as a boundary conditions
!------------------------------------------------------------------------------
      dtx=dt/dx
      dty=dt/dy
      !print*,u(1,1,1), u(2,1,1),u(3,1,1),u(4,1,1)

      if (choose_f .eq. 1) then
          call fluxes(nx,ny,neq,gamma,u,f,g,bound)
	      
	     
      do i=1,nx
        do j=1,ny
          up(:,i,j)=0.25*( u(:,i-1,j)+u(:,i+1,j)+u(:,i,j-1)+u(:,i,j+1) ) &
                         -dtx*0.5*(f(:,i+1,j)-f(:,i-1,j) ) &
                         -dty*0.5*(g(:,i,j+1)-g(:,i,j-1) )

           end do
      	 end do

      elseif(choose_f .eq. 2)then
        
        call fluxesHLL(fhll, ghll)
      
      	do i=1,nx
	        do j=1,ny

             up(:,i,j) = u(:,i,j) + 0.5*dtx*(fhll(:,i-1,j)-fhll(:,i+1,j))      &
                         + 0.5*dty*(ghll(:,i,j-1)-ghll(:,i,j+1))
        	end do
      	  end do

      endif

!------------------------------------------------------------------------------
!   Boundary conditions to the U^n+1
!------------------------------------------------------------------------------
      call boundaries(nx,ny,neq,up,time,gamma,bound)
!------------------------------------------------------------------------------
!   copy the up to the u
!------------------------------------------------------------------------------
      u(:,:,:)=up(:,:,:)
      !print*, dt, u(1,1,1), u(2,1,1),u(3,1,1),u(4,1,1)

      return
      end subroutine ulax
!------------------------------------------------------------------------------
!   end of the lax integration method module
!==============================================================================

!--------

!==============================================================================
! Calculation of the fluxes module
!------------------------------------------------------------------------------
      subroutine fluxes(nx,ny,neq,gamma,u,f,g,bound)
      use globals, only : qu, qpp,qp, relativity
      implicit none
      integer, intent(in) :: nx,ny,neq
      real, intent(in) :: gamma
      real, intent(in) :: bound
      real, intent(in) :: u(neq,0:nx+1,0:ny+1)
      real, intent(out) :: f(neq,0:nx+1,0:ny+1), g(neq,0:nx+1,0:ny+1)
      integer :: i, j
      real :: rho, vx, vy, P, lor, h

      if(relativity==0)then
          do i=1,nx
            do j=1,ny
              rho=u(1,i,j)
              vx=u(2,i,j)/rho
              vy=u(3,i,j)/rho
              P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)

              f(1,i,j)=rho*vx
              f(2,i,j)=rho*vx*vx+P
              f(3,i,j)=rho*vx*vy
              f(4,i,j)=vx*(u(4,i,j)+P)

              g(1,i,j)=rho*vy
              g(2,i,j)=rho*vx*vy
              g(3,i,j)=rho*vy*vy+P
              g(4,i,j)=vy*(u(4,i,j)+P)

            end do
          end do

        elseif(relativity==1)then
          
          do i=0,nx+1
            do j=0,ny+1
              

              
              qu(1) = u(1,i,j)
              qu(2) = u(2,i,j)
              qu(3) = u(3,i,j)
              qu(4) = u(4,i,j)

              call uprim(qu,qp)

              qpp(:,i,j)=qp

              rho=qpp(1,i,j)
              vx= qpp(2,i,j)
              vy= qpp(3,i,j)
              P= qpp(4,i,j)

              lor=1/sqrt(1-(vx**2+vy**2))
              h=1.+gamma/(gamma-1.)*P/rho
              
              f(1,i,j)=rho*vx*lor
              f(2,i,j)=rho*vx*vx*lor**2*h+P
              f(3,i,j)=rho*vx*vy*lor**2*h
              f(4,i,j)=rho*vx*lor**2*h

              g(1,i,j)=rho*vy*lor
              g(2,i,j)=rho*vx*vy*lor**2*h
              g(3,i,j)=rho*vy*vy*lor**2*h+P
              g(4,i,j)=rho*vy*lor**2*h

             ! print*,g

            end do
          end do

      endif

      return
      end subroutine fluxes
!------------------------------------------------------------------------------
! end of the fluxes module
!==============================================================================

!--------

!==============================================================================
! Boundary conditions module
!------------------------------------------------------------------------------
      subroutine boundaries(nx,ny,neq,u,time,gamma,bound)
      use globals, only : relativity
      implicit none
      integer, intent(in) :: nx,ny,neq
      real, intent(in) :: bound
      real, intent(out) ::u(neq,0:nx+1,0:ny+1)
      real, intent(in) :: time,gamma
      integer :: j , i
      real :: rhojet , vjet, Pjet, Taujet , rhoijet, theta, lorjet, hjet, vxjet, vyjet

!------------------------------------------------------------------------------
! Periodic or open boundaries, or your own boundary condition
!------------------------------------------------------------------------------

        ! u(:,0   ,:)=u(:,1 ,:)
        ! u(:,nx+1,:)=u(:,nx,:)
        ! u(:,:,0   )=u(:,:,1 )
        ! u(:,:,ny+1)=u(:,:,ny)

!---------------------------------------

        if(bound.eq.0.)then
          u(:,0   ,:)=u(:,1 ,:)
          u(:,nx+1,:)=u(:,nx,:)
          u(:,:,0   )=u(:,:,1 )
          u(:,:,ny+1)=u(:,:,ny)
        end if

!---------------------------------------

        if(bound.eq.1.)then
           u(1,0   ,:)=u(1,1 ,:)
           u(2,0   ,:)=-u(2,2 ,:)  !-u(2,1 ,:) si se quiere ref en x=0
           u(3,0   ,:)=u(3,1 ,:)
           u(4,0   ,:)=u(4,1 ,:)

         end if

!---------------------------------------

        if(bound.eq.2.)then
          u(1,nx+1,:)=u(1,nx,:)
          u(2,nx+1,:)=-u(2,nx,:)  !-u(2,nx,:) si se quiere ref en x=nx
          u(3,nx+1,:)=u(3,nx,:)
          u(4,nx+1,:)=u(4,nx,:)
         end if

!---------------------------------------

        if(bound.eq.3.)then
          u(1,:,0   )=u(1,:,1 )
          u(2,:,0   )=u(2,:,1 )
          u(3,:,0   )=-u(3,:,1 )  !-u(3,:,1 ) si se quiere ref en y=0
          u(4,:,0   )=u(4,:,1 )
         end if

!---------------------------------------

        if(bound.eq.4.)then
          u(1,:,ny+1)=u(1,:,ny)
          u(2,:,ny+1)=u(2,:,ny)
          u(3,:,ny+1)=-u(3,:,ny)  !-u(3,:,ny) si se quiere ref en y=ny
          u(4,:,ny+1)=u(4,:,ny)
         end if

!---------------------------------------

         if(bound.eq.5.)then
           u(1,0   ,:)=u(1,1 ,:)
           u(2,0   ,:)=-u(2,1 ,:)  !-u(2,1 ,:) si se quiere ref en x=0
           u(3,0   ,:)=u(3,1 ,:)
           u(4,0   ,:)=u(4,1 ,:)

          u(1,nx+1,:)=u(1,nx,:)
          u(2,nx+1,:)=-u(2,nx,:)  !-u(2,nx,:) si se quiere ref en x=nx
          u(3,nx+1,:)=u(3,nx,:)
          u(4,nx+1,:)=u(4,nx,:)
         end if

!---------------------------------------

        if(bound.eq.6.)then
          u(1,:,0   )=u(1,:,1 )
          u(2,:,0   )=u(2,:,1 )
          u(3,:,0   )=-u(3,:,1 )  !-u(3,:,1 ) si se quiere ref en y=0
          u(4,:,0   )=u(4,:,1 )

          u(1,:,ny+1)=u(1,:,ny)
          u(2,:,ny+1)=u(2,:,ny)
          u(3,:,ny+1)=-u(3,:,ny)  !-u(3,:,ny) si se quiere ref en y=ny
          u(4,:,ny+1)=u(4,:,ny)
         end if

!---------------------------------------

         if(bound.eq.7.)then
          u(1,0   ,:)=u(1,1 ,:)
          u(2,0   ,:)=-u(2,1,:)  !-u(2,1 ,:) si se quiere ref en x=0
          u(3,0   ,:)=u(3,1 ,:)
          u(4,0   ,:)=u(4,1 ,:)

          u(1,nx+1,:)=u(1,nx,:)
          u(2,nx+1,:)=-u(2,nx,:)  !-u(2,nx,:) si se quiere ref en x=nx
          u(3,nx+1,:)=u(3,nx,:)
          u(4,nx+1,:)=u(4,nx,:)

          u(1,:,0   )=u(1,:,1 )
          u(2,:,0   )=u(2,:,1 )
          u(3,:,0   )=-u(3,:,1 )  !-u(3,:,1 ) si se quiere ref en y=0
          u(4,:,0   )=u(4,:,1 )

          u(1,:,ny+1)=u(1,:,ny)
          u(2,:,ny+1)=u(2,:,ny)
          u(3,:,ny+1)=-u(3,:,ny)  !-u(3,:,ny) si se quiere ref en y=ny
          u(4,:,ny+1)=u(4,:,ny)
         end if

!---------------------------------------
!  Jet injection
!---------------------------------------
         if(bound.eq.10.)then
           rhojet=40.5
           Taujet=100.2
           vjet=0.5 !+ 0.3*sin(2*3.1416*time/Tauj)
           vxjet=0.99
           vyjet=0.10
           Pjet=0.1     !0.1
           theta=50


           u(:,0   ,:)=u(:,1 ,:)
           u(:,nx+1,:)=u(:,nx,:)
           u(:,:,0   )=u(:,:,1 )
           u(:,:,ny+1)=u(:,:,ny)

          if(relativity==0)then
             do j=1,ny
                if (abs(j-ny/2) <= ny/20)  then
                  u(1,0,j)=rhojet
                  u(2,0,j)=rhojet*vxjet
                  if((j-ny/2).ge.0)then
                    u(3,0,j)=-rhojet*vyjet
                  else
                    u(3,0,j)=rhojet*vyjet
                  end if
                  u(4,0,j)=0.5*rhojet*(vxjet**2+vyjet**2) + Pjet/(gamma-1.)
                end if
              end do

                !u(:,0   ,:)=u(:,1 ,:)
                !u(:,nx+1,:)=u(:,nx,:)
                !u(:,:,0   )=u(:,:,1 )
                !u(:,:,ny+1)=u(:,:,ny)

              ! do j=1,nx
              !    if (abs(j-nx/2) <= nx/20)  then
              !      u(1,j,0)=rhoj
                   
              !      if((j-nx/2).ge.0)then
                   
              !        u(2,j,0)= -rhoj*vxj     !*COS(45*3.14*time/180)
                   
              !      else
              !        u(2,j,0)=  rhoj*vxj   !*COS(45*3.14*time/180)
                   
              !      end if
                   

              !      u(3,j,0)= rhoj*vyj   !*sin(45*3.14*time/180)
                   
              !      u(4,j,0)=0.5*rhoj*(vxj**2+vyj**2) + Pj/(gamma-1.)
              !    end if
              !  end do

              !   u(:,nx+1,:)=u(:,nx,:)
                !u(:,:,0   )=u(:,:,1 )
                !u(:,:,ny+1)=u(:,:,ny)
          elseif(relativity==1)then

            do j=0,ny
           
            lorjet=1/sqrt(1-(vxjet**2+vyjet**2))
            hjet=1.+gamma/(gamma-1.)*Pjet/rhojet

           if(abs(j-ny/2) <= ny/20) then
            
            u(1,0,j)=rhojet*lorjet

            u(2,0,j)=rhojet*vxjet*lorjet**2*hjet!*abs(COS(2*3.14*time/Taujet)
              
              if((j-ny/2).ge.0)then
               
                u(3,0,j) =-rhojet*vyjet*lorjet**2*hjet!*COS(2*3.14*time/Taujet)
              else

                u(3,0,j) = rhojet*vyjet*lorjet**2*hjet!*COS(2*3.14*time/Taujet)
              end if
            

            u(4,0,j)=rhojet*lorjet**2*hjet-Pjet
           end if

          end do


          endif                                                                                                                                                                                                                                                                                                     


         end if
!---------------------------------------
!  End jet injection
!---------------------------------------

  return
end subroutine boundaries
!=======================================================================

subroutine RL(nx,ny, neq, gamma, u, rho_l,rho_d, rho_r,rho_u, vx_l,vx_d, vx_r,vx_u, &
  vy_l,vy_r,vy_u, vy_d,P_l,P_d, P_r,P_u)
  !Este modulo, digamos es para seleccionar, las densidades, velocidades y presiones tanto de izquierda como de derecha
  !Tambien sirve como desacoplamiento para mis primitivas
  use globals, only: relativity, qu, qpp, qp
  implicit none
  integer, intent(in) :: nx,ny,neq
  real, intent(in) :: gamma
  real, intent(in) :: u(neq,0:nx+1,0:ny+1)
  real, intent(out):: rho_l(0:nx+1,0:ny+1), rho_r(0:nx+1,0:ny+1),rho_u(0:nx+1,0:ny+1),rho_d(0:nx+1,0:ny+1)
  real, intent(out):: vx_l(0:nx+1,0:ny+1), vx_r(0:nx+1,0:ny+1), vx_u(0:nx+1,0:ny+1), vx_d(0:nx+1,0:ny+1)
  real, intent(out):: vy_l(0:nx+1,0:ny+1), vy_r(0:nx+1,0:ny+1), vy_u(0:nx+1,0:ny+1), vy_d(0:nx+1,0:ny+1)
  real, intent(out):: P_l(0:nx+1, 0:ny+1), P_r(0:nx+1, 0:ny+1), P_u(0:nx+1,0:ny+1), P_d(0:nx+1,0:ny+1)

  !real, parameter ::
  integer :: i,j


 
!Tuve que convertir mis variables primitivas en vectores para tener un mejor control de lo que toca derecha a izquirda,
! no se si es lo mej
! notar que corre de 0 a nx 
  
  if(relativity==0)then
    
    do i=1,nx
      do j=1,ny 
      
      rho_l(i,j)= u(1,i-1,j)
      vx_l(i,j) = u(2,i-1,j)/u(1,i-1,j)
      vy_l(i,j) = u(3,i-1, j)/u(1,i-1,j)
      P_l(i,j)  = (u(4,i-1,j)-0.5*rho_l(i,j)*(vx_l(i,j)**2+vy_l(i,j)**2))*(gamma-1.)

      rho_d(i,j)= u(1,i,j-1)
      vx_d(i,j) = u(2,i,j-1)/u(1,i,j-1)
      vy_d(i,j) = u(3,i,j-1)/u(1,i,j-1)
      P_d(i,j)  = (u(4,i,j-1)-0.5*rho_d(i,j)*(vx_d(i,j)**2+vy_d(i,j)**2))*(gamma-1.)

      

      rho_r(i,j)= u(1,i+1,j)
      vx_r(i,j) = u(2,i+1,j)/u(1,i+1,j)
      vy_r(i,j) = u(3,i+1,j)/u(1,i+1,j)
      P_r(i,j)  = (u(4,i+1,j)-0.5*rho_r(i,j)*(vx_r(i,j)**2+vy_r(i,j)**2))*(gamma-1.)

      rho_u(i,j)= u(1,i,j+1)
      vx_u(i,j) = u(2,i,j+1)/u(1,i,j+1)
      vy_u(i,j) = u(3,i, j+1)/u(1,i,j+1)
      P_u(i,j)  = (u(4,i,j+1)-0.5*rho_u(i,j)*(vx_u(i,j)**2+vy_u(i,j)**2))*(gamma-1.)
      
      end do
    end do

  elseif(relativity ==1)then


    do i=1,nx
      do j=1,ny

         qu(1) = u(1,i,j)
         qu(2) = u(2,i,j)
         qu(3) = u(3,i,j)
         qu(4) = u(4,i,j)

         call uprim(qu,qp)

         qpp(:,i,j)=qp

        !lor=1/sqrt(1-(vx**2+vy**2))

        rho_l(i,j)= qpp(1,i-1,j)
        vx_l(i,j) = qpp(2,i-1,j)
        vy_l(i,j) = qpp(3,i-1,j)
        P_l(i,j)  = qpp(4,i-1,j)

        rho_d(i,j)= qpp(1,i,j-1)
        vx_d(i,j) = qpp(2,i,j-1)
        vy_d(i,j) = qpp(3,i,j-1)
        P_d(i,j)  = qpp(4,i,j-1)

        

        rho_r(i,j)= qpp(1,i+1,j)
        vx_r(i,j) = qpp(2,i+1,j)
        vy_r(i,j) = qpp(3,i+1,j)
        P_r(i,j)  = qpp(4,i+1,j)

        rho_u(i,j)= qpp(1,i,j+1)
        vx_u(i,j) = qpp(2,i,j+1)
        vy_u(i,j) = qpp(3,i,j+1)
        P_u(i,j)  = qpp(4,i,j+1)

      end do
     end do





  endif
  

  
  return

end subroutine RL
!=======================================================================


 !calculate the wave speeds
   subroutine wavespeeds(s_l,s_d, s_r,s_u,rho_l,rho_d, rho_r,rho_u, vx_l,vx_d, vx_r,vx_u, &
  vy_l,vy_r,vy_u, vy_d,P_l,P_d, P_r,P_u )
    ! Este modulo calcula las velocidadesd de la onda
     use globals
     implicit none
!    integer, intent(in)::nx,neq  
    real ::s_l(0:nx+1,0:ny+1),s_r(0:nx+1,0:ny+1), s_d(0:nx+1,0:ny+1), s_u(0:nx+1,0:ny+1) !Regresa las velocidades de la onda (s_l, s_r)
    
    real :: rho_l(0:nx+1,0:ny+1), rho_r(0:nx+1,0:ny+1),rho_u(0:nx+1,0:ny+1),rho_d(0:nx+1,0:ny+1)
    real :: vx_l(0:nx+1,0:ny+1), vx_r(0:nx+1,0:ny+1), vx_u(0:nx+1,0:ny+1), vx_d(0:nx+1,0:ny+1)
    real :: vy_l(0:nx+1,0:ny+1), vy_r(0:nx+1,0:ny+1), vy_u(0:nx+1,0:ny+1), vy_d(0:nx+1,0:ny+1)
    real :: P_l(0:nx+1, 0:ny+1), P_r(0:nx+1, 0:ny+1), P_u(0:nx+1,0:ny+1), P_d(0:nx+1,0:ny+1)
    
    real:: cs_l(0:nx+1,0:ny+1), cs_r(0:nx+1,0:ny+1), cs_u(0:nx+1,0:ny+1), cs_d(0:nx+1,0:ny+1) !velocidades del sonido
    integer:: i,j
 
 ! bucle para crear las velocidades del sonido tanto de izquierda como de derecha
  do i=1,nx
    do j=1,ny
      cs_l(i,j)=sqrt(gamma*P_l(i,j)/rho_l(i,j))  
      cs_d(i,j)=sqrt(gamma*P_d(i,j)/rho_d(i,j))

      cs_r(i,j)=sqrt(gamma*P_r(i,j)/rho_r(i,j))
      cs_u(i,j)=sqrt(gamma*P_u(i,j)/rho_u(i,j))
    end do
  end do

  

 

  !Este bucle calcula las velocidades de las ondas 
  do i=1,nx
    do j=1,ny
      s_l(i,j)=vx_l(i,j)-cs_l(i,j) !velocidad izquierda - velocidad del sonido izquierda
      s_d(i,j)=vy_d(i,j)-cs_d(i,j)

      s_r(i,j)=vx_r(i,j)+cs_r(i,j) !velocidad derecha + velocidad del sonido derecha
      s_u(i,j)=vy_u(i,j)+cs_u(i,j)
    end do
  end do

  


  
  return
   end subroutine wavespeeds


!    ! ! !=======================================================================
   subroutine fluxesHLL(fhll, ghll)
!     !Este modulo calcula los flujos hll
      use globals
     implicit none
!     real, intent(in) :: time
    
     real, intent(out)::fhll(neq,0:nx+1,0:ny+1), ghll(neq,0:nx+1,0:ny+1)! esta es la variable que se obtiene y va para la subrutina ulax

     real :: rho_l(0:nx+1,0:ny+1), rho_r(0:nx+1,0:ny+1),rho_u(0:nx+1,0:ny+1),rho_d(0:nx+1,0:ny+1)
     real :: vx_l(0:nx+1,0:ny+1), vx_r(0:nx+1,0:ny+1), vx_u(0:nx+1,0:ny+1), vx_d(0:nx+1,0:ny+1)
     real :: vy_l(0:nx+1,0:ny+1), vy_r(0:nx+1,0:ny+1), vy_u(0:nx+1,0:ny+1), vy_d(0:nx+1,0:ny+1)
     real :: P_l(0:nx+1, 0:ny+1), P_r(0:nx+1, 0:ny+1), P_u(0:nx+1,0:ny+1), P_d(0:nx+1,0:ny+1)

     real :: lor_l(0:nx+1, 0:ny+1), lor_r(0:nx+1, 0:ny+1), lor_u(0:nx+1, 0:ny+1), lor_d(0:nx+1, 0:ny+1)
     real :: h_l(0:nx+1,0:ny+1), h_r(0:nx+1,0:ny+1), h_u(0:nx+1,0:ny+1), h_d(0:nx+1,0:ny+1)

     real ::s_l(0:nx+1,0:ny+1),s_r(0:nx+1,0:ny+1), s_d(0:nx+1,0:ny+1), s_u(0:nx+1,0:ny+1)
    
     

     real :: f_l(neq, 0:nx+1,0:ny+1), f_r(neq,0:nx+1,0:ny+1), u_l(neq,0:nx+1,0:ny+1), u_r(neq,0:nx+1,0:ny+1)
     real :: g_d(neq, 0:nx+1,0:ny+1), g_u(neq,0:nx+1,0:ny+1), u_d(neq,0:nx+1,0:ny+1), u_u(neq,0:nx+1,0:ny+1)

     integer :: i,j

     call RL(nx,ny, neq, gamma, u, rho_l,rho_d, rho_r,rho_u, vx_l,vx_d, vx_r,vx_u, &
  vy_l,vy_r,vy_u, vy_d,P_l,P_d, P_r,P_u)

     call wavespeeds(s_l,s_d, s_r,s_u,rho_l,rho_d, rho_r,rho_u, vx_l,vx_d, vx_r,vx_u, &
  vy_l,vy_r,vy_u, vy_d,P_l,P_d, P_r,P_u )

    
!       !igual calculamos los flujos y conservadas de lado derecho e izquierdo

  if(relativity==0)then

       do i=1,nx
        do j=1,ny


           u_l(1,i,j) = rho_l(i,j)
           u_l(2,i,j) = rho_l(i,j)*vx_l(i,j)
           u_l(3,i,j) = rho_l(i,j)*vy_l(i,j)
           u_l(4,i,j) = 0.5*rho_l(i,j)*(vx_l(i,j)**2+vy_l(i,j)**2)+P_l(i,j)/(gamma-1.)


           u_d(1,i,j) = rho_d(i,j)
           u_d(2,i,j) = rho_d(i,j)*vx_d(i,j)
           u_d(3,i,j) = rho_d(i,j)*vy_d(i,j)
           u_d(4,i,j) = 0.5*rho_d(i,j)*(vx_d(i,j)**2+vy_d(i,j)**2)+P_d(i,j)/(gamma-1.)


           u_r(1,i,j) = rho_r(i,j)
           u_r(2,i,j) = rho_r(i,j)*vx_r(i,j)
           u_r(3,i,j) = rho_r(i,j)*vy_r(i,j)
           u_r(4,i,j) = 0.5*rho_r(i,j)*(vx_r(i,j)**2+vy_r(i,j)**2)+P_r(i,j)/(gamma-1.)


           u_u(1,i,j) = rho_u(i,j)
           u_u(2,i,j) = rho_u(i,j)*vx_u(i,j)
           u_u(3,i,j) = rho_u(i,j)*vy_u(i,j)
           u_u(4,i,j) = 0.5*rho_u(i,j)*(vx_u(i,j)**2+vy_u(i,j)**2)+P_u(i,j)/(gamma-1.)

          


    !========================================================

           f_l(1,i,j) = rho_l(i,j)*vx_l(i,j)
           f_l(2,i,j) = rho_l(i,j)*vx_l(i,j)**2+P_l(i,j)
           f_l(3,i,j) = rho_l(i,j)*vx_l(i,j)*vy_l(i,j)
           f_l(4,i,j) = vx_l(i,j)*(u_l(4,i,j)+P_l(i,j))

           

           f_r(1,i,j) = rho_r(i,j)*vx_r(i,j)
           f_r(2,i,j) = rho_r(i,j)*vx_r(i,j)**2+P_r(i,j)
           f_r(3,i,j) = rho_r(i,j)*vx_r(i,j)*vy_r(i,j)
           f_r(4,i,j) = vx_r(i,j)*(u_r(4,i,j)+P_r(i,j))


           g_d(1,i,j) = rho_d(i,j)*vy_d(i,j)
           g_d(2,i,j) = rho_d(i,j)*vx_d(i,j)*vy_d(i,j)
           g_d(3,i,j) = rho_d(i,j)*vy_d(i,j)**2 + P_d(i,j)
           g_d(4,i,j) = vy_d(i,j)*(u_d(4,i,j)+P_d(i,j))


           g_u(1,i,j) = rho_u(i,j)*vy_u(i,j)
           g_u(2,i,j) = rho_u(i,j)*vx_u(i,j)*vy_u(i,j)
           g_u(3,i,j) = rho_u(i,j)*vy_u(i,j)**2 + P_u(i,j)
           g_u(4,i,j) = vy_u(i,j)*(u_u(4,i,j)+P_u(i,j))

    !      
    !       ! Estas conservadas nos sirven para calcular el fhll
           

    !       
        end do
      end do

elseif(relativity==1)then

      do i=1,nx
        do j=1,ny

           lor_l(i,j) = 1/sqrt(1-(vx_l(i,j)**2+vy_l(i,j)**2))
           lor_d(i,j) = 1/sqrt(1-(vx_d(i,j)**2+vy_d(i,j)**2))
           lor_r(i,j) = 1/sqrt(1-(vx_r(i,j)**2+vy_r(i,j)**2))
           lor_u(i,j) = 1/sqrt(1-(vx_u(i,j)**2+vy_u(i,j)**2))

           h_l(i,j)=1.+gamma/(gamma-1.)*P_l(i,j)/rho_l(i,j)
           h_d(i,j)=1.+gamma/(gamma-1.)*P_d(i,j)/rho_d(i,j)
           h_r(i,j)=1.+gamma/(gamma-1.)*P_r(i,j)/rho_r(i,j)
           h_u(i,j)=1.+gamma/(gamma-1.)*P_u(i,j)/rho_u(i,j)



           u_l(1,i,j) = rho_l(i,j)*lor_l(i,j)
           u_l(2,i,j) = rho_l(i,j)*vx_l(i,j)*lor_l(i,j)**2*h_l(i,j)
           u_l(3,i,j) = rho_l(i,j)*vy_l(i,j)*lor_l(i,j)**2*h_l(i,j)
           u_l(4,i,j) = rho_l(i,j)*lor_l(i,j)**2*h_l(i,j)-P_l(i,j)


           u_d(1,i,j) = rho_d(i,j)*lor_d(i,j)
           u_d(2,i,j) = rho_d(i,j)*vx_d(i,j)*lor_d(i,j)**2*h_d(i,j)
           u_d(3,i,j) = rho_d(i,j)*vy_d(i,j)*lor_d(i,j)**2*h_d(i,j)
           u_d(4,i,j) = rho_d(i,j)*lor_d(i,j)**2*h_d(i,j)-P_d(i,j)


           u_r(1,i,j) = rho_r(i,j)*lor_r(i,j)
           u_r(2,i,j) = rho_r(i,j)*vx_r(i,j)*lor_r(i,j)**2*h_r(i,j)
           u_r(3,i,j) = rho_r(i,j)*vy_r(i,j)*lor_r(i,j)**2*h_r(i,j)
           u_r(4,i,j) = rho_r(i,j)*lor_r(i,j)**2*h_r(i,j)-P_r(i,j)

           u_u(1,i,j) = rho_u(i,j)*lor_u(i,j)
           u_u(2,i,j) = rho_u(i,j)*vx_u(i,j)*lor_u(i,j)**2*h_u(i,j)
           u_u(3,i,j) = rho_u(i,j)*vy_u(i,j)*lor_u(i,j)**2*h_u(i,j)
           u_u(4,i,j) = rho_u(i,j)*lor_u(i,j)**2*h_u(i,j)-P_u(i,j)


!==========================================================================
 
          f_l(1,i,j) = rho_l(i,j)*vx_l(i,j)*lor_l(i,j)
          f_l(2,i,j) = rho_l(i,j)*vx_l(i,j)**2*lor_l(i,j)**2*h_l(i,j)+P_l(i,j)
          f_l(3,i,j) = rho_l(i,j)*vx_l(i,j)*vy_l(i,j)*lor_l(i,j)**2*h_l(i,j)
          f_l(4,i,j) = rho_l(i,j)*vx_l(i,j)*lor_l(i,j)**2*h_l(i,j)

          f_r(1,i,j) = rho_r(i,j)*vx_r(i,j)*lor_r(i,j)
          f_r(2,i,j) = rho_r(i,j)*vx_r(i,j)**2*lor_r(i,j)**2*h_r(i,j)+P_r(i,j)
          f_r(3,i,j) = rho_r(i,j)*vx_r(i,j)*vy_r(i,j)*lor_r(i,j)**2*h_r(i,j)
          f_r(4,i,j) = rho_r(i,j)*vx_r(i,j)*lor_r(i,j)**2*h_r(i,j)


          g_d(1,i,j) = rho_d(i,j)*vy_d(i,j)*lor_d(i,j)
          g_d(2,i,j) = rho_d(i,j)*vx_d(i,j)*vy_d(i,j)*lor_d(i,j)**2*h_d(i,j)
          g_d(3,i,j) = rho_d(i,j)*vy_d(i,j)**2*lor_d(i,j)**2*h_d(i,j)+P_d(i,j)
          g_d(4,i,j) = rho_d(i,j)*vy_d(i,j)*lor_d(i,j)**2*h_d(i,j)

          g_u(1,i,j) = rho_u(i,j)*vy_u(i,j)*lor_u(i,j)
          g_u(2,i,j) = rho_u(i,j)*vx_u(i,j)*vy_u(i,j)*lor_u(i,j)**2*h_u(i,j)
          g_u(3,i,j) = rho_u(i,j)*vy_u(i,j)**2*lor_u(i,j)**2*h_u(i,j)+P_u(i,j)
          g_u(4,i,j) = rho_u(i,j)*vy_u(i,j)*lor_u(i,j)**2*h_u(i,j)


        end do
      end do

           


endif


     do i=1,nx
        do j=1,ny

           
          if (0 .le. s_l(i,j)) then !less or equal 0<=sl
            fhll(:,i,j)=f_l(:, i,j)

        
      

          else if (s_l(i,j) .le. 0 .and. 0 .le. s_r(i,j)) then !sl<=0<=sr
        
            fhll(:,i,j)=(s_r(i,j)*f_l(:,i,j)-s_l(i,j)*f_r(:,i,j)+s_l(i,j)*s_r(i,j)*(u_r(:,i,j)-u_l(:,i,j)))/&
            (s_r(i,j)-s_l(i,j))

      
          else if (s_r(i,j) .le. 0) then !sr<=0
              fhll(:,i,j)=f_r(:,i,j)

          endif

          

          if(0 .le. s_d(i,j)) then
            ghll(:,i,j)= g_d(:,i,j)


           else if(s_d(i,j) .le. 0 .and. 0 .le. s_u(i,j)) then


            ghll(:,i,j)=(s_u(i,j)*g_d(:,i,j)-s_d(i,j)*g_u(:,i,j)+s_d(i,j)*s_u(i,j)*(u_u(:,i,j)-u_d(:,i,j)))/&
            (s_u(i,j)-s_d(i,j))


          else if(s_u(i,j) .le. 0) then
              ghll(:,i,j) = g_u(:,i,j)


          endif
        end do
     end do

    


   
   fhll(:,nx+1,:)=fhll(:,nx,:)
   fhll(:,0,:)=fhll(:,1,:)  
   fhll(:,:,ny+1)=fhll(:,:,ny)
   fhll(:,:,0)=fhll(:,:,1)

   ghll(:,nx+1,:)=ghll(:,nx,:)
   ghll(:,0,:)=ghll(:,1,:)
   ghll(:,:,ny+1)=ghll(:,:,ny)
   ghll(:,:,0)=ghll(:,:,1)
  


  
  
   return

   end subroutine fluxesHLL
!------------------------------------------------------------

   subroutine uprim(qu,qp)

!    use parameters, only : neq, gamma
     use globals, only : neq, gamma
    implicit none
    real :: qu(neq), qp(neq) ! qp Devuelve las primitivas
    real :: w, m2, lor
    real :: u2, alpha, chi

    m2 = sum(qu(2:3)**2)               ! v^2

    call newrap(qu, w, m2)

    alpha = m2 / w**2   ! alpha < 1 !
    u2  = alpha/(1.0-alpha)

    lor = sqrt(1.0 + u2)

    ! velocities
    qp(2:3) = qu(2:3) / w

    ! determination of the mass density
    qp(1) = qu(1)/lor

    ! thermal pressure
    chi = (w - qu(1)*(1.0+u2/(lor+1.0)))/(1.0+u2)

    qp(4)   = (gamma - 1.0)/gamma * chi

    qp(4) =max(qp(4),1d-10*qp(1))

    if(qp(1).lt.0.0) then
      print*,'negative density in uprim'
      stop
    end if

  end subroutine uprim
! !------------------------------------------------------------------------------!

! !------------------------------------------------------------------------------!
  subroutine newrap(qu, w, m2)

!    use parameters, only : neq, gamma
    use globals , only : neq, gamma
    implicit none
   ! real, intent(in) :: qu(neq)
    real, intent(in)  :: m2
    real, intent(out) :: w
    real, parameter :: eps = 1d-10
    real :: a, b, c, mu, alpha, u2, lor, chi, dpdchi, dpdrho
    real :: w0, dpdw, f, dfdw, pg, dv2dw, dchidw, drhodw, qu(neq)
    integer :: k

    !print*, qu(4)**2, m2+qu(1)**2, qu(4)

    if(qu(4)**2 .lt. m2+qu(1)**2 .or. qu(4).le. 0.0) then !.lt. -> '<'; .le. -> '<=='
      print*,'error in newrap'
      stop
    end if

    a = 3.0;   b = 2.0 * (-qu(4));   c = m2
    if(b**2-a*c .lt. 0.0) then
       print*,'b**2-a*c<0'
       stop
    end if
    w = ( - b + sqrt(b**2-a*c) ) / a   ! initial guess for w = rho * h * lor**2

    w0 = w
    mu = 1.0
  100 continue
    do k = 1, 40

      alpha = m2 / w**2   ! alpha < 1 !

      u2  = alpha/(1.0-alpha)

      if(u2 .lt. 0.0) then
        print*,'u2<0cc'
        print*,qu
        stop
      end if

      lor = sqrt(1.0 + u2)

      chi = (w - qu(1)*(1.0+u2/(lor+1.0)))/(1.0+u2)

      ! ideal gas case        
      pg     = (gamma - 1.0)/gamma * chi
      dpdchi = (gamma - 1.0)/gamma
      dpdrho = 0.0

      f = w - pg - qu(4)  ! f(w) = 0

      if(abs(f) .lt. eps) return

      dv2dw  = lor/w**3*m2
      dchidw = 1.0/lor**2 + dv2dw*(qu(1)+2.0*lor*chi)

      drhodw = dv2dw*qu(1)

      dpdw = dpdchi*dchidw + dpdrho*drhodw

      dfdw   = 1.0 - dpdw    ! df/dw

      w  = w0 - mu * f / dfdw                       ! Newton-raphson iteration

      if(abs(w-w0).lt.eps) return

      w0 = w

    end do

    if(mu .gt. 0.1) then
      mu = mu/2.0
      goto 100
    end if

  end subroutine newrap
!------------------------------------------------------------------------------!