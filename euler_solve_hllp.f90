!=======================================================================
!   This program solves the Euler equations with the Lax Method
!=======================================================================

!=======================================================================
!   This module contains global variables
module globals
  implicit none
  !
  !   This is the number of points used to discretize X
!----------------------
! esto si se modifica
!----------------------  
  real, parameter :: tmax= 0.5            ! maximumn integration time
  integer, parameter :: nx=1000  
  ! define courant
  integer, parameter :: sn = 5 !numero de capturas (snapshot) 
  
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=1.0, dx=xmax/float(nx)
  real, parameter :: gamma=5./3.

  !real, parameter :: Co=0.9 !Courant number

  real, parameter :: rhol=1.0
  real, parameter :: vl=0.75
  real, parameter :: pl=1.

  real, parameter :: rhor=0.125
  real, parameter :: vr=0.
  real, parameter :: pr=0.1
 
  
  real, parameter :: xc = 0.3

!----------------------
! esto no se modifica
!----------------------  
  real, parameter :: dtprint=tmax/float(sn)      ! interval between outputs
  !   This is a vector that contains u(x)
  integer, parameter :: nequ=3          !numero de ecuaciones a manejar
  real :: u(nequ,0:nx+1) !vector importante
  
  !
end module globals
!=======================================================================



!   main program
  program euler_lax
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real           :: time, dt             !  t, $\Delta t$
  real           :: tprint               ! time of next output
  integer        :: itprint              ! number of current output

  ! This subroutine generates the initial conditions
  call initconds(time, tprint, itprint)

  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)

  ! output at tprint intervals
    if(time.ge.tprint) then
      write(*,*) time,tmax,dt
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the $\Delta t$ allowed by the CFL criterium
    call courant(dt)
    !
    ! Integrate u fom t to t+dt
    call ulax(dt,time)
    ! time counter increases
    time=time+dt

  end do

  stop
end program euler_lax

!=======================================================================
! generates initial condition


subroutine initconds(time, tprint, itprint)
 use globals
  implicit none
  !real ::sl, sr
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here is the Sod tube test
  integer :: i
  real :: x

  !call wavespeeds(vl,vr, sl,sr, rhol, rhor, pl, pr, gamma)

  !  fill the vector u
  do i=0,nx+1
    x=float(i)*dx   ! obtain the position $x_i$
    if (x < xc ) then
      u(1,i)=rhol
      u(2,i)=u(1,i)*vl
      u(3,i)=(1.0/2.0)*u(2,i)*vl+pl/(gamma-1.)


    else
      u(1,i)=rhor
      u(2,i)=u(1,i)*vr
      u(3,i)=(1.0/2.0)*u(2,i)*vr+pr/(gamma-1.)
    end if
  end do

  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0

  return
end subroutine initconds
!=======================================================================
! output to file


subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  character (len=20) file1
  integer :: i
  real :: rho,vx,P
  ! open output file
  write(file1,'(a,i2.2,a)') 'euler_lax-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and rho, u(=vx) and P
  do i=1,nx
    rho=u(1,i)
    vx=u(2,i)/rho
    P=(u(3,i)-0.5*rho*vx**2)*(gamma-1.)
    write(10,*) float(i)*dx,rho, vx, P
  end do

  ! closes output file
  close(10)

  return
end subroutine output
!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine courant(dt)
  use globals
  implicit none
  real, intent(out) ::dt!out value
  ! Courant number =0.9
  real, parameter :: Co=0.8
  real :: rho, vx , P, cs
  integer :: i

  dt=1E30
  do i=0,nx+1
    rho=u(1,i)
    vx=u(2,i)/rho
    P=(u(3,i)-0.5*rho*vx**2)*(gamma-1.)
    cs=sqrt(gamma*P/rho) ! sound velocity
    dt=min( dt,Co*dx/(abs(vx)+cs) )
  end do

  return
end subroutine courant
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine ulax(dt,time)
  use globals
  implicit none
  real, intent(in) :: dt, time
  real :: up(nequ,0:nx+1), f(nequ,0:nx+1), fhll(nequ,0:nx+1)
  real :: dtx 
  integer :: i

  !  obtain the fluxes
 ! print*, 'dt',dt
  !call fluxesHLL(fhll, time) !flujos usando HLL
  call fluxesLF(nx,nequ,gamma,u,f) !flujos usando Lax–Friedrichs
  
  !   Here is the Lax method, notice that the values at the extremes can
  !   not be calculated, we need to enter then as boundary conditions
  
  dtx=dt/dx

!!  print*, 'u_antes_lax',u(1,30), u(1,29) , 'dtx',dtx
  do i=1,nx
     !up(:,i)=u(:,i)+dtx*(fhll(:,i-1)-fhll(:,i+1)) !esto es para HLL
     up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1))) !esto es para Lax–Friedrichs

    !print*,i, up(1,i)+dtx, dt, dtx  
  end do

  
!!  print*, 'fhll_ulax_d', fhll(1,31), 'fhll_ulax_a', fhll(1 ,29)

  !   Boundary conditions to the U^n+1
  call boundaries(nx,nequ,up)

  ! copy the up to the u
  u(:,:)=up(:,:)

  
  

  return

end subroutine ulax
!=======================================================================
!Obtain the fluxes F

subroutine fluxesLF(nx,nequ,gamma,u,f)
  !use globals
  implicit none
  integer, intent(in) :: nx,nequ
  real, intent(in) :: gamma
  real, intent(in) :: u(nequ,0:nx+1)
  real, intent(out):: f(nequ,0:nx+1)
  integer :: i
  real :: rho, vx , P

  do i=0,nx+1
    rho=u(1,i)
    vx= u(2,i)/rho
    P= (u(3,i)-0.5*rho*vx**2)*(gamma-1.)

    f(1,i)=rho*vx
    f(2,i)=rho*vx**2+P
    f(3,i)=vx*(u(3,i)+P)

  end do




  return
end subroutine fluxesLF


!=======================================================================
! Set boundary conditions
subroutine boundaries(nx,nequ,u)
  !use globals
  implicit none
  integer, intent(in) :: nx,nequ
  real,    intent(out):: u(nequ,0:nx+1)
  integer :: i

  !   Periodic boundary conditions
  !u(:,0 )=u(:,nx)
  !u(:,nx+1)=u(:,1)

  !  open (outlow) boundary conditions
  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)

  return
end subroutine boundaries
!=======================================================================

subroutine RL(nx, nequ, gamma, u, rho_l, rho_r, vx_l, vx_r, P_l, P_r, time)
  !Este modulo, digamos es para seleccionar, las densidades, velocidades y presiones tanto de izquierda como de derecha
  !Tambien sirve como desacoplamiento para mis primitivas
  implicit none
  integer, intent(in) :: nx,nequ
  real, intent(in) :: gamma
  real, intent(in) :: time
  real, intent(in) :: u(nequ,0:nx+1)
  real, intent(out):: rho_l(1,0:nx+1), rho_r(1,0:nx+1), vx_l(1,0:nx+1), vx_r(1,0:nx+1), P_l(1,0:nx+1), P_r(1,0:nx+1)
  real :: rho_promF, rho_promI, vx_promF, vx_promI, P_promF, P_promI
  integer :: i
  real :: rho, vx , P

 
!Tuve que convertir mis variables primitivas en vectores para tener un mejor control de lo que toca derecha a izquirda,
! no se si es lo mejor
! notar que corre de 0 a nx 
  do i=1,nx    
    rho_l(1,i)= u(1,i-1)
    vx_l(1,i) = u(2,i-1)/u(1,i-1)
    P_l(1,i)  = (u(3,i-1)-0.5*u(1,i-1)*(u(2,i-1)/u(1,i-1))*(u(2,i-1)/u(1,i-1)))*(gamma-1.)

    rho_r(1,i)= u(1,i+1)
    vx_r(1,i) = u(2,i+1)/u(1,i+1)
    P_r(1,i)  = (u(3,i+1)-0.5*u(1,i+1)*(u(2,i+1)/u(1,i+1))*(u(2,i+1)/u(1,i+1)))*(gamma-1.)
  end do
  
  ! Tuve que definir unas "condiciones de frontera" para los valores en los extremos
  ! ya que es dificil hacerlo en un bucle o bueno a lo mejor se puede hacer todo en  un bucle 
  ! pero todavia no se como se hace eso
  ! De todos modos siento yo que no tiene porque interferir

  !Aqui se esta definiendo outflow en los dos bordes de x, ie en x=0 y en x=nx+1

  rho_promI = 0.5*(rho_r(1,1)+rho_l(1,1))
  rho_promF = 0.5*(rho_r(1,nx)+rho_l(1,nx))

  vx_promI = 0.5*(vx_r(1,1)+vx_l(1,1))
  vx_promF = 0.5*(vx_r(1,nx)+vx_l(1,nx))

  P_promI = 0.5*(P_r(1,1)+P_l(1,1))
  P_promF = 0.5*(P_r(1,nx)+P_l(1,nx))

  rho_r(1,nx+1)=rho_promF
  rho_l(1,nx+1)=rho_promF
  
  rho_r(1,0)= rho_promI
  rho_l(1,0)= rho_promI 

  !ojo corregir  
  vx_r(1,nx+1)=vx_promF
  vx_l(1,nx+1)=vx_promF

  vx_r(1,0)= vx_promI
  vx_l(1,0)= vx_promI    

  !ojo corregir  
  P_r(1,nx+1)=P_promF
  P_l(1,nx+1)=P_promF

  P_r(1,0)= P_promI
  P_l(1,0)= P_promI      
!do i=0, 101

!! print*, '===================primitivas=================='
!! print*,'rho_r', rho_r(1,30),'rho_l', rho_l(1,30)
!! print*, 'vx_r', vx_r(1,30), 'vx_l', vx_l(1,30)
!! print*, 'P_r', P_r(1,30),  'P_l', P_l(1,30)
!! print*, '===================primitivas=================='
!end do
  return

end subroutine RL
!=======================================================================


 !calculate the wave speeds
   subroutine wavespeeds(s_l, s_r, rho_l, rho_r, vx_l, vx_r, P_l, P_r)
    ! Este modulo calcula las velocidadesd de la onda
     use globals
     implicit none
!    integer, intent(in)::nx,nequ
    real :: time    
    real, intent(out)::s_l(1,0:nx+1),s_r(1,0:nx+1) !Regresa las velocidades de la onda (s_l, s_r)
    real:: rho_l(1,0:nx+1), rho_r(1,0:nx+1), vx_l(1,0:nx+1), vx_r(1,0:nx+1), P_l(1,0:nx+1), P_r(1,0:nx+1)
    
    real:: cs_l(1,0:nx+1), cs_r(1,0:nx+1) !velocidades del sonido
    real:: cs_promF, cs_promI, s_promF, s_promI
    integer:: i
 
 ! bucle para crear las velocidades del sonido tanto de izquierda como de derecha
  do i=1,nx
   cs_l(1,i)=sqrt(gamma*P_l(1,i)/rho_l(1,i))    
   cs_r(1,i)=sqrt(gamma*P_r(1,i)/rho_r(1,i))
  end do

  cs_promI = 0.5*(cs_r(1,1)+cs_l(1,1))
  cs_promF = 0.5*(cs_r(1,nx)+cs_l(1,nx))
  
  !Igual una especie de condicienes de frontera
  cs_l(1,nx+1)=cs_promF
  cs_r(1,nx+1)=cs_promF
  
  cs_l(1,0)=cs_promI
  cs_r(1,0)=cs_promI

  !print*, 'cs_r', cs_r(1,100)


  !Este bucle calcula las velocidades de las ondas 
  do i=1,nx
   s_l(1,i)=vx_l(1,i)-cs_l(1,i) !velocidad izquierda - velocidad del sonido izquierda
   s_r(1,i)=vx_r(1,i)+cs_r(1,i) ! velocidad derecha + velocidad del sonido derecha
  end do

  s_promI = 0.5*(s_r(1,1)+s_l(1,1))
  s_promF = 0.5*(s_r(1,nx)+s_l(1,nx))

  !ojo corregir  

  !s_l(1,nx+1)=s_promF
  !s_r(1,nx+1)=s_promF

  s_l(1,nx+1)=s_l(1, nx)
  s_r(1,nx+1)=s_r(1, nx)
  
  !s_l(1,0)=s_promI
  !s_r(1,0)=s_promI

  s_l(1,0)=s_l(1, 0)
  s_r(1,0)=s_r(1, 0)
  

!!  print*, '==========================Sonido ===================='
!!  print*,'rho_r', rho_r(1,30),'rho_l', rho_l(1,30)
!!  print*, 'vx_r', vx_r(1,30), 'vx_l', vx_l(1,30)
!!  print*, 'P_r', P_r(1,30),  'P_l', P_l(1,30)
!!  print*, '                                                  '
!!  print*, 'cs_r', cs_r(1,30), 'cs_l', cs_l(1,30)
!!  print*, 's_r', s_r(1,30), 's_l', s_l(1,30)
!!  print*, '==========================Sonido ===================='


  !print*,'s_l' ,s_l(1,25),'s_r' ,s_r(1,25)

     return
   end subroutine wavespeeds


   ! ! !=======================================================================
  subroutine fluxesHLL(fhll, time)
    !Este modulo calcula los flujos hll
    use globals
    implicit none
    real, intent(in) :: time
    
    real, intent(out)::fhll(3,0:nx+1)! esta es la variable que se obtiene y va para la subrutina ulax
    real :: rho_l(1,0:nx+1), rho_r(1,0:nx+1), vx_l(1,0:nx+1), vx_r(1,0:nx+1), P_l(1,0:nx+1), P_r(1,0:nx+1) !variables primitivas tanto de derecha
    !como de izquierda
    real :: s_l(1,0:nx+1),s_r(1,0:nx+1) !velocidad de la onda

    real :: f_l(nequ, 0:nx+1), f_r(nequ,0:nx+1), u_l(nequ,0:nx+1), u_r(nequ,0:nx+1), dt
            !flujo de izquierda, flujo de derecha, conservads izquierda, conservadas derecha
    real:: f_promF1, f_promI1, u_promF1, u_promI1
    real:: f_promF2, f_promI2, u_promF2, u_promI2
    real:: f_promF3, f_promI3, u_promF3, u_promI3


    integer :: i

    call RL(nx, nequ, gamma, u, rho_l, rho_r, vx_l, vx_r ,  P_l ,  P_r,time) !llamamos a nuestras primitivas
    call wavespeeds(s_l, s_r, rho_l, rho_r, vx_l, vx_r, P_l, P_r) ! llamamos a nuestras velocidades del sonido
    
      !igual calculamos los flujos y conservadas de lado derecho e izquierdo
  do i=1,nx

      f_l(1,i) = rho_l(1,i)*vx_l(1,i)
      f_l(2,i) = rho_l(1,i)*vx_l(1,i)**2+P_l(1,i)
      f_l(3,i) = vx_l(1,i)*((1.0/2.0)*rho_l(1,i)*vx_l(1,i)**2+P_l(1,i)/(gamma-1)+P_l(1,i))

      f_r(1,i) = rho_r(1,i)*vx_r(1,i)
      f_r(2,i) = rho_r(1,i)*vx_r(1,i)**2+P_r(1,i)
      f_r(3,i) = vx_r(1,i)*((1.0/2.0)*rho_r(1,i)*vx_r(1,i)**2+P_r(1,i)/(gamma-1)+P_r(1,i))

      ! Estas conservadas nos sirven para calcular el fhll
      u_l(1,i) = rho_l(1,i)
      u_l(2,i) = rho_l(1,i)*vx_l(1,i)
      u_l(3,i) = (1.0/2.0)*rho_l(1,i)*vx_l(1,i)*vx_l(1,i)+P_l(1,i)/(gamma-1)

      u_r(1,i) = rho_r(1,i)
      u_r(2,i) = rho_r(1,i)*vx_r(1,i)
      u_r(3,i) = (1.0/2.0)*rho_r(1,i)*vx_r(1,i)*vx_r(1,i)+P_r(1,i)/(gamma-1)

  end do

!!      print*,'======================flujos=================='
 !!     print*,'f_r', f_r(1,30),'f_l', f_l(1,30)
 !!     print*,'u_r', u_r(1,30),'u_l', u_l(1,30)
 !!     print*,'======================flujos=================='
    !ojo corregir
    f_promI1 = 0.5*(f_r(1,1)+f_l(1,1))
    f_promF1 = 0.5*(f_r(1,nx)+f_l(1,nx))

    u_promI1 = 0.5*(u_r(1,1)+u_l(1,1))
    u_promF1 = 0.5*(u_r(1,nx)+u_l(1,nx))

    f_promI2 = 0.5*(f_r(2,1)+f_l(2,1))
    f_promF2 = 0.5*(f_r(2,nx)+f_l(2,nx))

    u_promI2 = 0.5*(u_r(2,1)+u_l(2,1))
    u_promF2 = 0.5*(u_r(2,nx)+u_l(2,nx))

    f_promI3 = 0.5*(f_r(3,1)+f_l(3,1))
    f_promF3 = 0.5*(f_r(3,nx)+f_l(3,nx))

    u_promI3 = 0.5*(u_r(3,1)+u_l(3,1))
    u_promF3 = 0.5*(u_r(3,nx)+u_l(3,nx))
    


    f_l(1,nx+1)=f_promF1
    f_r(1,nx+1)=f_promF1
  
    f_l(1,0)=f_promI1
    f_r(1,0)=f_promI1
 

    u_l(1,nx+1)=u_promF1
    u_r(1,nx+1)=u_promF1
  
    u_l(1,0)=u_promI1
    u_r(1,0)=u_promI1


    f_l(2,nx+1)=f_promF2
    f_r(2,nx+1)=f_promF2
  
    f_l(2,0)=f_promI2
    f_r(2,0)=f_promI2
 

    u_l(2,nx+1)=u_promF2
    u_r(2,nx+1)=u_promF2
  
    u_l(2,0)=u_promI2
    u_r(2,0)=u_promI2 


    f_l(3,nx+1)=f_promF3
    f_r(3,nx+1)=f_promF3
  
    f_l(3,0)=f_promI3
    f_r(3,0)=f_promI3
 

    u_l(3,nx+1)=u_promF3
    u_r(3,nx+1)=u_promF3
  
    u_l(3,0)=u_promI3
    u_r(3,0)=u_promI3

    ! print*,'================================================'
    ! print*, 'u_r', u_r(2,100),'f_r', f_r(2,101), 'f_l', f_l(2,101)
    ! print*, 'u_l', u_l(2,100), 's_l', s_l(1,101)
    ! print*, 's_r', s_r(1,101), 's_l', s_l(1,101)
    !print*, (s_r(1,101)*f_l(2,101)-s_l(1,101)*f_r(2,101)+s_l(1,101)*s_r(1,101)*(u_r(2,101)-u_l(2,101)))/(s_r(1,101)-s_l(1,101))
    !Este bucle calcula fhll
    do i=0,nx+1 !ojo aqui si se necesita defini i=0 y i=nx+1
      if (0 .le. s_l(1,i)) then !less or equal 0<=sl
        fhll(:,i)=f_l(:, i)
        !print*,'a'

      else if (s_l(1,i) .le. 0 .and. 0 .le. s_r(1,i)) then !sl<=0<=sr
        fhll(1,i)=(s_r(1,i)*f_l(1,i)-s_l(1,i)*f_r(1,i)+s_l(1,i)*s_r(1,i)*(u_r(1,i)-u_l(1,i)))/(s_r(1,i)-s_l(1,i))
        fhll(2,i)=(s_r(1,i)*f_l(2,i)-s_l(1,i)*f_r(2,i)+s_l(1,i)*s_r(1,i)*(u_r(2,i)-u_l(2,i)))/(s_r(1,i)-s_l(1,i))
        fhll(3,i)=(s_r(1,i)*f_l(3,i)-s_l(1,i)*f_r(3,i)+s_l(1,i)*s_r(1,i)*(u_r(3,i)-u_l(3,i)))/(s_r(1,i)-s_l(1,i))
        !print*,'b'
        !print*, 'fhll_0', fhll(1,i), i, 's_l', s_l(1,i)
        !print*,'hs',(s_r(1,100)*f_l(1,100)-s_l(1,100)*f_r(1,100)+s_l(1,100)*s_r(1,100)*(u_r(1,100)-u_l(1,100)))&
        !/(s_r(1,100)-s_l(1,100))
      
      else if (s_r(1,i) .le. 0) then !sr<=0
        fhll(:,i)=f_r(:,i)
! !        !write(*,*) i
          !print*,'c'
      endif
    end do
   

!!     print*,'======================fhll=================='
!!      print*,'f_hll31', fhll(1,31),'f_hll29', fhll(1,29), 'f_hll30', fhll(1,30)
!!      !print*,'u_r', u_r(1,31),'u_l', u_l(1,30)
!!      print*,'======================fhll=================='
  
  return

  end subroutine fluxesHLL