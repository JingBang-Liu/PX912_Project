! The purpose of this program is to use finite difference method to solve
! the thin film equations

MODULE Runge_Kutta_4

USE kinds

  IMPLICIT NONE

  CONTAINS
  ! Generates function
  FUNCTION fun(n,hh,mu)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(n+4), INTENT(IN) :: hh
    REAL(KIND=dbl), INTENT(IN) :: mu
    REAL(KIND=dbl), DIMENSION(n+4) :: fun
    INTEGER :: i

    fun(3) = -2.0_dbl*mu*((hh(4)**3)*(hh(6)-2.0_dbl*hh(5)+2.0_dbl*hh(3)-hh(2)) &
             - (hh(3)**3)*(hh(5)-2.0_dbl*hh(4)+2.0_dbl*hh(2)-hh(1)))
    fun(n+2) = fun(1)
    DO i=4,n+1
      fun(i) = -mu*((hh(i+1)**3)*(hh(i+3)-2.0_dbl*hh(i+2)+2.0_dbl*hh(i)-hh(i-1)) &
              - (hh(i-1)**3)*(hh(i+1)-2.0_dbl*hh(i)+2.0_dbl*hh(i-2)-hh(i-3)))
    END DO
    fun(1) = fun(n); fun(2) = fun(n+1)
    fun(n+3) = fun(4); fun(n+4) = fun(5)
  END FUNCTION fun

  FUNCTION RK4(n,hh,mu,dt)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2,n+4), INTENT(IN) :: hh
    REAL(KIND=dbl), INTENT(IN) :: mu, dt
    REAL(KIND=dbl), DIMENSION(n+4) :: k1,k2,k3,k4,RK4
    
    k1 = dt*fun(n,hh(1,:),mu)
    k2 = dt*fun(n,hh(1,:)+k1/2.0_dbl,mu)
    k3 = dt*fun(n,hh(1,:)+k2/2.0_dbl,mu)
    k4 = dt*fun(n,hh(1,:)+k3,mu)

    RK4 = hh(1,:) + (k1+2.0_dbl*k2+2.0_dbl*k3+k4)/6.0_dbl
  END FUNCTION RK4
  
END MODULE Runge_Kutta_4

PROGRAM MAIN

USE kinds
USE command_line
USE Runge_Kutta_4

  IMPLICIT NONE

  LOGICAL :: success, exists

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0
  INTEGER :: n! = 51
  INTEGER, PARAMETER :: m = 100000 
  REAL(KIND=dbl) :: dt! = 1e-3_dbl ! time step size
  REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi/w ! boundary
  REAL(KIND=dbl) :: dx != (upper_bnd-lower_bnd)/real(n-1) ! grid step size
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: x ! grid array
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: hh
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: h ! matrix to store solutions
  INTEGER :: k = 1 ! counter for what time we are on
  INTEGER :: i,j
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl
  REAL(KIND=dbl) :: mu

  CALL parse_args

  ! Try to grab n, the grid step size
  success = get_arg("n", n, exists=exists)

  ! Try to grab dt, the time step size
  success = get_arg("dt", dt, exists=exists)

  ALLOCATE(x(n))
  ALLOCATE(h(m,n+4))
  ALLOCATE(hh(2,n+4))
  
  dx = (upper_bnd-lower_bnd)/real(n-1)

  PRINT*, dx
  PRINT*, dt

  ! initialize grid
  do i=1,n
    x(i) = lower_bnd + dx * (real(i)-1.0_dbl)
  end do

  ! initial condition
  do i=3,n+2
    h(1,i) =1.0_dbl + 0.9_dbl * cos(x(i-2)*w)
  end do

!  PRINT*, x

  ! boundary condition
  h(1,1) = h(1,n)
  h(1,2) = h(1,n+1)
  h(1,n+4) = h(1,5)
  h(1,n+3) = h(1,4)

  do while ((k<m).and.(test==0))
    k = k + 1
    ! calculate mu
    mu = dt/(Ca*12.0_dbl*(dx**4))
    ! initiate hh
    hh = 0.0_dbl
    hh(1,:) = h(k-1,:)
    hh(2,:) = RK4(n,hh,mu,dt)
    h(k,:) = hh(2,:)
    do i=3,n+2
      if (h(k,i)<0) then
        test = 1
      end if
    end do
  end do

!  PRINT*, h(k,:)
  PRINT*, k

  open(9, file = 'data1.txt',form = 'formatted')
  23 FORMAT(3 (ES23.12E3))

  do i=1,n
    write(9,23) x(i), h(1,i+2), h(k-1,i+2)
  end do

  close(9)
  


  



END PROGRAM MAIN



































