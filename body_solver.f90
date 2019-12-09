PROGRAM MAIN

  USE kinds        
  USE lapack_precision
  USE lapack_interfaces
  USE generate_grid
  USE non_uniform_newton
  USE command_line

  IMPLICIT NONE

  LOGICAL :: success, exists

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0 
  INTEGER :: n ! input from command line
  INTEGER :: m = 100000
  REAL(KIND=dbl) :: dt ! input from command line
  REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi/w
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: x ! grid depend on n
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: dx ! dx depend on x
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl, A_bar = 1.0_dbl
  REAL(KIND=dbl) :: theta ! input from command line
  REAL(KIND=dbl), PARAMETER :: tol_newton = 1e-2_dbl
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: h
  INTEGER :: i, k
  

  CALL parse_args

  ! Try to grab n, the grid steps
  success = get_arg("n",n,exists=exists)

  ! Try to grab dt, the time step size
  success = get_arg("dt",dt,exists=exists)

  ! Try to grab theta, the implicit factor
  success = get_arg("theta",theta,exists=exists)

  !!!! create non-uniform grid
  CALL grid_square(x,dx,n,lower_bnd,upper_bnd)

  ALLOCATE(h(m,2*n+9))

  ! initial condition
  DO i=5,2*n+5
    h(1,i) = 1.0_dbl + 0.4_dbl * cos(x(i)*w)
  END DO

  ! initial condition
  h(1,1) = h(1,2*n+1); h(1,2) = h(1,2*n+2); h(1,3) = h(1,2*n+3); h(1,4) = h(1,2*n+4)
  h(1,2*n+6) = h(1,6); h(1,2*n+7) = h(1,7); h(1,2*n+8) = h(1,8); h(1,2*n+9) = h(1,9)

  DO WHILE ((k<m).and.(test==0))
    k = k + 1
    h(k,:) = non_uniform_2nd_order_newton(n,h(k-1,:),dx,theta,dt,A_bar,Ca,tol_newton)
    DO i=5,2*n+5
      IF ((h(k,i)<0).or.(h(k,i)/=h(k,i))) THEN
        test = 1
      END IF
    END DO
    IF (mod(k,1000)==0) THEN
      PRINT*, k
      PRINT*, minval(h(k,:))
    END IF
  END DO


  !PRINT*, ALLOCATED(x)
  !PRINT*, ALLOCATED(dx)
END PROGRAM

