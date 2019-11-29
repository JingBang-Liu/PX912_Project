! The purpose of this program is to use finite difference method to solve
! the thin film equations. Here we use Implicit method to solve the problem.

MODULE newton_iteration

USE kinds

  IMPLICIT NONE

  !INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(15, 307)

  CONTAINS
  ! Generates function
  FUNCTION fun(n,hh,c,mu,theta)
    INTEGER, INTENT(IN) :: n ! number of colums
    REAL(KIND=dbl), DIMENSION(2,n+4), INTENT(IN) :: hh ! h is the initial h
    REAL(KIND=dbl), DIMENSION(n), INTENT(IN) :: c ! input c
    REAL(KIND=dbl), INTENT(IN) :: mu, theta  ! input parameters
    REAL(KIND=dbl), DIMENSION(n+4) :: fun ! output
    INTEGER :: i

    ! got from boundary conditions
    fun(1) = hh(1,1) - hh(1,n)
    fun(2) = hh(1,2) - hh(1,n+1)
    fun(n+2) = hh(1,n+2) - hh(1,3)
    fun(n+4) = hh(1,n+4) - hh(1,5)
    fun(n+3) = hh(1,n+3) - hh(1,4)
    ! forward euler 
    !fun(3) = 0.0_dbl
    fun(3) = hh(1,3) - c(1) + theta*2.0_dbl*mu*((hh(1,6)-2.0_dbl*hh(1,5)+2.0_dbl*hh(1,3)-hh(1,2))*hh(1,4)**3 &
            - (hh(1,5)-2.0_dbl*hh(1,4)+2.0_dbl*hh(1,2)-hh(1,1))*hh(1,3)**3)
    ! backward euler
    !fun(n+2) = 0.0_dbl
    !fun(n+2) = hh(1,n+2) - c(n) + theta*2.0_dbl*mu*((hh(1,n+4)-2.0_dbl*hh(1,n+3)+2.0_dbl*hh(1,n+1)-hh(1,n))*hh(1,n+2)**3 &
    !          - (hh(1,n+3)-2.0_dbl*hh(1,n+2)+2.0_dbl*hh(1,n)-hh(1,n-1))*hh(1,n+1)**3)
    ! midpoint 
    do i=4,n+1
      fun(i) = hh(1,i) - c(i-2) + theta*mu*((hh(1,i+3)-2.0_dbl*hh(1,i+2)+2.0_dbl*hh(1,i)-hh(1,i-1))*hh(1,i+1)**3 &
              - (hh(1,i+1)-2.0_dbl*hh(1,i)+2.0_dbl*hh(1,i-2)-hh(1,i-3))*hh(1,i-1)**3)
    end do
  END FUNCTION fun

  ! Generates Jacobian
  FUNCTION Jac(n,hh,mu,theta)
    INTEGER, INTENT(IN) :: n ! number of colums
    REAL(KIND=dbl), DIMENSION(2,n+4), INTENT(IN) :: hh ! h is the inital h
    REAL(KIND=dbl), INTENT(IN) :: mu, theta ! input parameters
    REAL(KIND=dbl), DIMENSION(n+4,n+4) :: Jac ! output matrix
    INTEGER :: i

    ! initialize Jac
    Jac = 0.0_dbl
    ! got from boundary condtions
    Jac(1,1) = 1.0_dbl; Jac(1,n) = -1.0_dbl; Jac(2,2) = 1.0_dbl; Jac(2,n+1) = -1.0_dbl
    Jac(n+4,n+4) = 1.0_dbl; Jac(n+4,5) = -1.0_dbl; Jac(n+3,n+3) = 1.0_dbl; Jac(n+3,4) = -1.0_dbl
    Jac(n+2,n+2) = 1.0_dbl; Jac(n+2,3) = -1.0_dbl
    ! forward euler
    Jac(3,1) = theta*2.0_dbl*mu*hh(1,3)**3; Jac(3,2) = -theta*2.0_dbl*mu*(hh(1,4)**3+2.0_dbl*hh(1,3)**3)
    Jac(3,3) = 1.0_dbl+theta*2.0_dbl*mu*(2.0_dbl*hh(1,4)**3-3.0_dbl*(hh(1,5)-2.0_dbl*hh(1,4)+2.0_dbl*hh(1,2)-hh(1,1))*hh(1,3)**2)
    Jac(3,4) = theta*2.0_dbl*mu*(3.0_dbl*(hh(1,4)**2)*(hh(1,6)-2.0_dbl*hh(1,5)+2.0_dbl*hh(1,3)-hh(1,2))+2.0_dbl*hh(1,3)**3)
    Jac(3,5) = -theta*2.0_dbl*mu*(2.0_dbl*hh(1,4)**3+hh(1,3)**3); Jac(3,6) = theta*2.0_dbl*mu*hh(1,4)**3
    ! backward euler
    !Jac(n+2,n-1) = theta*2.0_dbl*mu*hh(1,n+1)**3; Jac(n+2,n) = -theta*2.0_dbl*mu*(hh(1,n+2)**3+2.0_dbl*hh(1,n+1)**3)
    !Jac(n+2,n+1) = -theta*2.0_dbl*mu*(3.0_dbl*(hh(1,n+1)**2)*(hh(1,n+3)-2.0_dbl*hh(1,n+2)+2.0_dbl*hh(1,n)-hh(1,n-1)) &
    !              - 2.0_dbl*hh(1,n+2)**3)
    !Jac(n+2,n+2) = 1.0_dbl+theta*2.0_dbl*mu*(3.0_dbl*(hh(1,n+2)**2)*(hh(1,n+4)-2.0_dbl*hh(1,n+3)+2.0_dbl*hh(1,n+1)-hh(1,n)) &
    !              + 2.0_dbl*hh(1,n+1)**3)
    !Jac(n+2,n+3) = -theta*2.0_dbl*mu*(2.0_dbl*hh(1,n+2)**3+hh(1,n+1)**3); Jac(n+2,n+4) = theta*2.0_dbl*mu*hh(1,n+2)**3
    ! midpoint
    do i=4,n+1
      Jac(i,i-3) = theta*mu*hh(1,i-1)**3; Jac(i,i-2) = -2.0_dbl*theta*mu*hh(1,i-1)**3
      Jac(i,i-1) = -theta*mu*(hh(1,i+1)**3+3.0_dbl*(hh(1,i-1)**2)*(hh(1,i+1)-2.0_dbl*hh(1,i)+2.0_dbl*hh(1,i-2)-hh(1,i-3)))
      Jac(i,i) = 1.0_dbl+theta*mu*2.0_dbl*(hh(1,i+1)**3+hh(1,i-1)**3)
      Jac(i,i+1) = theta*mu*(3.0_dbl*(hh(1,i+1)**2)*(hh(1,i+3)-2.0_dbl*hh(1,i+2)+2.0_dbl*hh(1,i)-hh(1,i-1))-hh(1,i-1)**3)
      Jac(i,i+2) = -theta*mu*2.0_dbl*hh(1,i+1)**3; Jac(i,i+3) = theta*mu*hh(1,i+1)**3
    end do 
  END FUNCTION Jac

END MODULE newton_iteration

PROGRAM MAIN

USE lapack_precision
USE lapack_interfaces
USE newton_iteration
USE kinds
USE command_line

  IMPLICIT NONE

  LOGICAL :: success, exists

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0
  INTEGER :: n! = 200
  INTEGER :: m = 10000000 
  REAL(KIND=dbl) :: dt != 1e-10_dbl ! time step size
  REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi/w ! boundary
  REAL(KIND=dbl) :: dx! = (upper_bnd-lower_bnd)/real(n-1) ! grid step size
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: x ! grid array
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: h ! matrix to store solutions
  INTEGER :: k = 1 ! counter for what time we are on
  INTEGER :: i,j
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl
  REAL(KIND=dbl) :: mu, theta! = 0.5_dbl
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: hh ! 2x(n+4) array to store results from newton iteration
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: J_0, J_inv, LU ! Jacobian matrix and it's inverse
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: c ! c array
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f ! array of function
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
  REAL(KIND=dbl), ALLOCATABLE :: work(:)
  INTEGER :: info, lwork
  REAL(KIND=dbl) :: error
  REAL(KIND=dbl), PARAMETER :: tol_newton = 1e-2_dbl
!  REAL(KIND=dbl) :: h3, hn2

  CALL parse_args

  ! Try to grab n, the grid step size
  success = get_arg("n", n, exists=exists)

  ! Try to grab dt, the time step size
  success = get_arg("dt", dt, exists=exists)

  ! Try to grab theta, the implicit factor
  success = get_arg("theta", theta, exists=exists)

  ALLOCATE(x(n))
  ALLOCATE(h(m,n+4))
  ALLOCATE(hh(2,n+4))
  ALLOCATE(J_0(n+4,n+4))
  ALLOCATE(J_inv(n+4,n+4))
  ALLOCATE(LU(n+4,n+4))
  ALLOCATE(c(n))
  ALLOCATE(f(n+4))
  ALLOCATE(ipiv(n+4))

  dx = (upper_bnd-lower_bnd)/real(n-1) 

  PRINT*, dx
  PRINT*, dt

  ! initialize grid
  do i=1,n
    x(i) = lower_bnd + dx * (real(i)-1.0_dbl)
  end do

  ! initial condition
  do i=3,n+2
    h(1,i) =1.0_dbl + 0.4_dbl * cos(x(i-2)*w)
  end do

!  PRINT*, x

  ! boundary condition
  h(1,1) = h(1,n)
  h(1,2) = h(1,n+1)
  h(1,n+4) = h(1,5)
  h(1,n+3) = h(1,4)
  !h3 = h(1,3)
  !hn2 = h(1,n+2)


  do while ((k<m).and.(test==0))
    k = k + 1
    ! calculate mu
    mu = dt/(Ca*12.0_dbl*(dx**4))
    ! calculate c
    c = 0.0_dbl
    c(1) = h(k-1,3) - 2.0_dbl*mu*(1.0_dbl-theta)*((h(k-1,6)-2.0_dbl*h(k-1,5)+2.0_dbl*h(k-1,3)-h(k-1,2))*h(k-1,4)**3 &
          - (h(k-1,5)-2.0_dbl*h(k-1,4)+2.0_dbl*h(k-1,2)-h(k-1,1))*h(k-1,3)**3)
    !c(n) = h(k-1,n+2) - 2.0_dbl*mu*(1.0_dbl-theta)*((h(k-1,n+4)-2.0_dbl*h(k-1,n+3)+2.0_dbl*h(k-1,n+1)-h(k-1,n))*h(k-1,n+2)**3 &
    !      - (h(k-1,n+3)-2.0_dbl*h(k-1,n+2)+2.0_dbl*h(k-1,n)-h(k-1,n-1))*h(k-1,n+1)**3)
    do i=4,n+1
      c(i-2) = h(k-1,i) - mu*(1.0_dbl-theta)*((h(k-1,i+3)-2.0_dbl*h(k-1,i+2)+2.0_dbl*h(k-1,i)-h(k-1,i-1))*h(k-1,i+1)**3 &
              - (h(k-1,i+1)-2.0_dbl*h(k-1,i)+2.0_dbl*h(k-1,i-2)-h(k-1,i-3))*h(k-1,i-1)**3)
    end do
    ! Newton iteration
    hh = 0.0_dbl ! initialize hh
    !do i=1,n+4
    !  hh(1,i) = h(k-1,i)
    !end do
    hh(1,:) = h(k-1,:)
    !hh(1,5) = hh(1,5)! - 1e-1_dbl
    error = 1.0_dbl
    do while (error>tol_newton)
      J_0 = Jac(n,hh,mu,theta)
  
      info = 0
      lwork = -1
      allocate(work(lwork))
      work = 0
      ipiv = 0
  
      LU = J_0
      CALL dgetrf(n+4,n+4,LU,n+4,ipiv,info) ! LU factorization
      J_inv = LU
      !PRINT*, info
      !J_inv = J_0
      CALL dgetri(n+4,J_inv,n+4,ipiv,work,lwork,info) ! inversion
      deallocate(work)
      f = fun(n,hh,c,mu,theta)
      hh(2,:) = hh(1,:) - matmul(J_inv,f)
      hh(1,:) = hh(2,:)
      f = fun(n,hh,c,mu,theta)
      error = norm2(f)
      !PRINT*, error
    end do
      h(k,:) = hh(1,:)
    do i=3,n+2
      if ((h(k,i)<0).or.(h(k,i)/=h(k,i))) then
        test = 1
      end if
    end do
    if (mod(k,1000)==0) then
      PRINT*, k
      PRINT*, minval(h(k,:))
    end if 
  end do


  open(9, file = 'data1.txt',form = 'formatted')
  23 FORMAT(3 (ES23.12E3))

  do i=1,n
    write(9,23) x(i), h(1,i+2), h(k,i+2)
  end do

  close(9)
  


  



END PROGRAM MAIN



































