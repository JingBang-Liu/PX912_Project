! The purpose of this program is to use finite difference method to solve
! the thin film equations

MODULE newton_iteration

  IMPLICIT NONE

  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(15, 307)
  INTEGER, PARAMETER :: int64 = SELECTED_INT_KIND(15)

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
    fun(n+4) = hh(1,n+4) - hh(1,5)
    fun(n+3) = hh(1,n+3) - hh(1,4)
    ! forward euler 
    fun(3) = c(1) - (1-theta)*2.0_dbl*mu*(hh(1,4)**3*(hh(1,6)-2.0_dbl*hh(1,5)+2.0_dbl*hh(1,3)-hh(1,2)) &
            - hh(1,3)**3*(hh(1,5)-2.0_dbl*hh(1,4)+2.0_dbl*hh(1,2)-hh(1,1))) - hh(1,3)
    ! backward euler
    fun(n+2) = c(n) - (1-theta)*2.0_dbl*mu*(hh(1,n+2)**3*(hh(1,n+4)-2.0_dbl*hh(1,n+3)+2.0_dbl*hh(1,n+1)-hh(1,n)) &
              - hh(1,n+1)**3*(hh(1,n+3)-2.0_dbl*hh(1,n+2)+2.0_dbl*hh(1,n)-hh(1,n-1))) - hh(1,n+2)
    ! midpoint 
    do i=4,n+1
      fun(i) = c(i-2) - (1-theta)*mu*(hh(1,i+1)**3*(hh(1,i+3)-2.0_dbl*hh(1,i+2)+2.0_dbl*hh(1,i)-hh(1,i-1)) &
              - hh(1,i-1)**3*(hh(1,i+1)-2.0_dbl*hh(1,i)+2.0_dbl*hh(1,i-2)-hh(1,i-3)))
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
    ! forward euler
    Jac(3,1) = -(1-theta)*2.0_dbl*mu*hh(1,3)**3; Jac(3,2) = (1-theta)*2.0_dbl*mu*(hh(1,4)**3+2.0_dbl*hh(1,3)**3)
    Jac(3,3) = (1-theta)*2.0_dbl*mu*(3.0_dbl*hh(1,3)**2*(hh(1,5)-2.0_dbl*hh(1,4)+2.0_dbl*hh(1,2)-hh(1,1))-2.0_dbl*hh(1,4)**3) &
              - 1.0_dbl
    Jac(3,4) = -(1-theta)*2.0_dbl*mu*(3.0_dbl*hh(1,4)**2*(hh(1,6)-2.0_dbl*hh(1,5)+2.0_dbl*hh(1,3)-hh(1,2))+2.0_dbl*hh(1,3)**3)
    Jac(3,5) = (1-theta)*2.0_dbl*mu*(2.0_dbl*hh(1,4)**3+hh(1,3)**3); Jac(3,6) = -(1-theta)*2.0_dbl*mu*hh(1,4)**3
    ! backward euler
    Jac(n+2,n-1) = -(1-theta)*2.0_dbl*mu*hh(1,n+1)**3; Jac(n+2,n) = (1+theta)*2.0_dbl*mu*(hh(1,n+2)**3+2.0_dbl*hh(1,n+1)**3)
    Jac(n+2,n+1) = (1-theta)*2.0_dbl*mu*(3.0_dbl*hh(1,n+1)*(hh(1,n+3)-2.0_dbl*hh(1,n+2)+2.0_dbl*hh(1,n)-hh(1,n-1)) &
                  - 2.0_dbl*hh(1,n+2)**3)
    Jac(n+2,n+2) = -(1-theta)*2.0_dbl*mu*(3.0_dbl*hh(1,n+2)**2*(hh(1,n+4)-2.0_dbl*hh(1,n+3)+2.0_dbl*hh(1,n+1)-hh(1,n)) &
                  + 2.0_dbl*hh(1,n+1)**3) - 1
    Jac(n+2,n+3) = (1-theta)*2.0_dbl*mu*(2.0_dbl*hh(1,n+2)**3+hh(1,n+1)**3); Jac(n+2,n+4) = -(1-theta)*2.0_dbl*mu*hh(1,n+2)**3
    ! midpoint
    do i=4,n+1
      Jac(i,i-3) = -(1-theta)*mu*hh(1,i-1)**3; Jac(i,i-2) = 2.0_dbl*(1-theta)*mu*hh(1,i-1)**3
      Jac(i,i-1) = (1-theta)*mu*(hh(1,i+1)**3+3.0_dbl*hh(1,i-1)**2*(hh(1,i+1)-2.0_dbl*hh(1,i)+2.0_dbl*hh(1,i-2)-hh(1,i-3)))
      Jac(i,i) = -(1-theta)*mu*2.0_dbl*(hh(1,i+1)**3+hh(1,i-1)**3)
      Jac(i,i+1) = -(1-theta)*mu*(3.0_dbl*hh(1,i+1)**2*(hh(1,i+3)-2.0_dbl*hh(1,i+2)+2.0_dbl*hh(1,i)-hh(1,i-1))-hh(1,i-1)**3)
      Jac(i,i+2) = (1-theta)*mu*2.0_dbl*hh(1,i+1)**3; Jac(i,i+3) = -(1-theta)*mu*hh(1,i+1)**3
    end do 
  END FUNCTION Jac

END MODULE newton_iteration

PROGRAM MAIN

USE lapack_precision
USE lapack_interfaces
USE newton_iteration

  IMPLICIT NONE

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0
  INTEGER, PARAMETER :: n = 200
  INTEGER(KIND=int64), PARAMETER :: m = 300 
  REAL(KIND=dbl) :: dt = 1e-10_dbl ! time step size
  REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi/w ! boundary
  REAL(KIND=dbl) :: dx = (upper_bnd-lower_bnd)/real(n-1) ! grid step size
  REAL(KIND=dbl), DIMENSION(n) :: x ! grid array
  REAL(KIND=dbl), DIMENSION(m,n+4) :: h ! matrix to store solutions
  INTEGER :: k = 1 ! counter for what time we are on
  INTEGER :: i,j
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl
  REAL(KIND=dbl) :: mu, theta = 0.5_dbl
  REAL(KIND=dbl), DIMENSION(2,n+4) :: hh ! 2x(n+4) array to store results from newton iteration
  REAL(KIND=dbl), DIMENSION(n+4,n+4) :: J_0, J_inv, LU ! Jacobian matrix and it's inverse
  REAL(KIND=dbl), DIMENSION(n) :: c ! c array
  REAL(KIND=dbl), DIMENSION(n+4) :: f ! array of function
  INTEGER, DIMENSION(n+4) :: ipiv
  REAL(KIND=dbl), ALLOCATABLE :: work(:)
  INTEGER :: info, lwork
  REAL(KIND=dbl) :: error
  REAL(KIND=dbl), PARAMETER :: tol_newton = 1e-3_dbl

  

  PRINT*, dx
  PRINT*, dt

  ! initialize grid
  do i=1,n
    x(i) = lower_bnd + dx * (real(i)-1.0_dbl)
  end do

  ! initial condition
  do i=3,n+2
    h(1,i) =1.0_dbl + 0.1_dbl * cos(x(i-2)*w)
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
    mu = dt/(12.0_dbl*dx**4*Ca)
    ! calculate c
    c(1) = h(k-1,3) - 2.0_dbl*mu*theta*(h(k-1,4)**3*(h(k-1,6)-2.0_dbl*h(k-1,5)+2.0_dbl*h(k-1,3)-h(k-1,2)) &
          - h(k-1,3)**3*(h(k-1,5)-2.0_dbl*h(k-1,4)+2.0_dbl*h(k-1,2)+h(k-1,1)))
    c(n) = h(k-1,n+2) - 2.0_dbl*mu*theta*(h(k-1,n+2)**3*(h(k-1,n+4)-2.0_dbl*h(k-1,n+3)+2.0_dbl*h(k-1,n+1)-h(k-1,n)) &
          - h(k-1,n+1)**3*(h(k-1,n+3)-2.0_dbl*h(k-1,n+2)+2.0_dbl*h(k-1,n)-h(k-1,n-1)))
    do i=4,n+1
      c(i-2) = h(k-1,i) - mu*theta*(h(k-1,i+1)**3*(h(k-1,i+3)-2.0_dbl*h(k-1,i+2)+2.0_dbl*h(k-1,i)-h(k-1,i-1)) &
              - h(k-1,i-1)**3*(h(k-1,i+1)-2.0_dbl*h(k-1,i)+2.0_dbl*h(k-1,i-2)+h(k-1,i-3)))
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
    PRINT*, info
    !J_inv = J_0
    CALL dgetri(n+4,J_inv,n+4,ipiv,work,lwork,info) ! inversion
    PRINT*, k
    deallocate(work)
    f = fun(n,hh,c,mu,theta)
    hh(2,:) = hh(1,:) - matmul(J_inv,f)
    hh(1,:) = hh(2,:)
    f = fun(n,hh,c,mu,theta)
    error = norm2(f)
    PRINT*, error
    end do
    h(k,:) = hh(1,:)
    do i=3,n+2
      if ((h(k,i)<0).or.(h(k,i)/=h(k,i))) then
        test = 1
      end if
    end do
  end do

!  PRINT*, h(k,:)
  PRINT*, k

  open(9, file = 'data1.txt',form = 'formatted')
  23 FORMAT(3 (ES23.12E3))

  do i=1,n
    write(9,23) x(i), h(1,i+2), h(k-2,i+2)
  end do

  close(9)
  


  



END PROGRAM MAIN



































