MODULE uniform_newton
 
  USE lapack_precision
  USE lapack_interfaces
  USE kinds
  USE generate_grid

  IMPLICIT NONE

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!! F !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION rhs_f_uni(h,dx,Ca,A_bar) RESULT(f)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h
  REAL(KIND=dbl), INTENT(IN) :: dx, Ca, A_bar
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f
  INTEGER :: i, n

  n = size(h)
  ALLOCATE(f(n))
  f = 0.0_dbl
  DO i=5,n-4
    f(i) = 1.0_dbl/2.0_dbl/dx*(h(i+1)**3/3.0_dbl/Ca*((h(i+3)-2.0_dbl*h(i+2)+2.0_dbl*h(i)-h(i-1))&
            /2.0_dbl/dx**3) + A_bar/h(i+1)*(h(i+2)-h(i))/2.0_dbl/dx - h(i-1)**3/3.0_dbl/Ca*(h(i+1) &
            -2.0_dbl*h(i)+2.0_dbl*h(i-2)-h(i-3))/2.0_dbl/dx**3 - A_bar/h(i-1)*(h(i)-h(i-2))&
            /2.0_dbl/dx)
  END DO
  END FUNCTION rhs_f_uni

!!!!!!!!!!!!!!!!!!!! function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION evaluate_f_uni(h_1,h_2,dx,Ca,A_bar,dt,theta) RESULT(f)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h_1, h_2
  REAL(KIND=dbl), INTENT(IN) :: dx, Ca, A_bar, dt, theta
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f, f_1, f_2
  INTEGER :: i, n 

  n = size(h_1)
  ALLOCATE(f(n))
  ALLOCATE(f_1(n))
  ALLOCATE(f_2(n))
  f_1 = rhs_f_uni(h_1,dx,Ca,A_bar)
  f_2 = rhs_f_uni(h_2,dx,Ca,A_bar)
  f = 0.0_dbl
  f(1) = h_2(1) - h_2(n-8); f(2) = h_2(2) - h_2(n-7); f(3) = h_2(3) - h_2(n-6); f(4) = h_2(4)&
          - h_2(n-5)
  f(n-3) = h_2(n-3) - h_2(6); f(n-2) = h_2(n-2) - h_2(7); f(n-1) = h_2(n-1) - h_2(8); f(n) = h_2(n)&
          - h_2(9)
  DO i=5,n-4
    f(i) = h_2(i) + dt*theta*f_2(i) - (h_1(i)-dt*(1.0_dbl-theta)*f_1(i))
  END DO
  END FUNCTION evaluate_f_uni

!!!!!!!!!!!!!!!!!!! Jacobian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Jacobian_uni(h,dx,Ca,A_bar,dt,theta) RESULT(J)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h
  REAL(KIND=dbl), INTENT(IN) :: dx, Ca, A_bar, dt, theta
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: J
  INTEGER :: i, n

  n = size(h)
  ALLOCATE(J(n,n))
  J = 0.0_dbl
  DO i=5,n-4
    J(i,i-3) = dt*theta/2.0_dbl/dx*h(i-1)**3/6.0_dbl/Ca/dx**3
    J(i,i-2) = dt*theta/2.0_dbl/dx*(-h(i-1)**3/3.0_dbl/Ca/dx**3+A_bar/2.0_dbl/dx/h(i-1))
    J(i,i-1) = dt*theta/2.0_dbl/dx*(-h(i+1)**3/6.0_dbl/Ca/dx**3-h(i-1)**2/Ca*(h(i+1)-2.0_dbl*h(i)&
            +2.0_dbl*h(i-2)-h(i-3))+A_bar/h(i-1)**2*(h(i)-h(i-2))/2.0_dbl/dx)
    J(i,i) = 1.0_dbl + dt*theta/2.0_dbl/dx*(h(i+1)**3/3.0_dbl/Ca/dx**3 - A_bar/2.0_dbl/h(i+1)/dx&
            +h(i-1)**3/3.0_dbl/Ca/dx**3 - A_bar/2.0_dbl/h(i-1)/dx)
    J(i,i+1) = dt*theta/2.0_dbl/dx*(h(i+1)**2/Ca*(h(i+3)-2.0_dbl*h(i+2)+2.0_dbl*h(i)-h(i-1))&
            /2.0_dbl/dx**3 - A_bar/h(i+1)**2*(h(i+2)-h(i))/2.0_dbl/dx - h(i-1)**3/6.0_dbl&
            /Ca/dx**3)
    J(i,i+2) = dt*theta/2.0_dbl/dx*(-h(i+1)**3/3.0_dbl/Ca/dx**3 + A_bar/2.0_dbl/h(i+1)/dx)
    J(i,i+3) = dt*theta/2.0_dbl/dx*h(i+1)**3/6.0_dbl/Ca/dx**3
  END DO
  J(1,1) = 1.0_dbl; J(1,n-8) = -1.0_dbl; J(2,2) = 1.0_dbl; J(2,n-7) = -1.0_dbl
  J(3,3) = 1.0_dbl; J(3,n-6) = -1.0_dbl; J(4,4) = 1.0_dbl; J(4,n-5) = -1.0_dbl
  J(n-3,n-3) = 1.0_dbl; J(n-3,6) = -1.0_dbl; J(n-2,n-2) = 1.0_dbl; J(n-2,7) = -1.0_dbl
  J(n-1,n-1) = 1.0_dbl; J(n-1,8) = -1.0_dbl; J(n,n) = 1.0_dbl; J(n,9) = -1.0_dbl
  END FUNCTION Jacobian_uni

!!!!!!!!!!!!!!!!!!!!!!!! explicit euler !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION explicit_euler_uni(h_1,dx,Ca,A_bar,dt) RESULT(h_2)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h_1
  REAL(KIND=dbl), INTENT(IN) :: dx, Ca, A_bar, dt
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: h_2, f
  INTEGER :: i, n

  n = size(h_1)
  ALLOCATE(h_2(n))
  h_2 = 0.0_dbl
  f = rhs_f_uni(h_1,dx,Ca,A_bar)
  DO i=5,n-4
    h_2(i) = h_1(i) - dt*f(i)
  END DO
  h_2(1) = h_2(n-8); h_2(2) = h_2(n-7); h_2(3) = h_2(n-6); h_2(4) = h_2(n-5)
  h_2(n-3) = h_2(6); h_2(n-2) = h_2(7); h_2(n-1) = h_2(8); h_2(n) = h_2(9)
  END FUNCTION explicit_euler_uni

!!!!!!!!!!!!!!!!!!!!!!! Newton !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION uniform_2nd_order_newton(h_1,dx,Ca,A_bar,dt,theta,tol) RESULT(h_2)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h_1
  REAL(KIND=dbl), INTENT(IN) :: dx, Ca, A_bar, dt, theta, tol
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: h_2, h_temp
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Jac, J_inv, LU
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f
  REAL(KIND=dbl) :: e
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: work(:)
  INTEGER :: info, lwork, n, k

  n = size(h_1)
  ALLOCATE(h_2(n))
  h_2 = explicit_euler_uni(h_1,dx,Ca,A_bar,dt)
  e = 1.0_dbl
  k = 0
  DO WHILE (e>tol)
    !k = k + 1
    !PRINT*, k
    e = 1.0_dbl
    ALLOCATE(Jac(n,n))
    Jac = Jacobian_uni(h_2,dx,Ca,A_bar,dt,theta)
    
    info = 0
    lwork = -1
    ALLOCATE(work(lwork))
    work = 0.0_dbl
    ALLOCATE(ipiv(n))
    ipiv = 0

    ALLOCATE(LU(n,n))
    LU = Jac
    CALL dgetrf(n,n,LU,n,ipiv,info)
    ALLOCATE(J_inv(n,n))
    CALL dgetri(n,J_inv,n,ipiv,work,lwork,info)
    DEALLOCATE(Jac)
    DEALLOCATE(LU)
    DEALLOCATE(work)
    DEALLOCATE(ipiv)
    ALLOCATE(f(n))
    f = evaluate_f_uni(h_1,h_2,dx,Ca,A_bar,dt,theta)
    ALLOCATE(h_temp(n))
    h_temp = h_2 - matmul(J_inv,f)
    h_2 = h_temp
    DEALLOCATE(h_temp)
    f = evaluate_f_uni(h_1,h_2,dx,Ca,A_bar,dt,theta)
    e = norm2(f)
    DEALLOCATE(f)
    DEALLOCATE(J_inv)
    !PRINT*, e
  END DO
  END FUNCTION uniform_2nd_order_newton

END MODULE
