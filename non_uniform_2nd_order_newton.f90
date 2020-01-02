MODULE non_uniform_newton

  USE lapack_precision
  USE lapack_interfaces
  USE kinds
  USE generate_grid
  USE Fornberg_coeff

  IMPLICIT NONE

  CONTAINS




  ! FUNCTION Jacobian_approx(n,h_1,h_2,a,b,c,theta,dt,dx,A_bar,Ca) RESULT(Jac)
  !   INTEGER, INTENT(IN) :: n
  !   REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h_1, h_2
  !   REAL(KIND=dbl), DIMENSION(2*n+9,6), INTENT(IN) :: a
  !   REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: b
  !   REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: c
  !   REAL(KIND=dbl), DIMENSION(2*n+8), INTENT(IN) :: dx
  !   REAL(KIND=dbl), INTENT(IN) :: theta, dt, A_bar, Ca
  !   REAL(KIND=dbl) :: Jac
  ! END FUNCTION Jacobian_approx

!!!!!!!!!!!!!!!!!! Create Coefficients using Fornberg scheme !!!!!!!!!!!!!!!!!!

  FUNCTION coef_3(x) RESULT(c)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(6) :: x_temp
  REAL(KIND=dbl), DIMENSION(3+1,6,6) :: c_temp
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: c
  INTEGER :: i, j, n

  n = size(x)
  ALLOCATE(c(n,6))
  x_temp = 0.0_dbl
  c = 0.0_dbl
  c_temp = 0.0_dbl

  DO i=4,n-3
    x_temp = 0.0_dbl
    DO j=1,3
      x_temp(j) = x(i-4+j)
      x_temp(7-j) = x(i+4-j)
    END DO
    c_temp = FB_coef(3,6,x(i),x_temp)
    DO j=1,6
      c(i,j) = c_temp(4,6,j)
    END DO
  END DO

  END FUNCTION

  FUNCTION coef_1(x) RESULT(c)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(3) :: x_temp
  REAL(KIND=dbl), DIMENSION(1+1,3,3) :: c_temp
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: c
  INTEGER :: i, j, n

  n = size(x)
  ALLOCATE(c(n,3))
  x_temp = 0.0_dbl
  c = 0.0_dbl
  c_temp = 0.0_dbl

  DO i=4,n-3
    x_temp = 0.0_dbl
    DO j=1,3
      x_temp(j) = x(i-2+j)
    END DO
    c_temp = FB_coef(1,3,x(i),x_temp)
    DO j=1,3
      c(i,j) = c_temp(2,3,j)
    END DO
  END DO
  
  END FUNCTION coef_1

  FUNCTION rhs_f_non(h,a,b,Ca,A_bar) RESULT(f)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h
  REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: a
  REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: b
  REAL(KIND=dbl), INTENT(IN) :: A_bar, Ca
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f
  INTEGER :: i, n

  n = size(h)
  ALLOCATE(f(n))
  f = 0.0_dbl
  DO i=5,n-4
    f(i) = b(i,1)*(h(i-1)*h(i-1)*h(i-1)/3.0_dbl/Ca*(a(i-1,1)*h(i-4)+a(i-1,2)*h(i-3)+a(i-1,3)*h(i-2) &
            +a(i-1,4)*h(i)+a(i-1,5)*h(i+1)+a(i-1,6)*h(i+2))+A_bar/h(i-1)*(b(i-1,1)*h(i-2) &
            +b(i-1,3)*h(i))+A_bar*b(i-1,2)) &
            +b(i,2)*(h(i)*h(i)*h(i)/3.0_dbl/Ca*(a(i,1)*h(i-3)+a(i,2)*h(i-2)+a(i,3)*h(i-1) &
            +a(i,4)*h(i+1)+a(i,5)*h(i+2)+a(i+1,6)*h(i+3))+A_bar/h(i)*(b(i,1)*h(i-1) &
            +b(i,3)*h(i+1)) + A_bar*b(i,2)) &
            +b(i,3)*(h(i+1)*h(i+1)*h(i+1)/3.0_dbl/Ca*(a(i+1,1)*h(i-2)+a(i+1,2)*h(i-1)+a(i+1,3)*h(i) &
            +a(i+1,4)*h(i+2)+a(i+1,5)*h(i+3)+a(i+1,6)*h(i+4))+A_bar/h(i+1)*(b(i+1,1)*h(i) &
            +b(i+1,3)*h(i+2)) + A_bar*b(i+1,2))
    END DO

  END FUNCTION rhs_f_non

  FUNCTION evaluate_f_non(h_1,h_2,a,b,Ca,A_bar,dt,theta) RESULT(f)
    REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h_1,h_2 ! h_1 is the j and h_2 is the j+1
    REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: a
    REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: b
    REAL(KIND=dbl), INTENT(IN) :: theta, dt, A_bar, Ca
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f_1, f_2, f
    INTEGER :: i, n
    
    n = size(h_1)
    ALLOCATE(f_1(n))
    ALLOCATE(f_2(n))
    ALLOCATE(f(n))
    f_1 = rhs_f_non(h_1,a,b,Ca,A_bar)
    f_2 = rhs_f_non(h_2,a,b,Ca,A_bar)
    DO i=5,n-4
    f(i) = h_2(i) + dt*theta*f_2(i) - (h_1(i)-dt*(1.0_dbl-theta)*f_1(i))
    END DO
    f(1) = h_2(1) - h_2(n-8)
    f(2) = h_2(2) - h_2(n-7)
    f(3) = h_2(3) - h_2(n-6)
    f(4) = h_2(4) - h_2(n-5)
    f(n-3) = h_2(n-3) - h_2(6)
    f(n-2) = h_2(n-2) - h_2(7)
    f(n-1) = h_2(n-1) - h_2(8)
    f(n) = h_2(n) - h_2(9)
  END FUNCTION evaluate_f_non

  FUNCTION Jacobian_non(h,a,b,Ca,A_bar,dt,theta) RESULT(J)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h
  REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: a
  REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: b
  REAL(KIND=dbl), INTENT(IN) :: theta, dt, A_bar, Ca
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: J
  INTEGER :: i, n

  n = size(h)
  ALLOCATE(J(n,n))
  J = 0.0_dbl
  J(1,1) = 1.0_dbl; J(1,n-8) = -1.0_dbl
  J(2,2) = 1.0_dbl; J(2,n-7) = -1.0_dbl
  J(3,3) = 1.0_dbl; J(3,n-6) = -1.0_dbl
  J(4,4) = 1.0_dbl; J(4,n-5) = -1.0_dbl
  J(n-3,6) = -1.0_dbl; J(n-3,n-3) = 1.0_dbl
  J(n-2,7) = -1.0_dbl; J(n-2,n-2) = 1.0_dbl
  J(n-1,8) = -1.0_dbl; J(n-1,n-1) = 1.0_dbl
  J(n,9) = -1.0_dbl; J(n,n) = 1.0_dbl
  DO i=5,n-4
    J(i,i-4) = dt*theta*(b(i,1)/(3.0_dbl*Ca)*a(i-1,1)*h(i-1)**3)
    J(i,i-3) = dt*theta*(b(i,1)/(3.0_dbl*Ca)*a(i-1,2)*h(i-1)**3 + b(i,2)/(3.0_dbl*Ca)*a(i,1)*h(i)**3)
    J(i,i-2) = dt*theta*(b(i,1)*(a(i-1,3)*h(i-1)**3/(3.0_dbl*Ca)+A_bar/h(i-1)*b(i-1,1)) &
              +b(i,2)/3.0_dbl/Ca*a(i,2)*h(i)**3 + b(i,3)/3.0_dbl/Ca*a(i+1,1)*h(i+1)**3)
    J(i,i-1) = dt*theta*(b(i,1)*(h(i-1)**2.0_dbl*(a(i-1,1)*h(i-4)+a(i-1,2)*h(i-3)+a(i-1,3)*h(i-2)&
              +a(i-1,4)*h(i)+a(i-1,5)*h(i+1)+a(i-1,6)*h(i+2))/Ca-A_bar*(b(i-1,1)*h(i-2)+b(i-1,3)&
              *h(i))/h(i-1)**2) + b(i,2)*(h(i)**3*a(i,3)/3.0_dbl/Ca+A_bar*b(i,1)/h(i)) + b(i,3)&
              *h(i+1)**3*a(i+1,2)/3.0_dbl/Ca)
    J(i,i) = 1.0_dbl + dt*theta*(b(i,1)*(a(i-1,4)/3.0_dbl/Ca*h(i-1)*h(i-1)*h(i-1)+A_bar/h(i-1)&
             *b(i-1,3))+b(i,2)*(1.0_dbl/Ca*h(i)*h(i)*(a(i,1)*h(i-3)+a(i,2)*h(i-2)+a(i,3)*h(i-1)& 
             +a(i,4)*h(i+1)+a(i,5)*h(i+2)+a(i,6)*h(i+3))-A_bar/h(i)/h(i)*(b(i,1)*h(i-1)&
             +b(i,3)*h(i+1))) &
             +b(i,3)*(a(i+1,3)/3.0_dbl/Ca*h(i+1)*h(i+1)*h(i+1)+A_bar/h(i+1)*b(i+1,1)))
    J(i,i+1) = dt*theta*(b(i,1)/3.0_dbl/Ca*a(i-1,5)*h(i-1)*h(i-1)*h(i-1) &
             +b(i,2)*(a(i,4)/3.0_dbl/Ca*h(i)*h(i)*h(i)+A_bar/h(i)*b(i,3)) &
             +b(i,3)*(1.0_dbl/Ca*h(i+1)*h(i+1)*(a(i+1,1)*h(i-2)+a(i+1,2)*h(i-1)+a(i+1,3)*h(i) &
             +a(i+1,4)*h(i+2)+a(i+1,5)*h(i+3)+a(i+1,6)*h(i+4))-A_bar/h(i+1)/h(i+1)*(b(i+1,1)&
             *h(i)+b(i+1,3)*h(i+2))))
    J(i,i+2) = (b(i,1)/3.0_dbl/Ca*a(i-1,6)*h(i-1)*h(i-1)*h(i-1) &
             +b(i,2)/3.0_dbl/Ca*a(i,5)*h(i)*h(i)*h(i)) &
             +b(i,3)*(a(i+1,4)/3.0_dbl/Ca*h(i+1)*h(i+1)*h(i+1)+A_bar/h(i+1)*b(i+1,3))*dt*theta
    J(i,i+3) = (b(i,2)/3.0_dbl/Ca*a(i,6)*h(i)*h(i)*h(i) + b(i,3)/3.0_dbl/Ca*a(i+1,5)*h(i+1)*h(i+1)*h(i+1))*dt*theta
    J(i,i+4) = (b(i,3)/3.0_dbl/Ca*a(i+1,6)*h(i+1)*h(i+1)*h(i+1))*dt*theta
  END DO


 END FUNCTION Jacobian_non 


!!!!!!!!!!!!!!!!!!! Explicit Euler for initial guess !!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION explicit_euler_non(h_1,a,b,Ca,A_bar,dt) RESULT(h_2)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h_1
  REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: a
  REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: b
  REAL(KIND=dbl), INTENT(IN) :: A_bar, Ca, dt
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: f,h_2
  INTEGER :: i, n

  n = size(h_1)
  ALLOCATE(f(n))
  ALLOCATE(h_2(n))
  f = rhs_f_non(h_1,a,b,Ca,A_bar) 
  DO i=5,n-4
    h_2(i) = h_1(i) - f(i)*dt
  END DO
  DEALLOCATE(f)
  h_2(1) = h_2(n-8); h_2(2) = h_2(n-7); h_2(3) = h_2(n-6); h_2(4) = h_2(n-5)
  h_2(n-3) = h_2(6); h_2(n-2) = h_2(7); h_2(n-1) = h_2(8); h_2(n) = h_2(9)
  END FUNCTION explicit_euler_non

  FUNCTION non_uniform_2nd_order_newton(h_1,a,b,Ca,A_bar,dt,theta,tol) RESULT(h_2)
    REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: h_1
    REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: a
    REAL(KIND=dbl), DIMENSION(:,:), INTENT(IN) :: b
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: h_2, f, h_temp
    REAL(KIND=dbl), INTENT(IN) :: A_bar, Ca, dt, tol, theta 
    REAL(KIND=dbl) :: e, e_temp
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: work(:)
    INTEGER :: info, lwork, i, j, n
    REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Jac, J_inv, LU

    n = size(h_1)
    ALLOCATE(h_2(n))
    h_2 = explicit_euler_non(h_1,a,b,Ca,A_bar,dt)
    e = 1.0_dbl
    e_temp = 1.0_dbl
    DO WHILE (e>tol)
      !e = 1.0_dbl
      ALLOCATE(Jac(n,n))
      Jac = Jacobian_non(h_2,a,b,Ca,A_bar,dt,theta) 
      DO i=1,n
        DO j=1,n
          IF (Jac(i,j)/=Jac(i,j)) THEN
            PRINT*, "JAC AT ", i, j, "IS NAN"
          END IF
        END DO
      END DO
      !PRINT*, "JAC AT 6, 10 IS ", Jac(6,10)
      !PRINT*, "JAC AT 7, 10 IS ", Jac(7,9)

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
      J_inv = LU
      CALL dgetri(n,J_inv,n,ipiv,work,lwork,info)
      DEALLOCATE(work)
      DEALLOCATE(Jac)
      DEALLOCATE(LU)
      DEALLOCATE(ipiv)
      ALLOCATE(f(n))
      f = evaluate_f_non(h_1,h_2,a,b,Ca,A_bar,dt,theta)
      ALLOCATE(h_temp(n))
      h_temp = h_2 - matmul(J_inv,f)
      h_2 = h_temp
      f = evaluate_f_non(h_1,h_2,a,b,Ca,A_bar,dt,theta)
      e = norm2(f)
      DEALLOCATE(h_temp)
      DEALLOCATE(f)
      DEALLOCATE(J_inv)
      IF (e_temp<e) THEN
        PRINT*, e
      END IF
      e_temp = e
    END DO
  END FUNCTION non_uniform_2nd_order_newton
END MODULE non_uniform_newton
















