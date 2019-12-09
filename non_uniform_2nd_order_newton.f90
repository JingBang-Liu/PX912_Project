MODULE non_uniform_newton

  USE lapack_precision
  USE lapack_interfaces
  USE kinds
  USE generate_grid

  IMPLICIT NONE

  CONTAINS



!!!!!!!!!!!!!!!!!!!!!!!!!! Preparation for non uniform newton !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION S1(x,dx)
  REAL(KIND=dbl), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  REAL(KIND=dbl) :: S1

  S1 = dx(5)*(dx(5)-dx(2))*(dx(6)+dx(2)-dx(1))+x*x*(dx(6)-dx(5)+dx(2)-dx(1)) &
      -(dx(6)-dx(1))*(dx(6)*dx(6)-dx(2)*dx(2)+dx(1)*dx(1)) + x*(dx(5)-dx(2)) &
      *(dx(5)-dx(6)-dx(2)+dx(1))
  END FUNCTION S1

  FUNCTION S2(x,dx)
  REAL(KIND=dbl), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  REAL(KIND=dbl) :: S2

  S2 = (dx(3)+dx(4))*(x+dx(5))*(x-dx(2))
  END FUNCTION S2

  FUNCTION S3(x,dx)
  REAL(KIND=dbl), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  REAL(KIND=dbl) :: S3

  S3 = dx(4)*(dx(4)+x)*(x-dx(6)+dx(1))+dx(3)*dx(3)*(dx(4)+x-dx(6)+dx(1)) &
      -dx(3)*(dx(4)+x)*(dx(4)+x-dx(6)+dx(1)) &
      +(dx(6)-dx(1))*(dx(6)*dx(6)-x*x+dx(1)*dx(1))
  END FUNCTION S3

  FUNCTION S4(x,dx)
  REAL(KIND=dbl), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  REAL(KIND=dbl) :: S4

  S4 = (dx(3)-x)*(dx(4)+x)*(dx(2)+dx(5))
  END FUNCTION S4


  FUNCTION coef_a(dx) RESULT(a)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: dx
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dxdx1
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: a
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: g1
  INTEGER :: i,j,n

  n = size(dx) + 1
  ALLOCATE(a(n,6))
  ALLOCATE(g1(n))
  ALLOCATE(dxdx1(n,6))
  a = 0.0_dbl
  g1 = 0.0_dbl
  dxdx1 = 0.0_dbl
  DO i=4,n-3
    dxdx1(i,1) = dx(i-3) + dx(i-2) + dx(i-1)
    dxdx1(i,2) = dx(i-2) + dx(i-1)
    dxdx1(i,3) = dx(i-1)
    dxdx1(i,4) = dx(i+1)
    dxdx1(i,5) = dx(i+1) + dx(i+2)
    dxdx1(i,6) = dx(i+1) + dx(i+2) + dx(i+3)
  END DO
  DO i=4,n-3
    g1(i) = (dxdx1(i,5)-dxdx1(i,6))*(dxdx1(i,6)-dxdx1(i,2))*(dxdx1(i,4)-dxdx1(i,6)+dxdx1(i,2)) &
          -dxdx1(i,1)*(dxdx1(i,1)+dxdx1(i,4))*(dxdx1(i,1)+dxdx1(i,5)-dxdx1(i,6)-dxdx1(i,2)) &
          +dxdx1(i,3)*(dxdx1(i,1)*dxdx1(i,1)-dxdx1(i,1)*(dxdx1(i,6)-dxdx1(i,5)+dxdx1(i,2))&
          +(dxdx1(i,6)+dxdx1(i,2))*(dxdx1(i,6)-dxdx1(i,5))+(dxdx1(i,4)*(dxdx1(i,5)-dxdx1(i,6)&
          -dxdx1(i,2)+dxdx1(i,1))))
  END DO
  DO i=4,n-3
    a(i,1) = (-dxdx1(i,3)+dxdx1(i,4)+dxdx1(i,5)-dxdx1(i,2))/(dxdx1(i,1)+dxdx1(i,6))/g1(i)
    a(i,2) = S3(dxdx1(i,5),dxdx1(i,:))/S4(dxdx1(i,5),dxdx1(i,:))/g1(i)
    a(i,3) = S1(dxdx1(i,4),dxdx1(i,:))/S2(dxdx1(i,3),dxdx1(i,:))/g1(i)
    a(i,4) = -S1(-dxdx1(i,3),dxdx1(i,:))/S2(-dxdx1(i,4),dxdx1(i,:))/g1(i)
    a(i,5) = -S3(-dxdx1(i,2),dxdx1(i,:))/S4(-dxdx1(i,2),dxdx1(i,:))/g1(i)
    a(i,6) = -a(i,1)
  END DO
  END FUNCTION coef_a 

  FUNCTION coef_b(dx) RESULT(b)
  REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: dx
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dxdx
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE ::b
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: g2
  INTEGER :: i,j,n

  n = size(dx) + 1
  ALLOCATE(b(n,3))
  ALLOCATE(g2(n))
  ALLOCATE(dxdx(n,3))
  b = 0.0_dbl
  g2 = 0.0_dbl
  dxdx = 0.0_dbl
  DO i=4,n-3
    dxdx(i,1) = dx(i-3) + dx(i-2) + dx(i-1)
    dxdx(i,2) = dx(i-2) + dx(i-1)
    dxdx(i,3) = dx(i-1)
    dxdx(i,4) = dx(i+1)
    dxdx(i,5) = dx(i+1) + dx(i+2)
    dxdx(i,6) = dx(i+1) + dx(i+2) + dx(i+3)
  END DO
  DO i=4,n-3
    g2(i) = dxdx(i,3)*dxdx(i,3)*dxdx(i,3) + dxdx(i,4)*dxdx(i,4)*dxdx(i,4)
  END DO
  DO i=4,n-3
    b(i,1) = -dxdx(i,3)*dxdx(i,3)/g2(i)
    b(i,2) = (dxdx(i,3)*dxdx(i,3) - dxdx(i,4)*dxdx(i,4))/g2(i)
    b(i,3) = dxdx(i,4)*dxdx(i,4) 
  END DO
  END FUNCTION coef_b

  FUNCTION Jacobian(n,h,a,b,theta,dt,A_bar,Ca) RESULT(J)
    INTEGER, INTENT(IN) ::n
    REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h
    REAL(KIND=dbl), DIMENSION(2*n+9,6), INTENT(IN) :: a
    REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: b
    REAL(KIND=dbl), INTENT(IN) :: theta, dt, A_bar, Ca
    REAL(KIND=dbl), DIMENSION(2*n+9,2*n+9) :: J
    INTEGER :: i

    J = 0.0_dbl
    J(1,1) = 1.0_dbl; J(1,2*n+1) = -1.0_dbl
    J(2,2) = 1.0_dbl; J(2,2*n+2) = -1.0_dbl
    J(3,3) = 1.0_dbl; J(3,2*n+3) = -1.0_dbl
    J(4,4) = 1.0_dbl; J(4,2*n+4) = -1.0_dbl
    J(2*n+6,6) = -1.0_dbl; J(2*n+6,2*n+6) = 1.0_dbl
    J(2*n+7,7) = -1.0_dbl; J(2*n+7,2*n+7) = 1.0_dbl
    J(2*n+8,8) = -1.0_dbl; J(2*n+8,2*n+8) = 1.0_dbl
    J(2*n+9,9) = -1.0_dbl; J(2*n+9,2*n+9) = 1.0_dbl
    DO i=5,2*n+5
      J(i,i-4) = (b(i,1)/3.0_dbl/Ca*a(i-1,1)*h(i-1)**3)*dt*theta
      J(i,i-3) = (b(i,1)/3.0_dbl/Ca*a(i-1,2)*h(i-1)**3 + b(i,2)/3.0_dbl/Ca*a(i,1)*h(i)**3)*dt*theta
      J(i,i-2) = (b(i,1)*(1.0_dbl/3.0_dbl/Ca*a(i-1,3)*h(i-1)**3+A_bar/h(i-1)*b(i-1,1)) &
                +b(i,2)/3.0_dbl/Ca*a(i,2)*h(i)**3 + b(i,3)/3.0_dbl/Ca*a(i+1,1)*h(i+1)**3)*dt*theta
      J(i,i-1) = (b(i,1)*(1.0_dbl/Ca*h(i-1)**2*(a(i-1,1)*h(i-4)+a(i-1,2)*h(i-3)+a(i-1,3)*h(i-2) &
                +a(i-1,4)*h(i)+a(i-1,5)*h(i+1)+a(i-1,6)*h(i+2))-A_bar/h(i-1)**2*(b(i-1,1)*h(i-2) &
                +b(i-1,3)*h(i))) &
                +b(i,2)*(a(i,3)/3.0_dbl/Ca*h(i-1)**3+A_bar/h(i)*b(i,1)) &
                +b(i,3)/3.0_dbl/Ca*a(i+1,2)*h(i+1)**3)*dt*theta
      J(i,i) = (b(i,1)*(a(i-1,4)/3.0_dbl/Ca*h(i-1)**3+A_bar/h(i-1)*b(i-1,3)) &
               +b(i,2)*(1.0_dbl/Ca*h(i)**2*(a(i,1)*h(i-3)+a(i,2)*h(i-2)+a(i,3)*h(i-1) &
               +a(i,4)*h(i+1)+a(i,5)*h(i+2)+a(i,6)*h(i+3))-A_bar/h(i)**2*(b(i,1)*h(i-1) &
               +b(i,3)*h(i+1))) &
               +b(i,3)*(a(i+1,3)/3.0_dbl/Ca*h(i+1)**3+A_bar/h(i)*b(i+1,1)))*dt*theta + 1.0_dbl
     J(i,i+1) = (b(i,1)/3.0_dbl/Ca*a(i-1,5)*h(i-1)**3 &
               +b(i,2)*(a(i,4)/3.0_dbl/Ca*h(i)**3+A_bar/h(i)*b(i,3)) &
               +b(i,3)*(1.0_dbl/Ca*h(i+1)**2*(a(i+1,1)*h(i-2)+a(i+1,2)*h(i-1)+a(i+1,3)*h(i) &
               +a(i+1,4)*h(i+2)+a(i+1,5)*h(i+3)+a(i+1,6)*h(i+4))-A_bar/h(i+1)**2*(b(i+1,1)*h(i) &
               +b(i+1,3)*h(i+2))))*dt*theta
     J(i,i+2) = (b(i,1)/3.0_dbl/Ca*a(i-1,6)*h(i-1)**3 &
               +b(i,2)/3.0_dbl/Ca*a(i,5)*h(i)**2) &
               +b(i,3)*(a(i+1,4)/3.0_dbl/Ca*h(i+1)**3+A_bar/h(i+1)*b(i+1,3))*dt*theta
     J(i,i+3) = (b(i,2)/3.0_dbl/Ca*a(i,6)*h(i)**3 + b(i,3)/3.0_dbl/Ca*a(i+1,5)*h(i+1)**3)*dt*theta
     J(i,i+4) = (b(i,3)/3.0_dbl/Ca*a(i+1,6)*h(i+1)**3)*dt*theta
    END DO
  END FUNCTION Jacobian 

  FUNCTION rhs_f(n,h,a,b,A_bar,Ca) RESULT(f)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h
    REAL(KIND=dbl), DIMENSION(2*n+9,6), INTENT(IN) :: a
    REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: b
    REAL(KIND=dbl), INTENT(IN) :: A_bar, Ca
    REAL(KIND=dbl), DIMENSION(2*n+9) :: f
    INTEGER :: i

    f = 0.0_dbl
    DO i=5,2*n+5
      f(i) = b(i,1)*(h(i-1)**3/3.0_dbl/Ca*(a(i-1,1)*h(i-4)+a(i-1,2)*h(i-3)+a(i-1,3)*h(i-2) &
            +a(i-1,4)*h(i)+a(i-1,5)*h(i+1)+a(i-1,6)*h(i+2))+A_bar/h(i-1)*(b(i-1,1)*h(i-2) &
            +b(i-1,3)*h(i))+A_bar*b(i-1,2)) &
            +b(i,2)*(h(i)**3/3.0_dbl/Ca*(a(i,1)*h(i-3)+a(i,2)*h(i-2)+a(i,3)*h(i-1) &
            +a(i,4)*h(i+1)+a(i,5)*h(i+2)+a(i+1,6)*h(i+3))+A_bar/h(i)*(b(i,1)*h(i-1) &
            +b(i,3)*h(i+1)) + A_bar*b(i,2)) &
            +b(i,3)*(h(i+1)**3/3.0_dbl/Ca*(a(i+1,1)*h(i-2)+a(i+1,2)*h(i-1)+a(i+1,3)*h(i) &
            +a(i+1,4)*h(i+2)+a(i+1,5)*h(i+3)+a(i+1,6)*h(i+4))+A_bar/h(i+1)*(b(i+1,1)*h(i) &
            +b(i+1,3)*h(i+2)) + A_bar*b(i+1,2))
    END DO
  END FUNCTION rhs_f

  FUNCTION evaluate_K(n,h,a,b,theta,dt,A_bar,Ca) RESULT(K)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h
    REAL(KIND=dbl), DIMENSION(2*n+9,6), INTENT(IN) :: a
    REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: b
    REAL(KIND=dbl), INTENT(IN) :: theta, dt, A_bar, Ca
    REAL(KIND=dbl), DIMENSION(2*n+9) :: f,K
    INTEGER :: i

    f = rhs_f(n,h,a,b,A_bar,Ca)
    K = 0.0_dbl
    DO i=5,2*n+5
      K(i) = h(i) - dt*(1.0_dbl-theta)*f(i)
    END DO
  END FUNCTION evaluate_K

  FUNCTION evaluate_f(n,h_1,h_2,a,b,theta,dt,A_bar,Ca) RESULT(ff)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h_1,h_2 ! h_1 is the j and h_2 is the j+1
    REAL(KIND=dbl), DIMENSION(2*n+9,6), INTENT(IN) :: a
    REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: b
    REAL(KIND=dbl), INTENT(IN) :: theta, dt, A_bar, Ca
    REAL(KIND=dbl), DIMENSION(2*n+9) :: ff,f,K
    INTEGER :: i
    
    f = rhs_f(n,h_2,a,b,A_bar,Ca)
    K = evaluate_K(n,h_1,a,b,theta,dt,A_bar,Ca)
    DO i=5,2*n+5
    ff(i) = h_2(i) + dt*theta*f(i) - K(i)
    END DO
    ff(1) = h_2(1) - h_2(2*n+1)
    ff(2) = h_2(2) - h_2(2*n+2)
    ff(3) = h_2(3) - h_2(2*n+3)
    ff(4) = h_2(4) - h_2(2*n+4)
    ff(2*n+6) = h_2(2*n+6) - h_2(6)
    ff(2*n+7) = h_2(2*n+7) - h_2(7)
    ff(2*n+8) = h_2(2*n+8) - h_2(8)
    ff(2*n+9) = h_2(2*n+9) - h_2(9)
  END FUNCTION evaluate_f

!!!!!!!!!!!!!!!!!!! Explicit Euler for initial guess !!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION explicit_euler(n,h_1,a,b,A_bar,Ca,dt) RESULT(h_2)
  INTEGER, INTENT(IN) :: n
  REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h_1
  REAL(KIND=dbl), DIMENSION(2*n+9,6), INTENT(IN) :: a
  REAL(KIND=dbl), DIMENSION(2*n+9,3), INTENT(IN) :: b
  REAL(KIND=dbl), INTENT(IN) :: A_bar, Ca, dt
  REAL(KIND=dbl), DIMENSION(2*n+9) :: f,h_2
  INTEGER :: i

  f = rhs_f(n,h_1,a,b,A_bar,Ca) 
  DO i=5,2*n+5
    h_2(i) = h_1(i) - f(i)*dt
  END DO
  h_2(1) = h_2(2*n+1); h_2(2) = h_2(2*n+2); h_2(3) = h_2(2*n+3); h_2(4) = h_2(2*n+4)
  h_2(2*n+6) = h_2(6); h_2(2*n+7) = h_2(7); h_2(2*n+8) = h_2(8); h_2(2*n+9) = h_2(9)
  END FUNCTION explicit_euler

  FUNCTION non_uniform_2nd_order_newton(n,h_1,dx,theta,dt,A_bar,Ca,tol) RESULT(h_2)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2*n+9), INTENT(IN) :: h_1
    REAL(KIND=dbl), DIMENSION(2*n+8), INTENT(IN) :: dx
    REAL(KIND=dbl), DIMENSION(2*n+9,6) :: a
    REAL(KIND=dbl), DIMENSION(2*n+9,3) :: b
    REAL(KIND=dbl), DIMENSION(2*n+9) :: h_2, f
    REAL(KIND=dbl), INTENT(IN) :: A_bar, Ca, dt, tol, theta 
    REAL(KIND=dbl) :: e
    INTEGER, DIMENSION(2*n+9) :: ipiv
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: work(:)
    INTEGER :: info, lwork
    REAL(KIND=dbl), DIMENSION(2*n+9,2*n+9) :: Jac, J_inv, LU

    a = coef_a(dx)
    b = coef_b(dx)  
    h_2 = explicit_euler(n,h_1,a,b,A_bar,Ca,dt)
    e = 1.0_dbl
    DO WHILE (e>tol)
      Jac = Jacobian(n,h_2,a,b,theta,dt,A_bar,Ca) 

      info = 0
      lwork = -1
      ALLOCATE(work(lwork))
      work = 0
      ipiv = 0

      LU = Jac
      CALL dgetrf(2*n+9,2*n+9,LU,2*n+9,ipiv,info)
      J_inv = LU
      CALL dgetri(2*n+9,J_inv,n+4,ipiv,work,lwork,info)
      DEALLOCATE(work)
      f = evaluate_f(n,h_1,h_2,a,b,theta,dt,A_bar,Ca)
      h_2 = h_2 - matmul(J_inv,f)
      f = evaluate_f(n,h_1,h_2,a,b,theta,dt,A_bar,Ca)
      e = norm2(f)
    END DO
  END FUNCTION non_uniform_2nd_order_newton
END MODULE non_uniform_newton
















