! The purpose of this program is to use finite difference method to solve
! the thin film equations

MODULE finite_difference_operators

  IMPLICIT NONE

  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(15, 307)
  INTEGER, PARAMETER :: int64 = SELECTED_INT_KIND(15)

  CONTAINS
  
  ! central difference to calculate first derivatives
  FUNCTION dhdx(dx,h2,h4)
    REAL(KIND=dbl), INTENT(IN) :: dx, h2, h4 ! input
    REAL(KIND=dbl) :: dhdx ! output
    
    dhdx = (h4-h2)/dx
  END FUNCTION dhdx

  ! central difference to calculate second derivatives
  FUNCTION d2hdx(dx,h2,h3,h4)
    REAL(KIND=dbl), INTENT(IN) :: dx, h2, h3, h4 ! input
    REAL(KIND=dbl) :: d2hdx ! output

    d2hdx = (h4-2.0_dbl*h3+h2)/(dx*dx)
  END FUNCTION d2hdx

  ! central difference to calculate third derivatives
  FUNCTION d3hdx(dx,h1,h2,h4,h5)
    REAL(KIND=dbl), INTENT(IN) :: dx, h1, h2, h4, h5 ! input
    REAL(KIND=dbl) :: d3hdx ! output

    d3hdx = (h5 -2*h4 +2*h2-h1)/2/(dx*dx*dx)
  END FUNCTION d3hdx

  ! central difference to calculate fourth derivatives
  FUNCTION d4hdx(dx,h1,h2,h3,h4,h5)
    REAL(KIND=dbl), INTENT(IN) :: dx, h1, h2, h3, h4, h5 ! input
    REAL(KIND=dbl) :: d4hdx ! output

    d4hdx = (h5-4*h4+6*h3-4*h2+h1)/(dx*dx*dx*dx)
  END FUNCTION d4hdx

  ! calculate f
  FUNCTION fun(dx,h1,h2,h3,h4,h5,Ca,A_f)
    REAL(KIND=dbl), INTENT(IN) :: dx, h1, h2, h3, h4, h5, Ca, A_f ! input
    REAL(KIND=dbl) :: fun ! output

    fun = -1.0_dbl/Ca*(h3*h3*dhdx(dx,h2,h4)*d3hdx(dx,h1,h2,h4,h5) &
            +h3*h3*h3*d4hdx(dx,h1,h2,h3,h4,h5)/3.0_dbl)&
            - A_f/h3*(d2hdx(dx,h2,h3,h4) - dhdx(dx,h2,h4)*dhdx(dx,h2,h4)/h3)
  END FUNCTION fun

END MODULE finite_difference_operators

PROGRAM MAIN

USE finite_difference_operators

  IMPLICIT NONE

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0
  INTEGER, PARAMETER :: n = 51
  INTEGER(KIND=int64), PARAMETER :: m = 300000 
  REAL(KIND=dbl) :: dt = 1e-4_dbl ! time step size
  REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi/w ! boundary
  REAL(KIND=dbl) :: dx = (upper_bnd-lower_bnd)/real(n-1) ! grid step size
  REAL(KIND=dbl), DIMENSION(n) :: x ! grid array
  REAL(KIND=dbl), DIMENSION(m,n+4) :: h ! matrix to store solutions
  INTEGER :: k = 1 ! counter for what time we are on
  INTEGER :: i,j
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl
  REAL(KIND=dbl), PARAMETER :: A_f = 1.0_dbl
  REAL(KIND=dbl), DIMENSION(n+4) :: k1, k2, k3, k4 ! for Runge-Kutta 4th order method
  

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
  h(1,1) = h(1,5)
  h(1,2) = h(1,4)
  h(1,n+4) = h(1,n)
  h(1,n+3) = h(1,n+1)

  do while ((k<m).and.(test==0))
    k = k + 1
    do i=3,n+2 ! calculate k1
      k1(i) = dt*fun(dx,h(k-1,i-2),h(k-1,i-1),h(k-1,i),h(k-1,i+1),h(k-1,i+2),Ca,A_f)
    end do
    k1(1) = k1(5)
    k1(2) = k1(4)
    k1(n+4) = k1(n)
    k1(n+3) = k1(n+1)
    do i=3,n+2 ! calculate k2
      k2(i) = dt*fun(dx,h(k-1,i-2)+k1(i-2)/2.0_dbl,h(k-1,i-1)+k1(i-1)/2.0_dbl, &
                 h(k-1,i)+k1(i)/2.0_dbl,h(k-1,i+1)+k1(i+1)/2.0_dbl,h(k-1,i+2)+k1(i+2)/2.0_dbl,Ca,A_f)
    end do
    k2(1) = k2(5)
    k2(2) = k2(4)
    k2(n+4) = k2(n)
    k2(n+3) = k2(n+1)
    do i=3,n+2 ! calcualte k3
      k3(i) = dt*fun(dx,h(k-1,i-2)+k2(i-2)/2.0_dbl,h(k-1,i-1)+k2(i-1)/2.0_dbl, &
                 h(k-1,i)+k2(i)/2.0_dbl,h(k-1,i+1)+k2(i+1)/2.0_dbl,h(k-1,i+2)+k2(i+2)/2.0_dbl,Ca,A_f)
    end do
    k3(1) = k3(5)
    k3(2) = k3(4)
    k3(n+4) = k3(n)
    k3(n+3) = k3(n+1)
    do i=3,n+2 ! calculate k4
      k4(i) = dt*fun(dx,h(k-1,i-2)+k3(i-2),h(k-1,i-1)+k3(i-1), &
                 h(k-1,i)+k3(i),h(k-1,i+1)+k3(i+1),h(k-1,i+2)+k3(i+2),Ca,A_f)
    end do
    k4(1) = k4(5)
    k4(2) = k4(4)
    k4(n+4) = k4(n)
    k4(n+3) = k4(n+1)
    do i=3,n+2 ! calculate h(k,:)
      h(k,i) = h(k-1,i) + (k1(i)+2.0_dbl*k2(i)+2.0_dbl*k3(i)+k4(i))/6.0_dbl
    end do
    h(k,1) = h(k,5)
    h(k,2) = h(k,4)
    h(k,n+4) = h(k,n)
    h(k,n+3) = h(k,n+1)
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



































