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

END MODULE finite_difference_operators

PROGRAM MAIN

USE finite_difference_operators

  IMPLICIT NONE

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0
  INTEGER, PARAMETER :: n = 401
  INTEGER(KIND=int64), PARAMETER :: m = 300000 
  REAL(KIND=dbl) :: dt = 1e-7_dbl ! time step size
  REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi/w ! boundary
  REAL(KIND=dbl) :: dx = (upper_bnd-lower_bnd)/real(n-1) ! grid step size
  REAL(KIND=dbl), DIMENSION(n) :: x ! grid array
  REAL(KIND=dbl), DIMENSION(m,n+4) :: h ! matrix to store solutions
  INTEGER :: k = 1 ! counter for what time we are on
  INTEGER :: i,j
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl
  

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
    do i=3,n+2
      h(k,i) = h(k-1,i) - &
      dt/3.0_dbl/Ca*(3.0_dbl*h(k-1,i)**2*dhdx(dx,h(k-1,i-1),h(k-1,i+1))*d3hdx(dx,h(k-1,i-2),h(k-1,i-1),h(k-1,i+1),h(k-1,i+2)) + &
      h(k-1,i)**3*d4hdx(dx,h(k-1,i-2),h(k-1,i-1),h(k-1,i),h(k-1,i+1),h(k-1,i+2)))  
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



































