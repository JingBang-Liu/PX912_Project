PROGRAM MAIN

  USE kinds
  USE Fornberg_coeff

  IMPLICIT NONE

  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: x
  REAL(KIND=dbl) :: x0
  INTEGER :: m, n_in
  REAL(KIND=dbl), DIMENSION(:,:,:), ALLOCATABLE :: c

  m = 1
  n_in = 3
  ALLOCATE(x(n_in))
  ALLOCATE(c(m+1,n_in,n_in))
  x = 0.0_dbl
  x0 = 0.0_dbl
  x(1) = -1.0_dbl
  x(2) = 0.0_dbl
  x(3) = 1.0_dbl

  c = FB_coef(m,n_in,x0,x)
  PRINT*, '1st derivative with 3 symmetrical points are', c(2,3,:)



END PROGRAM
