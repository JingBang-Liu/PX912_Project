MODULE Fornberg_coeff

  USE kinds

  IMPLICIT NONE

  CONTAINS
  !!!!!! function FB_coef takes in: 'm' the maximum order of derivative, 'n_in' the number of grid
  !!!!!! points taken into consideration, 'x0' the point of the derivatives, 'x' the grid points
  FUNCTION FB_coef(m,n_in,x0,x) RESULT(c)
  INTEGER, INTENT(IN) :: m
  INTEGER, INTENT(IN) :: n_in
  REAL(KIND=dbl), INTENT(IN) :: x0
  REAL(KIND=dbl), DIMENSION(n_in), INTENT(IN) :: x
  REAL(KIND=dbl), DIMENSION(m+1,n_in,n_in) :: c
  REAL(KIND=dbl) :: c1, c2, c3
  INTEGER :: i, j, k, n

  n = n_in - 1
  c = 0.0_dbl
  c(1,1,1) = 1.0_dbl
  c1 = 1.0_dbl
  c3 = 0.0_dbl
  DO i=1,n !! this is for n
    c2 = 1.0_dbl
    DO j=0,i-1 !! this is for nu
      c3 = x(i+1) - x(j+1)
      c2 = c2*c3
      DO k=1, min(i,m) !! this is for m
        c(k+1,i+1,j+1) = ((x(i+1) - x0)*c(k+1,i,j+1) - k*c(k,i,j+1))/c3
      END DO
    END DO
    DO k=1, min(i,m)
      c(k+1,i+1,i+1) = c1*(k*c(k,i,i) - (x(i) - x0)*c(k+1,i,i))
    END DO
    c1 = c2
  END DO
  END FUNCTION FB_coef

END MODULE Fornberg_coeff
