MODULE generate_grid

  USE kinds

  IMPLICIT NONE

  CONTAINS
  !Generate a uniform grid
  SUBROUTINE grid_uni(x,dx,n,lower_bnd,upper_bnd)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: x
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: dx
    REAL(KIND=dbl), INTENT(IN) :: lower_bnd, upper_bnd
    INTEGER :: i

    ALLOCATE(x(2*n+9))
    ALLOCATE(dx(1))
    dx = (upper_bnd-lower_bnd)/(real(2*n))

    DO i=5,2*n+5
      x(i) = lower_bnd + dx(1) * (real(i-4)-1.0_dbl)
    END DO
    DO i=1,4
      x(i) = lower_bnd - dx(1) * (real(5-i))
      x(2*n+5+i) = upper_bnd + dx(1) * real(i)
    END DO
  END SUBROUTINE grid_uni

  SUBROUTINE grid_square(x,dx,n,lower_bnd,upper_bnd)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: x
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: dx
    REAL(KIND=dbl), INTENT(IN) :: lower_bnd, upper_bnd
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: y ! uniform gird that project to a
                                        ! non-uniform gird x
    REAL(KIND=dbl) :: dy
    INTEGER :: i

    ALLOCATE(y(n+1))
    dy = 1.0_dbl/real(n)
    DO i=1,n+1
      y(i) = 1.0_dbl - dy*(real(i)-1.0_dbl)
    END DO    

    ALLOCATE(dx(2*n+8))
    ALLOCATE(x(2*n+9))

    DO i=1,n+1
      x(i+4) = lower_bnd + (upper_bnd-lower_bnd)/2.0_dbl*sqrt(1.0_dbl-y(i))
    END DO
    DO i=n+6,2*n+5
      x(i) = upper_bnd + lower_bnd - x(2*n+10-i)
    END DO
    DO i=1,4
      x(i) = 2.0_dbl*lower_bnd - x(10-i)
    END DO
    DO i=2*n+6,2*n+9
      x(i) = 2.0_dbl*upper_bnd - x(4*n+10-i)
    END DO
    DO i=1,2*n+8
      dx(i) = x(i+1) - x(i)
    END DO
  END SUBROUTINE grid_square
END MODULE generate_grid

