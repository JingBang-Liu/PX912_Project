MODULE generate_grid

  USE kinds

  IMPLICIT NONE

  CONTAINS
  !Generate a uniform grid
  SUBROUTINE gird_uni(x,n,lower_bnd,upper_bnd)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2*n+9), ALLOCATABLE, INTENT(INOUT) :: x
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: dx
    REAL(KIND=dbl), INTENT(IN) :: lower_bnd, upper_bnd
    INTEGER :: i

    ALLOCATE(dx(1))
    dx = (upper_bnd-lower_bnd)/(real(2*n))

    ALLOCATE(x(2*n+1))
    DO i=5,2*n+5
      x(i) = lower_bnd + dx * (real(i-4)-1.0_dbl)
    END DO
    DO i=1,4
      x(i) = lower_bnd - dx * (real(5-i))
      x(2*n+5+i) = upper_bnd + dx * real(i)
    END DO
  END SUBROUTINE grid_uni

  SUBROUTINE grid_square(x,n,lower_bnd,upper_bnd)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(2*n+9), ALLOCATABLE, INTENT(INOUT) :: x
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: dx
    REAL(KIND=dbl), INTENT(IN) : lower_bnd, upper_bnd
    REAL(KIND=dbl), DIMENSION(n+1) :: y ! uniform gird that project to a
                                        ! non-uniform gird x
    REAL(KIND=dbl) :: dy
    INTEGER :: i

    dy = 1.0_dbl/real(n)
    DO i=1,n+1
      y(i) = 1.0_dbl - dy*(real(i)-1.0_dbl)
    END DO    

    ALLOCATE(dx(2*n+8))
    ALLOCATE(x(2*n+9))

    DO i=1,n+1
      x(i+4) = lower_bnd + (upper_bnd-lower_bnd)/2.0_dbl*sqrt(1.0_dbl-y(i))
    END DO
    DO i=n+2,2*n+1
      x(i+4) = upper_bnd + lower_bnd - x(2*n+2-i)
    END DO
    DO i=1,4
      x(i) = x(2*n+i)
      x(2*n+5+i) = x(5+i)
    END DO
    DO i=1,2*n+8
      dx(i) = x(i+1) - x(i)
    END DO
  END SUBROUTINE grid_square
END MODULE generate_grid

