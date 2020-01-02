MODULE generate_grid

  USE kinds

  IMPLICIT NONE

  CONTAINS
  !Generate a uniform grid
  SUBROUTINE grid_uni(x,dx,lower_bnd,upper_bnd)
    REAL(KIND=dbl), DIMENSION(:), INTENT(INOUT) :: x
    REAL(KIND=dbl), INTENT(INOUT) :: dx
    REAL(KIND=dbl), INTENT(IN) :: lower_bnd, upper_bnd
    INTEGER :: i, n
    
    n = size(x)
    dx = (upper_bnd-lower_bnd)/(real(n-9))

    DO i=5,n-4
      x(i) = lower_bnd + dx * (real(i-4)-1.0_dbl)
    END DO
    DO i=1,4
      x(i) = lower_bnd - dx * (real(5-i))
      x(n-4+i) = upper_bnd + dx * real(i)
    END DO
  END SUBROUTINE grid_uni

  SUBROUTINE grid_square(n,x,lower_bnd,upper_bnd)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(:), INTENT(INOUT) :: x
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
  END SUBROUTINE grid_square

  SUBROUTINE grid_sin(n,x,lower_bnd, upper_bnd)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dbl), DIMENSION(:), INTENT(INOUT) :: x
    REAL(KIND=dbl), INTENT(IN) :: lower_bnd, upper_bnd
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: y
    REAL(KIND=dbl), PARAMETER :: pi=3.1415926535897932_dbl

    REAL(KIND=dbl) :: dy
    INTEGER :: i

    ALLOCATE(y(n+1))
    dy = 1.0_dbl/real(n)
    DO i=1,n+1
      y(i) = dy*(real(i-1))
    END DO
    DO i=1,n+1
      x(i+n+4) = (upper_bnd-lower_bnd)/2.0_dbl + (upper_bnd-lower_bnd)*asin(y(i))/pi
    END DO 
    DO i=5,n+4
      x(i) = upper_bnd + lower_bnd - x(2*n+10-i)
    END DO
    DO i=1,4
      x(i) = 2.0_dbl*lower_bnd - x(10-i)
    END DO
    DO i=2*n+6,2*n+9
      x(i) = 2.0_dbl*upper_bnd - x(4*n+10-i)
    END DO
  END SUBROUTINE grid_sin
END MODULE generate_grid

