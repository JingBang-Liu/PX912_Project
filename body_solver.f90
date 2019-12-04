PROGRAM MAIN

  USE kinds        
  USE lapack_precision
  USE lapack_interfaces
  USE generate_grid
  USE non_uniform_newton
  USE command_line

  IMPLICIT NONE

  LOGICAL :: success, exists

  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0 
  INTEGER :: n
  INTEGER :: m
  REAL(KIND=dbl)

END PROGRAM

