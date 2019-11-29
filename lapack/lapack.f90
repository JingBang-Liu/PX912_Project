PROGRAM matsolve

  USE lapack_interfaces
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  REAL(dp), DIMENSION(2, 2) :: A
  REAL(dp), DIMENSION(2) :: B
  INTEGER, DIMENSION(2) :: ipiv = [1,2]
  INTEGER :: info

  A(1,:) = [2, 1]
  A(2,:) = [0, 1]
  B = [1,3]

  CALL dgesv(2, 1, A, 2, ipiv, B, 2, info)

  PRINT *, B
  PRINT *, info

END PROGRAM matsolve
