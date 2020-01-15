!!!!!!!!!!!!!!!!!!!!!!!!!! Preparation for non uniform newton !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FUNCTION S1(x,dx)
  ! REAL(KIND=dbl), INTENT(IN) :: x
  ! REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  ! REAL(KIND=dbl) :: S1

  ! S1 = dx(5)*(dx(5)-dx(2))*(dx(6)+dx(2)-dx(1))+x*x*(dx(6)-dx(5)+dx(2)-dx(1)) &
  !     -(dx(6)-dx(1))*(dx(6)*dx(6)-dx(2)*dx(2)+dx(1)*dx(1)) - x*(dx(5)-dx(2)) &
  !     *(dx(5)-dx(6)-dx(2)+dx(1))
  ! END FUNCTION S1

  ! FUNCTION S2(x,dx)
  ! REAL(KIND=dbl), INTENT(IN) :: x
  ! REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  ! REAL(KIND=dbl) :: S2

  ! S2 = (dx(3)+dx(4))*(x+dx(5))*(x-dx(2))
  ! END FUNCTION S2

  ! FUNCTION S3(x,dx)
  ! REAL(KIND=dbl), INTENT(IN) :: x
  ! REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  ! REAL(KIND=dbl) :: S3

  ! S3 = dx(4)*(dx(4)+x)*(x-dx(6)+dx(1))+dx(3)*dx(3)*(dx(4)+x-dx(6)+dx(1)) &
  !     -dx(3)*(dx(4)+x)*(dx(4)+x-dx(6)+dx(1)) &
  !     +(dx(6)-dx(1))*(dx(6)*dx(6)-x*x+dx(1)*dx(1))
  ! END FUNCTION S3

  ! FUNCTION S4(x,dx)
  ! REAL(KIND=dbl), INTENT(IN) :: x
  ! REAL(KIND=dbl), DIMENSION(6), INTENT(IN) :: dx
  ! REAL(KIND=dbl) :: S4

  ! S4 = (dx(3)-x)*(dx(4)+x)*(dx(2)+dx(5))
  ! END FUNCTION S4

  ! FUNCTION S5(x,dx)
  ! REAL(KIND=dbl), INTENT(IN) :: x
  ! REAL(KIND=dbl), DIMENSION(3), INTENT(IN) :: dx
  ! REAL(KIND=dbl) :: S5

  ! S5 = 1.0_dbl + 2.0_dbl*x*x*dx(3) - 2.0_dbl*dx(3)*dx(3)*dx(3)
  ! END FUNCTION S5

  ! FUNCTION S6(x,dx)
  ! REAL(KIND=dbl), INTENT(IN) :: x
  ! REAL(KIND=dbl), DIMENSION(3), INTENT(IN) :: dx
  ! REAL(KIND=dbl) :: S6

  ! S6 = 2.0_dbl*x*(dx(1)-dx(2))*(dx(1)+dx(2))
  ! END FUNCTION S6

  ! FUNCTION coef_a(dx) RESULT(a)
  ! REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: dx
  ! REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dxdx1
  ! REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: a
  ! REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: g1
  ! INTEGER :: i,j,n

  ! n = size(dx) + 1
  ! ALLOCATE(a(n,6))
  ! ALLOCATE(g1(n))
  ! ALLOCATE(dxdx1(n,6))
  ! a = 0.0_dbl
  ! g1 = 0.0_dbl
  ! dxdx1 = 0.0_dbl
  ! DO i=4,n-3
  !   dxdx1(i,1) = dx(i-3) + dx(i-2) + dx(i-1)
  !   dxdx1(i,2) = dx(i-2) + dx(i-1)
  !   dxdx1(i,3) = dx(i-1)
  !   dxdx1(i,4) = dx(i)
  !   dxdx1(i,5) = dx(i) + dx(i+1)
  !   dxdx1(i,6) = dx(i) + dx(i+1) + dx(i+2)
  ! END DO
  ! DO i=4,n-3
  !   g1(i) = (dxdx1(i,5)-dxdx1(i,6))*(dxdx1(i,6)+dxdx1(i,2))*(dxdx1(i,4)-dxdx1(i,6)+dxdx1(i,1)-dxdx1(i,3)) &
  !          +dxdx1(i,1)*(dxdx1(i,4)-dxdx1(i,3))*(-dxdx1(i,5)+dxdx1(i,6)+dxdx1(i,2))-dxdx1(i,1)*dxdx1(i,1)*dxdx1(i,1) &
  !          +(-dxdx1(i,4)-dxdx1(i,5)+dxdx1(i,6)+dxdx1(i,2)+dxdx1(i,3))*dxdx1(i,1)*dxdx1(i,1) &
  !          +dxdx1(i,3)*dxdx1(i,4)*(dxdx1(i,5)-dxdx1(i,6)-dxdx1(i,2)+dxdx1(i,1)) 
  ! END DO
  ! DO i=4,n-3
  !   a(i,1) = (-dxdx1(i,3)+dxdx1(i,4)+dxdx1(i,5)-dxdx1(i,2))/(dxdx1(i,1)+dxdx1(i,6))/g1(i)
  !   a(i,2) = S3(dxdx1(i,5),dxdx1(i,:))/S4(dxdx1(i,2),dxdx1(i,:))/g1(i)
  !   a(i,3) = S1(dxdx1(i,4),dxdx1(i,:))/S2(dxdx1(i,3),dxdx1(i,:))/g1(i)
  !   a(i,4) = -S1(-dxdx1(i,3),dxdx1(i,:))/S2(-dxdx1(i,4),dxdx1(i,:))/g1(i)
  !   a(i,5) = -S3(-dxdx1(i,2),dxdx1(i,:))/S4(-dxdx1(i,5),dxdx1(i,:))/g1(i)
  !   a(i,6) = -a(i,1)
  !   DO j=1,6
  !     IF (a(i,j)/=a(i,j)) THEN
  !       PRINT*, "THE dx is", dxdx1(i,:)
  !     END IF
  !   END DO
  ! END DO
  ! END FUNCTION coef_a 

  ! FUNCTION coef_b(dx) RESULT(b)
  ! REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: dx
  ! REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dxdx
  ! REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: b
  ! REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: g2
  ! INTEGER :: i,j,n

  ! n = size(dx) + 1
  ! ALLOCATE(b(n,3))
  ! ALLOCATE(g2(n))
  ! ALLOCATE(dxdx(n,2))
  ! b = 0.0_dbl
  ! g2 = 0.0_dbl
  ! dxdx = 0.0_dbl
  ! DO i=4,n-3
  !   dxdx(i,1) = dx(i-1)
  !   dxdx(i,2) = dx(i)
  ! END DO
  ! DO i=4,n-3
  !   g2(i) = dxdx(i,1)*dxdx(i,2)*(dxdx(i,1)+dxdx(i,2))
  ! END DO
  ! DO i=4,n-3
  !   b(i,1) = -dxdx(i,2)*dxdx(i,2)/g2(i)
  !   b(i,2) = (dxdx(i,2)*dxdx(i,2)-dxdx(i,1)*dxdx(i,1))/g2(i)
  !   b(i,3) = dxdx(i,1)*dxdx(i,1)/g2(i) 
  ! END DO
  ! END FUNCTION coef_b

  ! FUNCTION coef_c(dx) RESULT(c)
  ! REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: dx
  ! REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dxdx
  ! REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: c
  ! INTEGER :: i, k, n

  ! n = size(dx) + 1
  ! ALLOCATE(c(n,3))
  ! ALLOCATE(dxdx(6,3))
  ! c = 0.0_dbl

  ! dxdx = 0.0_dbl
  ! dxdx(1,1) = dx(4); dxdx(1,2) = dx(4) + dx(5); dxdx(1,3) = dx(4) + dx(5) + dx(6)
  ! dxdx(2,1) = dx(6); dxdx(2,2) = dx(6)+dx(7); dxdx(2,3)=dx(6)+dx(7)+dx(8)
  ! dxdx(3,1) = dx(n+4); dxdx(3,2) = dx(n+4)+dx(n+5); dxdx(3,3) = dx(n+4)+dx(n+5)+dx(n+6)
  ! dxdx(4,1) = dx(n+6); dxdx(4,2) = dx(n+6)+dx(n+7); dxdx(4,3) = dx(n+6)+dx(n+7)+dx(n+8)
  ! dxdx(5,1) = dx(2*n+4); dxdx(4,2) = dx(2*n+4)+dx(2*n+5); dxdx(5,3) = dx(2*n+4)+dx(2*n+5)+dx(2*n+6)
  ! dxdx(6,1) = dx(2*n+6); dxdx(6,2) = dx(2*n+6)+dx(2*n+7); dxdx(6,3) = dx(2*n+6)+dx(2*n+7)+dx(2*n+8)
  
  ! DO i=1,3
  !   k = (i-1)*n+5
  !   c(k-1,1) = S5(dxdx(2*i-1,2),dxdx(2*i-1,:))/S6(dxdx(2*i-1,1),dxdx(2*i-1,:))
  !   c(k-1,2) = S5(dxdx(2*i-1,1),dxdx(2*i-1,:))/S6(-dxdx(2*i-1,2),dxdx(2*i-1,:))
  !   c(k-1,3) = 1.0_dbl
  !   c(k+1,1) = S5(dxdx(2*i,2),dxdx(2*i,:))/S6(dxdx(2*i,1),dxdx(2*i,:))
  !   c(k+1,2) = S5(dxdx(2*i,1),dxdx(2*i,:))/S6(-dxdx(2*i,2),dxdx(2*i,:))
  !   c(k+1,3) = 1.0_dbl
  ! END DO
  ! END FUNCTION coef_c
