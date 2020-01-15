PROGRAM MAIN

  USE kinds        
  USE lapack_precision
  USE lapack_interfaces
  USE generate_grid
  USE non_uniform_newton
  USE uniform_newton
  USE command_line
  USE netcdf_write

  IMPLICIT NONE

  !LOGICAL :: success, exists
  
  CHARACTER(LEN=20) :: GRID
  CHARACTER(LEN=20) :: SHAPE
  REAL(KIND=dbl), PARAMETER :: pi = 3.14159265359_dbl
  REAL(KIND=dbl), PARAMETER :: w = sqrt(0.5_dbl)
  INTEGER :: test = 0 
  INTEGER :: n ! input from command line
  INTEGER(KIND=INT64) :: m = 1e10
  REAL(KIND=dbl) :: dt ! input from command line
  !REAL(KIND=dbl), PARAMETER :: lower_bnd = 0.0_dbl, upper_bnd = 2.0_dbl*pi
  !REAL(KIND=dbl), PARAMETER :: lower_bnd = -pi/W, upper_bnd = pi/W
  REAL(KIND=dbl) :: lower_bnd, upper_bnd
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: x ! grid depend on n
  REAL(KIND=dbl) :: dx ! dx depend on x
  REAL(KIND=dbl) :: amp ! amplitude for pertubation
  REAL(KIND=dbl), PARAMETER :: Ca = 1.0_dbl, A_bar = 1.0_dbl
  REAL(KIND=dbl) :: theta ! input from command line
  REAL(KIND=dbl), PARAMETER :: tol_newton = 1e-2_dbl
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: h
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: h0
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: a, b
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: h_min, h_min_temp
  REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: time, time_temp
  REAL(KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: h_his, h_his_temp
  INTEGER :: time_gap, his_count
  TYPE(run_data) :: r_d
  INTEGER :: i, j, k

  !!! variables for Netcdf
  CHARACTER(LEN=25) :: filename="project3.nc"
  INTEGER :: ierr
  

  !CALL parse_args

  ! Try to grab n, the grid steps
  !success = get_arg("n",n,exists=exists)

  ! Try to grab dt, the time step size
  !success = get_arg("dt",dt,exists=exists)

  ! Try to grab theta, the implicit factor
  !success = get_arg("theta",theta,exists=exists)
  
  SHAPE = "1997"
  amp = 0.2
  ! options: uniform, non_uniform_square, non_uniform_sin
  GRID = "uniform"
  dt = 1e-6
  n = 100
  theta = 0.5_dbl
  time_gap = 100
  PRINT*, theta
  PRINT*, dt

  IF (SHAPE=="1999") THEN
    lower_bnd = 0.0_dbl
    upper_bnd = 2.0_dbl*pi
  ELSE IF (SHAPE=="1997") THEN
    lower_bnd = -pi/W
    upper_bnd = pi/W 
  END IF



  r_d%run_data_n = n
  r_d%run_data_dt = dt
  r_d%run_data_theta = theta
  r_d%run_data_grid = GRID

  !!!!!!!!!!!!!!!!!!!! for uniform grid !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (GRID == "uniform") THEN
    !dt = 1e-9
    ALLOCATE(x(2*n+9))
    !!!! create non-uniform grid
    CALL grid_uni(x,dx,lower_bnd,upper_bnd)
    !CALL grid_uni(x,dx,n,lower_bnd,upper_bnd)

    ALLOCATE(h(2,2*n+9))
    ALLOCATE(h0(2*n+9))

    h = 0.0_dbl
    ! initial condition
    DO i=1,2*n+9
      IF (SHAPE=="1999") THEN
        h(1,i) = 1.0_dbl + amp * cos(x(i))
      ELSE IF (SHAPE=="1997") THEN
        h(1,i) = 1.0_dbl - amp * cos(W*x(i))
      END IF
    END DO

    h0 = h(1,:)
    !PRINT*, "h0 is ", h(1,:)
    IF (h(1,2)/=h(1,2)) THEN
      PRINT*, "STUPID h0"
    END IF

    ! initial condition
    !h(1,1) = h(1,2*n+1); h(1,2) = h(1,2*n+2); h(1,3) = h(1,2*n+3); h(1,4) = h(1,2*n+4)
    !h(1,2*n+6) = h(1,6); h(1,2*n+7) = h(1,7); h(1,2*n+8) = h(1,8); h(1,2*n+9) = h(1,9)
    k = 1
    his_count = 1
    ALLOCATE(h_min(1))
    ALLOCATE(time(1))
    ALLOCATE(h_his(1,2*n+9))
    h_min(1) = minval(h(1,:))
    time(1) = 0.0_dbl
    h_his(1,:) = h(1,:)

    DO WHILE ((k<m).and.(test==0))
      k = k + 1
      h(2,:) = uniform_2nd_order_newton(h(1,:),dx,Ca,A_bar,dt,theta,tol_newton)
      DO i=5,2*n+5
        IF ((h(2,i)<0)) THEN
          test = 1
          PRINT*, "0", k
        ELSE IF (h(2,i)/=h(2,i)) THEN
          test = 1
        END IF
      END DO
      IF (mod(k,10)==0) THEN
        PRINT*, k
        PRINT*, minval(h(2,:))
      END IF
      h(1,:) = h(2,:)
      IF (mod(k,time_gap)==0) THEN
        his_count = his_count + 1
        ALLOCATE(h_his_temp(his_count-1,2*n+9))
        h_his_temp = h_his
        DEALLOCATE(h_his)
        ALLOCATE(h_his(his_count,2*n+9))
        h_his(1:his_count-1,:) = h_his_temp
        DEALLOCATE(h_his_temp)
        h_his(his_count,:) = h(1,:)
      END IF
      ALLOCATE(h_min_temp(k-1))
      ALLOCATE(time_temp(k-1))
      h_min_temp = h_min
      time_temp = time
      DEALLOCATE(h_min)
      DEALLOCATE(time)
      ALLOCATE(h_min(k))
      ALLOCATE(time(k))
      h_min(1:k-1) = h_min_temp
      time(1:k-1) = time_temp
      DEALLOCATE(h_min_temp)
      DEALLOCATE(time_temp)
      h_min(k) = minval(h(1,:))
      time(k) = dt*(k-1)
    END DO
  END IF

!!!!!!!!!!!!!!!!!!!!!!!!! for non_uniform_square gird !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (GRID == "non_uniform_square") THEN
    !dt = 2e-9
    ALLOCATE(x(2*n+9))
    CALL grid_square(n,x,lower_bnd,upper_bnd)

    ALLOCATE(h(2,2*n+9))
    ALLOCATE(h0(2*n+9))

    h = 0.0_dbl
    ! initial condition
    DO i=1,2*n+9
      IF (SHAPE=="1999") THEN
        h(1,i) = 1.0_dbl + amp * cos(x(i))
      ELSE IF (SHAPE=="1997") THEN
        h(1,i) = 1.0_dbl - amp * cos(W*x(i))
      END IF
    END DO

    h0 = h(1,:)
    !PRINT*, "h0 is ", h(1,:)
    IF (h(1,2)/=h(1,2)) THEN
      PRINT*, "STUPID h0"
    END IF
    
    ALLOCATE(a(2*n+9,6))
    ALLOCATE(b(2*n+9,3))
    a = coef_3(x)
    b = coef_1(x)

    k = 1
    his_count = 1

    ALLOCATE(h_min(1))
    ALLOCATE(time(1))
    ALLOCATE(h_his(1,2*n+9))
    h_min(1) = minval(h(1,:))
    time(1) = 0.0_dbl
    h_his(1,:) = h(1,:)

    DO WHILE ((k<m).and.(test==0))
      k = k + 1
      h(2,:) = non_uniform_2nd_order_newton(h(1,:),a,b,Ca,A_bar,dt,theta,tol_newton)
      DO i=5,2*n+5
        IF ((h(2,i)<0)) THEN
          test = 1
          PRINT*, "0", k
        ELSE IF (h(2,i)/=h(2,i)) THEN
          test = 1
          !PRINT*, "STUPID", k, i
        END IF
      END DO
      IF (mod(k,10)==0) THEN
        PRINT*, k
        PRINT*, minval(h(2,:))
      END IF
      h(1,:) = h(2,:)
      IF (mod(k,time_gap)==0) THEN
        his_count = his_count + 1
        ALLOCATE(h_his_temp(his_count-1,2*n+9))
        h_his_temp = h_his
        DEALLOCATE(h_his)
        ALLOCATE(h_his(his_count,2*n+9))
        h_his(1:his_count-1,:) = h_his_temp
        DEALLOCATE(h_his_temp)
        h_his(his_count,:) = h(1,:)
      END IF
      ALLOCATE(h_min_temp(k-1))
      ALLOCATE(time_temp(k-1))
      h_min_temp = h_min
      time_temp = time
      DEALLOCATE(h_min)
      DEALLOCATE(time)
      ALLOCATE(h_min(k))
      ALLOCATE(time(k))
      h_min(1:k-1) = h_min_temp
      time(1:k-1) = time_temp
      DEALLOCATE(h_min_temp)
      DEALLOCATE(time_temp)
      h_min(k) = minval(h(1,:))
      time(k) = dt*(k-1)
    END DO
  END IF

  !!!!!!!!!!!!!!!!!!!!!!! For non_uniform_sin grid !!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (GRID == "non_uniform_sin") THEN
    !dt = 2e-9
    ALLOCATE(x(2*n+9))
    CALL grid_sin(n,x,lower_bnd,upper_bnd)

    ALLOCATE(h(2,2*n+9))
    ALLOCATE(h0(2*n+9))

    h = 0.0_dbl
    ! initial condition
    DO i=1,2*n+9
      IF (SHAPE=="1999") THEN
        h(1,i) = 1.0_dbl + amp * cos(x(i))
      ELSE IF (SHAPE=="1997") THEN
        h(1,i) = 1.0_dbl - amp * cos(W*x(i))
      END IF
    END DO

    h0 = h(1,:)
    !PRINT*, "h0 is ", h(1,:)
    IF (h(1,2)/=h(1,2)) THEN
      PRINT*, "STUPID h0"
    END IF
    
    ALLOCATE(a(2*n+9,6))
    ALLOCATE(b(2*n+9,3))
    a = coef_3(x)
    b = coef_1(x)

    k = 1
    his_count = 1

    ALLOCATE(h_min(1))
    ALLOCATE(time(1))
    ALLOCATE(h_his(1,2*n+9))
    h_min(1) = minval(h(1,:))
    time(1) = 0.0_dbl
    h_his(1,:) = h(1,:)

    DO WHILE ((k<m).and.(test==0))
      k = k + 1
      h(2,:) = non_uniform_2nd_order_newton(h(1,:),a,b,Ca,A_bar,dt,theta,tol_newton)
      DO i=5,2*n+5
        IF ((h(2,i)<0)) THEN
          test = 1
          PRINT*, "0", k
        ELSE IF (h(2,i)/=h(2,i)) THEN
          test = 1
          !PRINT*, "STUPID", k, i
        END IF
      END DO
      IF (mod(k,10)==0) THEN
        PRINT*, k
        PRINT*, minval(h(2,:))
      END IF
      h(1,:) = h(2,:)
      IF (mod(k,time_gap)==0) THEN
        his_count = his_count + 1
        ALLOCATE(h_his_temp(his_count-1,2*n+9))
        h_his_temp = h_his
        DEALLOCATE(h_his)
        ALLOCATE(h_his(his_count,2*n+9))
        h_his(1:his_count-1,:) = h_his_temp
        DEALLOCATE(h_his_temp)
        h_his(his_count,:) = h(1,:)
      END IF
      ALLOCATE(h_min_temp(k-1))
      ALLOCATE(time_temp(k-1))
      h_min_temp = h_min
      time_temp = time
      DEALLOCATE(h_min)
      DEALLOCATE(time)
      ALLOCATE(h_min(k))
      ALLOCATE(time(k))
      h_min(1:k-1) = h_min_temp
      time(1:k-1) = time_temp
      DEALLOCATE(h_min_temp)
      DEALLOCATE(time_temp)
      h_min(k) = minval(h(1,:))
      time(k) = dt*(k-1)
    END DO
  END IF

  CALL write_project3(x,h0,h(2,:),h_min,time,h_his,r_d,filename,ierr)



  ! write data
  OPEN(9, FILE = 'data1.txt', FORM = 'formatted')
  23 FORMAT(3 (ES23.12E3))

  DO i=1,2*n+9
    WRITE(9,23) x(i), h0(i), h(2,i)
  END DO

  CLOSE(9)

  !PRINT*, ALLOCATED(x)
  !PRINT*, ALLOCATED(dx)
END PROGRAM

