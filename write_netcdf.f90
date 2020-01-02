!> @brief Module that writes charge density, potential, field strength
!>        particle trajectory, particle velocity, particle acceleration
!>        into a netCDF file.
!> @Author: JingBang Liu 
!>       
!> @note the is greatly inspired from the
!>       module write_netcdf written by CS Brady and H Ratcliffe
MODULE netcdf_write

  USE kinds
  USE netcdf

  IMPLICIT NONE 

  TYPE :: run_data
    INTEGER :: run_data_n
    REAL(KIND=dbl) :: run_data_dt, run_data_theta
    CHARACTER(LEN=25) :: run_data_grid
  END TYPE

  CONTAINS

  SUBROUTINE write_project3(x,h0,h,h_min,time,r_d,filename,ierr)
    !!!! Define input variables
    REAL(KIND=dbl), DIMENSION(:), INTENT(IN) :: x, h0, h, h_min, time
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(run_data) :: r_d
    INTEGER :: ierr
    !!!! Define dimensions
    INTEGER, PARAMETER :: ndims = 1
    !!!! Define dimensions for output variables
    CHARACTER(LEN=1), DIMENSION(ndims) :: dims_x=(/"x"/)
    CHARACTER(LEN=2), DIMENSION(ndims) :: dims_h0=(/"h0"/)
    CHARACTER(LEN=1), DIMENSION(ndims) :: dims_h=(/"h"/)
    CHARACTER(LEN=5), DIMENSION(ndims) :: dims_h_min=(/"h_min"/)
    CHARACTER(LEN=1), DIMENSION(ndims) :: dims_time=(/"t"/)
    !!!! Define the size ids for output variables
    INTEGER, DIMENSION(ndims) :: sizes_x, sizes_h0, sizes_h, sizes_h_min, sizes_time
    !!!! Define the dimension of ids for output variables
    INTEGER, DIMENSION(ndims) :: dim_ids_x, dim_ids_h0, dim_ids_h, dim_ids_h_min, dim_ids_time
    !!!! Define variable ids for output variables
    INTEGER :: var_id_x, var_id_h0, var_id_h, var_id_h_min, var_id_time
    !!!! Define file id
    INTEGER :: file_id
    INTEGER :: i

    !!!! Get sizes of variables
    sizes_x = size(x)
    sizes_h0 = size(h0)
    sizes_h = size(h)
    sizes_h_min = size(h_min)
    sizes_time = size(time)

    !!!! Create file
    ierr = nf90_create(filename,NF90_CLOBBER,file_id)
    !!!! Check if the file is created successfully
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Add run datas to global attribute
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"n",r_d%run_data_n)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"dt",r_d%run_data_dt)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"theta",r_d%run_data_theta)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
    ierr = nf90_put_att(file_id,NF90_GLOBAL,"gird",r_d%run_data_grid)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Define dim ids
    DO i=1,ndims
      ierr = nf90_def_dim(file_id,dims_x(i),sizes_x(i),dim_ids_x(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    DO i=1,ndims
      ierr = nf90_def_dim(file_id,dims_h0(i),sizes_h0(i),dim_ids_h0(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    DO i=1,ndims
      ierr = nf90_def_dim(file_id,dims_h(i),sizes_h(i),dim_ids_h(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    DO i=1,ndims
      ierr = nf90_def_dim(file_id,dims_h_min(i),sizes_h_min(i),dim_ids_h_min(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    DO i=1,ndims
      ierr = nf90_def_dim(file_id,dims_time(i),sizes_time(i),dim_ids_time(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    !!!! Define var ids
    ierr = nf90_def_var(file_id, "x", NF90_DOUBLE, dim_ids_x, var_id_x)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "h0", NF90_DOUBLE, dim_ids_h0, var_id_h0)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "h", NF90_DOUBLE, dim_ids_h, var_id_h)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "h_min", NF90_DOUBLE, dim_ids_h_min, var_id_h_min)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "time", NF90_DOUBLE, dim_ids_time, var_id_time)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Finish defining metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    !!!! Write x into file
    ierr = nf90_put_var(file_id, var_id_x, x) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "x"
      RETURN
    END IF

    !!!! Write h0 into file
    ierr = nf90_put_var(file_id, var_id_h0, h0) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "h0"
      RETURN
    END IF

    !!!! Write h into file
    ierr = nf90_put_var(file_id, var_id_h, h) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "h"
      RETURN
    END IF

    !!!! Write h_min into file
    ierr = nf90_put_var(file_id, var_id_h_min, h_min) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "h_min"
      RETURN
    END IF

    !!!! Write time into file
    ierr = nf90_put_var(file_id, var_id_time, time) 
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr)), "time"
      RETURN
    END IF

    !!!! Close the file
    ierr = nf90_close(file_id)
    IF (ierr /=nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF
  
  END SUBROUTINE write_project3

END MODULE netcdf_write

