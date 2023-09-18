
! Double precsion (real*8)  or  Single Precision (real*4)       (default: double precision)

module OCFD_precision
implicit none
include "mpif.h"
!------For Doubleprecision  (real*8)--------------------------------------------------------------------------
integer,parameter::OCFD_REAL_KIND=8,  OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION   ! double precison computing
! ------For Single precision (real*4)-------------------------------------------------------------------------
!     integer,parameter::OCFD_REAL_KIND=4,  OCFD_DATA_TYPE=MPI_REAL             !  single precision computing
 end module OCFD_precision