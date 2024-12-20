! -------------------------------------------------------------------------------------------------   
! OpenCFD-Comb Ver 3.x:   3D Chemical Reaction Naver-Stokes Solver by using Finite Difference method                                  
! Copyright by Li Xinliang, LHD, Institute of Mechanics.   Email: lixl@imech.ac.cn        
!  New feature in version 2:  
!   Character flux resconstruction is supported
!   WENO5-type boundary scheme
!   More simple: 3d Jocabian solver is used for all geometry
!  New feature in Version 3:
!   Particles are support  (Langrange-Euler sover)
!------------------------------------------------------------------------------------------------------  
!  OpenCFD-SC/OpenCFD-Comb use Double-Precision for default, If you want SINGLE-PRECISION computation,
!  you can change OCFD_REAL_KIND=8 to OCFD_REAL_KIND=4, and change 
!  OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION to OCFD_DATA_TYPE=MPI_REAL in OCFD_parameters.f90
!-------------------------------------------------------------------------------------------------------
! OpenCFD-Comb using DIMENSIONAL VALUES with SI UNIT  
! 本程序使用国际单位制有量纲计算  (化学反应部分子程序采用g-cm-s-mole 制单位， 程序接口处有单位转换）
!-------------------------------------------------------------------------------------------------------

! 2024-8-28, 为了节省内存，提升稳定性，不再使用MPI_Bsend(),改用MPI_Send(). 不使用MPI_BUFFER_ATTACH
     use flow_para 
	 implicit none 
!     integer,parameter ::IBUFFER_SIZE=10000000
!     real(kind=8):: BUFFER_MPI(IBUFFER_SIZE)
     integer:: ierr 
!----------------------------------------------	  
	   call mpi_init(ierr)                       ! initial of MPI 
	   call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)   ! get my_id 
!      call MPI_BUFFER_ATTACH(BUFFER_MPI,8*IBUFFER_SIZE,ierr)    ! attach buffer for MPI_Bsend( )
!----------------------------

	   call read_parameter      ! read parameters:  opencfd-comb2d.in (flow parameters) , chemical and transport parameters
 	   call partation3d_mpi     ! partation for MPI  (npx0*npy0*npz0  MPI precesses)
	   call set_parameters      ! Set parameters such as Scheme%Bound_Index
	   call NS_Solver	        ! Solve 3d N-S eq. 
	   call mpi_finalize(ierr)
	  end

