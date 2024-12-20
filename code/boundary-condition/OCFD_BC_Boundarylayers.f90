! boundary conditions： boundary-layers  (transition/ STBLI/ ...) 
! 2024-5-12 (ver 3.1)  Wall jet boundary condition is added
!---------------------------------------------------------------------
   module BC_data 
   use OCFD_precision
   implicit none 
   TYPE bc_type
   integer:: bc_init=0
   Real(kind=OCFD_REAL_KIND),allocatable:: flow_inlet0d(:),flow_upper0d(:), &
                         flow_inlet1d(:,:), flow_upper1d(:,:), &
                         flow_inlet2d(:,:,:), flow_upper2d(:,:,:) 
   Real(kind=OCFD_REAL_KIND),allocatable:: wall_pertb(:,:)      ! Wall blow-and-suction
   integer,allocatable:: NonReflect_upper(:,:)          
   
   end TYPE 
   TYPE (bc_type):: BC 
   end    
   
	
	subroutine ocfd_bc_BoundaryLayer
    use flow_data
	use BC_data  
	implicit none 
	integer i,j,k,m
	integer:: BcData_inlet, BcData_upper,bc_upper_nonref,bc_outlet,bc_dis_type,bc_dis_mt,bc_dis_mz, Iflag_wall_jet
	Real(kind=OCFD_REAL_KIND):: Tw,Wall_Xinit,bc_dis_A,bc_dis_Xbegin,bc_dis_Xend,bc_dis_ZL,bc_dis_freq
	Real(kind=OCFD_REAL_KIND)::  ht, tmp,A0
	Real(kind=OCFD_REAL_KIND):: di1(N_SPEC),di2(N_SPEC),di3(N_SPEC),p1,p2,p3,d0,d1
!---------------------------------------------------
     BcData_inlet=nint(Para%BC_para(1))      ! inlet: (0 free-stream;  1 1d data;  2 2d data)
	 BcData_upper=nint(Para%Bc_para(2))      ! upper boundary: (0 free-stream, 1 1d data ,  2 2d data, -1 None)
	 
	 bc_upper_nonref=nint(Para%Bc_para(3))   ! upper boundary: 0 auto,  1 Non-Reflection,  2 Dirichlet 
	 bc_outlet=nint(Para%Bc_para(4))         ! outlet: 0 Non-Reflection,   1 1st order extrapolation,  2 2nd order extrapolation
     Tw=Para%Bc_para(5)                      ! wall temperature
	 Wall_Xinit=Para%Bc_para(6)              ! x location of the wall leading 

! wall perturbation--------	 
	 bc_dis_type=nint(Para%Bc_para(7))     ! wall disturbance type (0 none ;  1 multi-wave blow-and-suction, Ref: Rai MM, AIAA 95-0583)
	 bc_dis_A=Para%Bc_para(8)              ! Amplitude of wall disturbance     
	 bc_dis_Xbegin=Para%Bc_para(9)	       ! Initial location of wall disturbance
	 bc_dis_Xend=Para%Bc_para(10)          ! End location of wall disturbance
	 bc_dis_mt=nint(Para%Bc_para(11))     ! multi-frequency
	 bc_dis_mz=nint(Para%Bc_para(12))     ! multi-wavenumber
	 bc_dis_ZL=Para%Bc_para(13)           ! Spanwise Length
     bc_dis_freq=Para%Bc_para(14)         ! base frequancy of disturbance  
!--------	 
	 Iflag_wall_jet=nint(Para%Bc_para(15))
	 
 !-------- j=1 wall ----------------		 
     if(BC%bc_init ==0 ) then 
	  call init_bcdata_boundary (BcData_inlet, BcData_upper)  ! read inlet & upper data 
	  call init_bc_nonRef (bc_upper_nonref)             ! upper boundary:  Non-Reflection or Dirichlet 
      call init_bc_wall_perturbation (bc_dis_type, bc_dis_A, bc_dis_Xbegin, bc_dis_Xend, bc_dis_mt, bc_dis_mz, bc_dis_ZL) 	   
      Bc%bc_init=1                     ! Run only for initial time 
	 endif 
	 
!----------------Inlet boundary -------------------------------------
  
 if(npx .eq. 0) then 
  if( BcData_inlet == 0) then 
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=BC%flow_inlet0d(1)
    u(1,j,k)=BC%flow_inlet0d(2)
    v(1,j,k)=BC%flow_inlet0d(3)
    w(1,j,k)=BC%flow_inlet0d(4) 
    T(1,j,k)=BC%flow_inlet0d(5)
    do m=1,N_SPEC 
    di(1,j,k,m)=BC%flow_inlet0d(5+m)	
	enddo 
	
   enddo 
   enddo 
  else if(BcData_inlet == 1) then   
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=BC%flow_inlet1d(j,1) 
    u(1,j,k)=BC%flow_inlet1d(j,2)  
    v(1,j,k)=BC%flow_inlet1d(j,3) 
    w(1,j,k)=BC%flow_inlet1d(j,4) 
    T(1,j,k)=BC%flow_inlet1d(j,5)
    do m=1,N_SPEC 
    di(1,j,k,m)=BC%flow_inlet1d(j,5+m)
    enddo 	
   enddo 
   enddo  	 
  else if (BcData_inlet == 2) then	 
   do k=1,nz 
   do j=1,ny 
    d(1,j,k)=BC%flow_inlet2d(j,k,1) 
    u(1,j,k)=BC%flow_inlet2d(j,k,2)  
    v(1,j,k)=BC%flow_inlet2d(j,k,3) 
    w(1,j,k)=BC%flow_inlet2d(j,k,4) 
    T(1,j,k)=BC%flow_inlet2d(j,k,5)  
    do m=1,N_SPEC 
    di(1,j,k,m)=BC%flow_inlet2d(j,k,5+m)
    enddo 	
   enddo 
   enddo  	   
  else 
   print*, "The inlet data is not supported !"
  endif 
 endif 
 
!----------------upper boundary ------------------------------------- 
  if(npy .eq. npy0-1) then 
   do k=1,nz 
   do i=1,nx 
   if(Bc%NonReflect_upper(i,k) .eq. 0) then   ! Dirichlet BC 
   if(BcData_upper .eq. 0) then 
    d(i,ny,k)=Bc%flow_upper0d(1) 
    u(i,ny,k)=Bc%flow_upper0d(2) 
    v(i,ny,k)=Bc%flow_upper0d(3) 
    w(i,ny,k)=Bc%flow_upper0d(4)  
    T(i,ny,k)=Bc%flow_upper0d(5) 
    do m=1,N_SPEC 
    di(i,ny,k,m)=	Bc%flow_upper0d(5+m)
	enddo 
	
   else if (BcData_upper .eq. 1) then
    d(i,ny,k)=BC%flow_upper1d(i,1) 
    u(i,ny,k)=BC%flow_upper1d(i,2)    
	v(i,ny,k)=BC%flow_upper1d(i,3) 
    w(i,ny,k)=BC%flow_upper1d(i,4)    
    T(i,ny,k)=BC%flow_upper1d(i,5)
	do m=1,N_SPEC 
	di(i,ny,k,m)=BC%flow_upper1d(i,5+m)
	enddo 
	
   else if(BcData_upper .eq. 2) then
    d(i,ny,k)=BC%flow_upper2d(i,k,1) 
    u(i,ny,k)=BC%flow_upper2d(i,k,2)    
	v(i,ny,k)=BC%flow_upper2d(i,k,3) 
    w(i,ny,k)=BC%flow_upper2d(i,k,4)    
    T(i,ny,k)=BC%flow_upper2d(i,k,5)  
	do m=1,N_SPEC 
	di(i,ny,k,m)=BC%flow_upper2d(i,k,5+m)
	enddo
   endif 
   endif 
   enddo 
   enddo
  endif 
!--------------------------------------------------------
!----wall----------------------------------------------- 
     call get_ht_multifrequancy(ht,tt,bc_dis_mt, bc_dis_freq)
     if(npy.eq.0) then
     do k=1,nz
     do i=1,nx
	 if(Axx(i,1,k) < Wall_Xinit ) then 
	  ! 壁面之前 （对称条件） 
	  d(i,1,k)=d(i,2,k) 
	  u(i,1,k)=u(i,2,k)
	  v(i,1,k)=v(i,2,k) 
	  w(i,1,k)=w(i,2,k) 
	  T(i,1,k)=T(i,2,k) 
	  do m=1,N_SPEC 
	  di(i,1,k,m)=di(i,2,k,m) 
	  enddo 
	  
	 else    ! 壁面
	   tmp=1.d0/(sqrt(Aix(i,1,k)**2+Aiy(i,1,k)**2+Aiz(i,1,k)**2)) 
	   A0=bc_dis_A*ht* BC%wall_pertb(i,k)
	   u(i,1,k)=A0*Aix(i,1,k)*tmp      !  Aix*tmp  normal vector 
	   v(i,1,k)=A0*Aiy(i,1,k)*tmp
	   w(i,1,k)=A0*Aiz(i,1,k)*tmp
       if(Tw.gt.0) then
        T(i,1,k)=Tw
       else
        T(i,1,k)=(4.d0*T(i,2,k)-T(i,3,k))/3.d0        ! dT/dy=0 
       endif
 
	    do m=1,N_SPEC 
		di2(m)=di(i,2,k,m)
		di3(m)=di(i,3,k,m)
		enddo 
		
		call comput_P(di2,T(i,2,k),p2)
		call comput_p(di3,T(i,3,k),p3)
		p1=(4.d0*p2-p3)/3.d0           !pw:  dp/dy=0
!-----------------------------
        if(Para%Wall_Catalysis==Non_Catalystic) then    ! 非催化壁 
         do m=1,N_SPEC 
         di1(m)=di2(m)/d(i,2,k)
		 enddo
        else 
         print*, "In this version, only Non Catalytic wall is supported !"
		 stop
		endif 
		
		call comput_d(di1,T(i,1,k),p1,d1)
		d(i,1,k)=d1 
		do m=1,N_SPEC 
		di(i,1,k,m)=di1(m)*d1 
		enddo 
	   endif	
	  enddo 
	  enddo 
	  
	 endif 
!-------outlet bounary ------------------------------------------
   if(npx .eq. npx0-1) then 
   if(bc_outlet==1) then  
	 do k=1,nz 
	 do j=1,ny 
	  d(nx,j,k)=d(nx-1,j,k)
	  u(nx,j,k)=u(nx-1,j,k)
	  v(nx,j,k)=v(nx-1,j,k)
	  w(nx,j,k)=w(nx-1,j,k)
	  T(nx,j,k)=T(nx-1,j,k)
	  do m=1,N_SPEC 
	  di(nx,j,k,m)=di(nx-1,j,k,m)
	  enddo 
	 enddo 
	 enddo 
    elseif (bc_outlet==2) then 
	 do k=1,nz 
	 do j=1,ny 
	  d(nx,j,k)=2.d0*d(nx-1,j,k)-d(nx-2,j,k)
	  u(nx,j,k)=2.d0*u(nx-1,j,k)-u(nx-2,j,k)
	  v(nx,j,k)=2.d0*v(nx-1,j,k)-v(nx-2,j,k)
	  w(nx,j,k)=2.d0*w(nx-1,j,k)-w(nx-2,j,k)
	  T(nx,j,k)=2.d0*T(nx-1,j,k)-T(nx-2,j,k)
	  do m=1,N_SPEC 
	  di(nx,j,k,m)=2.d0*di(nx-1,j,k,m)-di(nx-2,j,k,m)
      enddo 	  
	 enddo 
	 enddo 	
   endif 	
   endif 

   
!  Wall Jet boundary condition  壁面射流边界条件
   if(Iflag_wall_jet .eq. 1) then 
    call bc_wall_jet1               ! 圆孔壁面射流
   elseif( Iflag_wall_jet .eq. 2) then 
    call bc_wall_jet2               ! 任意形状壁面射流， 通过窗口文件进行描述
   endif 
   


    end 	
	
!---------------------------------------------------------------------	
! read  boundary data (inlet & upper) 
 subroutine init_bcdata_boundary (BcData_inlet, BcData_upper)
  use flow_data
  use BC_data 
  implicit none  
  integer::  BcData_inlet, BcData_upper
  integer::  i,j,k,m,i1,j1,k1,ierr
  Real(kind=OCFD_REAL_KIND):: tmp,d0 
  Real(kind=OCFD_REAL_KIND),allocatable:: flow1d0(:,:),flow2d0(:,:)

  
!  BcData_inlet=nint(Para%BC_para(1))      ! inlet type (1 free-stream;  2 1d data;  2 2d data)
!  BcData_upper=nint(Para%Bc_para(2))
  


!----------------read inlet data file  (1d formatted or 2d unformatted )  
 if(BcData_inlet == 0) then          ! 0d (uniform) inlet data file
  allocate(BC%flow_inlet0d(5+N_SPEC))
  if(my_id .eq. 0) then 
	open(99,file="flow0d-inlet-comb.dat")
	read(99,*) 
    read(99,*)
	read(99,*) (Bc%flow_inlet0d(k),k=1,5+N_SPEC)        ! d,u,v,w,T,di(1),di(2),....,di(N_SPEC) (di(:) 也可以是质量比分 ai(:))
    close(99)
    print*, "read 1d inlet data OK"
! Re-scaling of di   归一化，使得 di(1)+di(2)+......+di(N_SPEC)=d
    d0=0.d0 
	do m=1,N_SPEC 
    d0=d0+Bc%flow_inlet0d(5+m)
	enddo 
	do m=1,N_SPEC             
	Bc%flow_inlet0d(5+m)=Bc%flow_inlet0d(5+m)*Bc%flow_inlet0d(1)/d0
	enddo 
   endif 
   
   call MPI_bcast(BC%flow_inlet0d(1),5+N_SPEC,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)   

  
 else if(BcData_inlet == 1) then        ! 1d inlet data file  
   allocate(BC%flow_inlet1d(ny,5+N_SPEC), flow1d0(ny_global,5+N_SPEC))
   if(my_id .eq. 0) then 
	open(99,file="flow1d-inlet-comb.dat")
	read(99,*) 
	do j=1,ny_global
	read(99,*) tmp, (flow1d0(j,k),k=1,5+N_SPEC)        ! d,u,v,w,T,di(1),di(2)...,di(N_SPEC)
	enddo 
    close(99)
    print*, "read 1d inlet data OK"
	
 ! Rescaleing of di(:)   
    do j=1,ny_global
	d0=0.d0 
	do m=1,N_SPEC 
    d0=d0+flow1d0(j,5+m)
	enddo 
	do m=1,N_SPEC             
	flow1d0(j,5+m)=flow1d0(j,5+m)*flow1d0(j,1)/d0
	enddo 
    enddo 
	
   endif 
   
   call MPI_bcast(flow1d0,ny_global*(5+N_SPEC),OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)   
   do j=1,ny 
   j1=j_offset(npy)+j-1
   BC%flow_inlet1d(j,:)=flow1d0(j1,:) 
   enddo 
   deallocate(flow1d0) 
  
  else if(BcData_inlet == 2)  then   !2d inlet data file 
   allocate(BC%flow_inlet2d(ny,nz,5+N_SPEC), flow2d0(ny_global,nz_global))
   if(my_id .eq. 0)  open(99,file="flow2d-inlet-comb.dat",form="unformatted")
    do m=1,5+N_SPEC
    if(my_id .eq. 0) read(99) flow2d0 
	
    call MPI_bcast(flow2d0,ny_global*nz_global ,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 	
	do k=1,nz 
	do j=1,ny 
	 j1=j_offset(npy)+j-1
	 k1=k_offset(npz)+k-1 
	 BC%flow_inlet2d(j,k,m)=flow2d0(j1,k1) 
	enddo 
	enddo 
	enddo 
	if(my_id .eq. 0) then 
	  close(99)
	  print*, "read 2d inlet data OK"
	endif 
    deallocate(flow2d0) 
	
! Rescaling of di(:) ------------------------------------ 
   do k=1,nz
   do j=1,ny 
   d0=0.d0 
   do m=1,N_SPEC 
   d0=d0+BC%flow_inlet2d(j,k,5+m)
   enddo 
   do m=1,N_SPEC 
   BC%flow_inlet2d(j,k,5+m)=BC%flow_inlet2d(j,k,5+m)*BC%flow_inlet2d(j,k,1)/d0
   enddo       
   enddo
   enddo    
	
  endif 
  
 !----------------read upper-boundary data file  (1d formatted or 2d unformatted )  
 
  if(BcData_upper == 0) then          ! 0d (uniform) upper data file
  allocate(BC%flow_upper0d(5+N_SPEC))
  if(my_id .eq. 0) then 
	open(99,file="flow0d-upper-comb.dat")
	read(99,*) 
    read(99,*)
	read(99,*) (Bc%flow_upper0d(k),k=1,5+N_SPEC)        ! d,u,v,w,T,di(1),di(2),....,di(N_SPEC)
    close(99)
    print*, "read 1d inlet data OK"
! Re-scaling of di   归一化，使得 di(1)+di(2)+......+di(N_SPEC)=d
    d0=0.d0 
	do m=1,N_SPEC 
    d0=d0+Bc%flow_upper0d(5+m)
	enddo 
	do m=1,N_SPEC             
	Bc%flow_upper0d(5+m)=Bc%flow_upper0d(5+m)*Bc%flow_upper0d(1)/d0
	enddo 
   endif 
   
   call MPI_bcast(BC%flow_upper0d(1),5+N_SPEC,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)   
 
  else if(BcData_upper == 1) then        ! 1d upper-boundary data file  
    allocate(BC%flow_upper1d(nx,5+N_SPEC), flow1d0(nx_global,5+N_SPEC))
   if(my_id .eq. 0) then 
    open(99,file="flow1d-upper-comb.dat")
	read(99,*) 
	do i=1,nx_global
	read(99,*) tmp, (flow1d0(i,k),k=1,5+N_SPEC)        ! d,u,v,w,T,di(:)
	enddo 
    close(99)
	print*, "read 1d upper boundary data OK"
 
 ! Rescaleing of di(:)   
    do i=1,nx_global
	d0=0.d0 
	do m=1,N_SPEC 
    d0=d0+flow1d0(i,5+m)
	enddo 
	do m=1,N_SPEC             
	flow1d0(i,5+m)=flow1d0(i,5+m)*flow1d0(i,1)/d0
	enddo 
    enddo 	
	
	
   endif 
   
   call MPI_bcast(flow1d0,nx_global*(5+N_SPEC),OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)   
   do i=1,nx 
   i1=i_offset(npx)+i-1
   BC%flow_upper1d(i,:)=flow1d0(i1,:) 
   enddo 
   deallocate(flow1d0) 
 
 
  else if(BcData_upper == 2)  then   !2d boundary data file 
   allocate(BC%flow_upper2d(nx,nz,5+N_SPEC), flow2d0(nx_global,nz_global))
    if(my_id .eq. 0)  open(99,file="flow2d-upper-comb.dat",form="unformatted")
    do m=1,5+N_SPEC
    if(my_id .eq. 0) read(99) flow2d0 
	
    call MPI_bcast(flow2d0,nx_global*nz_global ,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 	
	do k=1,nz 
	do i=1,nx 
	 i1=i_offset(npx)+i-1
	 k1=k_offset(npz)+k-1 
	 BC%flow_upper2d(i,k,m)=flow2d0(i1,k1) 
	enddo 
	enddo 
	enddo 
	if(my_id .eq. 0) then 
	 close(99)
     print*, "read 2d upper boundary data OK"
	endif
	deallocate(flow2d0) 

! Rescaling of di(:) ------------------------------------ 
   do k=1,nz
   do i=1,nx 
   d0=0.d0 
   do m=1,N_SPEC 
   d0=d0+BC%flow_upper2d(i,k,5+m)
   enddo 
   do m=1,N_SPEC 
   BC%flow_upper2d(i,k,5+m)=BC%flow_upper2d(i,k,5+m)*BC%flow_upper2d(i,k,1)/d0
   enddo       
   enddo
   enddo    	
	
	
  endif  
 end 
 
 
 !-----------------------------------------------------
  subroutine init_bc_nonRef (Bc_NonReflect_upper)
  use flow_data
  use BC_data 
  implicit none  
  integer::  Bc_NonReflect_upper,i,j,k
  Real(kind=OCFD_REAL_KIND)::  sy ,sy0
  Real(kind=OCFD_REAL_KIND),parameter:: epsl=1.d-6 
  
!  Bc_NonReflect_upper=nint(Para%Bc_para(3))      ! 0 Auto,  1 Non-Reflection,  2 Dirichlet boudary
  
  sy0=tan(Para%AoA)
  allocate(BC%NonReflect_upper(nx,nz)) 
  Bc%NonReflect_upper=0
  
  
!  if(npy .eq. 0) then                           ! Bug removed 2021-5-4
  if(npy .eq. npy0-1) then 
    if(Bc_NonReflect_upper == 1) then 
	  BC%NonReflect_upper=1            ! Non-Reflection boundary 
   else if (Bc_NonReflect_upper == 2) then  
      BC%NonReflect_upper=0            ! Dirichlet boundary 
   else if (Bc_NonReflect_upper == 0) then      ! Auto type (Non-Reflection /Dirichlet)
    do k=1,nz 
	do i=1,nx 
	  if(i==1 .and. npx .eq. 0) then
	   sy=(Ayy(i+1,ny,k)-Ayy(i,ny,k))/(Axx(i+1,ny,k)-Axx(i,ny,k))
	  else 
	   sy=(Ayy(i,ny,k)-Ayy(i-1,ny,k))/(Axx(i,ny,k)-Axx(i-1,ny,k))
      endif 
	 if(sy >= sy0+ epsl) then 
	   Bc%NonReflect_upper(i,k)=0       ! Dirichlet BC 
	 else  
	  Bc%NonReflect_upper(i,k)=1        ! Non-Reflection BC
	 endif 
	enddo
	enddo 
   endif    
  endif 
  
end 	
	
   subroutine init_bc_wall_perturbation (bc_dis_type, bc_dis_A, bc_dis_Xbegin, bc_dis_Xend, bc_dis_mt, bc_dis_mz, bc_dis_ZL) 
! TM(k), Amplitude; faiz,fait: phase (random); mzmax, mtmax: wavenumbers
! consider wall jet
   use flow_data
   use BC_data
   implicit none  
   integer i,j,k,m,ierr
   real(kind=OCFD_REAL_KIND),allocatable :: faiz(:),Zl(:)
   real(kind=OCFD_REAL_KIND),parameter:: PI=3.14159265358979d0
   real(kind=OCFD_REAL_KIND):: ztmp,seta,fx,gz,rtmp
   integer:: bc_dis_type,bc_dis_mt,bc_dis_mz
   real(kind=OCFD_REAL_KIND):: bc_dis_A,bc_dis_Xbegin,bc_dis_Xend,bc_dis_ZL
   
   
!	 bc_dis_type=nint(Para%Bc_para(7))     ! wall disturbance type (0 none ;  1 multi-wave blow-and-suction, Ref: Rai MM, AIAA 95-0583)
!	 bc_dis_A=Para%Bc_para(8)              ! Amplitude of wall disturbance     
!	 bc_dis_Xbegin=Para%Bc_para(9)	       ! Initial location of wall disturbance
!	 bc_dis_Xend=Para%Bc_para(10)          ! End location of wall disturbance
!	 bc_dis_mt=nint(Para%Bc_para(11))     ! multi-frequency
!	 bc_dis_mz=nint(Para%Bc_para(12))     ! multi-wavenumber
!	 bc_dis_ZL=Para%Bc_para(13)           ! Spanwise Length
 
	 allocate(BC%wall_pertb(nx,nz))
     BC%wall_pertb=0.d0 
	 
   if(bc_dis_type == 0) then 
      BC%wall_pertb=0.d0       ! no disturbation
 !-----------------------------------------------
   else if (	bc_dis_type == 1) then 	 ! multi-wavenumber perturbation (Rai MM AIAA 95-0583) 
	 
	 if(bc_dis_mz >0) then  
	  allocate(faiz(bc_dis_mz),Zl(bc_dis_mz))
      ztmp=0.d0
      do k=1,bc_dis_mz
      call random_number(faiz(k))
      if(k.eq.1) then
       Zl(k)=1.d0
      else
       zl(k)=zl(k-1)/1.25d0
      endif
       ztmp=ztmp+Zl(k)
      enddo
      do k=1,bc_dis_mz
       zl(k)=zl(k)/ztmp
      enddo
      call MPI_bcast(faiz(1),bc_dis_mz,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     endif 

    do k=1,nz
	do i=1,nx 
	 if(Axx(i,1,k) >=bc_dis_Xbegin .and. Axx(i,1,k) <= bc_dis_Xend) then 
	   seta=2.d0*PI*(Axx(i,1,k)-bc_dis_Xbegin)/(bc_dis_Xend-bc_dis_Xbegin)
	   fx=4.d0/sqrt(27.d0)*sin(seta)*(1.d0-cos(seta))
	 else 
       fx=0.d0 
     endif 
   
     gz=0.d0
     seta=Azz(i,1,k)/bc_dis_ZL
     if(bc_dis_mz > 0) then
	  do m=1,bc_dis_mz
       gz=gz+Zl(m)*sin(2.d0*PI*m* (seta+faiz(m)) )
      enddo
     else if(bc_dis_mz == 0) then
      gz=1.d0
	 else if(bc_dis_mz < 0) then
      gz=sin(-2.d0*PI*bc_dis_mz*seta)
	 endif
	 BC%wall_pertb(i,k)=fx*gz
	enddo
    enddo 
    if(bc_dis_mz >0)  deallocate(faiz,Zl)
!----------------------------------------------- 
  else if (	bc_dis_type == 2) then 	       ! random disturbance 
	do k=1,nz 
	do i=1,nx 
	 if(Axx(i,1,k) >=bc_dis_Xbegin .and. Axx(i,1,k) <= bc_dis_Xend) then         ! a Bug removed (2021-6-21)
	  call random_number(rtmp)
	  BC%wall_pertb(i,k)=2.d0*(rtmp-0.5d0)          ! -1< wall_pertb < 1	  
	 else 
	  BC%wall_pertb(i,k)=0.d0 
	 endif  
    enddo 
	enddo 
  endif 
  end  

  
!--------------------------------	
! perturbation of Rai, see:  Rai MM, AIAA 95-0583
     subroutine get_ht_multifrequancy(ht,tt,mtmax,beta)
     use OCFD_constants
	 use Para_mpi
     implicit none
     integer mtmax,k,m,ierr
     integer,save:: Kflag=0
     real(kind=OCFD_REAL_KIND):: ht,tt,beta,Ttmp,rand_x
     real(kind=OCFD_REAL_KIND),parameter:: PI=3.14159265358979d0
     real(kind=OCFD_REAL_KIND),allocatable,save:: TM(:),fait(:)          ! Amplitute and random phase angle
     
     if(Kflag .eq. 0) then
 	 Kflag=1
     allocate(TM(mtmax),fait(mtmax))

	 Ttmp=0.d0
     do k=1,mtmax
      call random_number(rand_x)
      fait(k)=rand_x
     if(k.eq.1) then
      TM(k)=1.d0
     else
      TM(k)=TM(k-1)/1.25d0
     endif
      Ttmp=Ttmp+TM(k)
     enddo
      do k=1,mtmax
        TM(k)=TM(k)/Ttmp
      enddo
     call MPI_bcast(fait(1),mtmax,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
     endif

     ht=0.d0
     if(beta .gt. 0.d0) then      
      do m=1,mtmax
       ht=ht+TM(m)*sin(m*beta*tt*2.d0*PI+2.d0*PI*fait(m))             ! beta is base frequency (in Hz)
       enddo
      else
       ht=1.d0
      endif
      end
	  

	  
	  