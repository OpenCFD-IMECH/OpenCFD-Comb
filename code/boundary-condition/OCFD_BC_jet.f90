	
	subroutine ocfd_bc_jet
    use flow_data
	use BC_data  
	implicit none 
	integer i,j,k,m
	integer:: BcData_inlet, BcData_upper,Iflag_inlet
	Real(kind=OCFD_REAL_KIND):: d0,di1(N_SPEC),ai1(N_SPEC),c1,p1,d1,Et1,d1new
!---------------------------------------------------
!     BcData_inlet=nint(Para%BC_para(1))      ! inlet: (0 free-stream;  1 1d data;  2 2d data)
!	  BcData_upper=nint(Para%Bc_para(2))      ! upper boundary: (0 free-stream, 1 1d data ,  2 2d data, -1 None)
	  BcData_inlet=2      ! 2d data 
	  BcData_upper=-1     ! None 
     Iflag_inlet=nint(Para%BC_para(1))       ! 0 as supersonic inlet (Dirichlet);  1 ( Supersonic/subsonic)
     if(BC%bc_init ==0 ) then 
	  call init_bcdata_boundary(BcData_inlet, BcData_upper)  ! read inlet & upper data 
      Bc%bc_init=1                     ! Run only for initial time 
	 endif 	  
   
!------------Inlet --------------------------------------------------------
  if(npx .eq. 0) then    
!------------Supersonic or subsounci inlet -------------------------
   if(Iflag_inlet ==0) then 
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
   
	do k=1,nz 
    do j=1,ny
     u(1,j,k)=BC%flow_inlet2d(j,k,2)  
     v(1,j,k)=BC%flow_inlet2d(j,k,3) 
     w(1,j,k)=BC%flow_inlet2d(j,k,4) 
     T(1,j,k)=BC%flow_inlet2d(j,k,5) 
	 d1=BC%flow_inlet2d(j,k,1)
	 
	 d0=0.d0 
	 do m=1,N_SPEC 
	 di1(m)=BC%flow_inlet2d(j,k,5+m)
	 d0=d0+di1(m)
	 enddo 
	 ai1=di1/d0 !  density ratio        质量比分 
	 
	 
	 call comput_P(di1,T(1,j,k),p1)
     call comput_E(di1,T(1,j,k),Et1)	
     call comput_soundspeed(c1,p1,d1,Et1) 
	 
	 if(u(1,j,k) < c1  ) then  ! subsounic,  p1=p2 
	  p(1,j,k)=p(2,j,k) 
	  call  comput_d(ai1,T(1,j,k),p(1,j,k),d1new)
	  d(1,j,k)=d1new 
	  do m=1,N_SPEC 
      di(1,j,k,m)=d1new*ai1(m)
	  enddo 
 	 else                 ! supersonic ,   d=d_inf, di=di_inf 
	  d(1,j,k)=d1 
	  do m=1,N_SPEC 
	   di(1,j,k,m)=di1(m) 
	  enddo 
	endif 
	enddo 
    enddo 
   endif 
   endif 
!---------------------------------------
    if(npx .eq. npx0-1) then
	do k=1,nz
	do j=1,ny
     d(nx,j,k)=d(nx-1,j,k)
     u(nx,j,k)=u(nx-1,j,k)
     v(nx,j,k)=v(nx-1,j,k)
     w(nx,j,k)=w(nx-1,j,k)
     T(nx,j,k)=T(nx-1,j,k)
     p(nx,j,k)=p(nx-1,j,k)
     Et(nx,j,k)=Et(nx-1,j,k)
	 do m=1,N_Spec
	 di(nx,j,k,m)=di(nx-1,j,k,m)
     enddo
    enddo
	enddo

   endif   
	 
 
 end 
 
 
 