   subroutine init3d 
	  use flow_data
	  implicit none 
	  integer:: i,j,k,m,ierr
      real(kind=OCFD_REAL_KIND):: di1(N_SPEC),Et1
	 call init_mesh
	 call read_flow_data 
	 call Rescaling_di	
	 
     do k=1,nz 
     do j=1,ny
	 do i=1,nx
       f(i,j,k,1)=d(i,j,k)
	   f(i,j,k,2)=d(i,j,k)*u(i,j,k)
	   f(i,j,k,3)=d(i,j,k)*v(i,j,k)
	   f(i,j,k,4)=d(i,j,k)*w(i,j,k)
	   
	   do m=1,N_SPEC 
       di1(m)=di(i,j,k,m)
       enddo 
	   call comput_E(di1,T(i,j,k),Et1)
	   f(i,j,k,5)=Et1+d(i,j,k)* (u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))*0.5d0 
     enddo
	 enddo
     enddo 
   end 
   