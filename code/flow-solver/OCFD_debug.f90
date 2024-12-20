
!---------------------------------------------		
     subroutine debug_NaN 
      use flow_data
      implicit none
      logical:: Inan
	  integer:: i,j,k,m
	  
	  do k=1,nz
	  do j=1,ny 
	  do i=1,nx
      Inan=.FALSE. 
	  do m=1,5 
	  Inan= Inan .or. isnan(f(i,j,k,m))
	  enddo 
	  do m=1,N_SPEC
	  Inan=Inan .or. isnan(di(i,j,k,m))
	  enddo 
	  
	  if(Inan) then 
	   print*, "NaN found !!!, i0,j0,k0=", i_offset(npx)+i-1, j_offset(npy)+j-1, k_offset(npz)+k-1 
	   print*, "f(i,j,k,:)=", f(i,j,k,:)
	   print*, "di(i,j,k,:)=",di(i,j,k,:)
	   print*, "npx,npy,npz,i,j,k=",npx,npy,npz,i,j,k 
	   stop
	  endif
	  
	  enddo
	  enddo 
	  enddo 

      end 


