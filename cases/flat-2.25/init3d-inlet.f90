! -----  Init 3D ----
! Generate initial data for opencfd-sc, Ver 1.3     
	implicit none
	integer:: nx,ny,nz,j,k,m,N_SPEC
    real*8,allocatable,dimension(:,:):: d,u,v,w,T
	real*8,allocatable:: di(:,:,:)
	real*8:: tmp,d1,u1,v1,w1,T1,di1(100)
	
    print*, "Please input nx,ny,nz,N_SPEC ?"
	read(*,*) nx,ny,nz,N_SPEC
	allocate(d(nx,ny),u(nx,ny),v(nx,ny),w(nx,ny),T(nx,ny))
    allocate(di(nx,ny,N_SPEC))
!--------------------

		 open(99,file="flow1d-inlet-comb.dat")
		 read(99,*)
		 do j=1,ny 
		 read(99,*) tmp, d1,u1,v1,w1,T1,(di1(m),m=1,N_SPEC) 
         d(:,j)=d1 
		 u(:,j)=u1
		 v(:,j)=v1
		 w(:,j)=w1 
		 T(:,j)=T1
		 do m=1,N_SPEC 
		 di(:,j,m)=di1(m) 
		 enddo 
		 enddo 
         close(99)



	open(99,file="opencfd-comb0.dat",form="unformatted")
    write(99) 0,0.d0
	call write3d1(99,nx,ny,nz,d)
	call write3d1(99,nx,ny,nz,u)
	call write3d1(99,nx,ny,nz,v)
	call write3d1(99,nx,ny,nz,w)
	call write3d1(99,nx,ny,nz,T)
	do m=1,N_SPEC 
	call write3d1(99,nx,ny,nz,di(1,1,m))
    enddo 	
    close(99)
	deallocate(d,u,v,w,T,di)
 	end
	
	
	subroutine write3d1(no,nx,ny,nz,u)
    implicit none
 	integer:: no,nx,ny,nz,k
	real*8:: u(nx,ny)
	do k=1,nz
	write(no) u
	enddo
	end

