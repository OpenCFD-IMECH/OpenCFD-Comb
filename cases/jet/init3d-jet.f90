! -----  Init 3D ----
! Generate initial data for opencfd-sc, Ver 1.3     
	implicit none
	integer:: nx,ny,nz,i,j,k,m,N_SPEC
    real*8,allocatable,dimension(:,:):: d,u,v,w,T
	real*8,allocatable:: di(:,:,:)
	real*8:: tmp,d1,u1,v1,w1,T1,di1(100)
	
    print*, "Please input nx,ny,nz,N_SPEC ?"
	read(*,*) nx,ny,nz,N_SPEC
	allocate(d(ny,nz),u(ny,nz),v(ny,nz),w(ny,nz),T(ny,nz))
    allocate(di(ny,nz,N_SPEC))
!--------------------

		 open(99,file="flow2d-inlet-jet.dat")
		 read(99,*)
		 read(99,*)
		 do k=1,nz
		 do j=1,ny 
		 read(99,*) tmp, tmp,d(j,k),u(j,k),v(j,k),w(j,k),T(j,k),(di(j,k,m),m=1,N_SPEC) 
  		 enddo 
		 enddo 
         close(99)

    open(100,file="flow2d-inlet.dat",form="unformatted")
	write(100) d
	write(100) u	
	write(100) v	
	write(100) w	
	write(100) T
	do m=1,N_SPEC
	write(100) di(:,:,m)
    enddo 
    close(100)	

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
 	integer:: no,nx,ny,nz,k,i
	real*8:: u(ny,nz)
	real*8,allocatable:: f(:,:,:)
	allocate(f(nx,ny,nz))
	do i=1,nx 
	f(i,:,:)=u(:,:)
	enddo 
	
	do k=1,nz
	write(no) f(:,:,k)
	enddo
	deallocate(f)
	
	end

