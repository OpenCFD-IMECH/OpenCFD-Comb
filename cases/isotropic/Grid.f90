  ! Mesh and Initial condition for  flow over a Half-Cylinder
 	implicit none 
	integer,parameter:: nx=128, ny=nx, nz=nx
	real*8,parameter::  PI=3.1415926535897d0 
    integer:: i,j,k
	real*8:: xx(nx),yy(ny),zz(nz),hx
    real*8:: Lref
 !-------------------------------------
    Lref=1.d-3   ! 1mm
	hx=Lref*2.d0*PI/nx 

    do i=1,nx 
    xx(i)=(i-1.d0)*hx
    enddo
    yy=xx
	zz=xx
	
   
   open(99,file="OCFD-grid.dat",form="unformatted")
   write(99) xx 
   write(99) yy 
   write(99) zz 
   close(99) 
   end 

!----------------------------------------------
	 



