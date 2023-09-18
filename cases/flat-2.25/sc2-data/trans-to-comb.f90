    implicit none 
    real*8,dimension(:),allocatable:: x,y,z,z1
	real*8,dimension(:),allocatable:: d1,u1,v1,w1,T1
    integer:: nx,ny,nz,nz1,j,k
    real*8::  Linf,Tinf,Re,gamma,Ama,R0,a0,Tw,Uinf,dinf,Amu_inf,Pinf,hz1,Lz,tmp,beta
    real*8:: ci(2)	
!------------------------------------------------------
    nx=2193
	ny=72
	nz=128
	nz1=64
	
	Re=635000.d0
	Ama=2.25d0
	gamma=1.4d0	
    Linf=0.0254d0          ! 1 Inch 
    Tinf=169.44d0 	
    Lz=0.175d0 
	beta=7.845
	
    R0= 287.06d0   ! 空气的气体常数R
	ci(1)=0.233    ! ratio of O2
	ci(2)=0.767    ! ration of N2
	
	a0= sqrt(gamma*R0*Tinf)    ! 参考温度下的声速 
    Uinf=Ama*a0 	
    call Amu_sutherland(Tinf,Amu_inf)
    Dinf=Re*Amu_inf/(Uinf*Linf)
	Pinf=Dinf*R0*Tinf
	print*, "Transfor from non-dimension to SI unit value"
	
	print*, "Dinf, Uinf, Linf,Pinf,t_inf=",Dinf,Uinf,Linf,Pinf,Linf/Uinf  
    print*, "base frequency in Hz is ",  beta/(2.d0*3.1415926535d0)*Uinf/Linf 


    allocate(x(nx),y(ny),z(nz),z1(nz1))
	allocate(d1(ny),u1(ny),v1(ny),w1(ny),T1(ny))



	open(99,file="ocfd-grid.dat",form="unformatted")
	read(99) x
	read(99) y
	close(99) 

    x=x*Linf 
	y=y*Linf 
	
	do k=1,nz1 
	z1(k)=(k-1.d0)*Lz*Linf/nz1  
	enddo 
	
	
    open(99,file="flow1d-inlet.dat")	
    read(99,*)
	do j=1,ny
    read(99,*) tmp, d1(j),u1(j),v1(j),T1(j)
	enddo 
    close(99)
    w1=0.d0 


	d1=d1*Dinf 
	u1=u1*Uinf 
	v1=v1*Uinf
	w1=w1*Uinf 
	T1=T1*Tinf 
	
   open(100,file="grid-si-unit.dat",form="unformatted") 
   write(100) x 
   write(100) y 
   write(100) z1 
   close(100) 
   
   open(100,file="flow1d-inlet-comb.dat")
   write(100,*) "variables=y,d,u,v,w,T,d1,d2"
   do j=1,ny 
    write(100,"(8F16.8)") y(j), d1(j),u1(j),v1(j),w1(j),T1(j),d1(j)*ci(1),d1(j)*ci(2)  
   enddo 
   close(100)
   
   open(100,file="flow0d-upper-comb.dat")
   write(100,*) "flow0d upper for comb"
   write(100,*) "d,u,v,w,T,d1,d2"
   write(100,"(7F16.8)") d1(ny),u1(ny),v1(ny),w1(ny),T1(ny),d1(ny)*ci(1),d1(ny)*ci(2)
   close(100)
   
   end 







    subroutine Amu_sutherland(T,Amu)
    implicit none    
	real*8:: T, Amu
	Amu=1.789d-5*(T/288.15d0)**1.5d0*(288.15d0+110.4d0)/(T+110.4d0)
	end 


