!--------------------------------------------------------------------
! Mesh for jet combusion, copyright by Li Xinliang, lixl@imech.ac.cn
! Ver 1.0
! 以圆管直径无量纲化
!
! y (或r)方向网格设置：   0<y<Ya1, 均匀网格 （ny1个网格）;
!                   y> Ya1, 等比拉伸网格； 其中前ny2个网格拉伸比例为Yat1, 后ny3个网格的拉伸比例为Yat2
! x 方向网格设置：
!                   0<x<Xa, 均匀网格 （nx1个网格）
!                   x>Xa, 等比拉伸网格 (nx2个网格，拉伸比例Xet)   
!  单位： 国际单位制 (m)
!----------------------------------------------------------------------
   
   implicit none
    real*8,allocatable,dimension(:):: xx,yy,zz,y1
	integer:: nx,ny,nz,nx1,nx2,ny1,ny2,ny3,nyp,i,j,k
    real*8:: Xa,Xet,Ya,Yet1,Yet2
    
    open(80,file="mesh.in")
    read(80,*)
	read(80,*)  Xa, nx1, Xet, nx2
	read(80,*)
	read(80,*)  Ya, ny1, Yet1, ny2, Yet2, ny3
	read(80,*) 
    close(80)

	nx=nx1+nx2
    nyp=ny1+ny2+ny3
	ny=2*nyp-1
    nz=ny

 	print*, "nx,ny,nz=",nx,ny,nz
    allocate(xx(nx),yy(ny),zz(nz),y1(nyp))

	do i=1,nx1
	  xx(i)=Xa*(i-1.d0)/(nx1-1.d0)
	enddo
    do i=nx1+1,nx
	  xx(i)=xx(i-1)+Xet*(xx(i-1)-xx(i-2))
	enddo
   
	do j=1,ny1
	  y1(j)=Ya*(j-1.d0)/(ny1-1.d0)
	enddo
	
    do j=ny1+1,ny1+ny2
	 y1(j)=y1(j-1)+Yet1*(y1(j-1)-y1(j-2))
	enddo

	do j=ny1+ny2+1, ny1+ny2+ny3
	 y1(j)=y1(j-1)+Yet2*(y1(j-1)-y1(j-2))
    enddo

	do j=1,nyp
	  yy(j)=-y1(nyp+1-j)
	enddo
	
	do j=nyp+1,ny
	  yy(j)=y1(j-nyp+1)
	enddo
    
	do j=1,ny
	 zz(j)=yy(j)
	enddo


 
!-------------------------------------------------------------------------------

    open(99,file="x1d.dat")
	do i=1,nx
	write(99,*) i, xx(i)
	enddo
    close(99)

    open(99,file="y1d.dat")
	do i=1,ny
	write(99,*) i, yy(i)
	enddo
    close(99)

    open(99,file="OCFD-grid.dat",form="unformatted")
	write(99) xx 
	write(99) yy 
	write(99) zz 
	close(99)
	
	open(99,file="mesh2d.dat")
	write(99,*) "variables= x, y"
	write(99,*) "zone i= ", nx, " j= ", ny
	do j=1,ny
	do i=1,nx
	write(99,*) xx(i), yy(j)
	enddo
	enddo
	close(99)

		
  end
