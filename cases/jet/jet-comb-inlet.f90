!--------------------------------------------------------------------
! Mesh for jet combusion, copyright by Li Xinliang, lixl@imech.ac.cn
!  Ver 3.0
!  单位： 国际单位制 (m)
!  射流入口条件生成程序：  r<rd 为射流区， 温度T1;   r> rd 为伴流区，温度T2
!  V3.0 内外两层之间存在过渡区, 厚度 rbuf
!  初始时刻 x> xh , r< rh 为高温区， 温度为T3
!  Ver 4.1 入口添加扰动  (有Bug)
!  Ver 4.2 Bug removed,  x>Xh 时，压强仍为 p1 (密度不同)
!----------------------------------------------------------------------
   
   implicit none
    integer:: nx,ny,nz, N_Spec,i,j,k,m,i0
	  real*8:: rd, u1, T1,P1,d1, u2,T2,P2,d2,  rr, AT , xh, Th, rbuf, d0, u0, T0, p0,tmp
	  real*8,allocatable,dimension(:):: ci1,ci2,Mi,ci0,x1d,y1d,z1d
    real*8,allocatable,dimension(:,:):: y2d,z2d,u2d,v2d,w2d,d2d1,T2d1,d2d2,T2d2
    real*8,allocatable,dimension(:,:,:):: di2d1, di2d2, f3d
    real*8,parameter:: R0= 8.314d0 , PI=3.14159265358979d0               
    real*8:: epsl, rdist, omeg, fai(3) , r1, r2, fr, seta
   

   open(88,file="Reaction.in")
   read(88,*)
   read(88,*)
   read(88,*) N_spec
   close(88)

   allocate(ci1(N_spec),ci2(N_spec), Mi(N_spec),ci0(N_spec))
	  
   open(99,file="specie.in")
    read(99,*)
	  read(99,*)
    do m=1,N_spec
    do k=1,4
	  read(99,*)
	  enddo
	  read(99,*) Mi(m)
	enddo

!--------------------------------------------------------------------
! x<xh:  inlet ;  x> xh initial high temperature zone
    open(80,file="inlet.in")
    read(80,*)
    read(80,*)
	  read(80,*)  nx, ny, nz, rd, rbuf 
	  read(80,*)
	  read(80,*)  u1,T1,P1, (ci1(k),k=1,N_spec)
	  read(80,*) 
	  read(80,*)  u2,T2,P2, (ci2(k),k=1,N_spec)
    read(80,*) 
    read(80,*)  xh, Th            ! initial T=Th when x> xh and r< rh
    read(80,*)
    read(80,*)  epsl, rdist, omeg
    close(80)
!----------------------------------------------------------------
      print*, "nx,ny,nz=", nx,ny,nz
  	  print*, "N_spec=", N_spec

   allocate(x1d(nx),y1d(ny),z1d(nz))
   allocate(y2d(ny,nz),z2d(ny,nz))
   allocate(f3d(nx,ny,nz))

   open(99,file="OCFD-grid.dat",form="unformatted")
   read(99) x1d 
   read(99) y1d 
   read(99) z1d 
   close(99)
   print*, "read mesh ok"
  !-----find i0  (at x=xh) 
    do i=1,nx
       if(x1d(i) .gt. xh) exit
    enddo
    i0=i
    print*, "i0, x(i0)=", i0, x1d(i0)
   !---------------------------
   do k=1,nz
   do j=1,ny
    y2d(j,k)=y1d(j)
	z2d(j,k)=z1d(k)
   enddo
   enddo

   
   allocate(u2d(ny,nz),v2d(ny,nz),w2d(ny,nz),d2d1(ny,nz),T2d1(ny,nz),d2d2(ny,nz),T2d2(ny,nz) )
   allocate(di2d1(ny,nz,N_spec),di2d2(ny,nz,N_spec) )

    do k=1,3
     call random_number(fai(k) )
    enddo
    print*, "Fai=", fai
   
   open(99,file="test-seta.dat")
    do k=1, 360
     seta=(k-1.d0)/360.d0* 2.d0*PI
     fr=  sin(2.d0*PI*(omeg*seta+ fai(1))) & 
        + 0.5d0*sin(2.d0*PI*(2.d0*omeg*seta+ fai(2))) +  0.25d0*sin(2.d0*PI*(4.d0*omeg*seta+ fai(3)))     
     write(99,*) seta, fr
   enddo
   close(99)
   
!---------------------------------
!---inlet condition (x< xh)

   do k=1,nz
   do j=1,ny
    rr=sqrt(y2d(j,k)**2+z2d(j,k)**2)
    if(rr< rd-0.5*rbuf ) then
       T0=T1
       p0=p1
       u0=u1
       ci0(:)=ci1(:)
    else if(rr > rd + 0.5*rbuf) then
      T0=T2
      p0=p2
      u0=u2
      ci0(:)=ci2(:)
   else
     tmp=(rr-(rd-0.5d0*rbuf))/rbuf
     T0=T1+ (T2-T1) * tmp
     p0=p1+ (p2-p1)  * tmp
     u0 =u1+(u2-u1) * tmp
     ci0(:)=ci1(:)+(ci2(:)-ci1(:))*tmp
  endif
    
     
    
     d0=0.d0
     do m=1,N_spec
       di2d1(j,k,m)= p0*ci0(m)/(R0/Mi(m)*T0)
	   d0=d0+di2d1(j,k,m)
     enddo
   
     d2d1(j,k)=d0
     T2d1(j,k)= T0
     u2d(j,k)=u0
     v2d(j,k)=0.d0
     w2d(j,k)=0.d0
  
  !------ disturbance -----------------------------
     seta=acos(y2d(j,k)/rr) 
     if(z2d(j,k) .lt. 0 ) seta= 2.d0*PI- seta
     r1=rd-0.5*rdist ; r2=rd+0.5*rdist
     
	 if(rr .ge. r1  .and. rr .le. r2) then
      fr=  u0*epsl/1.75d0*   sin(PI*(rr-r1)/(r2-r1))  * ( sin(2.d0*PI*(omeg*seta+ fai(1))) & 
           + 0.5d0*sin(2.d0*PI*(2.d0*omeg*seta+ fai(2))) +  0.25d0*sin(2.d0*PI*(4.d0*omeg*seta+ fai(3)))     )
      v2d(j,k)= fr* y2d(j,k) /rr
      w2d(j,k)= fr* z2d(j,k) /rr  
      u2d(j,k)=u0+fr
   endif
    
  enddo
  enddo

    
!---------------------------------------------
!---inlet condition (x>xh)
   do k=1,nz
   do j=1,ny
    rr=sqrt(y2d(j,k)**2+z2d(j,k)**2)
    if(rr< rd-0.5*rbuf ) then
      T0=Th
       p0=p1
       u0=u1
       ci0(:)=ci1(:)
    else if(rr > rd + 0.5*rbuf) then
      T0=T2
      p0=p2
      u0=u2
      ci0(:)=ci2(:)
   else
     tmp=(rr-(rd-0.5d0*rbuf))/rbuf
     T0=Th+ (T2-Th) * tmp
     p0=p1+ (p2-p1)  * tmp
     u0 =u1+(u2-u1) * tmp
     ci0(:)=ci1(:)+(ci2(:)-ci1(:))*tmp
  endif
    
     d0=0.d0
     do m=1,N_spec
     di2d2(j,k,m)= p0*ci0(m)/(R0/Mi(m)*T0)
!	  d0=d0+di2d1(j,k,m)              ! Bug removed
	  d0=d0+di2d2(j,k,m)
    enddo
   
     d2d2(j,k)=d0
     T2d2(j,k)= T0
  enddo
  enddo

   open(100,file="flow2d-inlet-comb.dat",form="unformatted")
   write(100) d2d1 
   write(100) u2d 
   write(100) v2d 
   write(100) w2d 
   write(100) T2d1 
   do m=1,N_SPEC 
   write(100) di2d1(:,:,m)
   enddo 
   close(100)

   open(100,file="flow2d-inlet.dat")
   write(100,*) "variables=y,z,d,u,v,w,T, d1,d2,d3,d4,d5,d6,d7,d8,d9"
   write(100,*) "zone i= ", ny, " j= ", nz
   do k=1,nz
   do j=1,ny
    write(100,"(20F20.10)") y2d(j,k), z2d(j,k) , d2d1(j,k), u2d(j,k), v2d(j,k), w2d(j,k), T2d1(j,k), &
	          (di2d1(j,k,m), m=1,N_spec) 
   enddo
   enddo
   close(100)
   
   !-----------------------------------------------------
   print*,  "write 3d init file: opencfd-comb0.dat ..."
    
   open(100,file="opencfd-comb0.dat",form="unformatted")
   write(100) 0, 0.d0
   
    call write2(100,nx,ny,nz,i0,d2d1, d2d2,  f3d)
    call write2(100,nx,ny,nz,i0,u2d, u2d,     f3d)
    call write2(100,nx,ny,nz,i0,v2d, v2d,     f3d)
    call write2(100,nx,ny,nz,i0,w2d, w2d,   f3d)
    call write2(100,nx,ny,nz,i0,T2d1, T2d2, f3d)
 
    do m=1,N_spec
    call write2(100,nx,ny,nz,i0, di2d1(1,1,m),di2d2(1,1,m),f3d)
    enddo
   close(100)

   
   
    deallocate(x1d,y1d,z1d,y2d,z2d, u2d,v2d,w2d,d2d1,T2d1,d2d2,t2d2, di2d1,di2d2, f3d)
  end

   
   
   
   
   
   subroutine read3d(no,nx,ny,nz,f)
   implicit none
   integer:: no,nx,ny,nz,k
   real*8:: f(nx,ny,nz)
   do k=1,nz
     read(no) f(:,:,k)
   enddo
   end
    
	

! extent 2d  to 3d  
   subroutine write2(no,nx,ny,nz,i0, f2d1,f2d2, f3d)
   implicit none
   integer:: no,nx,ny,nz,i,j,k,i0
   real*8:: f2d1(ny,nz),f2d2(ny,nz), f3d(nx,ny,nz)
   print*, "write 3d file ..."
      
    do k=1,nz
    do j=1,ny
      do i=1,i0
       f3d(i,j,k)=f2d1(j,k)
      enddo
      do i=i0+1,nx
       f3d(i,j,k)=f2d2(j,k)
      enddo
    enddo
   enddo
   
   
   do k=1,nz
     write(no) f3d(:,:,k)
   enddo
   end
    
   
   
   
