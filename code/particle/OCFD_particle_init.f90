! 颗粒计算所使用的数据。  与流场计算的数据尽量分开，尽量不改变流场计算程序 
! 本版本仅支持直角网格
! 2024-3-6:  颗粒-流体 双向耦合
! 2024-8-14: ver 3.1a, Num_particle_Max, Np_del_Max, NP_send_Max 值可由用户设定
!---------flow data---------------------------------------------------------
module particle_data
 use OCFD_precision 
 implicit none
   TYPE particle  ! 颗粒信息
    real(kind=OCFD_REAL_KIND):: idp, dp, rhop,xp,yp,zp,up,vp,wp    ! 颗粒全局编号，颗粒直径，颗粒密度，颗粒位置， 颗粒速度
    integer::  ip,jp,kp  !  颗粒所在网格位置 （所在网格的i-,j-,k- 信息）  
   end TYPE  

 real(kind=OCFD_REAL_KIND),allocatable,dimension(:):: x1d,y1d,z1d   ! Coordinates  (1d)
 real(kind=OCFD_REAL_KIND):: xmin,xmax,ymin,ymax,zmin,zmax    ! 本MPI进程（局部）计算域的边界
 real(kind=OCFD_REAL_KIND):: xmin0,xmax0,ymin0,ymax0,zmin0,zmax0   !全局计算域边界
! integer,parameter:: Num_particle_Max=10000000, Np_del_Max=100000, NP_send_Max=10000    ! 每个区域（进程）最大的粒子数目； 流出粒子的最大数目
 integer:: Num_particle_Max=10000000, Np_del_Max=200000, Np_send_Max=20000    ! 每个区域（进程）最大的粒子数目； 流出粒子的最大数目
 integer,parameter:: NVar_P=9   ! 需要传递的变量数（9个： idp, d,rhop, Xp,Yp,Zp,Up,Vp,Wp  )
 integer:: Num_particle, Np_del   !粒子数目，流出计算域的粒子数目
 integer:: Nid_send(26), Nid_recv(26)       ! 26个相邻区域对应的MPI id （按顺序 发送及接收） 
 TYPE(particle),dimension(:),pointer:: ptcals	 ! 颗粒信息	   
 integer,dimension(:),pointer:: idp_del(:)   ! 流出粒子的（进程内）编号
!---------------------------------------------------- 
 integer,parameter:: OCFD_ANA_particle_Zplane=1001
  TYPE PARA_particle 
  integer:: Coupling_particle                        ! 耦合方式 （1 单向耦合， 2 双向耦合 ， 3,4 四向耦合）
! 颗粒对流体作用力反馈方式 （1 1阶精度插值， 2 2阶精度插值）, 颗粒初始速度（0 强制为0 1 跟随流体，2指定）  
  integer:: Force_feedback , Initial_velocity_Particle     
  integer:: IBC_particle,ANA_Number_particle 
  Real(kind=OCFD_REAL_KIND):: BC_Para_particle(100), ANA_Para_particle(100,10)
  end TYPE 
  TYPE(PARA_particle):: Para_P
end 

subroutine init_particle
 use flow_data
 use particle_data
 implicit none 
 call read_parameter_particle
 call allocate_mem_particle
 call Neighber_particle_MPI
 call init_particle_mesh
 call read_particle_data 
 if(Para_P%Initial_velocity_Particle .ne. 2) then 
   call set_particle_velocity           ! 设定颗粒初始速度 
 endif 
 end 
 

 subroutine allocate_mem_particle
 use flow_data
 use particle_data
 implicit none 
  allocate(x1d(1-LAP:nx+LAP),y1d(1-LAP:ny+LAP),z1d(1-LAP:nz+LAP))
  allocate(ptcals(NUM_particle_Max)) !颗粒信息
  allocate(idp_del(Np_del_Max))      !流出本计算域颗粒的编号
 end 


 subroutine read_particle_data
 use flow_data
 use particle_data
 implicit none 
 integer:: Num_particle0,i,j,k,k0,NP_read,nk,Npk,ierr
 real(kind=OCFD_REAL_KIND),allocatable:: pdata(:,:) 
 real(kind=OCFD_REAL_KIND):: xp,yp,zp
 integer,parameter:: NP_read_onetime=1000000     ! 每次读取颗粒的最大数目 （避免一次读取大量颗粒信息，内存无法承受）
 
!----------网格信息---------------------------
! 仅适用于直角网格  

! read data file of particles  
  if(my_id==0) then
   open(99,file="ocfd-particle.dat",form="unformatted")
   read(99) Num_particle0
   print*, "read ocfd-particle.dat, Num_particle=", Num_particle0
  endif
  call MPI_bcast(Num_particle0,1,MPI_INTEGER,0, MPI_COMM_WORLD,ierr) 
  allocate(pdata(NVar_P,NP_read_onetime))

!---------------------------------
  k0=0 
  Npk=int(Num_particle0/Np_read_onetime)+1
 do nk=1,Npk
  if(nk<Npk) then 
   Np_read=Np_read_onetime           ! 每次读取的数目
  else 
   Np_read= Num_particle0 - (Npk-1)* Np_read_onetime       ! 剩余的数目
  endif
  
  if(my_id ==0) then 
  do k=1, Np_read
  read(99) (pdata(j,k),j=1,NVar_P)  !idp, dp, rhop,Xp, YP, Zp, Up, Vp, Wp
  enddo 
  endif 
  
  call MPI_bcast(pdata,Np_read*NVar_P,OCFD_DATA_TYPE,0, MPI_COMM_WORLD,ierr) 
  do k=1,Np_read
   xp=pdata(4,k); yp=pdata(5,k); zp=pdata(6,k)
   ! 判断是否在本进程的计算域内
   if(xp>= xmin .and. xp < xmax .and. yp >= ymin .and. yp < ymax .and. zp>= zmin .and. zp < zmax) then  
   k0=k0+1
   call pmsg_cpout(ptcals(k0),pdata(1,k))
   call comput_ijk_particle(ptcals(k0)%xp, ptcals(k0)%yp,ptcals(k0)%zp, ptcals(k0)%ip, ptcals(k0)%jp,ptcals(k0)%kp)
   endif
  enddo
  enddo 

  Num_particle=k0       ! 本计算域内的粒子总数 
  
  if(my_id .eq. 0) then 
  close (99)  
  print*, "read ocfd-particle.dat OK"
  endif 
  
  deallocate(pdata)
   
 end 

!----------目前版本仅支持直角网格-------------------------- 
 subroutine init_particle_mesh
 use flow_data
 use particle_data
 implicit none 
 integer:: i,j,k,ierr
 real(kind=OCFD_REAL_KIND),dimension(3):: xmin3,xmin30,xmax3,xmax30
  do i=1-LAP,nx+LAP
  x1d(i)= Axx(i,1,1)
  enddo 
  do j=1-LAP,ny+LAP
  y1d(j)=Ayy(1,j,1)
  enddo 
  do k=1-LAP,nz+LAP 
  z1d(k)=Azz(1,1,k) 
  enddo 
  
 ! 计算域的边界坐标 
  xmin=x1d(1)
  ymin=y1d(1)
  zmin=z1d(1)
  if(npx ==npx0-1 .and. Para%Iperiodic_X ==0) then 
   xmax=x1d(nx)
  else 
   xmax=x1d(nx+1)
  endif 
  
  if(npy ==npy0-1 .and. Para%Iperiodic_Y ==0) then 
   ymax=y1d(ny)
  else 
   ymax=y1d(ny+1)
  endif 
  
  if(npz ==npz0-1 .and. Para%Iperiodic_Z ==0) then 
   zmax=z1d(nz)
  else 
   zmax=z1d(nz+1)
  endif 
  xmin3(1)=xmin; xmin3(2)=ymin; xmin3(3)=zmin
  xmax3(1)=xmax; xmax3(2)=ymax; xmax3(3)=zmax 
  call MPI_ALLREDUCE(xmin3,xmin30,3,OCFD_DATA_TYPE,MPI_MIN,MPI_Comm_World,ierr)  
  call MPI_ALLREDUCE(xmax3,xmax30,3,OCFD_DATA_TYPE,MPI_MAX,MPI_Comm_World,ierr)  
  xmin0=xmin30(1); ymin0=xmin30(2); zmin0=xmin30(3)     ! 全局计算域边界
  xmax0=xmax30(1); ymax0=xmax30(2); zmax0=xmax30(3)
  end 
  

!========================
 subroutine comput_ijk_particle(xp,yp,zp,ip,jp,kp)         ! 计算粒子所在的网格单元 （搜索算法，计算量较大，用于初值设置）
 use flow_data
 use particle_data
 implicit none 
 real(kind=OCFD_REAL_KIND):: xp,yp,zp
 integer:: ip,jp,kp,i
 if(xp< xmin .or. xp> xmax .or. yp< ymin .or. yp> ymax .or. zp<zmin .or. zp> zmax) then 
   print*, "A pratical is out of the domain !", xp,yp,zp, my_id 
   stop
 endif
 ip=nx; jp=ny; kp=nz  
 do i=1,nx 
 if(xp>= x1d(i) .and. xp < x1d(i+1)) then 
  ip=i
  exit
 endif 
 enddo 
 do i=1,ny 
 if(yp>= y1d(i) .and. yp < y1d(i+1)) then 
  jp=i
  exit
 endif 
 enddo 
 do i=1,nz 
  if(zp>= z1d(i) .and. zp < z1d(i+1)) then 
  kp=i
  exit
 endif
 enddo 
 end 
 
!------本进程（my_id)相邻的26个区域对应的MPI_id (注意发送和接收对应的次序不同）--------- 
! Nid_send(1:26) 26个相邻进程的id号 （发送消息时使用）
! Nid_recv(1:26) 26个相邻进程的id号 （接收消息时使用，次序于发送消息相对应）
!=================================== 
 subroutine Neighber_particle_MPI 
 use flow_data
 use particle_data
 implicit none  
 integer:: ks,i1,j1,k1,npx1,npy1,npz1,npx2,npy2,npz2
 do ks=1,26
  k1=int(ks/9)
  j1=int((ks-k1*9)/3)
  i1=mod(ks,3)
! 获得相邻进程的1维id； npx1: i1=0,1,2 对应 中、左、右侧进程； npx2: i1=0,1,2 对应 中、右、左侧进程 
  call comput_np1(npx1,i1,npx,npx0,Para%IPeriodic_X,1)
  call comput_np1(npy1,j1,npy,npy0,Para%IPeriodic_Y,1)
  call comput_np1(npz1,k1,npz,npz0,Para%IPeriodic_Z,1)
  if(npx1<0 .or. npy1 <0 .or. npz1 <0) then 
    Nid_send(ks)=MPI_PROC_NULL
  else 
    Nid_send(ks)=npx1+npy1*npx0+npz1*npx0*npy0 
  endif 
  
  call comput_np1(npx2,i1,npx,npx0,Para%IPeriodic_X,-1)
  call comput_np1(npy2,j1,npy,npy0,Para%IPeriodic_Y,-1)
  call comput_np1(npz2,k1,npz,npz0,Para%IPeriodic_Z,-1) 
  
  if(npx2<0 .or. npy2 <0 .or. npz2 <0) then 
    Nid_recv(ks)=MPI_PROC_NULL
  else 
    Nid_recv(ks)=npx2+npy2*npx0+npz2*npx0*npy0 
  endif
 enddo 
 end 

  subroutine comput_np1(np1,i1,np,np0,IPeriodic,Iflag)
  implicit none 
  integer:: np1,i1,np,np0,IPeriodic,Iflag,npx_nb
  if(i1==0) then 
   np1=np
  else if(i1==1) then 
   np1= npx_nb(np-Iflag,np0,IPeriodic)
  else 
   np1=npx_nb(np+Iflag,np0,IPeriodic)
  endif 
  end 
  
 ! 计算计算域的npx值， 考虑周期性  
  function npx_nb(npx1,npx0,IPeriodic)
  implicit none 
  integer:: npx_nb,npx1,npx0,Iperiodic
  npx_nb=npx1 
  if(npx1 < 0) then 
   if(Iperiodic == 1) then 
   npx_nb=npx1+npx0 
   else 
   npx_nb=-1
   endif 
  endif 
  
  if(npx1> npx0-1) then 
   if(Iperiodic==1) then 
   npx_nb=npx1-npx0 
   else 
   npx_nb=-1
   endif 
  endif 
  end  
 
 
 subroutine save_particle_data
 use flow_data
 use particle_data
 implicit none 
 integer:: j,k,kp,ks,num_particle0,NP,ierr,Status(MPI_status_Size)
 real(kind=OCFD_REAL_KIND),allocatable:: pdata(:,:)
 real(kind=OCFD_REAL_KIND):: pdata1(NVar_P)
 integer,allocatable:: Npk(:)
 character(len=100) filename1
 
 allocate(Npk(0:np_size-1))       ! 每个进程的粒子数
 call MPI_Gather(Num_particle,1,MPI_INTEGER,Npk,1,MPI_INTEGER,0,MPI_COMM_WORLD)  !得到每个进程的粒子数
 
 if(my_id .eq. 0) then 
  num_particle0=0        !计算粒子总数
  do k=0,np_size-1 
  num_particle0=num_particle0+Npk(k)
  enddo 
  
  write(filename1,"('ocfd-particle'I8.8'.dat')") Istep  
  print*, "write particle data ", filename1 
  print*, "Total number of particles is ", num_particle0
  open(88,file=filename1,form="unformatted")
  write(88) num_particle0 
! 将0进程的粒子信息写入数据文件  
  do ks=1,Num_particle
   call pmsg_cpin(ptcals(ks),pdata1)
   write(88) (pdata1(j),j=1,NVar_P)
  enddo 
 endif  

 if(my_id .ne. 0) then 
   allocate(pdata(NVar_P,Num_particle))
   do ks=1,Num_particle
    call pmsg_cpin(ptcals(ks),pdata(1,ks))
   enddo 
!   call MPI_BSend(pdata,NVar_P*Num_particle,OCFD_DATA_TYpe,0,999,MPI_COMM_WORLD,ierr)
    call MPI_Send(pdata,NVar_P*Num_particle,OCFD_DATA_TYpe,0,999,MPI_COMM_WORLD,ierr)    ! 2024-8-28
   deallocate(pdata)
 else  
  do kp=1,np_size-1   
   Np=Npk(kp) 
   allocate(pdata(NVar_P,NP))
   call MPI_Recv(pdata,NVar_P*Np,OCFD_DATA_TYPE,kp,999,MPI_COMM_WORLD,status,ierr)
   do ks=1,Np 
    write(88) (pdata(j,ks),j=1,NVar_P)
   enddo 
  deallocate(pdata)
 enddo  
 endif 
 
 if(my_id .eq. 0) then 
  print*, "Write particle data OK "
  close(88)
 endif 
 deallocate(Npk) 
 end 
 
 ! -------------------------------------------------------------------------------   
! read flow parameters  (Namelist type)      
	 subroutine read_parameter_particle
      use flow_para 
	  use particle_data
	  implicit none 
	  integer:: IBC_particle,ANA_Number_particle,Coupling_particle,Force_feedback,Initial_velocity_Particle,ierr
	  real(kind=OCFD_REAL_KIND):: BC_para_particle(100), ANA_Para_particle(100,10)
	  namelist /control_particle/ IBC_particle,ANA_Number_particle, &
			BC_para_particle, ANA_Para_particle,Coupling_particle,Force_feedback,Initial_velocity_Particle, &
			Num_particle_Max,Np_del_Max,Np_send_Max
			
	  integer:: nparameters(100)
	  real(kind=OCFD_REAL_KIND):: rparameters(100)
!--------Set defalut parameter ----------------
!--------particle para------	  
	  IBC_particle=0
	  ANA_Number_particle=0
	  BC_para_particle(:)=0.d0
	  ANA_Para_particle(:,:)=0.d0
	  Coupling_particle=2         ! 默认双向耦合
	  Force_feedback=2          ! 反馈力默认2阶精度 
	  Initial_velocity_Particle=1  ! 默认颗粒初始速度与当地流场相同
	  Num_particle_Max=10000000
	  Np_del_Max=200000
	  Np_send_Max=20000
!----------------------------------------------	  
	  
    if(my_id .eq. 0)  then
 	   open(99,file="opencfd-comb2.in")
	   read(99,nml=control_particle)
       close(99)
	  nparameters(1)=IBC_particle
	  nparameters(2)=ANA_Number_particle
	  nparameters(3)=Coupling_particle
	  nparameters(4)=Force_feedback
	  nparameters(5)=Initial_velocity_Particle
	  nparameters(6)=Num_particle_Max
	  nparameters(7)=Np_del_Max
	  nparameters(8)=Np_send_Max
	 endif   
  
	  call MPI_bcast(nparameters(1),100,MPI_INTEGER,0,  MPI_COMM_WORLD,ierr)  
	  call MPI_bcast(rparameters(1),100,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	
	 Para_P%IBC_particle= nparameters(1)
	 Para_P%ANA_Number_particle=nparameters(2)
	 Para_P%Coupling_particle=nparameters(3)
	 Para_P%Force_feedback=nparameters(4)
	 Para_P%Initial_velocity_Particle=nparameters(5) 
	 Num_particle_Max=nparameters(6)
	 Np_del_Max=nparameters(7)
	 Np_send_Max=nparameters(8)
	 
	 
     Para_P%BC_Para_particle(:)=BC_Para_particle(:)
     Para_P%ANA_Para_particle(:,:)=ANA_Para_particle(:,:)
     call MPI_bcast(Para_P%BC_Para_particle(1),100,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)	
	 call MPI_bcast(Para_P%ANA_Para_particle(1,1),1000,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)		 
	 
	 end 

 
    subroutine set_particle_velocity
    use flow_data
    use particle_data
    implicit none 
    integer:: k
	real(kind=OCFD_REAL_KIND):: d0,u0,v0,w0,mu0
    TYPE(particle),pointer:: Pt
    do k=1, Num_particle
    Pt=>ptcals(k)
    if(Para_P%Initial_velocity_Particle ==0) then 
     Pt%up=0.d0 
	 Pt%vp=0.d0 
	 Pt%wp=0.d0 
	else 
     call comput_u_particle(Pt%xp,Pt%yp,Pt%zp,Pt%ip,Pt%jp,Pt%kp,d0,u0,v0,w0,mu0)  ! 计算颗粒位置处的流体速度
     Pt%up=u0
	 Pt%vp=v0 
	 Pt%wp=w0 
	endif 
	enddo 
	end 
 
 
 
 
 
 
 
 
 
 
 