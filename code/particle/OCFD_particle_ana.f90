! code for particles, Code by Li Xinliang, lixl@imech.ac.cn------------------
! 后处理程序
! 本版本仅支持直角网格
!------------------------------------------------------------------
subroutine Ana_particle
 use flow_data
 use particle_data
 implicit none 
 integer m,Kana,Kstep_ana
 real(kind=OCFD_REAL_KIND):: ana_data(50)
 
  
  do m=1, Para_P%ANA_Number_particle  
   Kana=nint(Para_P%ANA_Para_particle(1,m))
   Kstep_ana=nint(Para_P%ANA_Para_particle(2,m))
   ana_data(:)=Para_P%ANA_Para_particle(3:52,m) 

 
 if(mod(Istep,Kstep_ana) .eq. 0) then 
   select case(Kana)    
    case (OCFD_ANA_particle_Zplane)
	 call ana_particle_Zplane(ana_data)
 	case default 
	 if(my_id .eq. 0) print*, "This analysis code is not supported!"
   end select 
 endif 
 enddo
 end 

!---------save particle in Zplane in the range  z1<z<z2 
subroutine Ana_particle_Zplane(ana_data)
 use flow_data
 use particle_data
 implicit none 
 real(kind=OCFD_REAL_KIND):: ana_data(50),z1,z2
 real(kind=OCFD_REAL_KIND),allocatable:: pdata(:,:)
! integer,parameter:: num_max_P=100000
 integer,parameter:: num_max_P=500000     ! 2024-8-26

 integer:: i,j,k,ks,ks0,kp,Np,ierr,Status(MPI_status_Size)
 integer,allocatable:: Npk(:)        !各进程的粒子数
 
 z1=ana_data(1)
 z2=ana_data(2)
 allocate(pdata(NVar_P,num_max_P),Npk(0:np_size-1))
 
 if(my_id .eq. 0) then
 print*, "wirte Zparticle_time.dat ..."
 open(99,file="Zparticle_time.dat",form="unformatted",position="append") 
 write(99) tt 
 endif 
 
  ks=0
  do k=1,Num_particle
  if(ptcals(k)%zp >=z1 .and. ptcals(k)%zp <=z2) then  
  ks=ks+1
  if(ks>num_max_P) then 
   print*, "Warning ! in Ana_particle_Zplane, ks>num_max_P, my_id=  ",my_id
   exit 
  endif 
  call pmsg_cpin(ptcals(k),pdata(1,ks))
  endif 
  enddo 
!==========================
  call MPI_Gather(ks,1,MPI_INTEGER,Npk,1,MPI_INTEGER,0,MPI_COMM_WORLD)  !得到每个进程的粒子数  
  ks0=sum(Npk)         ! 粒子数总和
  
  if(my_id .eq. 0) then 
  write(99) ks0         ! 颗粒总数 
  do k=1,ks             ! 0号进程的粒子
  write(99) (pdata(j,k),j=1,NVar_P)        
  enddo 
  endif 
  
  if(my_id .ne. 0) then 
   if(ks .ne. 0) then 
!    call MPI_BSend(pdata,NVar_P*ks,OCFD_DATA_TYpe,0,999,MPI_COMM_WORLD,ierr)   
     call MPI_Send(pdata,NVar_P*ks,OCFD_DATA_TYpe,0,999,MPI_COMM_WORLD,ierr)  ! 2024-8-28
   endif 
  else 
  do kp=1,np_size-1   
   Np=Npk(kp) 
   if(Np .ne. 0) then 
   call MPI_Recv(pdata,NVar_P*Np,OCFD_DATA_TYPE,kp,999,MPI_COMM_WORLD,status,ierr)
   do ks=1,Np 
    write(99) (pdata(j,ks),j=1,NVar_P)
   enddo 
   endif
  enddo 
  endif 
  
  if(my_id .eq. 0) then 
   close(99)
  endif 
  deallocate(pdata,Npk)
  end 
  