! 进程间重新分布颗粒（删除流出颗粒，录入流入颗粒）

 subroutine particle_redistribute 
 use flow_data
 use particle_data
 implicit none  
 integer:: i,j,k,i1,j1,k1,ks,id_send,id_recv,kptr,Np_in_total,kdel
 TYPE(particle),pointer:: Pt
 real(kind=OCFD_REAL_KIND),allocatable:: pdata_send(:,:,:),pdata_recv(:,:,:)
 integer:: Np_out(26),Np_in(26) ! 本区域周边的26个区域：流出粒子的数目，流入粒子的数目
 integer:: Status(MPI_status_Size),ierr

 allocate(pdata_send(NVar_P,NP_send_Max,26), pdata_recv(NVar_P,NP_send_Max,26))
 
!====================================================================================  
! 记录流出计算域的粒子，将编号记录在idp_del(:)数组中 
 kdel=0
  do k=1, Num_particle
  Pt=>ptcals(k)
  if(Pt%xp<xmin .or. Pt%xp>=xmax .or. Pt%yp<ymin .or. Pt%yp>=ymax .or. Pt%zp<zmin .or. Pt%zp>=zmax) then  
 ! 颗粒流出（本进程）计算域外
  kdel=kdel+1
  if(kdel > Np_del_Max ) print*, "Warning !!! kdel > Np_del_max", "my_id=", my_id  ! 出错信息
  idp_del(kdel)=k    !流出粒子的（本地）编号
  endif 
  enddo 
  Np_del=kdel  ! 流出粒子的数目

 !------------------------------------------
 !将流出的粒子分组 （26组），对应周围26个区域  （每个区域与周围26个区域相邻）
 Np_out(:)=0    ! 流出到（26个区域中）每个区域的粒子数
 Np_in(:)=0     ! 从周围26个区域的每个区域流入本进程的粒子数  
 
 do k=1,Np_del
  i1=0; j1=0; k1=0
  Pt=>ptcals(idp_del(k))
  if(Pt%xp<xmin) i1=1
  if(Pt%xp>=xmax) i1=2
  if(Pt%yp<ymin) j1=1
  if(pt%yp>=ymax) j1=2
  if(pt%zp<zmin) k1=1
  if(pt%zp>=zmax) k1=2
  ks=i1+j1*3+k1*9        ! 相邻计算域的编号
  Np_out(ks)=Np_out(ks)+1      ! 流出的粒子数
  call ptcal_perodic(ptcals(idp_del(k)))   ! 根据周期性调整粒子坐标 （穿越周期边界需要+/-周期）
  call pmsg_cpin( ptcals(idp_del(k)), pdata_send(1,Np_out(ks),ks))  ! copy into pdata_send 
  Pt%idp=-1.d0           !做标记（流出的粒子，需要删掉）
  enddo 
  
! 进程间交换信息，获得流入本计算域粒子的数目信息 
! 每个进程向周围26个进程 发送并接收 信息
  
! 交换流入/流出颗粒的数目------
  Np_in_total=0       ! 流入粒子的总数
  do ks=1,26
  call mpi_sendrecv(Np_out(ks),1,MPI_INTEGER,Nid_send(ks),8888,  &
        Np_in(ks),1,MPI_INTEGER,Nid_recv(ks),8888,MPI_COMM_WORLD,Status,ierr)
  Np_in_total=Np_in_total+Np_in(ks)
  enddo   
!--交换流入/流出的颗粒信息------ 
  do ks=1,26
  id_send=Nid_send(ks) 
  id_recv=Nid_recv(ks) 
  if(Np_out(ks) ==0) id_send= MPI_PROC_NULL   ! 如果流出颗粒数目为0， 则不发送数据
  if(Np_in(ks)==0)   id_recv= MPI_PROC_NULL
  call mpi_sendrecv(pdata_send(1,1,ks),Np_out(ks)*nVar_P, OCFD_DATA_TYPE ,id_send,9999,  &
        pdata_recv(1,1,ks),Np_in(ks)*nVar_P,OCFD_DATA_TYPE,id_recv,9999,MPI_COMM_WORLD,Status,ierr)
  enddo 
 
!------压缩主数据ptcals, 将流出的粒子信息删除，用队尾的粒子信息填充 
 kptr=Num_particle  ! 指针，指向队尾
 kdel=0  ! 处理过的流出粒子
 if(Num_particle-Np_del >= 1) then 
 do k1=1, Np_del
   do while( kptr>=1 )
    if(ptcals(kptr)%idp >=0 ) then 
      exit 
    else  
     kptr=kptr-1
     kdel=kdel+1
    endif 
   enddo 
   if(kdel >=Np_del) exit 
   ptcals(idp_del(k1))=ptcals(kptr)
   kptr=kptr-1
   kdel=kdel+1
   if(kdel >=Np_del) exit 
  enddo 
 else 
 kptr=0         ! 全部颗粒都流出计算域
 endif   
 
!---------将接收到的粒子信息装入主数据ptcals  
 do ks=1,26
 do k1=1,Np_in(ks)
 kptr=kptr+1
 call pmsg_cpout( ptcals(kptr), pdata_recv(1,k1,ks))
 Pt=>ptcals(kptr) 
 call comput_ijk_particle(Pt%xp,Pt%yp,Pt%zp,Pt%ip,Pt%jp,Pt%kp)        ! 搜索，建立颗粒所在的网格信息
 enddo 
 enddo 
 
  Num_particle=Num_particle-Np_del+Np_in_total  ! 更新后本计算域内的粒子总数     
 ! debug message  
   if(Num_particle .ne. kptr) then 
    print*, "Warning , may be error in particle_redistribute "
	print*, "Num_particle, kptr,my_id=", Num_particle, kptr,my_id
   endif 
!------------------------------------------------------------------------------
  
 deallocate(pdata_send,pdata_recv)
 end 
 
 
 !-----------------------------------------
 ! 处理周期边界条件： 穿越周期边界调整坐标
 subroutine  ptcal_perodic (Pt)
 use flow_data
 use particle_data
 implicit none  
 TYPE(particle):: Pt 
 if(Para%Iperiodic_X == 1) then 
  if(Pt%xp <xmin0 ) pt%xp=Pt%xp+Para%Periodic_ISpan(1) 
  if(Pt%xp >=xmax0 ) pt%xp=Pt%xp-Para%Periodic_ISpan(1)
 endif 
 if(Para%Iperiodic_Y == 1) then 
  if(Pt%yp <ymin0 ) pt%yp=Pt%yp+Para%Periodic_JSpan(2) 
  if(Pt%yp >=ymax0 ) pt%yp=Pt%yp-Para%Periodic_JSpan(2)
 endif 
 if(Para%Iperiodic_Z == 1) then 
  if(Pt%zp <zmin0 ) pt%zp=Pt%zp+Para%Periodic_KSpan(3) 
  if(Pt%zp >=zmax0 ) pt%zp=Pt%zp-Para%Periodic_KSpan(3)
 endif  
 end 
 
! copy into pdata 
subroutine pmsg_cpin(ptc,pdata)
 use particle_data
 implicit none 
 TYPE(particle):: ptc 
 real(kind=OCFD_REAL_KIND):: pdata(nVar_P)
 pdata(1)=ptc%idp 
 pdata(2)=ptc%dp
 pdata(3)=ptc%rhop 
 pdata(4)=ptc%xp 
 pdata(5)=ptc%yp
 pdata(6)=ptc%zp
 pdata(7)=ptc%up  
 pdata(8)=ptc%vp
 pdata(9)=ptc%wp
 end 
 
subroutine pmsg_cpout(ptc,pdata)
 use particle_data
 implicit none 
 TYPE(particle):: ptc 
 real(kind=OCFD_REAL_KIND):: pdata(nVar_P)
 ptc%idp =pdata(1) 
 ptc%dp  =pdata(2)
 ptc%rhop=pdata(3) 
 ptc%xp  =pdata(4) 
 ptc%yp  =pdata(5)
 ptc%zp  =pdata(6)
 ptc%up  =pdata(7)  
 ptc%vp  =pdata(8)
 ptc%wp  =pdata(9)
 end  
 
 