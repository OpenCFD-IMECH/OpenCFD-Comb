! boundary conditions： boundary-layers  (transition/ STBLI/ ...) 
! 2024-5-12 (ver 3.1)  Wall jet boundary condition is added
!---------------------------------------------------------------------
   module BC_wall_jet_data 
   use OCFD_precision
   implicit none 
   integer:: Init_walljet=0
!-------bc_wall_jet2---------------   
   integer:: ib_jet,kb_jet,ie_jet,ke_jet, nx_jet,nz_jet    ! 射流窗口的网格位置 （MPI 当地）
   integer:: ib_jet0,kb_jet0,ie_jet0,ke_jet0,nx_jet0,nz_jet0        ! 射流窗口的位置 （全局）
   real(kind=OCFD_REAL_KIND),allocatable:: flow2d_walljet(:,:,:) 
!-------- bc_wall_jet1-----------   
   real(kind=OCFD_REAL_KIND):: Pjet, Tjet, djet, ujet,vjet,wjet,xjet,zjet,Rjet
   real*8,allocatable:: di_jet(:)   ! 射流中各组分的质量比分
   end     

!-----------------------------------------------------
!  圆孔射流
   subroutine bc_Wall_jet1      
   use flow_data
   use BC_wall_jet_data
   implicit none 
   integer:: i,k,m
   real(kind=OCFD_REAL_KIND):: rr 

   if(Init_walljet==0) then 
   call Wall_jet_init1
   Init_walljet=1
   endif 
  
  if(npy .eq. 0) then 
   do k=1,nz 
   do i=1,nx 
    rr=sqrt((Axx(i,1,k)-xjet)**2+(Azz(i,1,k)-zjet)**2)
	if(rr<Rjet) then           ! 射流区， Dirichlet边界条件 
     d(i,1,k)=djet
     u(i,1,k)=ujet
     v(i,1,k)=vjet
     w(i,1,k)=wjet
	 T(i,1,k)=Tjet	 
    do m=1,N_SPEC 
	 di(i,1,k,m)=d(i,1,k)*di_jet(m)
	enddo
   endif 
   enddo 
   enddo 
   endif 
   end 




!----------------------------------------------------------------
!  利用射流窗口数据设置射流区壁面边界条件， 可处理任意形状的壁面射流	
	subroutine bc_Wall_jet2
    use flow_data
    use BC_wall_jet_data
	implicit none 
	integer i,j,k,m,i1,k1, ierr

!---------------------------------------------------
    if(Init_walljet .eq. 0) then 
	call Wall_jet_init2
	Init_walljet=1
	endif 

!----------Wall boudnary ---------------------------------------
 if(npy .eq. 0  .and. nx_jet .ne. 0 .and. nz_jet .ne. 0)  then
   do k=1,nz_jet
   do i=1,nx_jet
    k1=k+kb_jet-1
	i1=i+ib_jet-1 
     d(i1,1,k1)=flow2d_walljet(i,k,1)
     u(i1,1,k1)=flow2d_walljet(i,k,2)
     v(i1,1,k1)=flow2d_walljet(i,k,3)
     w(i1,1,k1)=flow2d_walljet(i,k,4)
	 T(i1,1,k1)=flow2d_walljet(i,k,5)

	 if(T(i1,1,k1) <=0 )  T(i1,1,k1)=T(i1,2,k1)                                 ! adabtic boundary   绝热壁
	 if(d(i1,1,k1) <=0)   d(i1,1,k1)=d(i1,2,k1)*T(i1,2,k1)/T(i1,1,k1)           ! dp/dy=0
    
	 do m=1,N_Spec
	  di(i1,1,k1,m)=d(i,1,k)*flow2d_walljet(i,k,5+m)                    ! flow2d_walljet(:,:,5+m)  质量比分ci(:)
      if(di(i1,1,k1,m) <0 )  di(i1,1,k1,m)=di(i1,2,k1,m)*d(i1,1,k1)/d(i1,2,k1)                ! dci/dy=0
	 enddo
  enddo
  enddo
 endif
 end 
 
!--------初始化，读入壁面射流数据 （适用于单个圆孔射流）
   subroutine Wall_jet_init1
    use flow_data
    use BC_wall_jet_data
	implicit none  
	integer:: m,ierr
	real(kind=OCFD_REAL_KIND):: para_real(10)
	
	allocate(di_jet(N_SPEC))
	if(my_id .eq. 0) then 
     open(100,file="BC_walljet.dat")
	 read(100,*)
	 read(100,*) 
	 read(100,*) Pjet, Tjet, ujet,vjet,wjet,xjet,zjet,Rjet
	 read(100,*)
	 read(100,*) (di_jet(m),m=1,N_SPEC)
	 close(100)
	endif 
	para_real(1)=Pjet
	para_real(2)=Tjet
	para_real(3)=ujet
	para_real(4)=vjet
	para_real(5)=wjet
	para_real(6)=xjet
	para_real(7)=zjet
	para_real(8)=Rjet	
	call MPI_bcast( para_real,10,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 	 
	Pjet=para_real(1)
	Tjet=para_real(2)
	ujet=para_real(3)
	vjet=para_real(4)
	wjet=para_real(5)
	xjet=para_real(6)
	zjet=para_real(7)
	Rjet=para_real(8)
	call MPI_bcast( di_jet,N_SPEC,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 	

	
! 根据压力、温度以及质量比分ai(:) ，计算密度 
    call comput_d(di_jet,Tjet,Pjet,djet)
 
 end 
	
   
   
   
   
   
   
   
   
!------------初始化， 读入窗口型射流数据---------------------------   
   subroutine Wall_jet_init2
    use flow_data
    use BC_wall_jet_data
	implicit none   
	integer:: para_int(10),NV
	integer:: npx1,npx2,i1,i2,npz1,npz2,k1,k2,ierr,i,k,m,i0,k0
	real(kind=OCFD_REAL_KIND),allocatable:: flow2d_walljet0(:,:,:)      ! 壁面射流（窗口）区域二维流场  （全局）
	
	NV=5+N_SPEC   ! 变量数 d,u,v,w,T, ci

!-----------读取壁面射流数据 （壁面射流出口物理量， nx_jet*nz_jet 网格的方形窗口）	
    if(my_id .eq. 0) then 
    open(100,file="BC_walljet.dat",form="unformatted") 
	read(100) nx_jet0,nz_jet0, ib_jet0,kb_jet0
    para_int(1)=nx_jet0
	para_int(2)=nz_jet0
	para_int(3)=ib_jet0
	para_int(4)=kb_jet0
	endif 
	
    call MPI_bcast(para_int,10,MPI_INTEGER,0,  MPI_COMM_WORLD,ierr) 
     nx_jet0=para_int(1)
	 nz_jet0=para_int(2)
	 ib_jet0=para_int(3)
	 kb_jet0=para_int(4)
	
	allocate(flow2d_walljet0(nx_jet0,nz_jet0,NV))
	
	if(my_id .eq. 0) then 
	read(100) flow2d_walljet0 
	close(100)
	endif 

	call MPI_bcast(flow2d_walljet0,nx_jet0*nz_jet0*NV,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr) 
	
!------------确定（射流）窗口的	
	ie_jet0=ib_jet0+nx_jet0-1
	ke_jet0=kb_jet0+nz_jet0-1
	
	call get_i_node(ib_jet0,npx1,i1)
	if(npx < npx1 ) then 
	  ib_jet=-1              ! 不在窗口区
	else if(npx ==npx1) then 
	  ib_jet=i1 
	else 
	  ib_jet=1
	endif 
	
	call get_i_node(ie_jet0,npx2,i2)	
	if(npx > npx2)then 
	 ie_jet=-1
	else if(npx ==npx2) then 
	 ie_jet=i2 
	else 
	 ie_jet=nx 
	endif 
	
	call get_k_node(kb_jet0,npz1,k1)
	if(npz < npz1 ) then 
	  kb_jet=-1              ! 不在窗口区
	else if(npz ==npz1) then 
	  kb_jet=k1 
	else 
	  kb_jet=1
	endif 
	
	call get_k_node(ke_jet0,npz2,k2)	
	if(npz > npz2)then 
	 ke_jet=-1
	else if(npz ==npz2) then 
	 ke_jet=k2 
	else 
	 ke_jet=nz 
	endif 	
    
	if(ib_jet ==-1 .or. ie_jet==-1) then 
	 nx_jet=0
	else 
	 nx_jet=ie_jet-ib_jet+1 
	endif 
	
	if(kb_jet ==-1 .or. ke_jet==-1) then 
	 nz_jet=0
	else 
	 nz_jet=ke_jet-kb_jet+1 
	endif 	
	
  if(npy .eq. 0) then 
   if(nx_jet .ne. 0 .and. nz_jet .ne. 0) then  
   allocate(flow2d_walljet(nx_jet,nz_jet,NV))
   
   do m=1,NV 
   do k=1,nz_jet 
   do i=1,nx_jet 
    k0=k_offset(npz)+ (k+kb_jet-1) -1      ! 全局index
	i0=i_offset(npx)+ (i+ib_jet-1) -1       
	k1=k0-kb_jet0+1                        ! (窗口中)局部index
	i1=i0-ib_jet0+1
	flow2d_walljet(i,k,m)=flow2d_walljet0(i1,k1,m) 
   enddo
   enddo
   enddo 
  endif 
  endif 
  
  deallocate(flow2d_walljet0)	
  end
	
	
	

