! code for particles, Code by Li Xinliang, lixl@imech.ac.cn------------------
! 颗粒流动计算程序， 颗粒相采用拉格朗日计算，追踪每个粒子
! 本程序被opencfd-sc 或opencfd-comb调用。 本程序不改变主流程的程序结构
! 本程序使用opencfd-sc/opencfd-comb的流场数据 （nx,ny,nz, 坐标数据、u,v,w,Amu等） 
! 本版本仅支持直角网格
!------------------------------------------------------------------
subroutine comput_particle
 use flow_data
 use particle_data
 implicit none 
 integer,save:: Kflag=0
 !--------------初始化（仅执行一次）--------------------
 if(Kflag ==0) then 
  Kflag=1
  call init_particle 
 endif 
 
!----------时间推进------------------------------- 
  call TimeAdv_particle       ! 时间推进， 1阶Euler方法
  call particle_redistribute  ! 根据颗粒的流向，重新分布颗粒：注销流出（本进程）计算域的颗粒，纳入流入（本进程）计算域的颗粒
  call Ana_particle           ! 后处理
!--------显示信息、输出数据------------ 	  
   if( mod(Istep,Para%Istep_show).eq.0  ) then 
     call show_particle_msg  ! particle message
   endif 
   if(mod(Istep, Para%Istep_Save).eq.0)   then 
     call save_particle_data  
   endif 
 end 
 

!颗粒方程时间推进
 subroutine TimeAdv_particle
 use flow_data
 use particle_data
 implicit none  
 real(kind=OCFD_REAL_KIND):: Fp(3),dtx,dty,dtz,dt,mp 
 integer:: k,kdel,k1,kp
 TYPE(particle),pointer:: Pt
 real(kind=OCFD_REAL_KIND),parameter:: PI=3.1415926535897932d0
 kdel=0
 dt=Para%dt
 do k=1, Num_particle
  Pt=>ptcals(k)
  call comput_ForceP(Fp,ptcals(k))
  
  if(Para_P%Coupling_particle > 1) then        ! 双向耦合 or 四向耦合
  ! 根据颗粒对流场的反作用力，更新流场（守恒变量）
   if(Para_P%Force_feedback ==1 ) then    
    call updata_flow_by_particle1(Fp,ptcals(k))   ! 1阶精度插值   
   else 
    call updata_flow_by_particle2(Fp,ptcals(k))   ! 2阶精度插值 
   endif
 endif 

! 更新颗粒速度及坐标
  mp=Pt%rhop*PI*(Pt%dp)**3/6.d0         ! 颗粒质量

 ! 时间推进,更新粒子的位置及速度信息， 1阶Euler 
  dtx=Pt%up*dt ; dty=Pt%vp*dt ; dtz=Pt%wp*dt
  Pt%xp=Pt%xp+dtx
  Pt%yp=Pt%yp+dty 
  Pt%zp=Pt%zp+dtz   
  Pt%up=Pt%up+Fp(1)/mp*dt
  Pt%vp=Pt%vp+Fp(2)/mp*dt 
  Pt%wp=Pt%wp+Fp(3)/mp*dt  
  
  call updata_Pijk(ptcals(k),dtx,dty,dtz) !更新颗粒所在网格位置（i-,j-,k-）
  enddo 

! 对于双向耦合（及4向耦合）， 由于颗粒改变了流场，因此需要更新物理量及边界条件  
  if( Para_P%Coupling_particle > 1) then  
   call comput_flow_variables           ! comput d,u,v,w,T,p,Et 
   call OCFD_bc             ! boundary condition
  endif 
 
 end 
 
 
 !颗粒在本计算域内，计算颗粒的新位置（颗粒所在网格的左上角点编号） 
 subroutine updata_Pijk(Pt,dtx,dty,dtz)
 use flow_data
 use particle_data
 implicit none 
 integer:: k1
 TYPE(particle) :: Pt
 real(kind=OCFD_REAL_KIND):: dtx,dty,dtz
 
 if(dtx >=0) then 
   do k1=Pt%ip,nx  
!   if(Pt%xp>= x1d(k1) .and. Pt%xp < x1d(k+1)) then   ! bug 
   if(Pt%xp>= x1d(k1) .and. Pt%xp < x1d(k1+1)) then    ! bug removed, 2024-3-21
   Pt%ip=k1
   exit
   endif 
   enddo 
  else 
   do k1=Pt%ip,1,-1  
   if(Pt%xp>= x1d(k1) .and. Pt%xp < x1d(k1+1)) then 
   Pt%ip=k1
   exit
  endif 
  enddo  
 endif 
 
 if(dty >=0) then 
   do k1=Pt%jp,ny  
   if(Pt%yp>= y1d(k1) .and. Pt%yp < y1d(k1+1)) then 
   Pt%jp=k1
   exit
   endif 
   enddo 
  else 
  do k1=Pt%jp,1,-1  
   if(Pt%yp>= y1d(k1) .and. Pt%yp < y1d(k1+1)) then 
   Pt%jp=k1
   exit
  endif 
  enddo  
 endif 
 
 if(dtz >=0) then 
   do k1=Pt%kp,nz  
   if(Pt%zp>= z1d(k1) .and. Pt%zp < z1d(k1+1)) then 
   Pt%kp=k1
   exit
   endif 
   enddo 
  else 
  do k1=Pt%kp,1,-1  
   if(Pt%zp>= z1d(k1) .and. Pt%zp < z1d(k1+1)) then 
   Pt%kp=k1
   exit
  endif 
  enddo   
 endif 
 
 end    
 



 
 ! 计算粒子的受力,d0,u0,mu0 是流体密度、速度、粘性系数； dp, rhop, up 是颗粒直径、密度、速度 
 subroutine comput_ForceP(Fp,Pt)
 use flow_data 
 use particle_data 
 implicit none 
 TYPE(particle):: Pt
 real(kind=OCFD_REAL_KIND):: Fp(3),d0,u0,v0,w0,mu0,vv, Rep,fd,Cd,Vp
 real(kind=OCFD_REAL_KIND),parameter:: PI=3.1415926535897932d0
 
!-------阻力（曳力）--------------------- 
 call comput_u_particle(Pt%xp,Pt%yp,Pt%zp,Pt%ip,Pt%jp,Pt%kp,d0,u0,v0,w0,mu0)
 
 vv=sqrt((u0-Pt%up)**2+(v0-Pt%vp)**2+(w0-Pt%wp)**2)        ! 颗粒与流体的相对速度 （滑移速度）
 Rep=d0*vv*Pt%dp/mu0   ! 颗粒Reynolds数
 if(Rep<1.E-3) then      ! 防止Cd分母为0  
 ! Cd=24.d0/Rep
 ! fd=d0*PI*Pt%dp*Pt%dp*vv*vv/8.d0* Cd 
  fd=3.d0*PI*Pt%dp*mu0
 else
 if(Rep>=1000.d0) then 
 Cd=0.44d0 
 else 
 Cd=24.d0/Rep*(1.d0+0.15d0*Rep**0.687d0)           ! 阻力系数
 endif 
  fd=d0*PI*Pt%dp*Pt%dp*Cd*vv/8.d0   ! 受力
 endif 
 
 Fp(1)=fd*(u0-Pt%up)            ! 颗粒受力
 Fp(2)=fd*(v0-Pt%vp)
 Fp(3)=fd*(w0-Pt%wp)
!--------压力梯度力-----------------
 Vp=PI*(Pt%dp)**3/6.d0      ! 颗粒体积
 Fp(1)=Fp(1)-Vp*(p(Pt%ip+1,Pt%jp,Pt%kp)-p(Pt%ip,Pt%jp,Pt%kp))/(x1d(Pt%ip+1)-x1d(pt%ip))             ! 压力梯度力
 Fp(2)=Fp(2)-Vp*(p(Pt%ip,Pt%jp+1,Pt%kp)-p(Pt%ip,Pt%jp,Pt%kp))/(y1d(Pt%jp+1)-y1d(pt%jp))             ! 压力梯度力
 Fp(3)=Fp(3)-Vp*(p(Pt%ip,Pt%jp,Pt%kp+1)-p(Pt%ip,Pt%jp,Pt%kp))/(z1d(Pt%kp+1)-z1d(pt%kp))             ! 压力梯度力

 end 
 
 
! 插值计算粒子位置（xp,yp,zp）处的流体密度、流体速度、流体粘性系数 （采用双线性插值计算）
 subroutine comput_u_particle(xp,yp,zp,ip,jp,kp,d0,u0,v0,w0,mu0)
 use flow_data
 use particle_data
 implicit none 
 real(kind=OCFD_REAL_KIND):: xp,yp,zp,d0,u0,v0,w0,mu0,dhx,dhy,dhz
 integer:: kflag,i,j,k,ip,jp,kp
 !---------------------
 i=ip
 j=jp
 k=kp
 dhx=(xp-x1d(i))/(x1d(i+1)-x1d(i))      ! dhx=dx/hx
 dhy=(yp-y1d(j))/(y1d(j+1)-y1d(j))
 dhz=(zp-z1d(k))/(z1d(k+1)-z1d(k))
 
 
! 判断是否为计算域 边或角区  (角区由于信息不足，采用1阶插值）
 kflag=0
 if(ip==nx) kflag=kflag+1
 if(jp==ny) kflag=kflag+1
 if(kp==nz) kflag=kflag+1
 if(kflag >=2 ) then 
  d0=d(i,j,k)
  u0=u(i,j,k)
  v0=v(i,j,k)
  w0=w(i,j,k)
 else           ! 2阶双线性插值 （Taylor展开）
  d0 =  d(i,j,k)+(d(i+1,j,k)-d(i,j,k))*dhx+(d(i,j+1,k)-d(i,j,k))*dhy+(d(i,j,k+1)-d(i,j,k))*dhz 
  u0=u(i,j,k)+(u(i+1,j,k)-u(i,j,k))*dhx+(u(i,j+1,k)-u(i,j,k))*dhy+(u(i,j,k+1)-u(i,j,k))*dhz 
  v0=v(i,j,k)+(v(i+1,j,k)-v(i,j,k))*dhx+(v(i,j+1,k)-v(i,j,k))*dhy+(v(i,j,k+1)-v(i,j,k))*dhz  
  w0=w(i,j,k)+(w(i+1,j,k)-w(i,j,k))*dhx+(w(i,j+1,k)-w(i,j,k))*dhy+(w(i,j,k+1)-w(i,j,k))*dhz
 endif 
 ! 粘性系数的插值 （粘性系数没有LAP层buffer区，边及角区采用1阶插值）
 if(kflag >= 1) then 
  mu0=Amu(i,j,k)
 else 
  mu0 = Amu(i,j,k)+(Amu(i+1,j,k)-Amu(i,j,k))*dhx+(Amu(i,j+1,k)-Amu(i,j,k))*dhy+(Amu(i,j,k+1)-Amu(i,j,k))*dhz 
 endif  
 end 
   
 
 
 
subroutine show_particle_msg  ! particle message
 use flow_data
 use particle_data
 implicit none 
 integer,parameter::Nmsg=2
 integer,dimension(Nmsg):: Npmsg1,Npmsg2,Npmsg   ! 粒子信息：粒子数目，流出数目 
 integer:: ierr
 Npmsg(1)= Num_particle  ! MPi进程内的粒子数目
 Npmsg(2)= Np_del        ! 流出MPI进程计算域的粒子数
 call MPI_Reduce(Npmsg,Npmsg1,Nmsg,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 call MPI_Reduce(Npmsg,Npmsg2,Nmsg,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
 if(my_id ==0) then 
  print*, "particle Number:",Npmsg1(1), "runout particle:",Npmsg1(2)
  print*, "Max local particle number:", Npmsg2(1), "Max local runout particle:", Npmsg2(2)
 endif 
 end 
  
! 考虑颗粒对流场的影响，更新流场（守恒变量），双向及四向耦合时使用  （插值到颗粒周围的8个网格点） 
 subroutine updata_flow_by_particle2(Fp,Pt) 
 use flow_data
 use particle_data
 implicit none 
 real(kind=OCFD_REAL_KIND):: d0,Fp(3),FF(3),wx(0:1),wy(0:1),wz(0:1),Fj
 integer:: i0,j0,k0,i,j,k,i1,j1,k1
 TYPE(particle):: Pt
  i0=Pt%ip
  j0=Pt%jp 
  k0=Pt%kp
! 三线性插值，将颗粒受力的反作用力，分配到网格的8个节点上  
  wx(0)=abs(x1d(i0+1)-Pt%xp)/(x1d(i0+1)-x1d(i0))  
  wy(0)=abs(y1d(j0+1)-pt%yp)/(y1d(j0+1)-y1d(j0))
  wz(0)=abs(z1d(k0+1)-pt%zp)/(z1d(k0+1)-z1d(k0))
  if(i0 == nx) wx(0)=1.d0 
  if(j0 == ny) wy(0)=1.d0 
  if(k0 == nz) wz(0)=1.d0 
  wx(1)=1.d0-wx(0)
  wy(1)=1.d0-wy(0)
  wz(1)=1.d0-wz(0) 
  i1=1; j1=1; k1=1
  if(i0==nx) i1=0
  if(j0==ny) j1=0
  if(k0==nz) k1=0
  Fj=Ajac(i0,j0,k0)/(hx*hy*hz)   ! 1/(V)   V 网格体积
! ----debug---------------------------  
  if(wx(0) <0 .or. wx(0) > 1) then 
    print*, "wx=",wx
    print*, "i0,j0,k0=",i0,j0,k0
	print*, "Pt%xp=", pt%xp, "x1d(i0),x1d(i0+1)=",x1d(i0),x1d(i0+1)
  endif 
  
  if(wy(0) <0 .or. wy(0) > 1) then 
    print*, "wy=",wy
    print*, "i0,j0,k0=",i0,j0,k0
	print*, "Pt%yp=", pt%yp, "y1d(j0),y1d(j0+1)=",y1d(j0),y1d(j0+1)
  endif 

  if(wz(0) <0 .or. wz(0) > 1) then 
    print*, "wz=",wz
    print*, "i0,j0,k0=",i0,j0,k0
	print*, "Pt%zp=", pt%zp, "z1d(k0),z1d(k0+1)=",z1d(k0),z1d(k0+1)
  endif 
  
! 考虑颗粒反作用力，更新守恒变量  
  do k=0,k1
  do j=0,j1
  do i=0,i1
   ff(:)=-Fp(:)*wx(i)*wy(j)*wz(k)*Fj    ! 网格节点的受力（体积力）
   f(i0+i,j0+j,k0+k,2)=f(i0+i,j0+j,k0+k,2)+ff(1)*Para%dt
   f(i0+i,j0+j,k0+k,3)=f(i0+i,j0+j,k0+k,3)+ff(2)*Para%dt
   f(i0+i,j0+j,k0+k,4)=f(i0+i,j0+j,k0+k,4)+ff(3)*Para%dt
   f(i0+i,j0+j,k0+k,5)=f(i0+i,j0+j,k0+k,5)+ &
    (ff(1)*u(i0+i,j0+j,k0+k)+ff(2)+v(i0+i,j0+j,k0+k)+ff(3)*w(i0+i,j0+j,k0+k))*Para%dt
  enddo 
  enddo 
  enddo 
  end 
  
 
! 考虑颗粒对流场的影响，更新流场（守恒变量），双向及四向耦合时使用  （施加到左上网格点，1阶精度） 
 subroutine updata_flow_by_particle1(Fp,Pt) 
 use flow_data
 use particle_data
 implicit none 
 real(kind=OCFD_REAL_KIND):: d0,Fp(3),FF(3),Fj
 integer:: i0,j0,k0
 TYPE(particle):: Pt
  i0=Pt%ip
  j0=Pt%jp 
  k0=Pt%kp

   Fj=Ajac(i0,j0,k0)/(hx*hy*hz)   ! 1/(V)   V 网格体积
   ff(:)=-Fp(:)*Fj 
   f(i0,j0,k0,2)=f(i0,j0,k0,2)+ff(1)*Para%dt
   f(i0,j0,k0,3)=f(i0,j0,k0,3)+ff(2)*Para%dt
   f(i0,j0,k0,4)=f(i0,j0,k0,4)+ff(3)*Para%dt
   f(i0,j0,k0,5)=f(i0,j0,k0,5)+ (ff(1)*u(i0,j0,k0)+ff(2)+v(i0,j0,k0)+ff(3)*w(i0,j0,k0))*Para%dt  

  end 
 
 
 
 
 
 
 
 
 
 
 
 