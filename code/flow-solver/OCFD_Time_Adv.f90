


! Time advance for flow 
! 流场时间推进（不考虑化学反应源项）

   subroutine flow_time_advance(dt1)
   use flow_data
   implicit none
   integer:: KRK,i,j,k,m,ierr 
   real(kind=OCFD_REAL_KIND):: dt1
        do m=1,5
		do k=1,nz 
        do j=1,ny
        do i=1,nx 
        fn(i,j,k,m)=f(i,j,k,m)
        enddo
        enddo
        enddo 
		enddo 
	
	    do m=1,N_SPEC
		do k=1,nz 
        do j=1,ny
        do i=1,nx 
        din(i,j,k,m)=di(i,j,k,m)
        enddo
        enddo
        enddo 
		enddo    
		

      do  KRK=1,3            ! 3-step Runge-Kutta 
!       call comput_flow_variables           ! comput d,u,v,w,T,p,Et 
       call comput_Amu_Amk_AMD
      
	     call exchange_boundary_xyz(d)
         call exchange_boundary_xyz(u)
         call exchange_boundary_xyz(v)
         call exchange_boundary_xyz(w)
         call exchange_boundary_xyz(T) 
         call exchange_boundary_xyz(p) 	  
         call exchange_boundary_xyz(Et) 
        
		if(N_SPEC ==1) then         ! 单组分
		 do k=1-LAP,nz+LAP 
		 do j=1-LAP,ny+LAP
		 do i=1-LAP,nx+LAP
		 di(i,j,k,1)=d(i,j,k)
         enddo 
         enddo 
         enddo 		 
		else 
	     do m=1,N_SPEC
 	     call exchange_boundary_xyz(di(1-LAP,1-LAP,1-LAP,m)) 
         enddo 
		endif 

!------------Scalar ------------------------	  
	  if(KRK .eq. 1 .and. Scheme%Scheme_Invis == OCFD_Scheme_Hybrid ) then 
       call comput_Rhybrid           ! comput only KRK=1 
      endif 	   
		

     if(N_SPEC .ne. 1) then         ! if N_SPEC ==1 单组分 
	  do m=1,N_SPEC 
       call  du_scaler_invis(di(1-LAP,1-LAP,1-LAP,m),du(1,1,1,1))	  
	   if(Para%IF_Viscous .eq. 1) then 
       call du_scaler_viscous(m,di(1-LAP,1-LAP,1-LAP,m),du(1,1,1,1))	   
	   endif  
	   call scalar_time_adv_RK3 (KRK,m,du(1,1,1,1),dt1) 
	  enddo 	 
	 endif 
	 
		 
      if(Para%IF_Scheme_Character .eq. 0) then 
        call du_inviscous    ! inviscous term  (non-character type)
      else 
        call du_inviscous_Character     ! inviscous term (character type)
      endif 
      
      if(Para%IF_Viscous .eq. 1)  then 
       call du_viscous
      endif 
      
      call OCFD_time_adv_RK3 (KRK,dt1)          ! time advance (3-step RK)
      
      call comput_flow_variables           ! comput d,u,v,w,T,p,Et 
      call OCFD_bc             ! boundary condition

      enddo
     end
	 


!-----------------------------------------------------------------------------------
! 3 steps 3rd order TVD type Runge-Kutta method by Jiang & Shu  
  subroutine OCFD_time_adv_RK3 (KRK,dt)
  use flow_data
  implicit none 
  integer  KRK,m,i,j,k
  real(kind=OCFD_REAL_KIND):: Ralfa(3),Rbeta(3),dt 
   
        Ralfa(1)=0.d0 ;       Rbeta(1)=1.d0
        Ralfa(2)=3.d0/4.d0 ;  Rbeta(2)=1.d0/4.d0
        Ralfa(3)=1.d0/3.d0 ;  Rbeta(3)=2.d0/3.d0

   if(KRK.eq.1) then
     do m=1,5
	 do k=1,nz
     do j=1,ny
     do i=1,nx
!       f(i,j,m)=Ralfa(KRK)*fn(i,j,m) +dt*du(i,j,m)*Rbeta(KRK)
        f(i,j,k,m)=fn(i,j,k,m) +dt*du(i,j,k,m)
     enddo
     enddo
     enddo
	 enddo 
   else
     do m=1,5
	 do k=1,nz
     do j=1,ny
     do i=1,nx
        f(i,j,k,m)=Ralfa(KRK)*fn(i,j,k,m)+Rbeta(KRK)* (f(i,j,k,m) + dt*du(i,j,k,m))
     enddo
     enddo
     enddo
	 enddo 
 endif
end

!-----------------------------------------------------------------------------------
! 3 steps 3rd order TVD type Runge-Kutta method by Jiang & Shu  
  subroutine Scalar_time_adv_RK3 (KRK,m,dus,dt)
  use flow_data
  implicit none 
  integer  KRK,m,i,j,k
  real(kind=OCFD_REAL_KIND):: Ralfa(3),Rbeta(3),dt 
  real(kind=OCFD_REAL_KIND):: dus(nx,ny,nz) 
        Ralfa(1)=0.d0 ;       Rbeta(1)=1.d0
        Ralfa(2)=3.d0/4.d0 ;  Rbeta(2)=1.d0/4.d0
        Ralfa(3)=1.d0/3.d0 ;  Rbeta(3)=2.d0/3.d0

   if(KRK.eq.1) then

	 do k=1,nz
     do j=1,ny
     do i=1,nx
!       f(i,j,m)=Ralfa(KRK)*fn(i,j,m) +dt*du(i,j,m)*Rbeta(KRK)
        di(i,j,k,m)=din(i,j,k,m) +dt*dus(i,j,k)
     enddo
     enddo
     enddo

   else

	 do k=1,nz
     do j=1,ny
     do i=1,nx
        di(i,j,k,m)=Ralfa(KRK)*din(i,j,k,m)+Rbeta(KRK)* (di(i,j,k,m) + dt*dus(i,j,k))
     enddo
     enddo
     enddo

 endif
end

!------化学反应时间推进， 推进方法由全局变量CHEM_TimeAdv 设定------------------------------------------------
    subroutine  Chemical_time_advance (dt2)
    use flow_data
	implicit none
	real(kind=OCFD_REAL_KIND):: Et1,Et0,di1(N_Spec),dt2
    integer:: i,j,k,m
     do k=1,nz
	 do j=1,ny
     do i=1,nx
     do m=1,N_SPEC
	  di1(m)=di(i,j,k,m)                                ! 组分密度
	 enddo
!       Et1=f(i,j,k,5)-0.5d0*(f(i,j,k,2)**2 + f(i,j,k,3)**2 +f(i,j,k,4)**2)/f(i,j,k,1)     ! 内能 （不含反应焓）= 总能-动能
       Et1=Et(i,j,k)         ! 内能     
	   Et0=Et1
! 化学反应时间推进	dt2 时间步;  支持1阶Euler, 2阶RK, 3阶RK, 1阶隐格式， 通过全局变量CHEM_TimeAdv设定 (1,2,3,-1)	
      call chemical_timeAdv(di1,Et1,dt2)             ! 利用OpenCFD-Comb 0维程序代码
 
  	  do m=1,N_SPEC
	   di(i,j,k,m)=di1(m)                                  ! 更新组分密度
	  enddo
       f(i,j,k,5)=f(i,j,k,5)+(Et1-Et0)                     ! 更新总能 （本程序中的内能不含反应焓）
       Et(i,j,k)=Et1                                       ! 更新内能 （不含反应焓)
	enddo
    enddo
    enddo


   end
!---------------------------------------------		

  
     subroutine show_flow_msg(wall_time) 
      use flow_data
      implicit none
      real*8:: wall_time, wtmp
	  real(kind=OCFD_REAL_KIND):: E0(2),E1(2),tmp0
 	  integer:: i,j,k,ierr
	  wtmp=wall_time
      wall_time=MPI_wtime()  

	  
	  E1(:)=0.d0
      E0(:)=0.d0 
	  do k=1,nz
	  do j=1,ny 
	  do i=1,nx 
	  tmp0=1.d0/Ajac(i,j,k)
	  E1(1)=E1(1)+f(i,j,k,5)*tmp0                ! Total Energy 
	  E1(2)=E1(2)+d(i,j,k)*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))*0.5d0*tmp0  !Kinetic Energy
	  enddo
	  enddo 
	  enddo 
	  
      call MPI_REDUCE(E1,E0,2,OCFD_DATA_TYPE,  MPI_SUM,0,MPI_COMM_WORLD,ierr)	  
	  E0=E0/(nx_global*ny_global*nz_global)	 

	  
	  if(my_id .eq. 0) then 
		print*, 'Istep=',Istep,'tt=',tt
        print*, 'CPU time per', Para%Istep_show,'step is:    ',wall_time-wtmp
		print*, "Averaged Total-energy E, Kinetic-energy K  are"
		write(*,"(3E25.15)") E0(1), E0(2) 
		
		open(66,file='opencfd.log',position='append')
        write(66,*) 'Istep=',Istep,'tt=',tt,'CPU time is',wall_time-wtmp
		write(66,*)  "Averaged Total-energy E, Kinetic-energy K  are"
		write(66,"(3E25.15)") E0(1), E0(2) 
        close(66)
      endif 

      end 


