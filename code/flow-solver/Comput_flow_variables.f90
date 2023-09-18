!----------------------------------------------

   subroutine comput_flow_variables
      Use flow_data
      implicit none
      integer:: Num_NegT,i,j,k,m
	  real(kind=OCFD_REAL_KIND):: di0(N_SPEC)
      Num_NegT=0
      do k=1,nz
      do j=1,ny
      do i=1,nx
        d(i,j,k)=f(i,j,k,1)
        u(i,j,k)=f(i,j,k,2)/f(i,j,k,1)
        v(i,j,k)=f(i,j,k,3)/f(i,j,k,1)
        w(i,j,k)=f(i,j,k,4)/f(i,j,k,1)
		Et(i,j,k)=f(i,j,k,5) -(f(i,j,k,2)*u(i,j,k) +f(i,j,k,3)*v(i,j,k) +f(i,j,k,4)*w(i,j,k))*0.5d0  ! 内能
       if(Et(i,j,k)  <= 0) then
	    call handle_NegativeEt(i,j,k,Num_NegT)
	   endif    
      enddo
	  enddo
	  enddo
!----------------------------------	  
      call Rescaling_di		! di(1)+....+di(N_SPEC)=d
!----------------------------------        
     call comput_T_P
 
      end
	
   subroutine comput_T_P
     Use flow_data
     implicit none 
	 integer:: i,j,k,m 
	 real(kind=OCFD_REAL_KIND):: di0(N_SPEC)	 
       do k=1,nz 
	   do j=1,ny 
	   do i=1,nx
        do m=1,N_SPEC
		di0(m)=di(i,j,k,m)
        enddo 
        call comput_T(di0,Et(i,j,k),T(i,j,k))
        call comput_P(di0,T(i,j,k),p(i,j,k))
       enddo
       enddo
       enddo
    end 
	
	
	
   subroutine comput_cc           ! 计算（近似）声速
      Use flow_data
      implicit none	 
	  integer:: i,j,k
	  real(kind=OCFD_REAL_KIND):: gamma
	  do k=1-LAP,nz+LAP 
	  do j=1-LAP,ny+LAP
	  do i=1-LAP,nx+LAP
	   gamma=p(i,j,k)/Et(i,j,k)+1.d0            ! 等效比热比
       cc(i,j,k)=sqrt(gamma*p(i,j,k)/d(i,j,k))  ! 近似声速 （通量分裂使用）
	  enddo
	  enddo
	  enddo 
	end 
	
	


! Transport coefficients
  subroutine comput_Amu_Amk_AMD
   use flow_data
   use CHEM, only: SPEC
!   use Transport
   implicit none
    real(kind=OCFD_REAL_KIND):: di0(N_SPEC),Xi(N_SPEC),AmDi(N_SPEC),X0,Cp 
	integer:: i,j,k,m 
	do k=1,nz
	do j=1,ny
	do i=1,nx
	  X0=0.d0 
      do m=1,N_SPEC 
	    di0(m)=di(i,j,k,m)
	    Xi(m)=di0(m)/SPEC(m)%Mi
		X0=X0+Xi(m)
	  enddo 
	  do m=1,N_SPEC
	   Xi(m)=Xi(m)/X0      ! mole ratio
	  enddo 
	  call comput_vis_coeff(Amu(i,j,k),T(i,j,k),Xi)
	  call comput_Cp(d(i,j,k),di0,T(i,j,k),Cp)
	  call comput_thermal_conductivity(Amk(i,j,k),Amu(i,j,k),Cp,T(i,j,k),Xi)
      call comput_diffusion_coeff(AmDi,Amu(i,j,k),T(i,j,k),d(i,j,k),p(i,j,k),Xi)
	  do m=1,N_SPEC 
	  AmD(i,j,k,m)=AmDi(m) 
	  enddo 
	enddo
	enddo
	enddo 
   end 	  	
 
! rescaling di so that di(1)+di(2)+...+di(N_SPEC)=d0 
     subroutine Rescaling_di
      Use flow_data
      implicit none 
	  integer:: i,j,k,m
	  real(kind=OCFD_REAL_KIND):: d0,tmp
	  do k=1,nz
	  do j=1,ny
	  do i=1,nx
	  d0=0.d0
	  do m=1,N_SPEC 
	  d0=d0+di(i,j,k,m)
	  enddo 
	  tmp=d(i,j,k)/d0 
	  do m=1,N_SPEC 
	  di(i,j,k,m)=di(i,j,k,m)*tmp
	  enddo 
	  enddo
	  enddo
	  enddo
	  end 