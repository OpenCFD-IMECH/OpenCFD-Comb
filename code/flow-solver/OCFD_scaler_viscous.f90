
  
!-------------Viscous term -----------------------------------------------------------------  
 subroutine du_scaler_viscous(K_SPEC,fs,Dus)
  Use flow_data
  Use viscous_data             ! 临时变量共用 粘性项的变量 （节省内存）
  implicit none
   real(kind=OCFD_REAL_KIND):: fs(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
   real(kind=OCFD_REAL_KIND),dimension(nx,ny,nz):: Dus       
   real(kind=OCFD_REAL_KIND):: cix,ciy,ciz,hi,E1,E2,E3  
   integer:: K_SPEC,i,j,k,m 
!------------------------------------------------------------------------
     if(K_SPEC == 1) then 
	 do k=1,nz 
	 do j=1,ny
	 do i=1,nx 
	 Qpx(i,j,k)=0.d0 
	 Qpy(i,j,k)=0.d0 
	 Qpz(i,j,k)=0.d0 
	 enddo
	 enddo
	 enddo 
	 endif 
	 


       do k=1-LAP,nz+LAP
	   do j=1-LAP,ny+LAP 
	   do i=1-LAP,nx+LAP 
	    Ev1(i,j,k,1)=fs(i,j,k)/d(i,j,k)               ! ci (质量比分）
	   enddo
	   enddo
	   enddo 
	  
		call OCFD_dx0(Ev1(1-LAP,1-LAP,1-LAP,1),uk,Scheme%Scheme_Vis)
        call OCFD_dy0(Ev1(1-LAP,1-LAP,1-LAP,1),ui,Scheme%Scheme_Vis)
        call OCFD_dz0(Ev1(1-LAP,1-LAP,1-LAP,1),us,Scheme%Scheme_Vis)

!-------------------------------------------------------------

      do k=1,nz
      do j=1,ny
      do i=1,nx
        call comput_hi(K_SPEC,T(i,j,k),hi)
		cix=uk(i,j,k)*Akx(i,j,k)+ui(i,j,k)*Aix(i,j,k)+us(i,j,k)*Asx(i,j,k)            ! dci/dx
        ciy=uk(i,j,k)*Aky(i,j,k)+ui(i,j,k)*Aiy(i,j,k)+us(i,j,k)*Asy(i,j,k)
        ciz=uk(i,j,k)*Akz(i,j,k)+ui(i,j,k)*Aiz(i,j,k)+us(i,j,k)*Asz(i,j,k)
        E1=d(i,j,k)*AmD(i,j,k,K_SPEC)*cix
        E2=d(i,j,k)*AmD(i,j,k,K_SPEC)*ciy
        E3=d(i,j,k)*AmD(i,j,k,K_SPEC)*ciz
	    
		Qpx(i,j,k)=Qpx(i,j,k)+E1*hi            ! 传质效应造成的能量输运项 （供能量方程使用）
	    Qpy(i,j,k)=Qpy(i,j,k)+E2*hi
	    Qpz(i,j,k)=Qpz(i,j,k)+E3*hi
		
		Ev1(i,j,k,1)=E1*Akx1(i,j,k)+E2*Aky1(i,j,k)+E3*Akz1(i,j,k) 
		Ev1(i,j,k,2)=E1*Aix1(i,j,k)+E2*Aiy1(i,j,k)+E3*Aiz1(i,j,k) 		
		Ev1(i,j,k,3)=E1*Asx1(i,j,k)+E2*Asy1(i,j,k)+E3*Asz1(i,j,k) 
	  enddo 
	  enddo 
	  enddo 
	  
       call exchange_boundary_x(Ev1(1-LAP,1-LAP,1-LAP,1))
       call exchange_boundary_y(Ev1(1-LAP,1-LAP,1-LAP,2))
       call exchange_boundary_z(Ev1(1-LAP,1-LAP,1-LAP,3))	
	   
	   call OCFD_dx0(Ev1(1-LAP,1-LAP,1-LAP,1),uk,Scheme%Scheme_Vis)
       call OCFD_dy0(Ev1(1-LAP,1-LAP,1-LAP,2),ui,Scheme%Scheme_Vis)
       call OCFD_dz0(Ev1(1-LAP,1-LAP,1-LAP,3),us,Scheme%Scheme_Vis)
	   
       do k=1,nz
       do j=1,ny
       do i=1,nx
        DuS(i,j,k)=Dus(i,j,k)+(uk(i,j,k)+ui(i,j,k)+us(i,j,k))*Ajac(i,j,k)
       enddo
       enddo
       enddo
	  
	  end
	  

