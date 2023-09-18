  subroutine du_scaler_invis(fs,Dus)
  use flow_data
  implicit none 
  real(kind=OCFD_REAL_KIND):: fs(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP)
  real(kind=OCFD_REAL_KIND),dimension(nx,ny,nz):: Dus       
  real(kind=OCFD_REAL_KIND):: fpx(1-LAP:nx+LAP),fmx(1-LAP:nx+LAP),fpy(1-LAP:ny+LAP),fmy(1-LAP:ny+LAP), &
                              fpz(1-LAP:nz+LAP),fmz(1-LAP:nz+LAP), &
							  hhx1(0:nx),hhx2(0:nx),hhy1(0:ny),hhy2(0:ny),hhz1(0:nz),hhz2(0:nz), &
	                          Scm_Hbx(0:nx),Scm_Hby(0:ny),Scm_Hbz(0:nz)						  
  
  real(kind=OCFD_REAL_KIND):: vs,hx_1, hy_1, hz_1,Rhb0
  integer:: i,j,k,m,set_hybrid_scheme
  
    hx_1=1.d0/hx ; 	hy_1=1.d0/hy ;   hz_1=1.d0/hz   

!---------Convertive terms --------------------------------  
!--------i-direction  ------------------------------- 
  do k=1,nz 
  do j=1,ny 
  
   if(Scheme%Scheme_Invis == 	OCFD_Scheme_Hybrid) then 
   do i=0,nx
	  Rhb0=max(Rhybrid(i,j,k),Rhybrid(i+1,j,k))
	  Scm_Hbx(i)=set_hybrid_scheme(Rhb0)
   enddo
   endif 
  
  do i=1-LAP,nx+LAP 
   vs=Akx1(i,j,k)*u(i,j,k)+Aky1(i,j,k)*v(i,j,k)+Akz1(i,j,k)*w(i,j,k)
   fpx(i)=0.5d0*(vs+abs(vs))*fs(i,j,k)                   ! S-W Splitting for scaler 
   fmx(i)=0.5d0*(vs-abs(vs))*fs(i,j,k)
  enddo 
   call OCFD2d_flux1(fpx(1-LAP),hhx1,nx,LAP,Scheme%Bound_index(1,1),Scm_Hbx, Scheme%Scheme_boundary(1:2))        
   call OCFD2d_flux2(fmx(1-LAP),hhx2,nx,LAP,Scheme%Bound_index(1,1),Scm_Hbx, Scheme%Scheme_boundary(1:2))    
   do i=1,nx   
   Dus(i,j,k)= (hhx1(i)-hhx1(i-1) + hhx2(i)-hhx2(i-1))*hx_1
   enddo 
  enddo 
  enddo 
  
 !--------j-direction ---------------------
 do k=1,nz 
 do i=1,nx 
  if(Scheme%Scheme_Invis == 	OCFD_Scheme_Hybrid) then 
   do j=0,ny
	  Rhb0=max(Rhybrid(i,j,k),Rhybrid(i,j+1,k))
	  Scm_Hby(j)=set_hybrid_scheme(Rhb0)
    enddo
   endif 
   do j=1-LAP,ny+LAP 
     vs=Aix1(i,j,k)*u(i,j,k)+Aiy1(i,j,k)*v(i,j,k)+Aiz1(i,j,k)*w(i,j,k)
     fpy(j)=0.5d0*(vs+abs(vs))*fs(i,j,k) 
     fmy(j)=0.5d0*(vs-abs(vs))*fs(i,j,k)
   enddo 
	 call OCFD2d_flux1(fpy(1-LAP),hhy1,ny,LAP,   Scheme%Bound_index(1,2),Scm_Hby, Scheme%Scheme_boundary(3:4) )      
     call OCFD2d_flux2(fmy(1-LAP),hhy2,ny,LAP,   Scheme%Bound_index(1,2),Scm_Hby, Scheme%Scheme_boundary(3:4) )

  do j=1,ny   
   Dus(i,j,k)= Dus(i,j,k)+(hhy1(j)-hhy1(j-1) + hhy2(j)-hhy2(j-1))*hy_1
  enddo 
 enddo 
 enddo  
!--------k-direction 
 do j=1,ny 
 do i=1,nx 
  if(Scheme%Scheme_Invis == 	OCFD_Scheme_Hybrid) then 
   do k=0,nz
	  Rhb0=max(Rhybrid(i,j,k),Rhybrid(i,j,k+1))
	  Scm_Hbz(k)=set_hybrid_scheme(Rhb0)	 
   enddo
  endif  
  do k=1-LAP,nz+LAP 
   vs=Asx1(i,j,k)*u(i,j,k)+Asy1(i,j,k)*v(i,j,k)+Asz1(i,j,k)*w(i,j,k)
   fpz(k)=0.5d0*(vs+abs(vs))*fs(i,j,k) 
   fmz(k)=0.5d0*(vs-abs(vs))*fs(i,j,k)
  enddo 
   call OCFD2d_flux1(fpz(1-LAP),hhz1,nz,LAP,   Scheme%Bound_index(1,3),Scm_Hbz, Scheme%Scheme_boundary(5:6) )    
   call OCFD2d_flux2(fmz(1-LAP),hhz2,nz,LAP,   Scheme%Bound_index(1,3),Scm_Hbz, Scheme%Scheme_boundary(5:6) )
  do k=1,nz 
   Dus(i,j,k)= -(Dus(i,j,k)+(hhz1(k)-hhz1(k-1) + hhz2(k)-hhz2(k-1))*hz_1)*Ajac(i,j,k)
  enddo
  enddo 
  enddo 
  end 
!---------------------------------------------------------------------
