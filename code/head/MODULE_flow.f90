
!---------flow data---------------------------------------------------------
module flow_data
 use flow_para 
 implicit none
 
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: Axx,Ayy,Azz,  &   ! Coordinates 
           Akx,Aky,Akz,Aix,Aiy,Aiz,Asx,Asy,Asz, Ajac , &       !  Jacobian coefficients
           Akx1,Aky1,Akz1,Aix1,Aiy1,Aiz1,Asx1,Asy1,Asz1          ! Akx1=Akx/Ajac 
		   
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:,:):: f,fn,du
 real(kind=OCFD_REAL_KIND),allocatable:: di(:,:,:,:),din(:,:,:,:)
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: d,u,v,w,T,p,cc,Et,Amu,AmK
 real(kind=OCFD_REAL_KIND),allocatable:: AmD(:,:,:,:)

end 


module inviscous_data                 ! data used by inviscous term
 use OCFD_precision
 implicit none 
! real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: cc          ! speed of sound
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:) ::   &
       Et1x,p1x,c1x,d1x,u1x,v1x,w1x,A1x,A2x,A3x,hhx1,hhx2, & 
	   Et1y,p1y,c1y,d1y,u1y,v1y,w1y,A1y,A2y,A3y,hhy1,hhy2, &
  	   Et1z,p1z,c1z,d1z,u1z,v1z,w1z,A1z,A2z,A3z,hhz1,hhz2
  real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:):: fpx,fmx,fpy,fmy,fpz,fmz, hhx,hhy,hhz
  integer,allocatable,dimension(:)::Scm_Hbx,Scm_Hby,Scm_Hbz        ! hybrid scheme index
end 



module  viscous_data                 ! data used by inviscous term
 use OCFD_precision
 implicit none 
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:,:)::  Ev1,Ev2,Ev3
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:)::  uk,ui,us,vk,vi,vs,wk,wi,ws,Tk,Ti,Ts
 real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:,:):: Qpx,Qpy,Qpz
end 
