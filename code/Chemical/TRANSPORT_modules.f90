!  计算输运参数：  粘性系数、导热系数、扩散系数
!  ver 1.5, 2024-12-18: 增加了多种形式的输运系数拟合公式，如Wilke, Gupta, Blottner 

!------ log-type fitting equation suggested by Chemkin 4.1, See Chemkin 4.1 theory manual (Chap 5.1);
! Gupta: NASA reference publication 1232, 1990

module Transport
 use Precision_EC
 implicit none 
 integer,parameter:: Amu_by_Sutherland=0,   Amk_by_Pr =0, AmD_by_Sc=0 , &
	                 Amu_by_Wilke=1,   Amk_by_Wilke =1, AmD_by_Wilke=1, &
                   Amu_by_Gupta = 2, Amk_by_Gupta = 2, AmD_by_Gupta = 2, &
                   Amu_by_Blottner = 3, Amk_by_Eucken = 3

 
 TYPE Transport_type 
 integer:: Flag_vis,Flag_cond,Flag_diff        ! 粘性系数、热导系数、扩散系数计算方式（0：Surtherland, Pr, Sc 计算；  1： 拟合+Wilke公式计算）
 real(PRE_EC) ::  Amu0, Pr, Sc , Amk0
 real(PRE_EC),dimension(:,:),pointer:: mu_an, k_an       ! 粘性系数， 热导系数  (见理论手册 2.2节）
 real(PRE_EC),dimension(:,:,:),pointer:: D_an            ! 输运系数
 end TYPE 
 TYPE (Transport_type):: Trans 
 
 
end 


subroutine init_transport_para
use CHEM, only: N_SPEC 
use Transport 
implicit none

! allocate(Trans%mu_an(4,N_SPEC), Trans%k_an(4,N_SPEC), Trans%D_an(4,N_SPEC,N_SPEC) )
! if( Trans%Flag_vis .eq. Amu_by_Wilke) then
!    call read_vis_coefficients
! endif

!  if (Trans%Flag_cond .eq. Amk_by_Wilke) then
!    call   read_conduct_coefficients
!  endif

!  if (Trans%Flag_diff .eq. AmD_by_Wilke) then
!     call read_diffusion_coefficients 
!  endif

   if (Trans%Flag_vis .eq. Amu_by_Wilke) then
      allocate (Trans%mu_an(4, N_SPEC))
      call read_vis_coefficients
   elseif (Trans%Flag_vis .eq. Amu_by_Gupta) then
      allocate (Trans%mu_an(5, N_SPEC))
      call read_vis_coefficients
   elseif (Trans%Flag_vis .eq. Amu_by_Blottner) then
      allocate (Trans%mu_an(3, N_SPEC))
      call read_vis_coefficients
   end if

   if (Trans%Flag_cond .eq. Amk_by_Wilke) then
      allocate (Trans%k_an(4, N_SPEC))
      call read_conduct_coefficients
   elseif (Trans%Flag_cond .eq. Amk_by_Gupta) then
      allocate (Trans%k_an(5, N_SPEC))
      call read_conduct_coefficients
   elseif (Trans%Flag_cond .eq. Amk_by_Eucken) then
      !无任何执行
   end if

   if (Trans%Flag_diff .eq. AmD_by_Wilke) then
      allocate (Trans%D_an(4, N_SPEC, N_SPEC))
      call read_diffusion_coefficients
   elseif (Trans%Flag_diff .eq. AmD_by_Gupta) then
      allocate (Trans%D_an(4, N_SPEC, N_SPEC))
      call read_diffusion_coefficients
   end if

end




!-------------------------------------------------
! 粘性系数 的拟合系数 (in g-cm-s )
 subroutine read_vis_coefficients
  use CHEM 
  use Transport  
  implicit none
  integer:: Nm,m,k
  real*8, allocatable:: at(:,:)  
  character(len=10),allocatable::  Name1(:)

 open(99,file="viscosity.in")
 read(99,*)
 read(99,*)
 read(99,*) Nm
!  allocate( at(4,Nm), Name1(Nm) )
!    do m=1,Nm
!    read(99,*) Name1(m)
!    read(99,*) (at(k,m),k=1,4)
!    enddo
   if (Trans%Flag_vis .eq. Amu_by_Wilke) then
      allocate (at(4, Nm), Name1(Nm))
      do m = 1, Nm
         read (99, *) Name1(m)
         read (99, *) (at(k, m), k=1, 4)
      end do
   elseif (Trans%Flag_vis .eq. Amu_by_Gupta) then
      allocate (at(5, Nm), Name1(Nm))
      do m = 1, Nm
         read (99, *) Name1(m)
         read (99, *) (at(k, m), k=1, 5)
      end do
   elseif (Trans%Flag_vis .eq. Amu_by_Blottner) then
      allocate (at(3, Nm), Name1(Nm))
      do m = 1, Nm
         read (99, *) Name1(m)
         read (99, *) (at(k, m), k=1, 3)
      end do
   end if
  close(99)

  do k=1,N_SPEC                  ! find spec(k)%name==Name1     寻找组分名称=Name1(m)的
     do m=1,Nm
     if(trim(Spec(k)%name) == trim (Name1(m)) ) then
	!  Trans%mu_an(1:4,k)=at(1:4,m)
  Trans%mu_an(:, k) = at(:, m)
     goto 120
     endif
   enddo
   print*, "Error ! Can not fine viscosity coeffcients of ", Spec(k)%name 
   stop
120 continue
 enddo

deallocate(at,Name1)

!  print*, "------------viscosity-------------------------"
! do k=1,N_SPEC
!  print*, Spec(k)%name
!  print*,mu_an(:,k)
! enddo


 end

!-------------------------------------------------
! 热导系数 的拟合系数
 subroutine read_conduct_coefficients
  use CHEM
  use Transport 
  implicit none
  integer:: Nm,m,k
  real*8, allocatable:: at(:,:)  
  character(len=10),allocatable::  Name1(:)

 open(99,file="conductivity.in")
 read(99,*)
 read(99,*)
 read(99,*) Nm
!  allocate( at(4,Nm), Name1(Nm) )
!    do m=1,Nm
!    read(99,*) Name1(m)
!    read(99,*) (at(k,m),k=1,4)
!    enddo
      if (Trans%Flag_cond .eq. Amk_by_Wilke) then
      allocate (at(4, Nm), Name1(Nm))
      do m = 1, Nm
         read (99, *) Name1(m)
         read (99, *) (at(k, m), k=1, 4)
      end do
   elseif (Trans%Flag_cond .eq. Amk_by_Gupta) then
      allocate (at(5, Nm), Name1(Nm))
      do m = 1, Nm
         read (99, *) Name1(m)
         read (99, *) (at(k, m), k=1, 5)
      end do
   end if
  close(99)

  do k=1,N_SPEC              ! find spec(k)%name==Name
     do m=1,Nm
     if(trim(Spec(k)%name) == trim (Name1(m)) ) then
	!  Trans%k_an(1:4,k)=at(1:4,m)
          Trans%k_an(:, k) = at(:, m)
     goto 130
     endif
   enddo
   print*, "Error ! Can not fine viscosity coeffcients of ", Spec(k)%name 
   stop
130 continue
 enddo

deallocate(at,Name1)

!  print*, "------------conductivity-------------------------"
! do k=1,N_SPEC
!  print*, Spec(k)%name
!  print*,ld_bn(:,k)
! enddo

 end


!------------------------------------------------------------
! 扩散系数的拟合系数
subroutine read_diffusion_coefficients 
use CHEM
use Transport 
implicit none
integer:: Np, i,j,k,m
real*8,allocatable:: dtp(:,:)
character(len=10),allocatable::  Name1(:), Name2(:)

open(99,file="diffusion_coefficients.in")
read(99,*)
read(99,*)
read(99,*)
read(99,*) 
read(99,*) Np
allocate(dtp(4,Np), Name1(Np), Name2(Np))
do m=1,Np
 read(99,*) Name1(m), Name2(m)
 read(99,*) (dtp(k,m), k=1,4)
enddo
close(99)


 do k=1,N_SPEC
 do j=1,N_SPEC
   do m=1,Np    
!   Find the diffusion coefficient;    symmetric of diffusion  dn(:,j,k)=dn(:,k,j)    
  if( (trim(Spec(j)%Name)==trim(Name1(m))  .and. trim(Spec(k)%Name)==trim(Name2(m)) )     &
   .or.  (trim(Spec(j)%Name)==trim(Name2(m))  .and. trim(Spec(k)%Name)==trim(Name1(m))) ) then
	Trans%D_an(1:4,j,k)=dtp(1:4,m)
    goto 140
    endif 
   enddo
   print*, "Error in read diffusion_coefficients.in, can not find ",  Spec(j)%Name, Spec(k)%Name
140 continue
   enddo
   enddo

deallocate(dtp,Name1,Name2)

! print*, "--------------diffusion----------------"
! do k=1,N_SPEC
! do j=1,N_SPEC
! print*, Spec(j)%Name, Spec(k)%Name
! print*, D_dn(:,j,k)
! enddo
! enddo

end 


! 计算粘性系数， 热传导系数， 扩散系数
!----------------------------------------------------

! 计算粘性系数（混合物的统一粘性系数Amu）
! Sutherland eq. or Wilke eq.
 subroutine comput_vis_coeff(Amu,T,Xi)        ! Xi 各组分的摩尔比分
  use CHEM
  use Transport  
  implicit none 
  real(PRE_EC):: Amu, T,Pk0,Pkj,mu
  real(PRE_EC):: mui(N_SPEC), Xi(N_SPEC) 
  real(PRE_EC),parameter:: Sutherland_TA=288.15d0, Sutherland_TC=110.4d0   ! Sutherland 公式中的常数
  integer:: m,j,k
  
   if (Trans%Flag_vis .eq. Amu_by_Sutherland) then          ! Sutherland eq.
      Amu=Trans%Amu0* sqrt( (T/Sutherland_TA)**3 ) * (Sutherland_TA+Sutherland_TC) /(Sutherland_TC+T)
   else 
      if(Trans%Flag_vis .eq. Amu_by_Wilke) then
         do  m=1,N_SPEC            ! 各组分粘性系数 (log fitting by CHEMKIN)
             mui(m)=exp(Trans%mu_an(1,m)+Trans%mu_an(2,m)*log(T)+Trans%mu_an(3,m)*log(T)**2    &        ! 注意单位制
                    +Trans%mu_an(4,m)*log(T)**3) *0.1d0                  ! 0.1  transform for g-cm-s unit to SI unit
         enddo
      elseif (Trans%Flag_vis .eq. Amu_by_Gupta) then
         if (T < 1000.) then
            Amu = Trans%Amu0*sqrt((T/Sutherland_TA)**3)*(Sutherland_TA + Sutherland_TC)/(Sutherland_TC + T)
            goto 150
         else
            do m = 1, N_SPEC            ! 各组分粘性系数 (log fitting by Gupta)
              !  mui(m) = exp(Trans%mu_an(1, m)*log(T)**4 + Trans%mu_an(2, m)*log(T)**3 + Trans%mu_an(3, m)*log(T)**2 &        ! ע�ⵥλ��
              !               + Trans%mu_an(4, m)*log(T) + Trans%mu_an(5, m))*0.1d0                  ! 0.1  transform for g-cm-s unit to SI unit
              mui(m) = 0.1d0*exp(Trans%mu_an(5, m))*T**(Trans%mu_an(1, m)*log(T)**3 + Trans%mu_an(2, m)*log(T)**2  &
                              + Trans%mu_an(3, m)*log(T) + Trans%mu_an(4, m))
            end do
         end if
      elseif (Trans%Flag_vis .eq. Amu_by_Blottner) then
         
         do m = 1, N_SPEC            ! 各组分粘性系数 (log fitting by Blottner)
                  ! 0.1  transform for g-cm-s unit to SI unit
            mui(m) = 0.1d0*exp((Trans%mu_an(1, m)*log(T) + Trans%mu_an(2, m))*log(T) + Trans%mu_an(3, m))
         end do
        
      end if
 
    ! Wilke Equation  (see: 理论手册 2.2节)
      mu=0.d0
      do k=1,N_Spec
         Pk0=0.d0
         do j=1,N_Spec
            Pkj= (1.d0+ sqrt(mui(k)/mui(j)* sqrt(Spec(j)%Mi/Spec(k)%Mi))  )**2/sqrt(8.d0* ( 1.d0+ Spec(k)%Mi/Spec(j)%Mi )  )
            Pk0=Pk0+Xi(j)*Pkj
         enddo
         mu=mu+Xi(k)*mui(k)/Pk0
      enddo
      Amu=mu
   endif    
   150 continue
end subroutine comput_vis_coeff
   
   
! 计算热传导系数k     
  subroutine comput_thermal_conductivity(Amk,Amu,Cp,T,Xi)    ! Amk 导热系数， Amu 粘性系数，Cp 定压比热，T温度， Xi 组分摩尔比分
  use CHEM, only: SPEC,N_SPEC
  use Transport  
  implicit none 
  real(PRE_EC):: Amk,Amu, Cp,T,Amk1,Amk2,Pk0,Pkj
  real(PRE_EC), parameter:: Sutherland_TA = 273.16d0, Sutherland_TC = 194.4d0   ! Sutherland
  
  real(PRE_EC):: ki(N_SPEC), Xi(N_SPEC) 
  real(PRE_EC):: Cv(N_SPEC)
  real(PRE_EC):: Amu_Blottner(N_SPEC)
  integer:: m, k, j
  
  if (Trans%Flag_cond .eq. Amk_by_Pr) then      ! by using Prandtl number 
      Amk=Amu*Cp/Trans%Pr
  elseif (Trans%Flag_cond .eq. Amk_by_Eucken) then
      call comput_Cv_EuckenFormula(Cv)
      do m = 1, N_SPEC            ! 各组分粘性系数 (log fitting by Blottner)
         Amu_Blottner(m) = 0.1d0*exp((Trans%mu_an(1, m)*log(T) + Trans%mu_an(2, m))*log(T) + Trans%mu_an(3, m))
         ki(m) = Amu_Blottner(m) * Cv(m)
      end do
      ! Wilke Equation 
      Amk1=0.d0
      do k=1,N_Spec
         Pk0=0.d0
         do j=1,N_Spec
            Pkj= (1.d0+ sqrt(ki(k)/ki(j)* sqrt(Spec(j)%Mi/Spec(k)%Mi))  )**2/sqrt(8.d0* ( 1.d0+ Spec(k)%Mi/Spec(j)%Mi )  )
            Pk0=Pk0+Xi(j)*Pkj
         enddo
         Amk1=Amk1+Xi(k)*ki(k)/Pk0
      enddo
      Amk=Amk1
  else
      if (Trans%Flag_cond .eq. Amk_by_Wilke) then
      do m=1,N_SPEC  
          ki(m)=exp(Trans%k_an(1,m)+Trans%k_an(2,m)*log(T)+Trans%k_an(3,m)*log(T)**2    & 
                +Trans%k_an(4,m)*log(T)**3) *1.d-5                   ! 1.d-5  transform for g-cm-s unit to SI unit  
      enddo  
      elseif (Trans%Flag_cond .eq. Amk_by_Gupta) then
          if (T < 1000.) then
              Amk = Trans%Amk0*sqrt((T/Sutherland_TA)**3)*(Sutherland_TA + Sutherland_TC)/(Sutherland_TC + T)
              goto 160
          else
              do m = 1, N_SPEC
                  ! ki(m) = exp(Trans%k_an(1, m)*log(T)**4 + Trans%k_an(2, m)*log(T)**3 + Trans%k_an(3, m)*log(T)**2 &
                  !          + Trans%k_an(4, m)*log(T) + Trans%k_an(5, m))*418.4                  ! 418.4  transform for cal/(cm*s*K) unit to SI unit
                  ki(m) = 418.4 * exp(Trans%k_an(5, m))*T**(Trans%k_an(1, m)*log(T)**3 + Trans%k_an(2, m)*log(T)**2  &
                              + Trans%k_an(3, m)*log(T) + Trans%k_an(4, m))
              end do
          end if
      end if
      Amk1=0.d0
      Amk2=0.d0
      do m=1,N_Spec
        Amk1=Amk1+Xi(m)*ki(m)
        Amk2=Amk2+Xi(m)/ki(m)
      enddo
      Amk=0.5d0*(Amk1+1.d0/Amk2)
  
  endif 
  160 continue
  end subroutine comput_thermal_conductivity

!-----------------------------------------------------------
  subroutine comput_diffusion_coeff(AmDi,Amu,T,d,p,Xi)
  use CHEM
  use Transport  
  implicit none 
  real(PRE_EC):: Amk,Amu,T,d,p,d0,ds
  real(PRE_EC):: AmDi(N_SPEC), Xi(N_SPEC),di(N_SPEC), D2(N_SPEC,N_SPEC) 
  integer:: m,j,k 
  real(PRE_EC):: epsl=1.d-20
  
  if (Trans%Flag_diff .eq. AmD_by_Sc) then      ! by using Prandtl number 
      d0=0.d0 
      do m=1,N_SPEC 
          d0=d0+Xi(m)*Spec(m)%Mi 
      enddo 
      do m=1,N_SPEC 
          AmDi(m)=Amu*(1.d0-Xi(m)*Spec(m)%Mi/d0)/((1.d0-Xi(m))*Trans%Sc*d +epsl)
      enddo 
  else 
      if (Trans%Flag_diff .eq. AmD_by_Wilke) then
        do  j=1,N_SPEC
	        do  k=1,j
              D2(k,j)=exp(Trans%D_an(1,k,j)+Trans%D_an(2,k,j)*log(T) & 
	                   +Trans%D_an(3,k,j)*log(T)**2+Trans%D_an(4,k,j)*log(T)**3) *1.d-4        ! 1.d-4:  g-cm-s to SI
          enddo
        enddo
      elseif (Trans%Flag_diff .eq. AmD_by_Gupta) then
         do j = 1, N_SPEC
            do k = 1, j
              !  D2(k, j) = exp(Trans%D_an(1, k, j)*log(T)**3 + Trans%D_an(2, k, j)*log(T)**2 &
              !                 + Trans%D_an(3, k, j)*log(T) + Trans%D_an(4, k, j))*1.d-4        ! 1.d-4:  g-cm-s to SI
              D2(k, j) = 1.d-4 * exp(Trans%D_an(4, k, j)) * T**(Trans%D_an(1, k, j)*log(T)**2 &
                                +Trans%D_an(2, k, j)*log(T) + Trans%D_an(3, k, j))
            end do
         end do
      end if
	    do j=1,N_SPEC
	        do k=j+1,N_SPEC
	          D2(k,j)=D2(j,k)
	        enddo
	    enddo

      do k=1,N_SPEC
        ds=0.d0
        do j=1,N_SPEC
          if(j .ne. k) ds=ds+Xi(j)/D2(k,j)
        enddo
        AmDi(k)=(1.d0-Xi(k))/(ds+1.d-20) *atm/p                         ! Di ~ 1/p
      enddo
  endif 
  end subroutine comput_diffusion_coeff

!---------------------------------------------------------------------------

