! ������ѧ�йصļ���
! Ver 0.36, 2016-10-16, ����Eʱ�����˿�Tct���� 
! Ver 0.37a, Newton�������¶�
! Ver 1.3a  2016-11-20  �����˸�д ����ʽ�޸ģ�����ѧ�����㲢δ�޸ģ�
! Ver 1.4, 2022-2-21 1) ����ѧ��֧�ֶ����ϣ�ԭ�汾ֻ֧��2����ϣ���  2�� �¶ȳ����������ʱ��Cp�趨Ϊ����
!----------------------------------------------------------------------------
!  read_Spec( ) ;  
!    comput_E(di,T,E) ;    comput_Cp(d,di,T,Cp) ;  comput_T(di,E,T) ;  comput_hi(k,T,hi) 
!    IF_RealGas ==  Perfect_GAS  or Real_GAS  (Ĭ��ΪReal_GAS)
!    update_Et(di,dinew,E,T)
!-------------------------------------------------------------------------------

  ! �����������ѧ����
  subroutine read_Spec
    use CHEM
    implicit none
    TYPE (SPECIE),pointer:: Sp
 
    integer:: i,j,k
    allocate(Spec(N_Spec))
 
    !�����������  
    open(99,file="specie.in")       
    read(99,*)
    read(99,*)
    read(99,*)  N_Tct, (Tct(k),k=1,N_Tct-1), Tct_max     ! N_Tct:����ѧ���ֶ���ϵĶ�����ͨ��Ϊ2���� ת���¶�(����1000K)�� �������¶ȣ�������ֵCpΪ��ֵ��
    read(99,*)
    do k=1,N_Spec
      SP=>Spec(k)
      read(99,*)
      read(99,*)  Sp%name                                                     
      read(99,*)  Sp%Mi                                                ! ���Ħ������
	  do j=1,N_Tct
      read(99,*) SP%Ai(j),SP%Bi(j),SP%Ci(j),SP%Di(j),SP%Ei(j),SP%Fi(j),SP%Gi(j)            ! �������ѧ�������ֶ���ϣ���j�Σ�
      enddo
      Sp%Ri=R0/Sp%Mi
	  
!     Sp%Ect=Sp%F2-Sp%F1   ! ��Tct�����ܣ����ʣ���
      Sp%Ect(1)=0.d0         ! ��1����ϣ�������
      do j=1,N_Tct-1      ! ����϶ε��ʲ�  �����䵽������Ϲ�ʽ�У� ʹ�����ʿ���϶�������   ����������Sp%F2-Sp%F1 ֵ��ͬ�������ȸ���Щ��
	  Sp%Ect(j+1)= -( (Sp%Ai(j+1)-Sp%Ai(j))*Tct(j) +(Sp%Bi(j+1)-Sp%Bi(j))*Tct(j)**2/2.d0 & 
	            + (Sp%Ci(j+1)-Sp%Ci(j))*Tct(j)**3/3.d0 +(Sp%Di(j+1)-Sp%Di(j))*Tct(j)**4/4.d0 &
				+ (Sp%Ei(j+1)-Sp%Ei(j))*Tct(j)**5/5.d0 ) + Sp%Ect(j) 
 	  enddo
	enddo
    close(99)
    ! print*, "read specie.in OK"
  end subroutine read_Spec

!---------------------------------------------------------------------------------
!--------״̬������ؼ��� ���ܶȡ� ѹ���� �¶� ������ �� ���� �ȣ� -------------------
 ! ��ȫ���� or ��ʵ���� ��Cp �ǳ�����
 
 ! �������� ��������Ӧ�ʣ�
  subroutine comput_E(di,T,E)           
    use CHEM
    implicit none
    real(PRE_EC):: T,E,di(N_SPEC) 

    if(IF_RealGas .eq. Perfect_GAS) then
      call   comput_E_pgas(di,T,E)
    else
       call   comput_E_rgas(di,T,E)
   endif
   end

  
! ���㶨ѹ����   
  subroutine comput_Cp(d,di,T,Cp)
    use CHEM
    implicit none
    real(PRE_EC):: T,Cp,d,di(N_SPEC) 

    if(IF_RealGas .eq. Perfect_GAS) then
      call  comput_Cp_pgas(d,di,T,Cp)
    else
      call  comput_Cp_rgas(d,di,T,Cp)
   endif
   end

! ��֪���ܼ����¶�
  subroutine comput_T(di,E,T)
    use CHEM
    implicit none
    real(PRE_EC):: E,T,di(N_SPEC)

    if(IF_RealGas .eq. Perfect_GAS) then
      call  comput_T_pgas(di,E,T)
    else
      call  comput_T_rgas(di,E,T)
   endif
   end

 !  �����
    subroutine comput_hi(k,T,hi)
    use CHEM
    implicit none
    integer k
    real(PRE_EC):: T,hi 
    if(IF_RealGas .eq. Perfect_GAS) then
      call comput_hi_pgas(k,T,hi)
    else
      call comput_hi_rgas(k,T,hi)
   endif
   end


!----------------��ȫ���� Cv=5/2 R, Cp=7/2 R --------------------

!  Perfect gas  ��ȫ�������  (Cv=5/2 R )
  subroutine comput_E_pgas(di,T,E)            ! ��ȫ����
    use CHEM
    implicit none
    integer k
    real(PRE_EC):: T,E,di(N_SPEC) 
     E=0.d0
     do k=1,N_SPEC
     E=E+2.5d0*SPEC(k)%Ri*di(k)*T     ! Cv=2.5R ,  E=rho*Cv*T
     enddo
   end

! ���㶨ѹ���� Cp   
  subroutine comput_Cp_pgas(d,di,T,Cp)
    use CHEM
    implicit none
    integer k
     real(PRE_EC):: T,Cp,d,di(N_SPEC) 
     Cp=0.d0
     do k=1,N_SPEC
      Cp=Cp+3.5d0*SPEC(k)%Ri*di(k)        ! Cp=3.5R 
     enddo
     Cp=Cp/d                 ! �ܶȼ�Ȩƽ��
   end

  subroutine comput_T_pgas(di,E,T)
    use CHEM
    implicit none
    real(PRE_EC):: E,T,di(N_SPEC),rCv
    integer k
    rCv=0.d0
     do k=1,N_SPEC
     rCv=rCv+2.5d0*SPEC(k)%Ri*di(k)      ! Cv=2.5R ,  E=rho*Cv*T
     enddo
     T=E/rCv
    end
 
  
! �����¶ȣ����������, ������Ӧ��  
  subroutine comput_hi_pgas(k,T,hi)
    use CHEM
    implicit none
    integer k
    real(PRE_EC):: T,hi 
      hi=3.5d0* SPEC(k)%Ri*T              ! Cp=7/2*R
   end

  

!---------��ѧ��Ӧ��������ܱ仯---------------------
!----------���ݷ�Ӧ�ʸ������� ���ѿ����˿�Tct������-------------
  subroutine update_Et(di,dinew,E,T)
    use CHEM
    implicit none
    real(PRE_EC) :: E,T,F0        ! ���ܣ�������Ӧ�ʣ��� �¶ȣ� F0 ��Ӧ��
    real(PRE_EC),dimension(N_Spec) :: di, dinew, bi      ! �ɵ�����ܶȡ��µ�����ܶ�
    integer :: i,j,k
 
    do i=1,N_SPEC
!      F0=Spec(i)%F1               ! �ѿ����˿�Tct����
      F0=Spec(i)%Fi(1)             ! �ѿ����˿�Tct����, ֻʹ�õ�1���������Fiֵ  ���������ε���ͨ��Ect�����������������У�
      E=E-(dinew(i)-di(i))*SPEC(i)%Ri*F0    ! �������� �����뻯ѧ��Ӧ�ͷ�/���յ�������
    end do
  end
  
   
! -------------��ʵ���� Cp= Cp (T) ----------------------------------
! �����¶ȼ�����ܶȣ��������ܣ�������ѧ�ܣ�  (E=Rho*e=Rho*Cv*T for perfect gas)
!        ���ܵ���Ϲ�ʽ Et=A*T+1/2*B*T**2+1/3*C*T**3+1/4*D*T**4+1/5*E*T**5  ;  
!        �����������ʣ� ����T> Tctʱ������ ��Ect=F2-F1)    
  subroutine comput_E_rgas(di,T,E)
    use CHEM
    implicit none
    integer k,j,jc
    real(PRE_EC):: T,E,R1,di(N_SPEC) 
    real(PRE_EC):: A1,B1,C1,D1,E1,Ect, Cv2
    A1=0.d0; B1=0.d0; C1=0.d0; D1=0.d0; E1=0.d0 ;  Ect=0.d0
	jc=1
	do j=1,N_Tct-1
     if(T>Tct(j)) jc=j+1   !jc �¶�T����������������
	enddo
	do k=1,N_SPEC
      R1=SPEC(k)%Ri
      A1=A1+(SPEC(k)%Ai(jc)-1.d0)*R1*di(k)          ! ϵ��������ּ�Ȩ��� (Ai-1) for Cv
      B1=B1+SPEC(k)%Bi(jc)*R1*di(k)
      C1=C1+SPEC(k)%Ci(jc)*R1*di(k)
      D1=D1+SPEC(k)%Di(jc)*R1*di(k)
      E1=E1+SPEC(k)%Ei(jc)*R1*di(k)
      Ect=Ect+SPEC(k)%Ect(jc)*R1*di(k)
	 end do

	 if(T<=Tct_max) then
	  E=A1*T+B1*T**2/2.d0+C1*T**3/3.d0+D1*T**4/4.d0+E1*T**5/5.d0 + Ect     ! Ect: ��Tct����
     else 
      Cv2=A1+B1*Tct_max+C1*Tct_max**2+D1*Tct_max**3+E1*Tct_max**4               ! Cv at Tct_max
	  E=A1*Tct_max+B1*Tct_max**2/2.d0+C1*Tct_max**3/3.d0+D1*Tct_max**4/4.d0+E1*Tct_max**5/5.d0 &
 	            + Ect +Cv2*(T-Tct_max)     ! ���� ��T> Tct_max), ����Tct����
	 endif
  end
!--------------------------------------------------------------------------  
! �����¶ȼ�����ܶȣ��������Cp   (�������ƽ�����ȣ�
  subroutine comput_Cp_rgas(d,di,T,Cp)
    use CHEM
    implicit none
    integer k,j,jc
    real(PRE_EC):: T,Cp,R1,d,di(N_SPEC) 
    real(PRE_EC):: A1,B1,C1,D1,E1
    
    A1=0.d0; B1=0.d0; C1=0.d0; D1=0.d0; E1=0.d0
	jc=1
	do j=1,N_Tct-1
     if(T>Tct(j)) jc=j+1
	enddo
    do k=1,N_SPEC
      R1=SPEC(k)%Ri
      A1=A1+(SPEC(k)%Ai(jc)-1.d0)*R1*di(k)          ! ϵ��������ּ�Ȩ��� (Ai-1) for Cv
      B1=B1+SPEC(k)%Bi(jc)*R1*di(k)
      C1=C1+SPEC(k)%Ci(jc)*R1*di(k)
      D1=D1+SPEC(k)%Di(jc)*R1*di(k)
      E1=E1+SPEC(k)%Ei(jc)*R1*di(k)
	 end do
    
    if(T <= Tct_max ) then
	  Cp=(A1+B1*T+C1*T**2+D1*T**3+E1*T**4)/d 
    else 
	  Cp=(A1+B1*Tct_max+C1*Tct_max**2+D1*Tct_max**3+E1*Tct_max**4)/d
    endif

  end

!---------------------------------------------------------------------------------
! �����¶ȣ���������� �����ʣ�, ������Ӧ��  
!  ���ʵ���Ϲ�ʽ ht=A*T+1/2*B*T**2+1/3*C*T**3+1/4*D*T**4+1/5*E*T**5     
	subroutine comput_hi_rgas(k,T,hi)
    use CHEM
    implicit none
    integer k,j,jc
    real(PRE_EC):: T,hi,Cp_max 
	TYPE (SPECIE),pointer:: Sp
    Sp=>SPEC(k)    ! �� i�����
	jc=1
	do j=1,N_Tct-1  ! ����T�������¶�����
     if(T>Tct(j)) jc=j+1
	enddo
    
	if(T<=Tct_max) then
      hi=(Sp%Ai(jc)*T+Sp%Bi(jc)*T**2/2.d0+Sp%Ci(jc)*T**3/3.d0+   &
	      Sp%Di(jc)*T**4/4.d0+Sp%Ei(jc)*T**5/5.d0 + Sp%Ect(jc))*Sp%Ri       ! Ect ��Tct����
    else 
	  Cp_max=Sp%Ai(jc)+Sp%Bi(jc)*Tct_max+Sp%Ci(jc)*Tct_max**2+Sp%Di(jc)*Tct_max**3 &
	         +Sp%Ei(jc)*Tct_max**4         ! Cp at Tct_max
      hi=(Sp%Ai(jc)*Tct_max+Sp%Bi(jc)*Tct_max**2/2.d0+Sp%Ci(jc)*Tct_max**3/3.d0+   &
	      Sp%Di(jc)*Tct_max**4/4.d0+Sp%Ei(jc)*Tct_max**5/5.d0 + Sp%Ect(jc) + Cp_max*(T-Tct_max) )*Sp%Ri 			 
   
    endif 
  end
  
  ! �������ֵ�Gibbs ������ (��������������ʹ�ÿ�Ect������ G = H-TS   ;       ���㻯ѧ��Ӧƽ�ⳣ��ʱʹ��
  subroutine  comput_gibbs(k,T,gi)         ! gi=Gi/Ri
    use CHEM
    implicit none
	TYPE (SPECIE),pointer:: Sp
    integer:: k,j,jc
    real(PRE_EC):: T, gi, Cp_max,h_max,S_max 
     SP=>Spec(k)    ! �� i�����
	 jc=1
	 do j=1,N_Tct-1       ! ����T�������¶�����
     if(T>Tct(j)) jc=j+1
	 enddo	 
     if(T<=Tct_max) then
       gi=Sp%Ai(jc)*(1.d0-log(T))*T-SP%Bi(jc)*T**2/2.d0-SP%Ci(jc)*T**3/6.d0-SP%Di(jc)*T**4/12.d0 & 
	     -SP%Ei(jc)*T**5/20.d0+SP%Fi(jc)-SP%Gi(jc)*T
      else
	  Cp_max=Sp%Ai(jc)+Sp%Bi(jc)*Tct_max+Sp%Ci(jc)*Tct_max**2+Sp%Di(jc)*Tct_max**3+Sp%Ei(jc)*Tct_max**4         ! Cp at Tct_max	  
      h_max= Sp%Ai(jc)*Tct_max+Sp%Bi(jc)*Tct_max**2/2.d0+Sp%Ci(jc)*Tct_max**3/3.d0    &
	         +Sp%Di(jc)*Tct_max**4/4.d0+Sp%Ei(jc)*Tct_max**5/5.d0 + Sp%Fi(jc)          ! �����Tct���� ������+�����ʣ�
	  S_max=Sp%Ai(jc)*log(Tct_max)+SP%Bi(jc)*Tct_max+SP%Ci(jc)*Tct_max**2/2.d0+SP%Di(jc)*Tct_max**3/3.d0 &
	         +SP%Ei(jc)*Tct_max**4/4.d0+SP%Gi(jc)
	  gi=h_max+Cp_max*(T-Tct_max)-T*(S_max+Cp_max*log(T/Tct_max))  ! g=h-TS
	  end if
    
	end       
!--------------------------------------------------------------------------  

!  ��������E ��������ѧ�ܣ�������ܶ�di,������¶�
!   Newton�� �����¶� 
  
!  di : ����ܶȣ� E ���� �����������ʣ��� T �¶� 
!  �����ټ��㣬 ���ʵ�λ�ƣ�
  subroutine comput_T_rgas(di,E,T)
    use CHEM
    implicit none
    real(PRE_EC),parameter:: T_Cr=1.d-3        ! ��������
    integer,parameter:: Kstep_lmt=30           ! ������������
    real(PRE_EC):: E,T,p,di(N_SPEC)
    real(PRE_EC),dimension(NTct_Max):: A1,B1,C1,D1,E1,Ect
    real(PRE_EC):: R1,T1,T2,Et1,Ex,error
    integer:: k,ks,j,jc
 
    A1(1:N_Tct)=0.d0; B1(1:N_Tct)=0.d0; C1(1:N_Tct)=0.d0; D1(1:N_Tct)=0.d0; E1(1:N_Tct)=0.d0
    Ect(1:N_Tct)=0.d0
   do j=1,N_Tct   
   do k=1,N_SPEC
    !      ���ܵ���Ϲ�ʽ (Et=A*T+B*T**2+C*T**3+D*T**4+E*T**5) �� �����������ʣ�
      R1=SPEC(k)%Ri
      A1(j)=A1(j)+(SPEC(k)%Ai(j)-1.d0)*R1*di(k)          ! ϵ��������ּ�Ȩ��ͣ���������
      B1(j)=B1(j)+SPEC(k)%Bi(j)*R1*di(k)
      C1(j)=C1(j)+SPEC(k)%Ci(j)*R1*di(k)
      D1(j)=D1(j)+SPEC(k)%Di(j)*R1*di(k)
      E1(j)=E1(j)+SPEC(k)%Ei(j)*R1*di(k)
      Ect(j)=Ect(j)+SPEC(k)%Ect(j)*R1*di(k)           ! ��-�����ʲ� (������Ϲ�ʽ֮��Ĳ��죬 �Ա��ԽTctʱ����������
   enddo
   enddo
      ! newton  iteration  
  
    T1=E/A1(1)            ! ��ֵ
  
     ks=0 ! ��������

     do while (ks <=Kstep_lmt)

	 jc=1
	 do j=1,N_Tct-1       ! ����T1�������¶�����
     if(T1>Tct(j)) jc=j+1
	 enddo	 
	 

    if(T1<=Tct_max) then                               ! ������ (>1000K)
      Et1=A1(jc)*T1+B1(jc)*T1**2/2.0+C1(jc)*T1**3/3.d0+D1(jc)*T1**4/4.d0+E1(jc)*T1**5/5.d0 + Ect(jc)              ! ���� (���� Tct����) 
      Ex=A1(jc)+B1(jc)*T1+C1(jc)*T1**2+D1(jc)*T1**3+E1(jc)*T1**4       ! ���� (Cv)
    else               ! T> Tct_max 
	  Ex=A1(jc)+B1(jc)*Tct_max+C1(jc)*Tct_max**2+D1(jc)*Tct_max**3+E1(jc)*Tct_max**4               ! Cv at Tct_max  
	  Et1=A1(jc)*Tct_max+B1(jc)*Tct_max**2/2.d0+C1(jc)*Tct_max**3/3.d0+D1(jc)*Tct_max**4/4.d0  &
	     +E1(jc)*Tct_max**5/5.d0 + Ect(jc) +Ex*(T1-Tct_max)     ! ���� (���� Tct����) 
	end if
	
     T2=T1-(Et1-E)/Ex           !  Newton iteration
     ks=ks+1
     error=abs(T2-T1)
	 T1=T2    ! ����
    if( error < T_Cr  ) exit
    enddo 
   T=T2
   
  end

! �����¶ȡ������ѹ�� ������ѹ֮�ͣ��� di(k) ����ֵ��ܶ�
  subroutine comput_P(di,T,p)
    use CHEM
    implicit none
    real(PRE_EC):: T,p,di(N_SPEC)
    integer:: k

    p=0.d0 
    do k=1,N_SPEC
      p=p+SPEC(k)%Ri*di(k)*T          ! ѹ������ѹ֮��
    end do
  end

! ����ѹ�����¶��Լ������ȷ�ai(:) �������ܶ� ��2016-5-11��
 
  subroutine comput_d(ai,T,p,d)
    use CHEM
    implicit none
    real(PRE_EC):: T,p,d,ai(N_SPEC),R1
	integer:: k
    R1=0.d0
	do k=1,N_SPEC
	 R1=R1+Spec(k)%Ri*ai(k)
	enddo
	d=p/(T*R1)
   end

! ------------Comput sound speed (approximate), 2022-4-24 --------------------------------------
  subroutine comput_soundspeed(c,p,d,Et)
    use CHEM
    implicit none
    real(PRE_EC):: c,p,d,Et,gamma  
    gamma=p/Et+1.d0            ! ��Ч���ȱ�
    c=sqrt(gamma*p/d)          ! �������� ��ͨ������ʹ�ã�
 end


