! �뻯ѧ��Ӧ�йص��ӳ���
! Ver 1.3, 2016-11-20  ���¸�д;   �����˸�ѹ��Ӧ�� ˫Arrhenius��ʽ��Ӧ��  Aʹ��mol/cm3��Ϊ��λ
! Ver 1.3a, 2017-2-28, �޸���  Reaction_rate( ) �е�Bug  (��������Ũ�� Ctrb ��
! Ver 1.5, 2022-2-24, ֧������ʽ
! Ver 1.5a, 2023-8-15, ��ѹ��Ӧ�����޸���Bug;  log(Fc) ==> log10(Fc)
! Ver 1.6, 2024-12-18,   ֧�� Park (1990) �淴Ӧ����ģ�ͣ�ʹ��5�������ƽ�ⳣ����Ref: Passiatore ��2022�� JFM 941, A21, doi:10.1017/jfm.2022.283�� 
!          ReReac=0 ������Gibbs�����ܵ�ƽ�ⳣ���� ReReac=1: ʹ��Arrhenius��ʽ�� ReReac=2: ʹ��Park (1990) �������ϵ�ƽ�ⳣ��ģ��
! --------------------------------------------------------
!  ���뻯ѧ��Ӧ�����йصĵ�λ���� mol-cm�� A*b**T*exp(-E/T)�е�A
!  ������Ľӿ��Բ��ù��ʵ�λ�� :  kg-m-s-K-mol
! --------------- using mol-cm  units-------------------   


! ��ȡ��ѧ��Ӧ���� �� A ʹ�� mol-cm ��λ;   b, E ����KΪ��λ��  K=A*T**b*exp(-E/T)��
  subroutine read_Reac
    use CHEM
    implicit none
    integer:: i,j
    TYPE (REACTION), pointer:: RC

!    print*, "--------------Read Reaction.in --------------------------------------"
    
    open(100,file="Reaction.in")                          ! ��ѧ��Ӧ����
    read(100,*)
    read(100,*) 
! N_spec �����Ŀ��  N_REAC ��ѧ��Ӧ��Ŀ;  ReReac �淴Ӧ�����㷨��0 ���û���Gibbs�����ܵ�ƽ�ⳣ��, 1 ����Arrhenius��ʽ��2 Park 1990 ƽ�ⳣ��ģ�� ��
! CHEM_TimeAdv ��ѧ��Ӧʱ���ƽ�������1�� 1��Euler; 2�� 2��RK�� 3�� 3��RK (Ŀǰ�汾�в�֧��)�� -1�� 1������ʽ

    read(100,*) N_Spec, N_REAC , ReReac     
    allocate(REAC(N_REAC))
    allocate(Vf(N_SPEC,N_REAC), Vr(N_SPEC,N_REAC) )          ! ��Ӧϵ������ ��ά���������*��Ӧ����
    read(100,*)
    read(100,*)                                       
	read(100,*)
    read(100,*)
    !-------------------------------------------------------------------------------                                           
  
    do j=1,N_REAC                                         ! ���뷴Ӧ��Ϣ
    !  print*, "-------Reaction --", j
      Rc=>REAC(j)
      read(100,*)
      read(100,*) 
	  read(100,*)  Rc%Af_type,  Rc%TrReac , Rc%Fc_troe      ! ����Ӧϵ�����ͣ�  �Ƿ����巴Ӧ ; Troeϵ�� ����ѹ��Ӧʹ�ã�
      
	  read(100,*) Rc%Af(1),Rc%Bf(1),Rc%Ef(1)                 ! ����Ӧϵ��
      if(Rc%Af_type .ne.  Af_nomal )     read(100,*)  Rc%Af(2), Rc%Bf(2), Rc%Ef(2)       ! ��2������Ӧϵ�� �����ڸ�ѹ��Ӧ �� ˫Arrheniusϵ��������
	   if( ReReac==0) read(100,*)   ! ��ƽ�ⳣ�������淴Ӧ����ʹ�ò����� ��ȡ���У�����Reaction.in��ͨ����
	   if( ReReac==1)  read(100,*) Rc%Ar, Rc%Br, Rc%Er                     ! Arrhenius ��ʽ���淴Ӧ����ϵ�� 
       if( ReReac==2)  read(100,*) Rc%A1r, Rc%A2r, Rc%A3r, Rc%A4r, Rc%A5r  ! Park ƽ�ⳣ��ģ���е�5�������� �������淴Ӧ����
	  
      read(100,*) (vf(i,j), i=1,N_spec)    ! ����Ӧϵ������
      read(100,*) (vr(i,j), i=1,N_spec)    ! �淴Ӧϵ������
      
	  if(Rc%TrReac .eq. 1) then
        allocate(Rc%TrEff(N_SPEC))                    ! ����Ӱ��ϵ��
        read(100,*) (Rc%TrEff(i), i=1,N_Spec)
      end if
      
	  Rc%sgm=0
      do i=1,N_spec
        Rc%sgm=Rc%sgm+vr(i,j)-vf(i,j)      ! ��-�� ��Ӧ�ļ���
      end do
    end do
    close(100)

    ! print*, "read Reaction.in OK"
    ! print*, "N_spec=", N_spec, "N_Reac=", N_Reac    

  end subroutine read_Reac





! ���㻯ѧƽ�ⳣ���� TΪ�¶�
  subroutine comput_KjX(T,Kjx)
    use CHEM
    implicit none
    TYPE (SPECIE),pointer:: Sp
    integer:: i,j
    real(PRE_EC) :: T,Kjx(N_Reac),dgi,gi(N_spec)
	real(PRE_EC),parameter:: Rs=R0/atm * 1.d6                          ! Rs,  ��λ cm3/(mol.K)  ���ο��ֲ� 1.22 ʽ

    ! �������ֵ�Gibbs ������

    do i=1,N_spec
      call comput_gibbs(i,T,gi(i))
!      SP=>Spec(i)    ! �� i�����
!      if(T<Tct) then
!         gi(i)=Sp%A1*(1.d0-log(T))*T-SP%B1*T**2/2.d0-SP%C1*T**3/6.d0-SP%D1*T**4/12.d0-SP%E1*T**5/20.d0+SP%F1-SP%G1*T
!       else
!        gi(i)=Sp%A2*(1.d0-log(T))*T-SP%B2*T**2/2.d0-SP%C2*T**3/6.d0-SP%D2*T**4/12.d0-SP%E2*T**5/20.d0+SP%F2-SP%G2*T
!      end if
    end do

    do j=1,N_Reac
      dgi=0.d0
      do i=1,N_Spec
        dgi=dgi+gi(i)*(Vr(i,j)-Vf(i,j))       ! �����ʲ�
       end do
      Kjx(j)=exp(-dgi/T) * (Rs*T)**(-Reac(j)%sgm)                    ! �������ֲ� 1.22 ʽ ; sgm ��-����Ӧ�ļ���
    end do
  end

!-----------------------------------------------------------------------
function Arrhenius(A,b,E,T)          ! Arrhenius ��ʽ
  	use Precision_EC	
	implicit none
  	real(PRE_EC)::  A,b,E,T,Arrhenius
  	Arrhenius=A*T**b*exp(-E/T)          ! Eʹ�� K ��Ϊ��λ��  
end

function ParkCurveFit(A1,A2,A3,A4,A5,T)       ! Park 5������Ϲ�ʽ�� ����ƽ�ⳣ��Kjx
    use Precision_EC
	implicit none
    real(PRE_EC)::  A1,A2,A3,A4,A5,T,ParkCurveFit,z
	z = 10000.d0/T
    ParkCurveFit = exp(A1 + A2*z + A3*(z**2) + A4*(z**3) + A5*(z**4))            
end
!---------------------------------------------------------------------




!   ���㻯ѧ��Ӧ����ϵ��  (��λ��mol-cm-K-s ��)
!   ע�⣬ �������ӳ���λ��ͬ�������� cm��
!   ����ϵ���а����������巴Ӧ�йص���
!   ��ѹ��Ӧ�е�log(Fc) �޸�Ϊ log10(Fc)  (�����е� log() ָ����log10(), ������ ln())
!   ����� Park �淴Ӧ����ģ�� ��2024-12-19��
	subroutine Reaction_rate_coef(ci,wf,wr,T)  ! wf, wr ����Ӧ���淴Ӧ������ϵ��
    use CHEM
	
 	implicit none
	real(PRE_EC),dimension(N_Reac):: wf,wr, Kjx  ! ��Ӧ����ϵ��;  ��ѧ��Ӧƽ�ⳣ�� (mol-cm-K unit)
	real(PRE_EC),dimension(N_Spec):: ci     ! Ħ��Ũ�� (mol/cm3) 
    real(PRE_EC):: T,Ctrb,Arrhenius ,ParkCurveFit, wf1,wf2, Prs, F,Fc, ac,an       !��ѧ��Ӧ�����й���
    TYPE (REACTION), pointer:: RC

	integer:: i,j,k
	! ����Gibbs�����ܣ����㻯ѧ��Ӧƽ�ⳣ��Kjx ��ReReac=1 ʹ��Arrhenius��ʽ�����淴Ӧ����,����ʹ��kjx;  ReReac=2 ����Park ��ʽ����kjx��
        if(ReReac ==0) then        
        	call comput_KjX(T,Kjx)
        endif

! ----------  ���㻯ѧ��Ӧ���ʳ��� wf, wr ----------------
		do j=1,N_REAC
	 		Rc=>REAC(j)
	 
			if(Rc%TrReac .eq. 1) then            ! Three-body reaction
        		Ctrb=0.d0
				do i=1,N_SPEC
        	   		Ctrb=Ctrb+ci(i)*Rc%TrEff(i)             ! ����Ũ�� ������ּ�Ȩ�ͣ�Rc%TrEff(:)ΪӰ������ ��
			 	enddo
        	endif

!     ����Ӧ���ʳ���
		  	wf1=Arrhenius(Rc%Af(1), Rc%bf(1), Rc%Ef(1), T)            ! ����Ӧ�ٶ�ϵ��
			if(  Rc%Af_type .eq.  Af_nomal)  then
				wf(j)=wf1                                                                        ! ���淴Ӧ
        	else if (Rc%Af_type .eq. Af_DualArrhenius) then                ! ˫Arrhenius��ʽ����
				wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
				wf(j)=wf1+wf2
        	else if (Rc%Af_type .eq. Af_Highpress) then                        !  ��ѹ��Ӧ��  Fall-off �����巴Ӧ�� �� H+O2 (+M) = HO2 (+M)
				wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
				Prs=wf1/wf2*Ctrb                          ! Ctrb=[M]
        		Fc=Rc%Fc_troe                              ! ��ѹ��Ӧ��Troe��ʽ�� ��Chemkin�����ֲ�
!				ac=-0.4d0-0.67*log(Fc)
!        		an=0.75d0-1.27*log(Fc)
!        		F=exp(log(Fc)/( 1.d0+ ( (log(Prs)+ac)/(an-0.14d0*(log(Prs)+ac) ) )**2) )
				ac=-0.4d0-0.67*log10(Fc)
        		an=0.75d0-1.27*log10(Fc)
        		F=10.d0**(log10(Fc)/( 1.d0+ ( (log10(Prs)+ac)/(an-0.14d0*(log10(Prs)+ac) ) )**2) )

!               wf(j)=wf2*Prs/(1.d0+Prs)                 ! Lindemann type
	    		wf(j)=wf2*Prs/(1.d0+Prs)*F               ! Troe type
        	endif

!     �����淴Ӧ���� ��3��ģ�ͣ�
!          2024-12-19, ����Park 1990 ģ��
            if(ReReac ==0 ) then 
			   wr(j)=wf(j)/kjx(j)           ! ����Gibbs������ƽ�ⳣ��ģ�� 
			else if (ReReac ==1) then 
			   wr(j)=Arrhenius(Rc%Ar,Rc%br,Rc%Er,T)            ! ����Arrhenius��ʽ���淴Ӧ����ģ��
            else if (ReReac ==2) then
			  kjx(j)=ParkCurveFit(Rc%A1r,Rc%A2r,Rc%A3r,Rc%A4r,Rc%A5r,T)   ! Park  ģ��, ����ƽ�ⳣ�� 
	          if(kjx(j) < 1.d-30) kjx(j)=1.d-30                ! �趨���޷�ֹ��ĸΪ0��1.d-30 ����Ϊ�ԡ� 2024-12-19  		  
			  wr(j)=wf(j)/kjx(j)
!              wr(j)=ParkCurveFit(Rc%A1r,Rc%A2r,Rc%A3r,Rc%A4r,Rc%A5r,T)     ! wrong test 
			endif 

		 
	    	if(Rc%TrReac .eq. 1  .and. Rc%Af_type .ne. Af_Highpress ) then       !  ���巴Ӧ ����Fall-off�ͣ� Fall-off�͵����ʳ����а���������Ũ��Ctrb)
        		wf(j)=wf(j)*Ctrb
				wr(j)=wr(j)*Ctrb
        	endif
	  	enddo 
	end subroutine Reaction_rate_coef
	  


!     ���㷴Ӧ��������  (mol/(cm3.s))
!     ��λ mol-cm-K   ��ע�⣬ �������ӳ���λ��ͬ�������� cm��)	    
	subroutine Reaction_rate(ci,ri,T)  
    use CHEM
 	implicit none
	real(PRE_EC),dimension(N_Reac):: wf,wr  ! ��Ӧ����ϵ�� (mol-cm-K unit)
	real(PRE_EC),dimension(N_Spec):: ci,ri     ! Ħ��Ũ�� (mol/cm3) , ��Ӧ�������� (mol/(cm3.s))
    real(PRE_EC):: T,wf1,wr1
	integer:: i,j
!   ���㷴Ӧ����ϵ�� wf, wr   ( ��λ mol-s-K-cm �� )
    	call Reaction_rate_coef(ci,wf,wr,T)
!   ���㻯ѧ��Ӧ���� wf1, wr1  	  ��λ�� mol/(cm3-s)
	 	ri(:)=0.d0 
	 	do j=1,N_REAC
	 	  	wf1=wf(j)        ! ����Ӧ���� ����j����Ӧ��
	 	  	wr1=wr(j)        ! �淴Ӧ���ʣ���j����Ӧ��
			do i=1,N_SPEC
				if(vf(i,j) .ne. 0) wf1=wf1*ci(i)**vf(i,j)
				if(vr(i,j) .ne. 0) wr1=wr1*ci(i)**vr(i,j)
			enddo
	 	    do i=1,N_SPEC         ! ���
     	    	ri(i)=ri(i)+ (wf1-wr1) *(vr(i,j)-vf(i,j))   ! ��i����ֵ������� (Ħ��/(��������.��))        
			enddo
     	enddo
    end subroutine Reaction_rate

!---��ѧ��Ӧ����,���䵼������ Arate(i,j)=  d(ri)/d(cj)          ( Ħ��Ũ�������ʶ�Ħ��Ũ�ȵĵ����� û�г���Mi/Mj)
  !     ���㷴Ӧ��������  (mol/(cm3.s))
!     ��λ mol-cm-K   ��ע�⣬ �������ӳ���λ��ͬ�������� cm��)	    
	   subroutine Reaction_rate_and_Matrix(ci,ri,Arate,T)  
       use CHEM
 	   implicit none
	   real(PRE_EC),dimension(N_Reac):: wf,wr  ! ��Ӧ����ϵ�� (mol-cm-K unit)
	   real(PRE_EC),dimension(N_Spec):: ci,ri       ! ci Ħ��Ũ�� (mol/cm3); ri Ħ��Ũ��������   
	   real(PRE_EC),dimension(N_Spec,N_Spec):: Arate        ! ��������
       real(PRE_EC):: T,wf1,wr1,Awf,Awr  !Awf=d(eta)/(cj)  ����Ӧ���ʶ�cj�ĵ����� Awr �淴Ӧ���ʶ�cj���� 
	   integer:: i,j,k
	   real(PRE_EC),parameter::  Aepsl=1.d-20
	   
!   ���㷴Ӧ����ϵ�� wf, wr   ( ��λ mol-s-K-cm �� )
      call Reaction_rate_coef(ci,wf,wr,T)
	  
!   ���㻯ѧ��Ӧ���� wf1, wr1  	  ��λ�� mol/(cm3-s)
	   ri(:)=0.d0 	
       Arate(:,:)=0.d0 
	   
	 do k=1,N_REAC
 	  
	   wf1=wf(k)
	   wr1=wr(k)
	   do i=1,N_SPEC
		  if(vf(i,k) .ne. 0) wf1=wf1*ci(i)**vf(i,k)
		  if(vr(i,k) .ne. 0) wr1=wr1*ci(i)**vr(i,k)
	   enddo
	   
  !  ��Ӧ�������� 
	   do i=1,N_SPEC         ! ���
        ri(i)=ri(i)+ (wf1-wr1) *(vr(i,k)-vf(i,k))   ! ��i����ֵ������� (Ħ��/(��������.��))        
       enddo
		
 ! ��������		
	   do j=1,N_SPEC                 ! d(eta_i)/d(c_j) 	 
        Awf=vf(j,k)*wf1/(ci(j)+Aepsl)          ! +Aepsl ��ֹ��ĸΪ0
        Awr=vr(j,k)*wr1/(ci(j)+Aepsl)          ! +Aepsl ��ֹ��ĸΪ0
       do i=1,N_SPEC 
	    Arate(i,j)=Arate(i,j)+(Awf-Awr) *(vr(i,k)-vf(i,k)) 
       enddo 
	  enddo
	enddo
!---------------------------------------------------------
  
   end  


!----------���Ի�ѧ��Ӧ����---------------------------------------------------------------------------------------------
	   subroutine test_Reaction_rate
       use CHEM
 	   implicit none
	   real(PRE_EC),dimension(N_Reac):: wf,wr, Kjx  ! ��Ӧ����ϵ��;  ��ѧ��Ӧƽ�ⳣ�� (mol-cm-K unit)
	   real(PRE_EC),dimension(N_Spec):: gi       
       real(PRE_EC):: p1, T,Ctrb,Arrhenius , wf1,wf2, Prs , Fc,F, ac,an      !��ѧ��Ӧ�����й���
       TYPE (REACTION), pointer:: RC
	   TYPE (SPECIE),pointer:: Sp
	   integer:: i,j,k

      print*, " ------test chemical reation rates ------"
	  print*, "please input p (in atm;  used for fall-off reactions) "
	  read(*,*) p1


  open(99,file="gibbs.dat")
  do k=300, 2000
    T=k*1.d0
    do i=1,N_spec
	call comput_gibbs(i,T,gi(i))          ! ����Gibbs���� g=h-TS
!      SP=>Spec(i)    ! �� i�����
!      if(T<Tct) then
!        gi(i)=(Sp%A1*(1.d0-log(T))*T-SP%B1*T**2/2.d0-SP%C1*T**3/6.d0-SP%D1*T**4/12.d0-SP%E1*T**5/20.d0+SP%F1-SP%G1*T)
!      else
!        gi(i)=(Sp%A2*(1.d0-log(T))*T-SP%B2*T**2/2.d0-SP%C2*T**3/6.d0-SP%D2*T**4/12.d0-SP%E2*T**5/20.d0+SP%F2-SP%G2*T)
!      end if
    end do

    write(99,"(30E16.8)") T, (gi(j), j=1,N_Spec)     ! gi ��λΪ R0
   enddo


	  open(99,file="Kjx.dat")
	  open(100,file="Reaction-wf.dat")
	  open(101,file="Reaction-wr.dat")

     do k=300, 2000
		   T=k*1.d0
		   ctrb= p1*atm/(R0*T) * 1.d-6     ! [M] in mol/cm3

           call comput_KjX(T,Kjx)
		   write(99,"(30E16.8)") T, (KjX(j), j=1,N_Spec)     ! gi ��λΪ R0
      

! ----------  ���㻯ѧ��Ӧ���ʳ��� wf, wr ----------------
	do j=1,N_REAC
	 Rc=>REAC(j)

!     ����Ӧ���ʳ���
		  wf1=Arrhenius(Rc%Af(1), Rc%bf(1), Rc%Ef(1), T)            ! ����Ӧ�ٶ�ϵ��
		  if(  Rc%Af_type .eq.  Af_nomal)  then
		     wf(j)=wf1                                                                        ! ���淴Ӧ
          else if (Rc%Af_type .eq. Af_DualArrhenius) then                ! ˫Arrhenius��ʽ����
		    wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
			wf(j)=wf1+wf2
         else if (Rc%Af_type .eq. Af_Highpress) then                        !  ��ѹ��Ӧ��  Fall-off �����巴Ӧ�� �� H+O2 (+M) = HO2 (+M)
		    wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
			Prs=wf1/wf2*Ctrb                          ! Ctrb=[M]
!			wf(j)=wf2*Prs/(1.d0+Prs)               ! Lindemann type
!  	         ac=-0.4d0-0.67*log(Fc)
!            an=0.75d0-1.27*log(Fc)
!            F=exp(log(Fc)/( 1.d0+ ( (log(Prs)+ac)/(an-0.14d0*(log(Prs)+ac) ) )**2) )
			 ac=-0.4d0-0.67*log10(Fc)
        	 an=0.75d0-1.27*log10(Fc)
        	 F=10.d0**(log10(Fc)/( 1.d0+ ( (log10(Prs)+ac)/(an-0.14d0*(log10(Prs)+ac) ) )**2) )			
	        wf(j)=wf2*Prs/(1.d0+Prs)*F               ! Troe type

	    endif
!     �淴Ӧ���ʳ���
		 if(ReReac ==1) then
		   wr(j)=Arrhenius(Rc%Ar,Rc%br,Rc%Er,T)            ! ����Ӧ�ٶ�ϵ��
         else
           wr(j)=wf(j)/kjx(j)      
         endif
 enddo

     write(100,"(30E16.8)") T, (wf(j), j=1,N_Spec)     ! gi ��λΪ R0
	 write(101,"(30E16.8)") T, (wr(j), j=1,N_Spec)     ! gi ��λΪ R0

enddo
 
 end
 
!--------------------------------------------------------------------
!   Time Advance for Chemical reaction terms  ( 1st Euler; 2nd RK;  1st Euler-implicit) 
!   ��ѧ��Ӧ���ʱ���ƽ�


!-------------------------------------------------------------------
! ���㻯ѧ��Ӧ��ƽ� dt ʱ�䲽 
!    ���������di ������ܶ� (Kg/m3)�� E ���ܣ�������Ӧ�ʣ�
!    �����������ѧ��Ӧdtʱ�䲽���di (Kg/m3),  E
!    ���������ټ��㣬 �ӿ�Ϊ���ʵ�λ��  ���ڲ����㻯ѧ��Ӧ����ʱ�� ʹ��mol-cm�ƣ�    
!    di  kg/m3;  E: J/m3, dt: s;
!-----------------------------------------------------------------------------------------
!  �ӿڵ�λ�� ������ ;    di : kg/m3, E: J,  dt : s

! CHEM_Time_Euler=1, CHEM_Time_RK2=2, CHEM_Time_RK3=3, CHEM_Time_implicit1=-1
     subroutine chemical_timeAdv(di,E,dt)
      use CHEM
 	  implicit none
	  real(PRE_EC):: E,dt,di(N_Spec)
       select case (CHEM_TimeAdv)
	   case(CHEM_Time_Euler)
	    call chemical_Euler(di,E,dt)           ! 1���Ը�ʽ
	   case(CHEM_Time_RK2)
	    call chemical_RK2(di,E,dt)             ! 2��Runge-Kutta �Ը�ʽ
	   case(CHEM_Time_RK3)
	    call chemical_RK3(di,E,dt)             ! 2��Runge-Kutta �Ը�ʽ	   
	   case (CHEM_Time_implicit1)
	    call chemical_Euler_implicit(di,E,dt)  ! 1������ʽ
	   case default
	    print*, "This time advance method (in chemical reaction) is not supported !"
	   end select 
	 end 
	  


      
	  subroutine chemical_Euler(di,E,dt)
      use CHEM
 	  implicit none
	  real(PRE_EC):: E,T,dt
	  real(PRE_EC),dimension(N_Spec)::di, ci, dinew, ri     ! �ܶȣ�Ħ��Ũ��, tʱ�䲽�������ܶ�
  !        di ����ܶȣ� ��λkg/m3; dinew  ��Ӧ������ܶ�
  !        ci Ħ��Ũ�ȣ� ��λ  mol/cm3,   ri ��Ӧ���� ��λ mol/(cm3.s)         �� �ڲ�ʹ�� mol-cm ֵ��λ  
	  integer:: i,j,k
!------------------------------------------------------

        call comput_T(di,E,T)    ! �����ܼ����¶�
	
	    do i=1,N_SPEC
		   ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! ���Ħ��Ũ��  (��λ mol/cm3;    mol/m3 --> mol/cm3)
        enddo

        call  Reaction_rate(ci,ri,T)    ! ��֪Ħ��Ũ��ci,�¶�T,  ����i ��ֵ�������ri ����λ�� mol/ (cm3.s) )
	  
	     do i=1,N_Spec
	       dinew(i)=di(i)+ri(i) *SPEC(i)%Mi *1.d6 * dt             	! ʱ���ƽ� ��1��Euler) , ��λ kg/m3  ( kg/cm3--> kg/m3)   
        enddo
   
 	   call update_Et(di,dinew,E,T)   ! �������� (���ڻ�ѧ��Ӧ���ӵ�����)

	  do i=1,N_Spec
	   di(i)=dinew(i)
	  enddo

	  end 



!-----����2nd Runge-Kutta ���� 
! di in kg/m3 ;   ci in mol/cm3
  subroutine chemical_RK2(di,E,dt)
    use CHEM
    implicit none
    real(PRE_EC):: E,T,dt
    real(PRE_EC),dimension(N_Spec)::di, di1, dinew, ci, ri      ! �ܶȣ�Ħ��Ũ��, tʱ�䲽�������ܶ�
    integer:: i,j,k
    !------------------------------------------------------

    !--------step 1-------------------------- 
    call comput_T(di,E,T)    ! �����ܼ����¶�

    do i=1,N_SPEC
	  ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! ���Ħ��Ũ��  (��λ mol/cm3;    mol/m3 --> mol/cm3)
    end do
    call  Reaction_rate(ci,ri,T)    ! ��֪Ħ��Ũ��ci,�¶�T,  ����i ��ֵ�������ri ����λ�� mol/ (cm3.s) )

    do i=1,N_Spec
      di1(i)=di(i) + ri(i) *SPEC(i)%Mi *1.d6 * dt       ! ʱ���ƽ� (step 1) 
    end do
   
    call update_Et(di,di1,E,T )   ! ��������
   
    !--------step 2-------------------------- 
    call comput_T(di,E,T)       ! �����¶�

    do i=1,N_SPEC
      ci(i)=di1(i)/SPEC(i)%Mi  * 1.d-6          ! ���Ħ��Ũ�� (mol/m3 --> mol/cm3)
    end do

    call Reaction_rate(ci,ri,T)    ! ��֪Ħ��Ũ��ci,�¶�T, ���㷴Ӧ������������ ri
     
    do i=1,N_Spec
      dinew(i)=0.5d0*di(i)+0.5d0*di1(i)+0.5d0*  ri(i) *SPEC(i)%Mi *1.d6 * dt      ! ʱ���ƽ� (step 2) 
    end do
 
    call update_ET(di1,dinew,E,T )   ! ��������
    !--------------------------------------------------------------------
    do i=1,N_Spec             ! ��������ܶ�
      di(i)=dinew(i)
    end do
  end 

!------------------------------------------------------

	  subroutine chemical_Euler_implicit(di,E,dt)
      use CHEM
 	  implicit none
	  real(PRE_EC):: E,T,dt
	  real(PRE_EC),dimension(N_Spec)::di, ci, dinew, ri, di1,dd     ! �ܶȣ�Ħ��Ũ��, tʱ�䲽�������ܶ�
  !        di ����ܶȣ� ��λkg/m3; dinew  ��Ӧ������ܶ�
  !        ci Ħ��Ũ�ȣ� ��λ  mol/cm3,   ri ��Ӧ���� ��λ mol/(cm3.s)         �� �ڲ�ʹ�� mol-cm ֵ��λ  
	  integer:: i,j,k
	   real(PRE_EC),dimension(N_Spec,N_Spec):: Arate ,A1  ! d(ri)/d(cj)
!------------------------------------------------------

        call comput_T(di,E,T)    ! �����ܼ����¶�
	
	    do i=1,N_SPEC
		   ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! ���Ħ��Ũ��  (��λ mol/cm3;    mol/m3 --> mol/cm3)
        enddo
        call Reaction_rate_and_Matrix(ci,ri,Arate,T)
		do j=1,N_SPEC
		do i=1,N_SPEC
          if(i==j) then		
		    A1(i,j)=1.d0
		  else
		    A1(i,j)=0.d0
		  endif 
		  A1(i,j)=A1(i,j)-Arate(i,j)*SPEC(i)%Mi/SPEC(j)%Mi*dt    ! Matrix
		enddo
		enddo 
		do i=1,N_SPEC 
		di1(i)=ri(i)*SPEC(i)%Mi*1.d6*dt    ! ��ʽ���֣�  ��λ kg/m3  (1.d6: ��λת�� kg/cm3--> kg/m3 ) 
		enddo 
		
		call Gauss(N_SPEC,A1,di1,dd)  ! ������Դ���������(A1*dd=di1)  (dd=dinew-di)
	  
	     do i=1,N_Spec
	       dinew(i)=di(i)+dd(i)           	! ʱ���ƽ� ��1��Euler) , ��λ kg/m3  ( kg/cm3--> kg/m3)   
        enddo
   
 	   call update_Et(di,dinew,E,T)   ! �������� (���ڻ�ѧ��Ӧ���ӵ�����)

	  do i=1,N_Spec
	   di(i)=dinew(i)
	  enddo

	  end 
! ==========================================================================
!  Ax=b  Solver                                                            =
!  Gauss  Elimination with Maximal Column Pivoting                         =            
!  Gauss  ����Ԫ��ȥ���� ���Ax=b                                          =
!  ������Ԫ���̣������ A��b ��ֵ�ᷢ���ı�                                =
! ========================================================================== 
	 subroutine Gauss(n,A,b,x)
	 use Precision_EC
	 implicit none
	 integer:: n,i,j,k,kk
	 real(PRE_EC):: A(n,n),b(n),x(n)
	 real(PRE_EC):: fmax,t,f
	do k=1,n-1
! ---- find main element ---        
	  fmax=abs(a(k,k))
	  kk=k
	  do i=k+1,n
	  if( abs(a(i,k)).gt.fmax )  then
	   fmax=abs(a(i,k))
	   kk=i
	  endif
	  enddo
!c ----- exchange line------         
	  do j=k,n
	  t=a(k,j)
	  a(k,j)=a(kk,j)
	  a(kk,j)=t
	  enddo
	  t=b(k)
	  b(k)=b(kk)
	  b(kk)=t
!c    -----Elimination --
	 do i=k+1,n
	  f=-a(i,k)/a(k,k)
	  do j=1,n
	  a(i,j)=a(i,j)+f*a(k,j)
	  enddo
	  b(i)=b(i)+f*b(k)
	 enddo
	enddo
!c  ----- get x ---------       
	x(n)=b(n)/a(n,n)
	do i=n-1,1,-1
	x(i)=b(i)
	do j=i+1,n
	x(i)=x(i)-a(i,j)*x(j)
	enddo
	x(i)=x(i)/a(i,i)
    enddo

 end
 
 
!-----����3rd Runge-Kutta ���� 
! di in kg/m3 ;   ci in mol/cm3
  subroutine chemical_RK3(di,E,dt)
    use CHEM
    implicit none
    real(PRE_EC):: E,T,dt
    real(PRE_EC),dimension(N_Spec)::di,din, dinew, ci, ri      ! �ܶȣ�Ħ��Ũ��, tʱ�䲽�������ܶ�
    integer:: i,j,k,KRK
    !------------------------------------------------------
    real(kind=OCFD_REAL_KIND):: Ralfa(3),Rbeta(3) 
   
        Ralfa(1)=0.d0 ;       Rbeta(1)=1.d0
        Ralfa(2)=3.d0/4.d0 ;  Rbeta(2)=1.d0/4.d0
        Ralfa(3)=1.d0/3.d0 ;  Rbeta(3)=2.d0/3.d0
		
	din(:)=di(:)	
   do KRK=1,3		
    !--------step 1-------------------------- 
    call comput_T(di,E,T)    ! �����ܼ����¶�

    do i=1,N_SPEC
	  ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! ���Ħ��Ũ��  (��λ mol/cm3;    mol/m3 --> mol/cm3)
    enddo
    call  Reaction_rate(ci,ri,T)    ! ��֪Ħ��Ũ��ci,�¶�T,  ����i ��ֵ�������ri ����λ�� mol/ (cm3.s) )

     do i=1,N_Spec
      dinew(i)= Ralfa(KRK)*din(i)+ Rbeta(KRK)*(di(i) +    ri(i) *SPEC(i)%Mi *1.d6 * dt)       ! ʱ���ƽ� (step K) 
     enddo
	 
    call update_Et(di,dinew,E,T )   ! ��������
    di=dinew 
   enddo  
 end 

!------------------------------------------------------
 
 