! OpenCFD-Comb �� Ver 1.3
! �뻯ѧ��Ӧ������ѧ�����йص�������
! ����������������ƣ� �����������������������ù��ʱ�׼��λ�� 
! ��ѧ��Ӧ����Arrhenius=A*T**b*exp(-E/T) A��λΪmol/cm3
! Ver1.3, 2016-11-20,  ���汾��֧�ֵ�����ʽ ������Ҫ��ʹ�����ڰ汾)
! Ver 1.4, 2022-2-16: Cp ����ʽ��Ϲ�ʽ��,��T> Tmaxʱ�� Cp=Cp_max (���ⳬ���¶ȷ�Χ����Cp�쳣��
! Ver 1.5, 2024-12-18: TYPE REACTION �������5 ���淴Ӧϵ��. ��ʽ��Passiatore 2022 (JFM 941, A21) or Park 1990
!---------------------------------------------------------------------------------------------------  
  module Precision_EC
  include "mpif.h"
  integer,parameter:: PRE_EC=8 , OCFD_REAL_KIND=8          ! Double Precision
  integer,parameter:: OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION
  end module  Precision_EC 

! ����
!-------------------------------------------------------------------- 
 module const_chem
  use Precision_EC
  implicit none
  real(PRE_EC),parameter::  R0=8.314d0             ! ͨ�����峣�� (J/mol-K)
  real(PRE_EC),parameter::  J2K=4.184d0            ! �ȹ�����  (Kal/J)              
  real(PRE_EC),parameter::  atm=1.013d5            ! ����ѹ (Pa)
  integer,parameter:: Real_GAS=1, Perfect_GAS=0
  integer,parameter:: Af_nomal=0, Af_Highpress=1, Af_DualArrhenius=2     ! ���棬��ѹ��Ӧ��˫Arrhenius������Ӧ
  integer,parameter:: NTct_Max=10                  ! ת���¶ȵ������Ŀ
  integer,parameter:: CHEM_Time_Euler=1, CHEM_Time_RK2=2, CHEM_Time_RK3=3, CHEM_Time_implicit1=-1 
  end module const_chem

!------------------------------------------------------------------------------------
! ����ѧ����ѧ��Ӧ����
 module CHEM
  use Const_chem
  implicit none
  TYPE SPECIE 
      character*10:: name        ! �������
	  real(PRE_EC):: Mi,Ri                 ! ��ֵ�Ħ���������������� kg/mol��, ���峣��(=R0/Mi)  
!	  real(PRE_EC):: A1,B1,C1,D1,E1,F1,G1   ! Cp���¶���ϲ��� ��������, <Tct ���ã�; F,G���㻯ѧ�ʡ���ʱʹ��
!	  real(PRE_EC):: A2,B2,C2,D2,E2,F2,G2   ! Cp���¶���ϲ��� �������� > Tct ���ã�
	  real(PRE_EC),dimension(NTct_Max):: Ai,Bi,Ci,Di,Ei,Fi,Gi   ! Cp���¶���ϲ��� ���ֶ���ϣ��������NTct_Max�Σ�; Fi, Gi���㻯ѧ�ʡ���ʱʹ��
	  real(PRE_EC):: Ect(NTct_Max)          ! Ect=F2-F1 (��ת���¶�Tct�����ܲ�, �Ա������ʿ������)
      integer:: SpecFlag   ! 1 ��ԭ��, 2 ˫ԭ��
 END TYPE SPECIE
  
!------------------------��ѧ��Ӧ��������-----------------------------------
  TYPE REACTION
	 real(PRE_EC):: Af(2),Bf(2),Ef(2)       ! ��ѧ��Ӧ����Arrhenius��ʽ�е�ϵ�� ������Ӧ��
	                                        !  ͨ��ֻʹ��1��ϵ���� �����ڸ�ѹ��Ӧ ��˫Arrhenius��ʽ�ķ�Ӧʹ��2��ϵ��
	 real(PRE_EC):: Ar,Br,Er                ! ��ѧ��Ӧ����Arrhenius��ʽ�е�ϵ�� ���淴Ӧ�� ����÷�Ӧƽ�ⳣ�����㣬��ɺ��ԣ�
     real(PRE_EC):: A1r, A2r, A3r, A4r, A5r  ! Park 1990 ģ�ͼ���ƽ�ⳣ��Keq,r ʹ�õ�5�����ϵ�� (Passiatore 2022 JFM; Park 1990)
	 integer:: sgm, TrReac, Af_type        ! sgm :  �淴Ӧ������Ӧ֮��ļ���  (�ο��ֲ�1.20ʽ)   ! TrReac �Ƿ����巴Ӧ (0/1)
                                                              ! Af_type =0, ���棻 1 ��-��ѹ��ӦFalloff�ͣ����巴Ӧ����  2 ��������Arrheniusϵ��
     real(PRE_EC):: Fc_troe                  ! Troe ϵ�� ����ѹ fall-off�ͷ�Ӧʹ�ã� 
     real(PRE_EC),dimension(:),pointer:: TrEff ! ���巴Ӧϵ��
  ENDTYPE
!-----------------------------------------------------------------------------------------
  integer:: CHEM_TimeAdv    ! ��ѧ��Ӧ����ʱ���ƽ�����
  
  integer:: N_Spec,N_Reac   ! ���������ѧ��Ӧ��
  integer:: ReReac          ! �淴Ӧ�����㷨(0 ���÷�Ӧƽ�ⳣ������; 1 Arrhenius��ʽ;  2 Park ��5�������ƽ�ⳣ����ʽ)
  integer:: IF_RealGas = Real_GAS ! ����-�¶ȼ��㷽�� ����ȫ���� Cp =7/2 R; ��ʵ���壺Cp ��Ϲ�ʽ; Ĭ��Real_GAS��

  TYPE (SPECIE), dimension(:),pointer:: SPEC            ! ���
  TYPE (REACTION),dimension(:),pointer:: REAC           ! ��Ӧ
  integer:: N_Tct      ! ����ѧ���ֶ���Ϲ�ʽ�� �ֶ���  ������Ϊ2���ֵ��¼��������Σ�
  real(PRE_EC):: Tct(NTct_Max)                                   ! �ֶ���Ϲ�ʽ��ת���¶�  (����CP, ����1000K)
  real(PRE_EC):: Tct_max              ! ����¶� ���������¶���CpΪ��ֵ��
  integer,dimension(:,:),pointer:: Vf, Vr    ! ��/�淴Ӧϵ������ ��Vf(N_spec,N_Reac) 

 END module CHEM
!------------------------------------------------------------------------------------------
