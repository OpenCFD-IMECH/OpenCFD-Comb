! OpenCFD-Comb ， Ver 1.3
! 与化学反应、热力学特性有关的物理量
! 本计算采用有量纲制； 如无特殊声明，物理量采用国际标准单位； 
! 化学反应速率Arrhenius=A*T**b*exp(-E/T) A单位为mol/cm3
! Ver1.3, 2016-11-20,  本版本不支持点隐格式 （如需要请使用早期版本)
! Ver 1.4, 2022-2-16: Cp 多项式拟合公式中,当T> Tmax时， Cp=Cp_max (以免超过温度范围出现Cp异常）
! Ver 1.5, 2024-12-18: TYPE REACTION 中添加了5 个逆反应系数. 格式见Passiatore 2022 (JFM 941, A21) or Park 1990
!---------------------------------------------------------------------------------------------------  
  module Precision_EC
  include "mpif.h"
  integer,parameter:: PRE_EC=8 , OCFD_REAL_KIND=8          ! Double Precision
  integer,parameter:: OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION
  end module  Precision_EC 

! 常数
!-------------------------------------------------------------------- 
 module const_chem
  use Precision_EC
  implicit none
  real(PRE_EC),parameter::  R0=8.314d0             ! 通用气体常数 (J/mol-K)
  real(PRE_EC),parameter::  J2K=4.184d0            ! 热功当量  (Kal/J)              
  real(PRE_EC),parameter::  atm=1.013d5            ! 大气压 (Pa)
  integer,parameter:: Real_GAS=1, Perfect_GAS=0
  integer,parameter:: Af_nomal=0, Af_Highpress=1, Af_DualArrhenius=2     ! 常规，高压反应，双Arrhenius描述反应
  integer,parameter:: NTct_Max=10                  ! 转换温度的最大数目
  integer,parameter:: CHEM_Time_Euler=1, CHEM_Time_RK2=2, CHEM_Time_RK3=3, CHEM_Time_implicit1=-1 
  end module const_chem

!------------------------------------------------------------------------------------
! 热力学、化学反应特性
 module CHEM
  use Const_chem
  implicit none
  TYPE SPECIE 
      character*10:: name        ! 组分名称
	  real(PRE_EC):: Mi,Ri                 ! 组分的摩尔质量（分子量， kg/mol）, 气体常数(=R0/Mi)  
!	  real(PRE_EC):: A1,B1,C1,D1,E1,F1,G1   ! Cp的温度拟合参数 （低温区, <Tct 适用）; F,G计算化学焓、熵时使用
!	  real(PRE_EC):: A2,B2,C2,D2,E2,F2,G2   ! Cp的温度拟合参数 （高温区 > Tct 适用）
	  real(PRE_EC),dimension(NTct_Max):: Ai,Bi,Ci,Di,Ei,Fi,Gi   ! Cp的温度拟合参数 （分段拟合，最多允许NTct_Max段）; Fi, Gi计算化学焓、熵时使用
	  real(PRE_EC):: Ect(NTct_Max)          ! Ect=F2-F1 (跨转换温度Tct的内能差, 以保障热焓跨段连续)
      integer:: SpecFlag   ! 1 单原子, 2 双原子
 END TYPE SPECIE
  
!------------------------化学反应速率描述-----------------------------------
  TYPE REACTION
	 real(PRE_EC):: Af(2),Bf(2),Ef(2)       ! 化学反应速率Arrhenius公式中的系数 （正反应）
	                                        !  通常只使用1套系数， 但对于高压反应 或双Arrhenius公式的反应使用2套系数
	 real(PRE_EC):: Ar,Br,Er                ! 化学反应速率Arrhenius公式中的系数 （逆反应， 如采用反应平衡常数计算，则可忽略）
     real(PRE_EC):: A1r, A2r, A3r, A4r, A5r  ! Park 1990 模型计算平衡常数Keq,r 使用的5个拟合系数 (Passiatore 2022 JFM; Park 1990)
	 integer:: sgm, TrReac, Af_type        ! sgm :  逆反应与正反应之间的级差  (参考手册1.20式)   ! TrReac 是否三体反应 (0/1)
                                                              ! Af_type =0, 常规； 1 高-低压反应Falloff型（三体反应）；  2 采用两套Arrhenius系数
     real(PRE_EC):: Fc_troe                  ! Troe 系数 （高压 fall-off型反应使用） 
     real(PRE_EC),dimension(:),pointer:: TrEff ! 三体反应系数
  ENDTYPE
!-----------------------------------------------------------------------------------------
  integer:: CHEM_TimeAdv    ! 化学反应计算时间推进方法
  
  integer:: N_Spec,N_Reac   ! 组分数，化学反应数
  integer:: ReReac          ! 逆反应速率算法(0 采用反应平衡常数计算; 1 Arrhenius公式;  2 Park 的5参数拟合平衡常数公式)
  integer:: IF_RealGas = Real_GAS ! 内能-温度计算方法 （完全气体 Cp =7/2 R; 真实气体：Cp 拟合公式; 默认Real_GAS）

  TYPE (SPECIE), dimension(:),pointer:: SPEC            ! 组分
  TYPE (REACTION),dimension(:),pointer:: REAC           ! 反应
  integer:: N_Tct      ! 热力学量分段拟合公式的 分段数  （常见为2，分低温及高温两段）
  real(PRE_EC):: Tct(NTct_Max)                                   ! 分段拟合公式的转换温度  (计算CP, 例如1000K)
  real(PRE_EC):: Tct_max              ! 最高温度 （超过该温度则Cp为定值）
  integer,dimension(:,:),pointer:: Vf, Vr    ! 正/逆反应系数矩阵 （Vf(N_spec,N_Reac) 

 END module CHEM
!------------------------------------------------------------------------------------------
