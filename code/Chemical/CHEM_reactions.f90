! 与化学反应有关的子程序
! Ver 1.3, 2016-11-20  重新改写;   考虑了高压反应， 双Arrhenius公式反应；  A使用mol/cm3作为单位
! Ver 1.3a, 2017-2-28, 修改了  Reaction_rate( ) 中的Bug  (关于三体浓度 Ctrb ）
! Ver 1.5, 2022-2-24, 支持隐格式
! Ver 1.5a, 2023-8-15, 高压反应部分修改了Bug;  log(Fc) ==> log10(Fc)
! Ver 1.6, 2024-12-18,   支持 Park (1990) 逆反应速率模型，使用5参数拟合平衡常数（Ref: Passiatore （2022） JFM 941, A21, doi:10.1017/jfm.2022.283） 
!          ReReac=0 ：基于Gibbs自由能的平衡常数； ReReac=1: 使用Arrhenius公式； ReReac=2: 使用Park (1990) 五参数拟合的平衡常数模型
! --------------------------------------------------------
!  仅与化学反应速率有关的单位采用 mol-cm制 A*b**T*exp(-E/T)中的A
!  本程序的接口仍采用国际单位制 :  kg-m-s-K-mol
! --------------- using mol-cm  units-------------------   


! 读取化学反应数据 （ A 使用 mol-cm 单位;   b, E 采用K为单位；  K=A*T**b*exp(-E/T)）
  subroutine read_Reac
    use CHEM
    implicit none
    integer:: i,j
    TYPE (REACTION), pointer:: RC

!    print*, "--------------Read Reaction.in --------------------------------------"
    
    open(100,file="Reaction.in")                          ! 化学反应特征
    read(100,*)
    read(100,*) 
! N_spec 组分数目；  N_REAC 化学反应数目;  ReReac 逆反应速率算法（0 采用基于Gibbs自由能的平衡常数, 1 采用Arrhenius公式，2 Park 1990 平衡常数模型 ）
! CHEM_TimeAdv 化学反应时间推进方法：1： 1阶Euler; 2： 2阶RK； 3： 3阶RK (目前版本尚不支持)； -1： 1阶隐格式

    read(100,*) N_Spec, N_REAC , ReReac     
    allocate(REAC(N_REAC))
    allocate(Vf(N_SPEC,N_REAC), Vr(N_SPEC,N_REAC) )          ! 反应系数矩阵 （维数：组分数*反应数）
    read(100,*)
    read(100,*)                                       
	read(100,*)
    read(100,*)
    !-------------------------------------------------------------------------------                                           
  
    do j=1,N_REAC                                         ! 读入反应信息
    !  print*, "-------Reaction --", j
      Rc=>REAC(j)
      read(100,*)
      read(100,*) 
	  read(100,*)  Rc%Af_type,  Rc%TrReac , Rc%Fc_troe      ! 正反应系数类型；  是否三体反应 ; Troe系数 （高压反应使用）
      
	  read(100,*) Rc%Af(1),Rc%Bf(1),Rc%Ef(1)                 ! 正反应系数
      if(Rc%Af_type .ne.  Af_nomal )     read(100,*)  Rc%Af(2), Rc%Bf(2), Rc%Ef(2)       ! 第2套正反应系数 （用于高压反应 或 双Arrhenius系数描述）
	   if( ReReac==0) read(100,*)   ! 用平衡常数计算逆反应，不使用参数； 读取空行，保持Reaction.in的通用性
	   if( ReReac==1)  read(100,*) Rc%Ar, Rc%Br, Rc%Er                     ! Arrhenius 格式的逆反应速率系数 
       if( ReReac==2)  read(100,*) Rc%A1r, Rc%A2r, Rc%A3r, Rc%A4r, Rc%A5r  ! Park 平衡常数模型中的5个参数， 用于算逆反应速率
	  
      read(100,*) (vf(i,j), i=1,N_spec)    ! 正反应系数矩阵
      read(100,*) (vr(i,j), i=1,N_spec)    ! 逆反应系数矩阵
      
	  if(Rc%TrReac .eq. 1) then
        allocate(Rc%TrEff(N_SPEC))                    ! 三体影响系数
        read(100,*) (Rc%TrEff(i), i=1,N_Spec)
      end if
      
	  Rc%sgm=0
      do i=1,N_spec
        Rc%sgm=Rc%sgm+vr(i,j)-vf(i,j)      ! 逆-正 反应的级差
      end do
    end do
    close(100)

    ! print*, "read Reaction.in OK"
    ! print*, "N_spec=", N_spec, "N_Reac=", N_Reac    

  end subroutine read_Reac





! 计算化学平衡常数， T为温度
  subroutine comput_KjX(T,Kjx)
    use CHEM
    implicit none
    TYPE (SPECIE),pointer:: Sp
    integer:: i,j
    real(PRE_EC) :: T,Kjx(N_Reac),dgi,gi(N_spec)
	real(PRE_EC),parameter:: Rs=R0/atm * 1.d6                          ! Rs,  单位 cm3/(mol.K)  见参考手册 1.22 式

    ! 计算各组分的Gibbs 自由焓

    do i=1,N_spec
      call comput_gibbs(i,T,gi(i))
!      SP=>Spec(i)    ! 第 i个组分
!      if(T<Tct) then
!         gi(i)=Sp%A1*(1.d0-log(T))*T-SP%B1*T**2/2.d0-SP%C1*T**3/6.d0-SP%D1*T**4/12.d0-SP%E1*T**5/20.d0+SP%F1-SP%G1*T
!       else
!        gi(i)=Sp%A2*(1.d0-log(T))*T-SP%B2*T**2/2.d0-SP%C2*T**3/6.d0-SP%D2*T**4/12.d0-SP%E2*T**5/20.d0+SP%F2-SP%G2*T
!      end if
    end do

    do j=1,N_Reac
      dgi=0.d0
      do i=1,N_Spec
        dgi=dgi+gi(i)*(Vr(i,j)-Vf(i,j))       ! 自由焓差
       end do
      Kjx(j)=exp(-dgi/T) * (Rs*T)**(-Reac(j)%sgm)                    ! 见理论手册 1.22 式 ; sgm 逆-正反应的级差
    end do
  end

!-----------------------------------------------------------------------
function Arrhenius(A,b,E,T)          ! Arrhenius 公式
  	use Precision_EC	
	implicit none
  	real(PRE_EC)::  A,b,E,T,Arrhenius
  	Arrhenius=A*T**b*exp(-E/T)          ! E使用 K 作为单位；  
end

function ParkCurveFit(A1,A2,A3,A4,A5,T)       ! Park 5参数拟合公式， 计算平衡常数Kjx
    use Precision_EC
	implicit none
    real(PRE_EC)::  A1,A2,A3,A4,A5,T,ParkCurveFit,z
	z = 10000.d0/T
    ParkCurveFit = exp(A1 + A2*z + A3*(z**2) + A4*(z**3) + A5*(z**4))            
end
!---------------------------------------------------------------------




!   计算化学反应速率系数  (单位：mol-cm-K-s 制)
!   注意， 与其他子程序单位不同，尤其是 cm制
!   速率系数中包含了与三体反应有关的量
!   高压反应中的log(Fc) 修改为 log10(Fc)  (文献中的 log() 指的是log10(), 而不是 ln())
!   添加了 Park 逆反应速率模型 （2024-12-19）
	subroutine Reaction_rate_coef(ci,wf,wr,T)  ! wf, wr 正反应及逆反应的速率系数
    use CHEM
	
 	implicit none
	real(PRE_EC),dimension(N_Reac):: wf,wr, Kjx  ! 反应速率系数;  化学反应平衡常数 (mol-cm-K unit)
	real(PRE_EC),dimension(N_Spec):: ci     ! 摩尔浓度 (mol/cm3) 
    real(PRE_EC):: T,Ctrb,Arrhenius ,ParkCurveFit, wf1,wf2, Prs, F,Fc, ac,an       !化学反应速率有关量
    TYPE (REACTION), pointer:: RC

	integer:: i,j,k
	! 基于Gibbs自由能，计算化学反应平衡常数Kjx （ReReac=1 使用Arrhenius公式计算逆反应速率,无需使用kjx;  ReReac=2 采用Park 公式计算kjx）
        if(ReReac ==0) then        
        	call comput_KjX(T,Kjx)
        endif

! ----------  计算化学反应速率常数 wf, wr ----------------
		do j=1,N_REAC
	 		Rc=>REAC(j)
	 
			if(Rc%TrReac .eq. 1) then            ! Three-body reaction
        		Ctrb=0.d0
				do i=1,N_SPEC
        	   		Ctrb=Ctrb+ci(i)*Rc%TrEff(i)             ! 三体浓度 （各组分加权和；Rc%TrEff(:)为影响因子 ）
			 	enddo
        	endif

!     正反应速率常数
		  	wf1=Arrhenius(Rc%Af(1), Rc%bf(1), Rc%Ef(1), T)            ! 正反应速度系数
			if(  Rc%Af_type .eq.  Af_nomal)  then
				wf(j)=wf1                                                                        ! 常规反应
        	else if (Rc%Af_type .eq. Af_DualArrhenius) then                ! 双Arrhenius公式描述
				wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
				wf(j)=wf1+wf2
        	else if (Rc%Af_type .eq. Af_Highpress) then                        !  高压反应；  Fall-off 型三体反应， 如 H+O2 (+M) = HO2 (+M)
				wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
				Prs=wf1/wf2*Ctrb                          ! Ctrb=[M]
        		Fc=Rc%Fc_troe                              ! 高压反应的Troe形式， 见Chemkin理论手册
!				ac=-0.4d0-0.67*log(Fc)
!        		an=0.75d0-1.27*log(Fc)
!        		F=exp(log(Fc)/( 1.d0+ ( (log(Prs)+ac)/(an-0.14d0*(log(Prs)+ac) ) )**2) )
				ac=-0.4d0-0.67*log10(Fc)
        		an=0.75d0-1.27*log10(Fc)
        		F=10.d0**(log10(Fc)/( 1.d0+ ( (log10(Prs)+ac)/(an-0.14d0*(log10(Prs)+ac) ) )**2) )

!               wf(j)=wf2*Prs/(1.d0+Prs)                 ! Lindemann type
	    		wf(j)=wf2*Prs/(1.d0+Prs)*F               ! Troe type
        	endif

!     计算逆反应速率 （3种模型）
!          2024-12-19, 引入Park 1990 模型
            if(ReReac ==0 ) then 
			   wr(j)=wf(j)/kjx(j)           ! 基于Gibbs函数的平衡常数模型 
			else if (ReReac ==1) then 
			   wr(j)=Arrhenius(Rc%Ar,Rc%br,Rc%Er,T)            ! 基于Arrhenius公式的逆反应速率模型
            else if (ReReac ==2) then
			  kjx(j)=ParkCurveFit(Rc%A1r,Rc%A2r,Rc%A3r,Rc%A4r,Rc%A5r,T)   ! Park  模型, 计算平衡常数 
	          if(kjx(j) < 1.d-30) kjx(j)=1.d-30                ! 设定下限防止分母为0，1.d-30 有人为性。 2024-12-19  		  
			  wr(j)=wf(j)/kjx(j)
!              wr(j)=ParkCurveFit(Rc%A1r,Rc%A2r,Rc%A3r,Rc%A4r,Rc%A5r,T)     ! wrong test 
			endif 

		 
	    	if(Rc%TrReac .eq. 1  .and. Rc%Af_type .ne. Af_Highpress ) then       !  三体反应 （非Fall-off型； Fall-off型的速率常数中包含了三体浓度Ctrb)
        		wf(j)=wf(j)*Ctrb
				wr(j)=wr(j)*Ctrb
        	endif
	  	enddo 
	end subroutine Reaction_rate_coef
	  


!     计算反应物生成率  (mol/(cm3.s))
!     单位 mol-cm-K   （注意， 与其他子程序单位不同，尤其是 cm制)	    
	subroutine Reaction_rate(ci,ri,T)  
    use CHEM
 	implicit none
	real(PRE_EC),dimension(N_Reac):: wf,wr  ! 反应速率系数 (mol-cm-K unit)
	real(PRE_EC),dimension(N_Spec):: ci,ri     ! 摩尔浓度 (mol/cm3) , 反应物生成率 (mol/(cm3.s))
    real(PRE_EC):: T,wf1,wr1
	integer:: i,j
!   计算反应速率系数 wf, wr   ( 单位 mol-s-K-cm 制 )
    	call Reaction_rate_coef(ci,wf,wr,T)
!   计算化学反应速率 wf1, wr1  	  单位： mol/(cm3-s)
	 	ri(:)=0.d0 
	 	do j=1,N_REAC
	 	  	wf1=wf(j)        ! 正反应速率 （第j个反应）
	 	  	wr1=wr(j)        ! 逆反应速率（第j个反应）
			do i=1,N_SPEC
				if(vf(i,j) .ne. 0) wf1=wf1*ci(i)**vf(i,j)
				if(vr(i,j) .ne. 0) wr1=wr1*ci(i)**vr(i,j)
			enddo
	 	    do i=1,N_SPEC         ! 组分
     	    	ri(i)=ri(i)+ (wf1-wr1) *(vr(i,j)-vf(i,j))   ! 第i个组分的生成率 (摩尔/(立方厘米.秒))        
			enddo
     	enddo
    end subroutine Reaction_rate

!---化学反应速率,及其导数矩阵 Arate(i,j)=  d(ri)/d(cj)          ( 摩尔浓度生成率对摩尔浓度的导数； 没有乘以Mi/Mj)
  !     计算反应物生成率  (mol/(cm3.s))
!     单位 mol-cm-K   （注意， 与其他子程序单位不同，尤其是 cm制)	    
	   subroutine Reaction_rate_and_Matrix(ci,ri,Arate,T)  
       use CHEM
 	   implicit none
	   real(PRE_EC),dimension(N_Reac):: wf,wr  ! 反应速率系数 (mol-cm-K unit)
	   real(PRE_EC),dimension(N_Spec):: ci,ri       ! ci 摩尔浓度 (mol/cm3); ri 摩尔浓度生成率   
	   real(PRE_EC),dimension(N_Spec,N_Spec):: Arate        ! 导数矩阵
       real(PRE_EC):: T,wf1,wr1,Awf,Awr  !Awf=d(eta)/(cj)  正反应速率对cj的导数， Awr 逆反应速率对cj导数 
	   integer:: i,j,k
	   real(PRE_EC),parameter::  Aepsl=1.d-20
	   
!   计算反应速率系数 wf, wr   ( 单位 mol-s-K-cm 制 )
      call Reaction_rate_coef(ci,wf,wr,T)
	  
!   计算化学反应速率 wf1, wr1  	  单位： mol/(cm3-s)
	   ri(:)=0.d0 	
       Arate(:,:)=0.d0 
	   
	 do k=1,N_REAC
 	  
	   wf1=wf(k)
	   wr1=wr(k)
	   do i=1,N_SPEC
		  if(vf(i,k) .ne. 0) wf1=wf1*ci(i)**vf(i,k)
		  if(vr(i,k) .ne. 0) wr1=wr1*ci(i)**vr(i,k)
	   enddo
	   
  !  反应物生成率 
	   do i=1,N_SPEC         ! 组分
        ri(i)=ri(i)+ (wf1-wr1) *(vr(i,k)-vf(i,k))   ! 第i个组分的生成率 (摩尔/(立方厘米.秒))        
       enddo
		
 ! 导数矩阵		
	   do j=1,N_SPEC                 ! d(eta_i)/d(c_j) 	 
        Awf=vf(j,k)*wf1/(ci(j)+Aepsl)          ! +Aepsl 防止分母为0
        Awr=vr(j,k)*wr1/(ci(j)+Aepsl)          ! +Aepsl 防止分母为0
       do i=1,N_SPEC 
	    Arate(i,j)=Arate(i,j)+(Awf-Awr) *(vr(i,k)-vf(i,k)) 
       enddo 
	  enddo
	enddo
!---------------------------------------------------------
  
   end  


!----------测试化学反应速率---------------------------------------------------------------------------------------------
	   subroutine test_Reaction_rate
       use CHEM
 	   implicit none
	   real(PRE_EC),dimension(N_Reac):: wf,wr, Kjx  ! 反应速率系数;  化学反应平衡常数 (mol-cm-K unit)
	   real(PRE_EC),dimension(N_Spec):: gi       
       real(PRE_EC):: p1, T,Ctrb,Arrhenius , wf1,wf2, Prs , Fc,F, ac,an      !化学反应速率有关量
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
	call comput_gibbs(i,T,gi(i))          ! 计算Gibbs函数 g=h-TS
!      SP=>Spec(i)    ! 第 i个组分
!      if(T<Tct) then
!        gi(i)=(Sp%A1*(1.d0-log(T))*T-SP%B1*T**2/2.d0-SP%C1*T**3/6.d0-SP%D1*T**4/12.d0-SP%E1*T**5/20.d0+SP%F1-SP%G1*T)
!      else
!        gi(i)=(Sp%A2*(1.d0-log(T))*T-SP%B2*T**2/2.d0-SP%C2*T**3/6.d0-SP%D2*T**4/12.d0-SP%E2*T**5/20.d0+SP%F2-SP%G2*T)
!      end if
    end do

    write(99,"(30E16.8)") T, (gi(j), j=1,N_Spec)     ! gi 单位为 R0
   enddo


	  open(99,file="Kjx.dat")
	  open(100,file="Reaction-wf.dat")
	  open(101,file="Reaction-wr.dat")

     do k=300, 2000
		   T=k*1.d0
		   ctrb= p1*atm/(R0*T) * 1.d-6     ! [M] in mol/cm3

           call comput_KjX(T,Kjx)
		   write(99,"(30E16.8)") T, (KjX(j), j=1,N_Spec)     ! gi 单位为 R0
      

! ----------  计算化学反应速率常数 wf, wr ----------------
	do j=1,N_REAC
	 Rc=>REAC(j)

!     正反应速率常数
		  wf1=Arrhenius(Rc%Af(1), Rc%bf(1), Rc%Ef(1), T)            ! 正反应速度系数
		  if(  Rc%Af_type .eq.  Af_nomal)  then
		     wf(j)=wf1                                                                        ! 常规反应
          else if (Rc%Af_type .eq. Af_DualArrhenius) then                ! 双Arrhenius公式描述
		    wf2=Arrhenius(Rc%Af(2), Rc%bf(2), Rc%Ef(2), T) 
			wf(j)=wf1+wf2
         else if (Rc%Af_type .eq. Af_Highpress) then                        !  高压反应；  Fall-off 型三体反应， 如 H+O2 (+M) = HO2 (+M)
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
!     逆反应速率常数
		 if(ReReac ==1) then
		   wr(j)=Arrhenius(Rc%Ar,Rc%br,Rc%Er,T)            ! 负反应速度系数
         else
           wr(j)=wf(j)/kjx(j)      
         endif
 enddo

     write(100,"(30E16.8)") T, (wf(j), j=1,N_Spec)     ! gi 单位为 R0
	 write(101,"(30E16.8)") T, (wr(j), j=1,N_Spec)     ! gi 单位为 R0

enddo
 
 end
 
!--------------------------------------------------------------------
!   Time Advance for Chemical reaction terms  ( 1st Euler; 2nd RK;  1st Euler-implicit) 
!   化学反应项的时间推进


!-------------------------------------------------------------------
! 计算化学反应项，推进 dt 时间步 
!    输入变量：di 各组分密度 (Kg/m3)； E 内能（不含反应焓）
!    输出变量：化学反应dt时间步后的di (Kg/m3),  E
!    采用有量纲计算， 接口为国际单位制  （内部计算化学反应速率时， 使用mol-cm制）    
!    di  kg/m3;  E: J/m3, dt: s;
!-----------------------------------------------------------------------------------------
!  接口单位： 国际制 ;    di : kg/m3, E: J,  dt : s

! CHEM_Time_Euler=1, CHEM_Time_RK2=2, CHEM_Time_RK3=3, CHEM_Time_implicit1=-1
     subroutine chemical_timeAdv(di,E,dt)
      use CHEM
 	  implicit none
	  real(PRE_EC):: E,dt,di(N_Spec)
       select case (CHEM_TimeAdv)
	   case(CHEM_Time_Euler)
	    call chemical_Euler(di,E,dt)           ! 1阶显格式
	   case(CHEM_Time_RK2)
	    call chemical_RK2(di,E,dt)             ! 2阶Runge-Kutta 显格式
	   case(CHEM_Time_RK3)
	    call chemical_RK3(di,E,dt)             ! 2阶Runge-Kutta 显格式	   
	   case (CHEM_Time_implicit1)
	    call chemical_Euler_implicit(di,E,dt)  ! 1阶隐格式
	   case default
	    print*, "This time advance method (in chemical reaction) is not supported !"
	   end select 
	 end 
	  


      
	  subroutine chemical_Euler(di,E,dt)
      use CHEM
 	  implicit none
	  real(PRE_EC):: E,T,dt
	  real(PRE_EC),dimension(N_Spec)::di, ci, dinew, ri     ! 密度，摩尔浓度, t时间步后的组分密度
  !        di 组分密度， 单位kg/m3; dinew  反应后组分密度
  !        ci 摩尔浓度， 单位  mol/cm3,   ri 反应速率 单位 mol/(cm3.s)         ！ 内部使用 mol-cm 值单位  
	  integer:: i,j,k
!------------------------------------------------------

        call comput_T(di,E,T)    ! 由内能计算温度
	
	    do i=1,N_SPEC
		   ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! 组分摩尔浓度  (单位 mol/cm3;    mol/m3 --> mol/cm3)
        enddo

        call  Reaction_rate(ci,ri,T)    ! 已知摩尔浓度ci,温度T,  计算i 组分的生成率ri （单位： mol/ (cm3.s) )
	  
	     do i=1,N_Spec
	       dinew(i)=di(i)+ri(i) *SPEC(i)%Mi *1.d6 * dt             	! 时间推进 （1阶Euler) , 单位 kg/m3  ( kg/cm3--> kg/m3)   
        enddo
   
 	   call update_Et(di,dinew,E,T)   ! 更新内能 (由于化学反应增加的内能)

	  do i=1,N_Spec
	   di(i)=dinew(i)
	  enddo

	  end 



!-----采用2nd Runge-Kutta 计算 
! di in kg/m3 ;   ci in mol/cm3
  subroutine chemical_RK2(di,E,dt)
    use CHEM
    implicit none
    real(PRE_EC):: E,T,dt
    real(PRE_EC),dimension(N_Spec)::di, di1, dinew, ci, ri      ! 密度，摩尔浓度, t时间步后的组分密度
    integer:: i,j,k
    !------------------------------------------------------

    !--------step 1-------------------------- 
    call comput_T(di,E,T)    ! 由内能计算温度

    do i=1,N_SPEC
	  ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! 组分摩尔浓度  (单位 mol/cm3;    mol/m3 --> mol/cm3)
    end do
    call  Reaction_rate(ci,ri,T)    ! 已知摩尔浓度ci,温度T,  计算i 组分的生成率ri （单位： mol/ (cm3.s) )

    do i=1,N_Spec
      di1(i)=di(i) + ri(i) *SPEC(i)%Mi *1.d6 * dt       ! 时间推进 (step 1) 
    end do
   
    call update_Et(di,di1,E,T )   ! 更新内能
   
    !--------step 2-------------------------- 
    call comput_T(di,E,T)       ! 计算温度

    do i=1,N_SPEC
      ci(i)=di1(i)/SPEC(i)%Mi  * 1.d-6          ! 组分摩尔浓度 (mol/m3 --> mol/cm3)
    end do

    call Reaction_rate(ci,ri,T)    ! 已知摩尔浓度ci,温度T, 计算反应物质量生成率 ri
     
    do i=1,N_Spec
      dinew(i)=0.5d0*di(i)+0.5d0*di1(i)+0.5d0*  ri(i) *SPEC(i)%Mi *1.d6 * dt      ! 时间推进 (step 2) 
    end do
 
    call update_ET(di1,dinew,E,T )   ! 更新内能
    !--------------------------------------------------------------------
    do i=1,N_Spec             ! 更新组分密度
      di(i)=dinew(i)
    end do
  end 

!------------------------------------------------------

	  subroutine chemical_Euler_implicit(di,E,dt)
      use CHEM
 	  implicit none
	  real(PRE_EC):: E,T,dt
	  real(PRE_EC),dimension(N_Spec)::di, ci, dinew, ri, di1,dd     ! 密度，摩尔浓度, t时间步后的组分密度
  !        di 组分密度， 单位kg/m3; dinew  反应后组分密度
  !        ci 摩尔浓度， 单位  mol/cm3,   ri 反应速率 单位 mol/(cm3.s)         ！ 内部使用 mol-cm 值单位  
	  integer:: i,j,k
	   real(PRE_EC),dimension(N_Spec,N_Spec):: Arate ,A1  ! d(ri)/d(cj)
!------------------------------------------------------

        call comput_T(di,E,T)    ! 由内能计算温度
	
	    do i=1,N_SPEC
		   ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! 组分摩尔浓度  (单位 mol/cm3;    mol/m3 --> mol/cm3)
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
		di1(i)=ri(i)*SPEC(i)%Mi*1.d6*dt    ! 显式部分；  单位 kg/m3  (1.d6: 单位转换 kg/cm3--> kg/m3 ) 
		enddo 
		
		call Gauss(N_SPEC,A1,di1,dd)  ! 求解线性代数方程组(A1*dd=di1)  (dd=dinew-di)
	  
	     do i=1,N_Spec
	       dinew(i)=di(i)+dd(i)           	! 时间推进 （1阶Euler) , 单位 kg/m3  ( kg/cm3--> kg/m3)   
        enddo
   
 	   call update_Et(di,dinew,E,T)   ! 更新内能 (由于化学反应增加的内能)

	  do i=1,N_Spec
	   di(i)=dinew(i)
	  enddo

	  end 
! ==========================================================================
!  Ax=b  Solver                                                            =
!  Gauss  Elimination with Maximal Column Pivoting                         =            
!  Gauss  列主元消去法， 求解Ax=b                                          =
!  由于消元过程，计算后 A和b 的值会发生改变                                =
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
 
 
!-----采用3rd Runge-Kutta 计算 
! di in kg/m3 ;   ci in mol/cm3
  subroutine chemical_RK3(di,E,dt)
    use CHEM
    implicit none
    real(PRE_EC):: E,T,dt
    real(PRE_EC),dimension(N_Spec)::di,din, dinew, ci, ri      ! 密度，摩尔浓度, t时间步后的组分密度
    integer:: i,j,k,KRK
    !------------------------------------------------------
    real(kind=OCFD_REAL_KIND):: Ralfa(3),Rbeta(3) 
   
        Ralfa(1)=0.d0 ;       Rbeta(1)=1.d0
        Ralfa(2)=3.d0/4.d0 ;  Rbeta(2)=1.d0/4.d0
        Ralfa(3)=1.d0/3.d0 ;  Rbeta(3)=2.d0/3.d0
		
	din(:)=di(:)	
   do KRK=1,3		
    !--------step 1-------------------------- 
    call comput_T(di,E,T)    ! 由内能计算温度

    do i=1,N_SPEC
	  ci(i)=di(i)/SPEC(i)%Mi *1.d-6           ! 组分摩尔浓度  (单位 mol/cm3;    mol/m3 --> mol/cm3)
    enddo
    call  Reaction_rate(ci,ri,T)    ! 已知摩尔浓度ci,温度T,  计算i 组分的生成率ri （单位： mol/ (cm3.s) )

     do i=1,N_Spec
      dinew(i)= Ralfa(KRK)*din(i)+ Rbeta(KRK)*(di(i) +    ri(i) *SPEC(i)%Mi *1.d6 * dt)       ! 时间推进 (step K) 
     enddo
	 
    call update_Et(di,dinew,E,T )   ! 更新内能
    di=dinew 
   enddo  
 end 

!------------------------------------------------------
 
 