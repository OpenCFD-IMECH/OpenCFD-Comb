! 与热力学有关的计算
! Ver 0.36, 2016-10-16, 计算E时考虑了跨Tct修正 
! Ver 0.37a, Newton法计算温度
! Ver 1.3a  2016-11-20  进行了改写 （版式修改，热力学量计算并未修改）
! Ver 1.4, 2022-2-21 1) 热力学量支持多段拟合（原版本只支持2段拟合）；  2） 温度超过拟合上限时，Cp设定为常数
!----------------------------------------------------------------------------
!  read_Spec( ) ;  
!    comput_E(di,T,E) ;    comput_Cp(d,di,T,Cp) ;  comput_T(di,E,T) ;  comput_hi(k,T,hi) 
!    IF_RealGas ==  Perfect_GAS  or Real_GAS  (默认为Real_GAS)
!    update_Et(di,dinew,E,T)
!-------------------------------------------------------------------------------

  ! 读入组分热力学特性
  subroutine read_Spec
    use CHEM
    implicit none
    TYPE (SPECIE),pointer:: Sp
 
    integer:: i,j,k
    allocate(Spec(N_Spec))
 
    !读入组分特性  
    open(99,file="specie.in")       
    read(99,*)
    read(99,*)
    read(99,*)  N_Tct, (Tct(k),k=1,N_Tct-1), Tct_max     ! N_Tct:热力学量分段拟合的段数（通常为2）； 转换温度(例如1000K)； 拟合最高温度（超过该值Cp为定值）
    read(99,*)
    do k=1,N_Spec
      SP=>Spec(k)
      read(99,*)
      read(99,*)  Sp%name                                                     
      read(99,*)  Sp%Mi                                                ! 组分摩尔质量
	  do j=1,N_Tct
      read(99,*) SP%Ai(j),SP%Bi(j),SP%Ci(j),SP%Di(j),SP%Ei(j),SP%Fi(j),SP%Gi(j)            ! 组分热力学参数（分段拟合，第j段）
      enddo
      Sp%Ri=R0/Sp%Mi
	  
!     Sp%Ect=Sp%F2-Sp%F1   ! 跨Tct的内能（或焓）差
      Sp%Ect(1)=0.d0         ! 第1段拟合，不修正
      do j=1,N_Tct-1      ! 跨拟合段的焓差  （补充到热焓拟合公式中， 使得热焓跨拟合段连续）   （理论上与Sp%F2-Sp%F1 值相同，但精度更高些）
	  Sp%Ect(j+1)= -( (Sp%Ai(j+1)-Sp%Ai(j))*Tct(j) +(Sp%Bi(j+1)-Sp%Bi(j))*Tct(j)**2/2.d0 & 
	            + (Sp%Ci(j+1)-Sp%Ci(j))*Tct(j)**3/3.d0 +(Sp%Di(j+1)-Sp%Di(j))*Tct(j)**4/4.d0 &
				+ (Sp%Ei(j+1)-Sp%Ei(j))*Tct(j)**5/5.d0 ) + Sp%Ect(j) 
 	  enddo
	enddo
    close(99)
    ! print*, "read specie.in OK"
  end subroutine read_Spec

!---------------------------------------------------------------------------------
!--------状态方程相关计算 （密度、 压力、 温度 、内能 、 比热 等） -------------------
 ! 完全气体 or 真实气体 （Cp 非常数）
 
 ! 计算内能 （不含反应焓）
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

  
! 计算定压比热   
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

! 已知内能计算温度
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

 !  组分焓
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


!----------------完全气体 Cv=5/2 R, Cp=7/2 R --------------------

!  Perfect gas  完全气体计算  (Cv=5/2 R )
  subroutine comput_E_pgas(di,T,E)            ! 完全气体
    use CHEM
    implicit none
    integer k
    real(PRE_EC):: T,E,di(N_SPEC) 
     E=0.d0
     do k=1,N_SPEC
     E=E+2.5d0*SPEC(k)%Ri*di(k)*T     ! Cv=2.5R ,  E=rho*Cv*T
     enddo
   end

! 计算定压比热 Cp   
  subroutine comput_Cp_pgas(d,di,T,Cp)
    use CHEM
    implicit none
    integer k
     real(PRE_EC):: T,Cp,d,di(N_SPEC) 
     Cp=0.d0
     do k=1,N_SPEC
      Cp=Cp+3.5d0*SPEC(k)%Ri*di(k)        ! Cp=3.5R 
     enddo
     Cp=Cp/d                 ! 密度加权平均
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
 
  
! 根据温度，计算组分焓, 不含反应焓  
  subroutine comput_hi_pgas(k,T,hi)
    use CHEM
    implicit none
    integer k
    real(PRE_EC):: T,hi 
      hi=3.5d0* SPEC(k)%Ri*T              ! Cp=7/2*R
   end

  

!---------化学反应引起的内能变化---------------------
!----------根据反应焓更新内能 （已考虑了跨Tct修正）-------------
  subroutine update_Et(di,dinew,E,T)
    use CHEM
    implicit none
    real(PRE_EC) :: E,T,F0        ! 内能（不含反应焓）， 温度， F0 反应焓
    real(PRE_EC),dimension(N_Spec) :: di, dinew, bi      ! 旧的组分密度、新的组分密度
    integer :: i,j,k
 
    do i=1,N_SPEC
!      F0=Spec(i)%F1               ! 已考虑了跨Tct修正
      F0=Spec(i)%Fi(1)             ! 已考虑了跨Tct修正, 只使用第1段拟合区的Fi值  （其他区段的已通过Ect修正，包含在热焓中）
      E=E-(dinew(i)-di(i))*SPEC(i)%Ri*F0    ! 更新内能 （加入化学反应释放/吸收的能量）
    end do
  end
  
   
! -------------真实气体 Cp= Cp (T) ----------------------------------
! 根据温度及组分密度，计算内能（不含化学能）  (E=Rho*e=Rho*Cv*T for perfect gas)
!        内能的拟合公式 Et=A*T+1/2*B*T**2+1/3*C*T**3+1/4*D*T**4+1/5*E*T**5  ;  
!        不考虑生成焓， 考虑T> Tct时的修正 （Ect=F2-F1)    
  subroutine comput_E_rgas(di,T,E)
    use CHEM
    implicit none
    integer k,j,jc
    real(PRE_EC):: T,E,R1,di(N_SPEC) 
    real(PRE_EC):: A1,B1,C1,D1,E1,Ect, Cv2
    A1=0.d0; B1=0.d0; C1=0.d0; D1=0.d0; E1=0.d0 ;  Ect=0.d0
	jc=1
	do j=1,N_Tct-1
     if(T>Tct(j)) jc=j+1   !jc 温度T所处的拟合区间段数
	enddo
	do k=1,N_SPEC
      R1=SPEC(k)%Ri
      A1=A1+(SPEC(k)%Ai(jc)-1.d0)*R1*di(k)          ! 系数，各组分加权求和 (Ai-1) for Cv
      B1=B1+SPEC(k)%Bi(jc)*R1*di(k)
      C1=C1+SPEC(k)%Ci(jc)*R1*di(k)
      D1=D1+SPEC(k)%Di(jc)*R1*di(k)
      E1=E1+SPEC(k)%Ei(jc)*R1*di(k)
      Ect=Ect+SPEC(k)%Ect(jc)*R1*di(k)
	 end do

	 if(T<=Tct_max) then
	  E=A1*T+B1*T**2/2.d0+C1*T**3/3.d0+D1*T**4/4.d0+E1*T**5/5.d0 + Ect     ! Ect: 跨Tct修正
     else 
      Cv2=A1+B1*Tct_max+C1*Tct_max**2+D1*Tct_max**3+E1*Tct_max**4               ! Cv at Tct_max
	  E=A1*Tct_max+B1*Tct_max**2/2.d0+C1*Tct_max**3/3.d0+D1*Tct_max**4/4.d0+E1*Tct_max**5/5.d0 &
 	            + Ect +Cv2*(T-Tct_max)     ! 内能 （T> Tct_max), 含跨Tct修正
	 endif
  end
!--------------------------------------------------------------------------  
! 根据温度及组分密度，计算比热Cp   (混合气体平均比热）
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
      A1=A1+(SPEC(k)%Ai(jc)-1.d0)*R1*di(k)          ! 系数，各组分加权求和 (Ai-1) for Cv
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
! 根据温度，计算组分焓 （热焓）, 不含反应焓  
!  热焓的拟合公式 ht=A*T+1/2*B*T**2+1/3*C*T**3+1/4*D*T**4+1/5*E*T**5     
	subroutine comput_hi_rgas(k,T,hi)
    use CHEM
    implicit none
    integer k,j,jc
    real(PRE_EC):: T,hi,Cp_max 
	TYPE (SPECIE),pointer:: Sp
    Sp=>SPEC(k)    ! 第 i个组分
	jc=1
	do j=1,N_Tct-1  ! 查找T所处的温度区段
     if(T>Tct(j)) jc=j+1
	enddo
    
	if(T<=Tct_max) then
      hi=(Sp%Ai(jc)*T+Sp%Bi(jc)*T**2/2.d0+Sp%Ci(jc)*T**3/3.d0+   &
	      Sp%Di(jc)*T**4/4.d0+Sp%Ei(jc)*T**5/5.d0 + Sp%Ect(jc))*Sp%Ri       ! Ect 跨Tct修正
    else 
	  Cp_max=Sp%Ai(jc)+Sp%Bi(jc)*Tct_max+Sp%Ci(jc)*Tct_max**2+Sp%Di(jc)*Tct_max**3 &
	         +Sp%Ei(jc)*Tct_max**4         ! Cp at Tct_max
      hi=(Sp%Ai(jc)*Tct_max+Sp%Bi(jc)*Tct_max**2/2.d0+Sp%Ci(jc)*Tct_max**3/3.d0+   &
	      Sp%Di(jc)*Tct_max**4/4.d0+Sp%Ei(jc)*Tct_max**5/5.d0 + Sp%Ect(jc) + Cp_max*(T-Tct_max) )*Sp%Ri 			 
   
    endif 
  end
  
  ! 计算各组分的Gibbs 自由焓 (函数连续，无需使用跨Ect修正） G = H-TS   ;       计算化学反应平衡常数时使用
  subroutine  comput_gibbs(k,T,gi)         ! gi=Gi/Ri
    use CHEM
    implicit none
	TYPE (SPECIE),pointer:: Sp
    integer:: k,j,jc
    real(PRE_EC):: T, gi, Cp_max,h_max,S_max 
     SP=>Spec(k)    ! 第 i个组分
	 jc=1
	 do j=1,N_Tct-1       ! 查找T所处的温度区段
     if(T>Tct(j)) jc=j+1
	 enddo	 
     if(T<=Tct_max) then
       gi=Sp%Ai(jc)*(1.d0-log(T))*T-SP%Bi(jc)*T**2/2.d0-SP%Ci(jc)*T**3/6.d0-SP%Di(jc)*T**4/12.d0 & 
	     -SP%Ei(jc)*T**5/20.d0+SP%Fi(jc)-SP%Gi(jc)*T
      else
	  Cp_max=Sp%Ai(jc)+Sp%Bi(jc)*Tct_max+Sp%Ci(jc)*Tct_max**2+Sp%Di(jc)*Tct_max**3+Sp%Ei(jc)*Tct_max**4         ! Cp at Tct_max	  
      h_max= Sp%Ai(jc)*Tct_max+Sp%Bi(jc)*Tct_max**2/2.d0+Sp%Ci(jc)*Tct_max**3/3.d0    &
	         +Sp%Di(jc)*Tct_max**4/4.d0+Sp%Ei(jc)*Tct_max**5/5.d0 + Sp%Fi(jc)          ! 无需跨Tct修正 （热焓+生成焓）
	  S_max=Sp%Ai(jc)*log(Tct_max)+SP%Bi(jc)*Tct_max+SP%Ci(jc)*Tct_max**2/2.d0+SP%Di(jc)*Tct_max**3/3.d0 &
	         +SP%Ei(jc)*Tct_max**4/4.d0+SP%Gi(jc)
	  gi=h_max+Cp_max*(T-Tct_max)-T*(S_max+Cp_max*log(T/Tct_max))  ! g=h-TS
	  end if
    
	end       
!--------------------------------------------------------------------------  

!  根据内能E （不含化学能）及组分密度di,计算出温度
!   Newton法 计算温度 
  
!  di : 组分密度； E 内能 （不含生成焓）； T 温度 
!  有量纲计算， 国际单位制；
  subroutine comput_T_rgas(di,E,T)
    use CHEM
    implicit none
    real(PRE_EC),parameter:: T_Cr=1.d-3        ! 收敛精度
    integer,parameter:: Kstep_lmt=30           ! 迭代次数限制
    real(PRE_EC):: E,T,p,di(N_SPEC)
    real(PRE_EC),dimension(NTct_Max):: A1,B1,C1,D1,E1,Ect
    real(PRE_EC):: R1,T1,T2,Et1,Ex,error
    integer:: k,ks,j,jc
 
    A1(1:N_Tct)=0.d0; B1(1:N_Tct)=0.d0; C1(1:N_Tct)=0.d0; D1(1:N_Tct)=0.d0; E1(1:N_Tct)=0.d0
    Ect(1:N_Tct)=0.d0
   do j=1,N_Tct   
   do k=1,N_SPEC
    !      内能的拟合公式 (Et=A*T+B*T**2+C*T**3+D*T**4+E*T**5) 。 不包括生成焓！
      R1=SPEC(k)%Ri
      A1(j)=A1(j)+(SPEC(k)%Ai(j)-1.d0)*R1*di(k)          ! 系数，各组分加权求和（低温区）
      B1(j)=B1(j)+SPEC(k)%Bi(j)*R1*di(k)
      C1(j)=C1(j)+SPEC(k)%Ci(j)*R1*di(k)
      D1(j)=D1(j)+SPEC(k)%Di(j)*R1*di(k)
      E1(j)=E1(j)+SPEC(k)%Ei(j)*R1*di(k)
      Ect(j)=Ect(j)+SPEC(k)%Ect(j)*R1*di(k)           ! 高-低温焓差 (两套拟合公式之间的差异， 以便跨越Tct时内能连续）
   enddo
   enddo
      ! newton  iteration  
  
    T1=E/A1(1)            ! 初值
  
     ks=0 ! 迭代次数

     do while (ks <=Kstep_lmt)

	 jc=1
	 do j=1,N_Tct-1       ! 查找T1所处的温度区段
     if(T1>Tct(j)) jc=j+1
	 enddo	 
	 

    if(T1<=Tct_max) then                               ! 高温区 (>1000K)
      Et1=A1(jc)*T1+B1(jc)*T1**2/2.0+C1(jc)*T1**3/3.d0+D1(jc)*T1**4/4.d0+E1(jc)*T1**5/5.d0 + Ect(jc)              ! 内能 (含跨 Tct修正) 
      Ex=A1(jc)+B1(jc)*T1+C1(jc)*T1**2+D1(jc)*T1**3+E1(jc)*T1**4       ! 导数 (Cv)
    else               ! T> Tct_max 
	  Ex=A1(jc)+B1(jc)*Tct_max+C1(jc)*Tct_max**2+D1(jc)*Tct_max**3+E1(jc)*Tct_max**4               ! Cv at Tct_max  
	  Et1=A1(jc)*Tct_max+B1(jc)*Tct_max**2/2.d0+C1(jc)*Tct_max**3/3.d0+D1(jc)*Tct_max**4/4.d0  &
	     +E1(jc)*Tct_max**5/5.d0 + Ect(jc) +Ex*(T1-Tct_max)     ! 内能 (含跨 Tct修正) 
	end if
	
     T2=T1-(Et1-E)/Ex           !  Newton iteration
     ks=ks+1
     error=abs(T2-T1)
	 T1=T2    ! 更新
    if( error < T_Cr  ) exit
    enddo 
   T=T2
   
  end

! 根据温度、计算出压力 （各分压之和）； di(k) 各组分的密度
  subroutine comput_P(di,T,p)
    use CHEM
    implicit none
    real(PRE_EC):: T,p,di(N_SPEC)
    integer:: k

    p=0.d0 
    do k=1,N_SPEC
      p=p+SPEC(k)%Ri*di(k)*T          ! 压力，分压之和
    end do
  end

! 根据压力、温度以及质量比分ai(:) ，计算密度 （2016-5-11）
 
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
    gamma=p/Et+1.d0            ! 等效比热比
    c=sqrt(gamma*p/d)          ! 近似声速 （通量分裂使用）
 end


