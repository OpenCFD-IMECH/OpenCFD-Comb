! 3D Compressible N-S Finite Difference Solver 
! CopyRight by Li Xinliang  Email: lixl@imech.ac.cn
! Version 2.x   2021-2
!---------------------------------------------------------------------
   subroutine   NS_solver
   use flow_data
   implicit none
   real*8:: wall_time 
   integer:: KRK,i,j,k,m,ierr 
   real(kind=OCFD_REAL_KIND):: dt1,dt2             ! dt1=dt/2;  dt2=dt/M_step
!-----------------------------------------------------------------------
   call   allocate_flow_data        ! f, fn, d,u,v,w,T, Axx,Ayy,....,Ajac 
   call   allocate_inviscous_data   ! work data for inviscous terms 
   call   allocate_vicous_data      ! work data for viscous terms
   call   init3d
   
   call   comput_flow_variables 
   call   OCFD_bc          ! boundary condition 
   
    dt1=0.5*Para%dt
	dt2=Para%dt/Para%Chemical_sub_step
	
   if(my_id.eq.0) print*, 'init ok'
      wall_time=MPI_wtime() 
!c-----------------------------------------------------------------------        
     do while (tt <   Para%End_time-Para%dt*1.d-4 )       ! -dt*1.d-4 , considering rounding error    

	   call flow_time_advance(dt1)
	   
	   if(Para%IF_Reaction == 1) then 
	   do m=1,Para%Chemical_sub_step
	   call Chemical_time_advance(dt2)
	   enddo
	   call comput_T_P         ! renew T and P
	   endif

	   call  flow_time_advance(dt1)  
!c ---------------4. Loop of t ------------------

	    Istep=Istep+1
        tt=tt+Para%dt

	   call filtering	
      
	  if( mod(Istep,Para%Istep_show).eq.0  ) then 
       call show_flow_msg(wall_time)  ! CPU time, Total energy, Kinetic energy
      endif 
		
 
        call OCFD_analysis
        if(mod(Istep, Para%Istep_Save).eq.0)   then 
		  call save_flow_data  
        endif 

!       call MPI_barrier(MPI_COMM_WORLD,ierr)

       enddo          
!c---------------------------------------------------------------------    

       call save_flow_data 
       if(my_id.eq.0) print*, 'OK The END of opencfd'

      call deallocate_flow_data
	  call deallocate_inviscous_data     
	  call deallocate_vicous_data
	 end

!-----------------------------------------------------------------------------------

