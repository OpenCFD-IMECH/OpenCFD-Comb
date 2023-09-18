  subroutine handle_NegativeEt(i,j,k,Num_NegT)
  use flow_data 
  implicit none 
  integer:: i,j,k,i1,j1,k1,Num_NegT
  integer,parameter:: Max_NegT=10
  real(kind=OCFD_REAL_KIND):: Et0
  ! input you code here !!!
  i1=i_offset(npx)+i-1
  j1=j_offset(npy)+j-1
  k1=k_offset(npz)+k-1 
  print*, "Negative T found !",  i1,j1,k1 
  Et0=Et(i,j,k)
  Et(i,j,k)=(Et(i+1,j,k)+Et(i-1,j,k)+Et(i,j+1,k)+Et(i,j-1,k)+ Et(i,j,k+1)+Et(i,j,k-1))/6.d0
  f(i,j,k,5)=f(i,j,k,5)+(Et(i,j,k)-Et0)
  Num_NegT=Num_NegT+1
  if(Num_NegT > 10) then 
   print*, "Number of Negative Temperature Points >  Limit, STOP !"  
   stop  
  endif 
  end 
  
