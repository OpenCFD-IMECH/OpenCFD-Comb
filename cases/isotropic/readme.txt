运行方法
1） 编译运行Grid.f90 生成网格
2） 编译运行 init-isotropic.f90 生成初值
  ln -s opencfd-comb0.dat opencfd-comb.dat
3) 运行opencfd-comb 2.x
4) 统计结果在statistics.dat 中
----------------------------------------------
注： OpenCFD-Comb 采用有量纲计算， 参考量为：

Reference Values:
 Lref=1mm
 Uref=100 m/s 
 Dref=0.0878271 Kg/m3
 Tref=288.15 K 
 Amu0 (at Tref)= 1.789E-5 Pa.s
Re_ref=490.92845
Ma_ref=0.29323717
Re_lamda=72.000
