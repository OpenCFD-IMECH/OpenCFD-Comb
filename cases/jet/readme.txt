1） 运行jet-comb-mesh.f90 生成网格
2)  运行jet-comb-inlet.f90 生成入口数据文件flow2d-inlet-comb.dat 及初值文件opencfd-comb0.dat
    ln -s opencfd-comb0.dat opencfd-comb.dat
3) 运行opencfd-comb 
4) 后处理， 可以使用../../util/comb-post3d-1.0.f90
