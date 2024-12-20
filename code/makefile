f77=mpif90
#f77 = ifort
opt= -O3   
#opt=-fast
srs=./head/OCFD_precision.f90  ./Chemical/CHEM_modules.f90 ./Chemical/CHEM_reactions.f90 \
   ./Chemical/CHEM_thermdyn.f90  ./Chemical/TRANSPORT_modules.f90   \
   ./head/MODULE_parameter.f90 ./head/MODULE_flow.f90 ./head/OCFD_allocate_data.f90 \
   $(wildcard ./boundary-condition/*.f90 ./flow-solver/*.f90  ./mpi-io/*.f90 ./others/*.f90 ./schemes/*.f90) \
  ./Opencfd-comb.f90 \
  ./particle/OCFD_particle_init.f90 ./particle/OCFD_particle_time_adv.f90 \
  ./particle/OCFD_particle_redistribute.f90 ./particle/OCFD_particle_ana.f90
  
	 
OBJS=$(srs:.f90=.o)

%.o:%.f90
	$(f77) $(opt) -c -o $@ $<

default: $(OBJS)
	$(f77) -O3  -o opencfd-comb-3.2.out $(OBJS)

clean:
	rm -f *.out *.o *.mod ./head/*.o ./head/*.mod ./Chemical/*.o ./Chemical/*.mod ./boundary-condition/*.o ./flow-solver/*.o ./mpi-io/*.o ./others/*.o ./schemes/*.o
