install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/flat-2.25)

# If testing active add test for flat-2.25 case

# data generation
set(case "flat-2.25-datagen")
set(case_dir "${test_dir}/flat-2.25")
file(MAKE_DIRECTORY ${case_dir})

add_executable(flat-2.25-initial-gen init3d-inlet.f90)
target_link_libraries(flat-2.25-initial-gen)
add_test(NAME ${case} COMMAND flat-2.25-initial-gen WORKING_DIRECTORY ${case_dir})

set(case "flat-2.25")
set(case_dir "${test_dir}/${case}")

file(COPY          opencfd-comb2.in DESTINATION ${case_dir})
file(COPY                 specie.in DESTINATION ${case_dir})
file(COPY              viscosity.in DESTINATION ${case_dir})
file(COPY           conductivity.in DESTINATION ${case_dir})
file(COPY diffusion_coefficients.in DESTINATION ${case_dir})
file(COPY     OCFD-grid.dat DESTINATION ${case_dir})
file(COPY     flow0d-upper-comb.dat DESTINATION ${case_dir})
file(COPY     flow1d-inlet-comb.dat DESTINATION ${case_dir})

set(MPIEXEC_FLAG "--oversubscribe")
set(NUMPROCS "8")

add_test(NAME ${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${MPIEXEC_FLAG}  $<TARGET_FILE:opencfd-comb> WORKING_DIRECTORY ${case_dir})
