install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/isotropic)

# If testing active add test for isotropic case

# data generation

# grid generation
set(case "isotropic-gridgen")
set(case_dir "${test_dir}/isotropic")
add_executable(gridgen-isotropic Grid.f90)
target_link_libraries(gridgen-isotropic)
add_test(NAME ${case} COMMAND gridgen-isotropic WORKING_DIRECTORY ${case_dir})

# initial field generation
set(case "isotropic-initgen")
set(case_dir "${test_dir}/isotropic")
add_executable(init-isotropic init-isotropic.f90)
target_link_libraries(init-isotropic)
add_test(NAME ${case} COMMAND init-isotropic WORKING_DIRECTORY ${case_dir})

# run opencfd-comb
set(case "isotropic")
set(case_dir "${test_dir}/${case}")

file(COPY          opencfd-comb2.in DESTINATION ${case_dir})
file(COPY           conductivity.in DESTINATION ${case_dir})
file(COPY diffusion_coefficients.in DESTINATION ${case_dir})
file(COPY                 specie.in DESTINATION ${case_dir})
file(COPY              viscosity.in DESTINATION ${case_dir})

set(MPIEXEC_FLAG "--oversubscribe")
set(NUMPROCS "8")

message(STATUS "MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message(STATUS "MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}")
message(STATUS "MPIEXEC_FLAG: ${MPIEXEC_FLAG}")
message(STATUS "NUMPROCS: ${NUMPROCS}")

add_test(NAME ${case} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROCS} ${MPIEXEC_FLAG}  $<TARGET_FILE:opencfd-comb> WORKING_DIRECTORY ${case_dir})
