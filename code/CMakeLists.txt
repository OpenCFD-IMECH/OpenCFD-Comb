message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

# Set the root directory for the source files
set(CODE_DIR ${CMAKE_SOURCE_DIR}/code)

# Add directories containing source files
set(SOURCE_DIRECTORIES
    ${CODE_DIR}/boundary-condition
    ${CODE_DIR}/Chemical
    ${CODE_DIR}/flow-solver
    ${CODE_DIR}/head
    ${CODE_DIR}/mpi-io
    ${CODE_DIR}/others
    ${CODE_DIR}/schemes
)

message(STATUS "SOURCE_DIRECTORIES ${SOURCE_DIRECTORIES}")

# Collect all Fortran source files recursively
file(GLOB_RECURSE Fortran_SOURCES
    ${CODE_DIR}/*.f90
    ${CODE_DIR}/*.F90
    ${CODE_DIR}/*.f
    ${CODE_DIR}/*.F 
)

# Include directories for Fortran modules
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(opencfd-comb ${Fortran_SOURCES})

target_link_libraries(opencfd-comb)
if (MPI_FOUND)
  target_link_libraries(opencfd-comb PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS opencfd-comb
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)