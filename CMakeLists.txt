cmake_minimum_required(VERSION 3.18)
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ruiqing Yang and Baofeng Ma

project(OpenGLSA VERSION 1.0 LANGUAGES C Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
list(PREPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/.deps")

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(arpackng CONFIG REQUIRED)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported LANGUAGES Fortran)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ${ipo_supported})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>"
        "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-Ofast;-march=native>")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-fp-model;precise>"
        "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3;-xHost>")
endif()

# 2D solver
set(SRC_2D src/2D)
set(OPENGLSA_2D_SOURCES
    ${SRC_2D}/Params.f90
    ${SRC_2D}/Centrals.f90
    ${SRC_2D}/WENOs.f90
    ${SRC_2D}/Flux_Splitting.f90
    ${SRC_2D}/Flux_Inviscid.f90
    ${SRC_2D}/Flux_Inviscid_d.f90
    ${SRC_2D}/Flux_Inviscid_Character.f90
    ${SRC_2D}/Flux_Inviscid_Character_PP.f90
    ${SRC_2D}/Flux_Viscous.f90
    ${SRC_2D}/Flux_Viscous_d.f90
    ${SRC_2D}/Flux_Viscous_cyl_d.f90
    ${SRC_2D}/SFD.f90
    ${SRC_2D}/IO.f90
    ${SRC_2D}/Alloc_Mem.f90
    ${SRC_2D}/Init.f90
    ${SRC_2D}/Read_Mesh.f90
    ${SRC_2D}/Read_Para.f90
    ${SRC_2D}/Partation2d_mpi.f90
    ${SRC_2D}/TimeStepping.f90
    ${SRC_2D}/NS2d_Solver.f90
    ${SRC_2D}/LNS2d_Solver.f90
    ${SRC_2D}/IRAM_Solver.f90
    ${SRC_2D}/OpenGLSA.f90
    ${SRC_2D}/Apply_BC.f90
    ${SRC_2D}/BoundaryConditions.f90
)

# 3D solver
set(SRC_3D src/3D)
set(OPENGLSA_3D_SOURCES
    ${SRC_3D}/Params.f90
    ${SRC_3D}/Centrals.f90
    ${SRC_3D}/WENOs.f90
    ${SRC_3D}/Flux_Splitting.f90
    ${SRC_3D}/Flux_Inviscid.f90
    ${SRC_3D}/Flux_Inviscid_d.f90
    ${SRC_3D}/Flux_Inviscid_Character.f90
    ${SRC_3D}/Flux_Inviscid_Character_PP.f90
    ${SRC_3D}/Flux_Viscous.f90
    ${SRC_3D}/Flux_Viscous_d.f90
    ${SRC_3D}/SFD.f90
    ${SRC_3D}/IO.f90
    ${SRC_3D}/Alloc_Mem.f90
    ${SRC_3D}/Init.f90
    ${SRC_3D}/Read_Mesh.f90
    ${SRC_3D}/Read_Para.f90
    ${SRC_3D}/Partation_mpi.f90
    ${SRC_3D}/TimeStepping.f90
    ${SRC_3D}/NS_Solver.f90
    ${SRC_3D}/LNS_Solver.f90
    ${SRC_3D}/IRAM_Solver.f90
    ${SRC_3D}/OpenGLSA3d.f90
    ${SRC_3D}/Apply_BC.f90
    ${SRC_3D}/BoundaryConditions.f90
)

# Both trees define modules with identical names, so each target gets its own module directory.
add_executable(OpenGLSA ${OPENGLSA_2D_SOURCES})
add_executable(OpenGLSA3D ${OPENGLSA_3D_SOURCES})
set_target_properties(OpenGLSA PROPERTIES OUTPUT_NAME openglsa
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/2D")
set_target_properties(OpenGLSA3D PROPERTIES OUTPUT_NAME openglsa3d
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/3D")
foreach(tgt IN ITEMS OpenGLSA OpenGLSA3D)
    target_link_libraries(${tgt} PRIVATE
        PARPACK::PARPACK MPI::MPI_Fortran LAPACK::LAPACK BLAS::BLAS)
endforeach()

foreach(util IN ITEMS x2op reset_step op2tec x2op3d op2tec3d)
    add_executable(${util} utils/${util}.f90)
endforeach()

option(BUILD_TECIO_UTILS "Build the binary Tecplot converters" ON)
if(BUILD_TECIO_UTILS)
    enable_language(CXX)
    find_package(Threads REQUIRED)
    find_path(TECIO_INCLUDE_DIR tecio.f90 PATH_SUFFIXES tecio)
    find_library(TECIO_LIBRARY NAMES tecio)
    if(NOT TECIO_INCLUDE_DIR OR NOT TECIO_LIBRARY)
        message(FATAL_ERROR "TecIO not found. Run bash install_deps.sh {gnu|intel}, or configure with -DBUILD_TECIO_UTILS=OFF.")
    endif()
    add_library(TecIO::TecIO STATIC IMPORTED)
    set_target_properties(TecIO::TecIO PROPERTIES
        IMPORTED_LOCATION "${TECIO_LIBRARY}"
        IMPORTED_LINK_INTERFACE_LANGUAGES CXX
        INTERFACE_INCLUDE_DIRECTORIES "${TECIO_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES Threads::Threads)
    foreach(util IN ITEMS op2plt op2plt3d)
        add_executable(${util} utils/${util}.f90)
        # The Fortran driver supplies the startup code for the Fortran main program.
        set_target_properties(${util} PROPERTIES LINKER_LANGUAGE Fortran)
        target_link_libraries(${util} PRIVATE TecIO::TecIO)
        if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
            target_link_options(${util} PRIVATE -cxxlib)
        endif()
    endforeach()
endif()
