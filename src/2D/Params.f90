! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module precision
    use mpi
    use iso_fortran_env, only: real64
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: OCFD_DATA_TYPE = MPI_DOUBLE_PRECISION
end module precision

module Constants
    use precision
    implicit none
    integer, parameter :: MAX1D = 20000, USER_DATA_LONG = 20000

    integer, parameter :: OCFD_Scheme_LINEAR5 = 1, OCFD_Scheme_WENOZ5 = 2, OCFD_Scheme_WENOZQ5 = 3, &
                          OCFD_Scheme_WENOMR5 = 4 
    integer, parameter :: del_right = 1, del_left = 2, del_none = 0
    integer, parameter :: OCFD_Scheme_CD4 = 10, OCFD_Scheme_CD6 = 11
    integer, parameter :: BC_None = 0, BC_Flat_Plate = 1, BC_Lid_Driven_Cavity = 2, &
                          BC_Sphere = 3, BC_Cylinder = 4, BC_User_Def = 99
end module Constants

module Para_mpi
    implicit none
    integer :: my_id, npx, npy, npx0, npy0, np_size
    integer :: MPI_COMM_X, MPI_COMM_Y, ID_XM1, ID_XP1, ID_YM1, ID_YP1
    integer :: i_offset(0:1024), j_offset(0:1024), i_nn(0:1024), j_nn(0:1024)
    integer :: LAP
end module Para_mpi

module Scheme_Para
    implicit none
    type TYPE_Scheme
        integer :: Scheme_Invis, Scheme_Vis
        integer :: Ka1, Kb1, Ka2, Kb2
    end type
    type(TYPE_Scheme) :: Scheme
end module Scheme_Para

module flow_para
    use precision
    use Constants
    use Para_mpi
    use Scheme_Para
    implicit none

    type TYPE_Para
        integer :: RunMode
        integer :: IF_init
        integer :: Iperiodic_X, Iperiodic_Y
        integer :: Isymmetric_X(2), Isymmetric_Y(2)
        integer :: IF_Cylindrical_Coordinate, IF_SFD, IF_Character, &
                   IF_Viscous, &
                   IBC, pp_limiter, &
                   Istep_Show, Istep_Save
        real(dp) :: Periodic_ISpan(2), Periodic_JSpan(2)
        real(dp) :: Re, Ma, gamma, Pr, T_inf, Cp, Cv
        real(dp) :: dt, End_time
        real(dp) :: BC_para(100)  ! parameters for boundary conditions
        real(dp) :: Chi, Delta    ! parameters for selective frequency damping
        !> parameters for eigenvalue solver
        integer :: nev, ncv, max_iter
        real(dp) :: tol, wavenumber
    end type
    !
    type(TYPE_Para) :: Para
    integer :: Istep
    real(dp) :: tt
    integer, allocatable, dimension(:, :) :: BcI, BcJ
    !
end module flow_para

module mesh_data
    use flow_para
    implicit none
    integer :: nx_global, ny_global, nx, ny
    real(dp), allocatable, dimension(:, :) :: Axx, Ayy
    real(dp), allocatable, dimension(:, :) :: Akx, Aky, Aix, Aiy, Ajac, Akx1, Aky1, Aix1, Aiy1
    real(dp), allocatable, dimension(:, :) :: Akxk, Akyk, Aixk, Aiyk, Akxi, Akyi, Aixi, Aiyi
    real(dp), allocatable, dimension(:, :) :: dydk, dydi ! for cylindrical coordinate viscous terms
end module mesh_data

module flow_data
    use flow_para
    implicit none
    real(dp), allocatable :: f(:, :, :), fn(:, :, :), du(:, :, :), fsfd(:, :, :)
    real(dp), allocatable, dimension(:, :) :: d, u, v, T, p, c, rhoe, mu, Hc
    !> Inviscid fluxes
    real(dp), allocatable, dimension(:, :) :: fpx, fmx, fpy, fmy, hhx, hhy
    real(dp), allocatable, dimension(:) :: Hc1x, c1x, d1x, u1x, v1x, A1x, A2x, hhx1, hhx2
    real(dp), allocatable, dimension(:) :: Hc1y, c1y, d1y, u1y, v1y, A1y, A2y, hhy1, hhy2
    real(dp), allocatable, dimension(:, :, :) :: Source_Invis
    !> Viscous fluxes
    real(dp), allocatable, dimension(:, :) :: dmudk, dmudi
    real(dp), allocatable, dimension(:, :) :: dudk, dvdk, dTdk, dudi, dvdi, dTdi
    real(dp), allocatable, dimension(:, :) :: d2udkdi, d2vdkdi, d2Tdkdi
    real(dp), allocatable, dimension(:, :) :: d2udk2, d2vdk2, d2Tdk2, d2udi2, d2vdi2, d2Tdi2
    real(dp), allocatable, dimension(:, :, :) :: Source_Vis
    !> Limiter statistics
    integer :: limiter_count_rho, limiter_count_pressure
    integer :: limiter_count_total
    !> Positivity-preserving limiter activation map (last RK sub-stage)
    integer, allocatable, dimension(:, :) :: pp_mark_rho_x, pp_mark_p_x
    integer, allocatable, dimension(:, :) :: pp_mark_rho_y, pp_mark_p_y
end module flow_data

module fluc_data
    use flow_para
    implicit none
    real(dp), allocatable, dimension(:, :) :: d_fluc, u_fluc, v_fluc, w_fluc, T_fluc, p_fluc, rhoe_fluc, mu_fluc
    !> Inviscid fluxes
    real(dp), allocatable, dimension(:) :: p2x, d2x, u2x, v2x, w2x
    real(dp), allocatable, dimension(:) :: p2y, d2y, u2y, v2y, w2y

    !> Viscous fluxes
    real(dp), allocatable, dimension(:, :) :: dmufdk, dmufdi
    real(dp), allocatable, dimension(:, :) :: dufdk, dvfdk, dwfdk, dTfdk, dufdi, dvfdi, dwfdi, dTfdi
    real(dp), allocatable, dimension(:, :) :: d2ufdkdi, d2vfdkdi, d2wfdkdi, d2Tfdkdi
    real(dp), allocatable, dimension(:, :) :: d2ufdk2, d2vfdk2, d2wfdk2, d2Tfdk2, d2ufdi2, d2vfdi2, d2wfdi2, d2Tfdi2
end module fluc_data

module ARNOLDI
    use precision
    implicit none
    include "debugF90.h"
    include "statF90.h"
    integer :: IPARAM(11), IPNTR(14)
    logical, allocatable :: select(:)
    real(dp), allocatable, dimension(:) :: RESID, WORKD, WORKL, WORKEV
    real(dp), allocatable, dimension(:, :) :: EIVEC, EIVAL
    character(LEN=1) :: BMAT
    character(LEN=2) :: WHICH
    integer :: IDO, NEV, NCV, LWORKL, INFO, NLOC, NCONV, LDV
    real(dp) :: TOL, SIGMAR, SIGMAI
    logical :: FIRST, RVEC
    real(dp), external :: DLAPY2, PDNORM2
end module ARNOLDI
