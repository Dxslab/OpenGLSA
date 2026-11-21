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
                          BC_User_Def = 99
end module Constants

module Para_mpi
    implicit none
    integer :: my_id, npx, npy, npz, npx0, npy0, npz0, np_size
    integer :: MPI_COMM_X, MPI_COMM_Y, MPI_COMM_Z, &
               MPI_COMM_XY, MPI_COMM_YZ, MPI_COMM_XZ, &
               ID_XM1, ID_XP1, ID_YM1, ID_YP1, ID_ZM1, ID_ZP1
    integer, allocatable, dimension(:) :: i_offset, j_offset, k_offset, i_nn, j_nn, k_nn
    integer :: LAP
    integer :: TYPE_LAPX1, TYPE_LAPY1, TYPE_LAPZ1, TYPE_LAPX2, TYPE_LAPY2, TYPE_LAPZ2
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
        integer :: Iperiodic_X, Iperiodic_Y, Iperiodic_Z
        integer :: Isymmetric_X(2), Isymmetric_Y(2), Isymmetric_Z(2)
        integer :: IF_SFD, IF_Character, IF_Viscous, &
                   IBC, pp_limiter, &
                   Istep_Show, Istep_Save
        real(dp) :: Periodic_ISpan(3), Periodic_JSpan(3), Periodic_KSpan(3)
        real(dp) :: Re, Ma, gamma, Pr, T_inf, Cp, Cv
        real(dp) :: dt, End_time
        real(dp) :: BC_para(100)  ! parameters for boundary conditions
        real(dp) :: Chi, Delta    ! parameters for selective frequency damping
        !> parameters for eigenvalue solver
        integer :: nev, ncv, max_iter
        real(dp) :: tol
    end type
    !
    type(TYPE_Para) :: Para
    integer :: Istep
    real(dp) :: tt
    integer, allocatable, dimension(:, :, :) :: BcI, BcJ, BcK
    !
end module flow_para

module mesh_data
    use flow_para
    implicit none
    integer :: nx_global, ny_global, nz_global, nx, ny, nz
    real(dp), allocatable, dimension(:, :, :) :: Axx, Ayy, Azz
    real(dp), allocatable, dimension(:, :, :) :: Akx, Aky, Akz
    real(dp), allocatable, dimension(:, :, :) :: Aix, Aiy, Aiz
    real(dp), allocatable, dimension(:, :, :) :: Asx, Asy, Asz
    real(dp), allocatable, dimension(:, :, :) :: Ajac
    real(dp), allocatable, dimension(:, :, :) :: Akx1, Aky1, Akz1
    real(dp), allocatable, dimension(:, :, :) :: Aix1, Aiy1, Aiz1
    real(dp), allocatable, dimension(:, :, :) :: Asx1, Asy1, Asz1
    !> for viscous fluxes
    real(dp), allocatable, dimension(:, :, :) :: Akxk, Akyk, Akzk
    real(dp), allocatable, dimension(:, :, :) :: Aixk, Aiyk, Aizk
    real(dp), allocatable, dimension(:, :, :) :: Asxk, Asyk, Aszk
    real(dp), allocatable, dimension(:, :, :) :: Akxi, Akyi, Akzi
    real(dp), allocatable, dimension(:, :, :) :: Aixi, Aiyi, Aizi
    real(dp), allocatable, dimension(:, :, :) :: Asxi, Asyi, Aszi
    real(dp), allocatable, dimension(:, :, :) :: Akxs, Akys, Akzs
    real(dp), allocatable, dimension(:, :, :) :: Aixs, Aiys, Aizs
    real(dp), allocatable, dimension(:, :, :) :: Asxs, Asys, Aszs
end module mesh_data

module flow_data
    use flow_para
    implicit none
    real(dp), allocatable :: f(:, :, :, :), fn(:, :, :, :), du(:, :, :, :), fsfd(:, :, :, :)
    real(dp), allocatable, dimension(:, :, :) :: d, u, v, w, T, p, c, rhoe, mu, Hc
    !> Inviscid fluxes
    real(dp), allocatable, dimension(:, :) :: fpx, fmx, fpy, fmy, fpz, fmz, hhx, hhy, hhz
    real(dp), allocatable, dimension(:) :: Hc1x, c1x, d1x, u1x, v1x, w1x, A1x, A2x, A3x, hhx1, hhx2
    real(dp), allocatable, dimension(:) :: Hc1y, c1y, d1y, u1y, v1y, w1y, A1y, A2y, A3y, hhy1, hhy2
    real(dp), allocatable, dimension(:) :: Hc1z, c1z, d1z, u1z, v1z, w1z, A1z, A2z, A3z, hhz1, hhz2
    !> Viscous fluxes
    real(dp), allocatable, dimension(:, :, :) :: dmudk, dmudi, dmuds
    real(dp), allocatable, dimension(:, :, :) :: dudk, dvdk, dwdk, dTdk
    real(dp), allocatable, dimension(:, :, :) :: dudi, dvdi, dwdi, dTdi
    real(dp), allocatable, dimension(:, :, :) :: duds, dvds, dwds, dTds
    real(dp), allocatable, dimension(:, :, :) :: d2udkdi, d2vdkdi, d2wdkdi, d2Tdkdi
    real(dp), allocatable, dimension(:, :, :) :: d2udkds, d2vdkds, d2wdkds, d2Tdkds
    real(dp), allocatable, dimension(:, :, :) :: d2udids, d2vdids, d2wdids, d2Tdids
    real(dp), allocatable, dimension(:, :, :) :: d2udk2, d2vdk2, d2wdk2, d2Tdk2
    real(dp), allocatable, dimension(:, :, :) :: d2udi2, d2vdi2, d2wdi2, d2Tdi2
    real(dp), allocatable, dimension(:, :, :) :: d2uds2, d2vds2, d2wds2, d2Tds2
    !> Positivity-preserving limiter statistics
    integer :: limiter_count_rho, limiter_count_pressure
    integer :: limiter_count_total
    !> Activation map of the positivity-preserving limiter (per direction)
    integer, allocatable, dimension(:, :, :) :: pp_mark_rho_x, pp_mark_p_x
    integer, allocatable, dimension(:, :, :) :: pp_mark_rho_y, pp_mark_p_y
    integer, allocatable, dimension(:, :, :) :: pp_mark_rho_z, pp_mark_p_z
end module flow_data

module fluc_data
    use flow_para
    implicit none
    real(dp), allocatable, dimension(:, :, :) :: d_fluc, u_fluc, v_fluc, w_fluc, T_fluc, p_fluc, rhoe_fluc, mu_fluc
    !> Inviscid fluxes
    real(dp), allocatable, dimension(:) :: p2x, d2x, u2x, v2x, w2x
    real(dp), allocatable, dimension(:) :: p2y, d2y, u2y, v2y, w2y
    real(dp), allocatable, dimension(:) :: p2z, d2z, u2z, v2z, w2z
    !> Viscous fluxes
    real(dp), allocatable, dimension(:, :, :) :: dmufdk, dmufdi, dmufds
    real(dp), allocatable, dimension(:, :, :) :: dufdk, dvfdk, dwfdk, dTfdk
    real(dp), allocatable, dimension(:, :, :) :: dufdi, dvfdi, dwfdi, dTfdi
    real(dp), allocatable, dimension(:, :, :) :: dufds, dvfds, dwfds, dTfds
    real(dp), allocatable, dimension(:, :, :) :: d2ufdkdi, d2vfdkdi, d2wfdkdi, d2Tfdkdi
    real(dp), allocatable, dimension(:, :, :) :: d2ufdkds, d2vfdkds, d2wfdkds, d2Tfdkds
    real(dp), allocatable, dimension(:, :, :) :: d2ufdids, d2vfdids, d2wfdids, d2Tfdids
    real(dp), allocatable, dimension(:, :, :) :: d2ufdk2, d2vfdk2, d2wfdk2, d2Tfdk2
    real(dp), allocatable, dimension(:, :, :) :: d2ufdi2, d2vfdi2, d2wfdi2, d2Tfdi2
    real(dp), allocatable, dimension(:, :, :) :: d2ufds2, d2vfds2, d2wfds2, d2Tfds2
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
