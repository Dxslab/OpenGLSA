! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
!OpenCFD ver 2, CopyRight by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module mpi_halo3d
    implicit none
    private
    public :: mpi_exchange_faces3d
contains
    subroutine mpi_exchange_faces3d(f, exchange_x, exchange_y, exchange_z)
        use mesh_data
        implicit none

        logical, intent(IN) :: exchange_x, exchange_y, exchange_z
        real(dp), asynchronous, intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        integer, parameter :: TAG_XM = 9000, TAG_XP = 8000, TAG_YM = 7000, TAG_YP = 6000, &
                              TAG_ZM = 5000, TAG_ZP = 4000
        integer :: requests(12), statuses(MPI_STATUS_SIZE, 12), ierr
        logical :: mpi_x, mpi_y, mpi_z

        mpi_x = exchange_x .and. ID_XM1 /= my_id
        mpi_y = exchange_y .and. ID_YM1 /= my_id
        mpi_z = exchange_z .and. ID_ZM1 /= my_id

        requests = MPI_REQUEST_NULL
        if (mpi_x) then
            call MPI_IRECV(f(1-LAP, 1, 1), 1, TYPE_LAPX2, ID_XM1, TAG_XP, &
                           MPI_COMM_WORLD, requests(1), ierr)
            call MPI_IRECV(f(nx+1, 1, 1), 1, TYPE_LAPX2, ID_XP1, TAG_XM, &
                           MPI_COMM_WORLD, requests(2), ierr)
        end if
        if (mpi_y) then
            call MPI_IRECV(f(1, 1-LAP, 1), 1, TYPE_LAPY2, ID_YM1, TAG_YP, &
                           MPI_COMM_WORLD, requests(3), ierr)
            call MPI_IRECV(f(1, ny+1, 1), 1, TYPE_LAPY2, ID_YP1, TAG_YM, &
                           MPI_COMM_WORLD, requests(4), ierr)
        end if
        if (mpi_z) then
            call MPI_IRECV(f(1, 1, 1-LAP), 1, TYPE_LAPZ2, ID_ZM1, TAG_ZP, &
                           MPI_COMM_WORLD, requests(5), ierr)
            call MPI_IRECV(f(1, 1, nz+1), 1, TYPE_LAPZ2, ID_ZP1, TAG_ZM, &
                           MPI_COMM_WORLD, requests(6), ierr)
        end if
        if (mpi_x) then
            call MPI_ISEND(f(1, 1, 1), 1, TYPE_LAPX2, ID_XM1, TAG_XM, &
                           MPI_COMM_WORLD, requests(7), ierr)
            call MPI_ISEND(f(nx-LAP+1, 1, 1), 1, TYPE_LAPX2, ID_XP1, TAG_XP, &
                           MPI_COMM_WORLD, requests(8), ierr)
        end if
        if (mpi_y) then
            call MPI_ISEND(f(1, 1, 1), 1, TYPE_LAPY2, ID_YM1, TAG_YM, &
                           MPI_COMM_WORLD, requests(9), ierr)
            call MPI_ISEND(f(1, ny-LAP+1, 1), 1, TYPE_LAPY2, ID_YP1, TAG_YP, &
                           MPI_COMM_WORLD, requests(10), ierr)
        end if
        if (mpi_z) then
            call MPI_ISEND(f(1, 1, 1), 1, TYPE_LAPZ2, ID_ZM1, TAG_ZM, &
                           MPI_COMM_WORLD, requests(11), ierr)
            call MPI_ISEND(f(1, 1, nz-LAP+1), 1, TYPE_LAPZ2, ID_ZP1, TAG_ZP, &
                           MPI_COMM_WORLD, requests(12), ierr)
        end if
        call MPI_WAITALL(12, requests, statuses, ierr)

        if (exchange_x .and. .not. mpi_x) then
            f(1-LAP:0, 1:ny, 1:nz) = f(nx-LAP+1:nx, 1:ny, 1:nz)
            f(nx+1:nx+LAP, 1:ny, 1:nz) = f(1:LAP, 1:ny, 1:nz)
        end if
        if (exchange_y .and. .not. mpi_y) then
            f(1:nx, 1-LAP:0, 1:nz) = f(1:nx, ny-LAP+1:ny, 1:nz)
            f(1:nx, ny+1:ny+LAP, 1:nz) = f(1:nx, 1:LAP, 1:nz)
        end if
        if (exchange_z .and. .not. mpi_z) then
            f(1:nx, 1:ny, 1-LAP:0) = f(1:nx, 1:ny, nz-LAP+1:nz)
            f(1:nx, 1:ny, nz+1:nz+LAP) = f(1:nx, 1:ny, 1:LAP)
        end if
    end subroutine mpi_exchange_faces3d
end module mpi_halo3d

subroutine partition_mpi
    use flow_para
    use mesh_data
    implicit none

    integer :: k, ka
    integer :: npx1, npx2, npy1, npy2, npz1, npz2, my_mod1
    integer :: color_x, color_y, color_z, key_x, key_y, key_z
    integer :: color_xy, color_yz, color_xz, key_xy, key_yz, key_xz
    integer :: IERR

    call compute_partition

    npx = mod(my_id, npx0)                  ! Block index in i-direction
    npy = mod(my_id/npx0, npy0)             ! Block index in j-direction
    npz = my_id/(npx0*npy0)                 ! Block index in k-direction

    !> Create communicators for coordinate lines.
    color_x = npy+npz*npy0
    color_y = npx+npz*npx0
    color_z = npx+npy*npx0
    key_x = npx
    key_y = npy
    key_z = npz
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_x, key_x, MPI_COMM_X, IERR)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_y, key_y, MPI_COMM_Y, IERR)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_z, key_z, MPI_COMM_Z, IERR)

    !> Create communicators for coordinate planes.
    color_xy = npz
    color_yz = npx
    color_xz = npy
    key_xy = npy*npx0+npx
    key_yz = npz*npy0+npy
    key_xz = npz*npx0+npx
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_xy, key_xy, MPI_COMM_XY, IERR)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_yz, key_yz, MPI_COMM_YZ, IERR)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_xz, key_xz, MPI_COMM_XZ, IERR)

    nx = nx_global/npx0
    ny = ny_global/npy0
    nz = nz_global/npz0
    if (npx < mod(nx_global, npx0)) nx = nx+1
    if (npy < mod(ny_global, npy0)) ny = ny+1
    if (npz < mod(nz_global, npz0)) nz = nz+1

    if (allocated(i_offset)) deallocate (i_offset)
    if (allocated(j_offset)) deallocate (j_offset)
    if (allocated(k_offset)) deallocate (k_offset)
    if (allocated(i_nn)) deallocate (i_nn)
    if (allocated(j_nn)) deallocate (j_nn)
    if (allocated(k_nn)) deallocate (k_nn)
    allocate (i_offset(0:npx0-1), j_offset(0:npy0-1), k_offset(0:npz0-1))
    allocate (i_nn(0:npx0-1), j_nn(0:npy0-1), k_nn(0:npz0-1))

    do k = 0, npx0-1
        ka = min(k, mod(nx_global, npx0))
        i_offset(k) = int(nx_global/npx0)*k+ka+1
        i_nn(k) = nx_global/npx0
        if (k < mod(nx_global, npx0)) i_nn(k) = i_nn(k)+1
    end do

    do k = 0, npy0-1
        ka = min(k, mod(ny_global, npy0))
        j_offset(k) = int(ny_global/npy0)*k+ka+1
        j_nn(k) = ny_global/npy0
        if (k < mod(ny_global, npy0)) j_nn(k) = j_nn(k)+1
    end do

    do k = 0, npz0-1
        ka = min(k, mod(nz_global, npz0))
        k_offset(k) = int(nz_global/npz0)*k+ka+1
        k_nn(k) = nz_global/npz0
        if (k < mod(nz_global, npz0)) k_nn(k) = k_nn(k)+1
    end do

    call New_MPI_datatype

    npx1 = my_mod1(npx-1, npx0)
    npx2 = my_mod1(npx+1, npx0)
    ID_XM1 = (npz*npy0+npy)*npx0+npx1
    ID_XP1 = (npz*npy0+npy)*npx0+npx2
    if (Para%Iperiodic_X == 0 .and. npx == 0) ID_XM1 = MPI_PROC_NULL
    if (Para%Iperiodic_X == 0 .and. npx == npx0-1) ID_XP1 = MPI_PROC_NULL

    npy1 = my_mod1(npy-1, npy0)
    npy2 = my_mod1(npy+1, npy0)
    ID_YM1 = (npz*npy0+npy1)*npx0+npx
    ID_YP1 = (npz*npy0+npy2)*npx0+npx
    if (Para%Iperiodic_Y == 0 .and. npy == 0) ID_YM1 = MPI_PROC_NULL
    if (Para%Iperiodic_Y == 0 .and. npy == npy0-1) ID_YP1 = MPI_PROC_NULL

    npz1 = my_mod1(npz-1, npz0)
    npz2 = my_mod1(npz+1, npz0)
    ID_ZM1 = (npz1*npy0+npy)*npx0+npx
    ID_ZP1 = (npz2*npy0+npy)*npx0+npx
    if (Para%Iperiodic_Z == 0 .and. npz == 0) ID_ZM1 = MPI_PROC_NULL
    if (Para%Iperiodic_Z == 0 .and. npz == npz0-1) ID_ZP1 = MPI_PROC_NULL

end subroutine partition_mpi

subroutine compute_partition
    use flow_para
    use mesh_data
    implicit none

    integer :: cand_npx0, cand_npy0, cand_npz0
    integer :: load_x, load_y, load_z, load_max
    integer :: best_npx0, best_npy0, best_npz0, best_load_max
    real(dp) :: ratio_xy, ratio_xz, cand_ratio_xy, cand_ratio_xz
    real(dp) :: ratio_diff_xy, ratio_diff_xz, shape_err, best_shape_err
    integer :: min_block, IERR

    min_block = max(LAP, 5)

    ratio_xy = real(nx_global, dp)/real(max(1, ny_global), dp)
    ratio_xz = real(nx_global, dp)/real(max(1, nz_global), dp)
    best_npx0 = 0
    best_npy0 = 0
    best_npz0 = 0
    best_load_max = huge(best_load_max)
    best_shape_err = huge(best_shape_err)

    do cand_npx0 = 1, np_size
        if (mod(np_size, cand_npx0) /= 0) cycle

        do cand_npy0 = 1, np_size/cand_npx0
            if (mod(np_size/cand_npx0, cand_npy0) /= 0) cycle

            cand_npz0 = np_size/(cand_npx0*cand_npy0)

            !> reject decompositions whose smallest block is too thin
            if (nx_global/cand_npx0 < min_block) cycle
            if (ny_global/cand_npy0 < min_block) cycle
            if (nz_global/cand_npz0 < min_block) cycle

            load_x = (nx_global+cand_npx0-1)/cand_npx0
            load_y = (ny_global+cand_npy0-1)/cand_npy0
            load_z = (nz_global+cand_npz0-1)/cand_npz0
            load_max = load_x*load_y*load_z

            cand_ratio_xy = real(cand_npx0, dp)/real(cand_npy0, dp)
            cand_ratio_xz = real(cand_npx0, dp)/real(cand_npz0, dp)
            ratio_diff_xy = cand_ratio_xy/ratio_xy
            ratio_diff_xz = cand_ratio_xz/ratio_xz
            shape_err = max(max(ratio_diff_xy, 1.0_dp/ratio_diff_xy), &
                            max(ratio_diff_xz, 1.0_dp/ratio_diff_xz))

            if (load_max < best_load_max) then
                best_load_max = load_max
                best_shape_err = shape_err
                best_npx0 = cand_npx0
                best_npy0 = cand_npy0
                best_npz0 = cand_npz0
            else if (load_max == best_load_max .and. shape_err < best_shape_err) then
                best_shape_err = shape_err
                best_npx0 = cand_npx0
                best_npy0 = cand_npy0
                best_npz0 = cand_npz0
            end if
        end do
    end do

    if (best_npx0 == 0) then
        if (my_id == 0) then
            print *, "*** Error: cannot decompose ", nx_global, "x", ny_global, &
                "x", nz_global, " over ", np_size, " processes ***"
            print *, "Each local block must be at least ", min_block, &
                " points wide in every direction (LAP = ", LAP, ")."
            print *, "Use at most ", (nx_global/min_block)*(ny_global/min_block) &
                *(nz_global/min_block), " processes for this grid."
        end if
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end if

    npx0 = best_npx0
    npy0 = best_npy0
    npz0 = best_npz0

    if (my_id == 0) then
        print *, "MPI decomposition: npx0 = ", npx0, ", npy0 = ", npy0, &
            ", npz0 = ", npz0
    end if

end subroutine compute_partition

function my_mod1(i, n)
    implicit none

    integer :: my_mod1, i, n
    if (i < 0) then
        my_mod1 = i+n
    else if (i > n-1) then
        my_mod1 = i-n
    else
        my_mod1 = i
    end if

end function my_mod1

subroutine New_MPI_datatype
    use flow_para
    use mesh_data
    implicit none

    integer :: IERR, TYPE_tmp, real_size
    integer(KIND=MPI_ADDRESS_KIND) :: stride_z

    call MPI_TYPE_SIZE(OCFD_DATA_TYPE, real_size, IERR)
    stride_z = int((nx+2*LAP)*(ny+2*LAP), KIND=MPI_ADDRESS_KIND)* &
               int(real_size, KIND=MPI_ADDRESS_KIND)

    call MPI_TYPE_VECTOR(ny, LAP, nx+2*LAP, OCFD_DATA_TYPE, TYPE_LAPX1, IERR)
    call MPI_TYPE_VECTOR(LAP, nx, nx+2*LAP, OCFD_DATA_TYPE, TYPE_LAPY1, IERR)
    call MPI_TYPE_VECTOR(LAP, nx, (nx+2*LAP)*(ny+2*LAP), OCFD_DATA_TYPE, TYPE_LAPZ1, IERR)

    call MPI_TYPE_VECTOR(ny, nx, nx+2*LAP, OCFD_DATA_TYPE, TYPE_tmp, IERR)
    call MPI_TYPE_CREATE_HVECTOR(nz, 1, stride_z, TYPE_LAPX1, TYPE_LAPX2, IERR)
    call MPI_TYPE_CREATE_HVECTOR(nz, 1, stride_z, TYPE_LAPY1, TYPE_LAPY2, IERR)
    call MPI_TYPE_CREATE_HVECTOR(LAP, 1, stride_z, TYPE_tmp, TYPE_LAPZ2, IERR)

    call MPI_TYPE_COMMIT(TYPE_LAPX1, IERR)
    call MPI_TYPE_COMMIT(TYPE_LAPY1, IERR)
    call MPI_TYPE_COMMIT(TYPE_LAPZ1, IERR)
    call MPI_TYPE_COMMIT(TYPE_LAPX2, IERR)
    call MPI_TYPE_COMMIT(TYPE_LAPY2, IERR)
    call MPI_TYPE_COMMIT(TYPE_LAPZ2, IERR)
    call MPI_TYPE_FREE(TYPE_tmp, IERR)

    call MPI_BARRIER(MPI_COMM_WORLD, IERR)

end subroutine New_MPI_datatype

subroutine mpi_exchange_xyz3d(f)
    use mpi_halo3d, only: mpi_exchange_faces3d
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

    call mpi_exchange_faces3d(f, .true., .true., .true.)
end subroutine mpi_exchange_xyz3d

subroutine mpi_exchange_yz3d(f)
    use mpi_halo3d, only: mpi_exchange_faces3d
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

    call mpi_exchange_faces3d(f, .false., .true., .true.)
end subroutine mpi_exchange_yz3d

subroutine mpi_exchange_x3d(f)
    use mpi_halo3d, only: mpi_exchange_faces3d
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

    call mpi_exchange_faces3d(f, .true., .false., .false.)
end subroutine mpi_exchange_x3d

subroutine mpi_exchange_y3d(f)
    use mpi_halo3d, only: mpi_exchange_faces3d
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

    call mpi_exchange_faces3d(f, .false., .true., .false.)
end subroutine mpi_exchange_y3d

subroutine mpi_exchange_z3d(f)
    use mpi_halo3d, only: mpi_exchange_faces3d
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)

    call mpi_exchange_faces3d(f, .false., .false., .true.)
end subroutine mpi_exchange_z3d

subroutine get_i_node(i_global, node_i, i_local)
    use Para_mpi
    implicit none

    integer :: i_global, node_i, i_local, ia

    node_i = npx0-1
    do ia = 0, npx0-2
        if (i_global >= i_offset(ia) .and. i_global < i_offset(ia+1)) node_i = ia
    end do
    i_local = i_global-i_offset(node_i)+1

end subroutine get_i_node

subroutine get_j_node(j_global, node_j, j_local)
    use Para_mpi
    implicit none

    integer :: j_global, node_j, j_local, ja

    node_j = npy0-1
    do ja = 0, npy0-2
        if (j_global >= j_offset(ja) .and. j_global < j_offset(ja+1)) node_j = ja
    end do
    j_local = j_global-j_offset(node_j)+1

end subroutine get_j_node

subroutine get_k_node(k_global, node_k, k_local)
    use Para_mpi
    implicit none

    integer :: k_global, node_k, k_local, ka

    node_k = npz0-1
    do ka = 0, npz0-2
        if (k_global >= k_offset(ka) .and. k_global < k_offset(ka+1)) node_k = ka
    end do
    k_local = k_global-k_offset(node_k)+1

end subroutine get_k_node

function get_id(npx1, npy1, npz1)
    use Para_mpi
    implicit none

    integer :: get_id, npx1, npy1, npz1

    get_id = (npz1*npy0+npy1)*npx0+npx1

end function get_id
