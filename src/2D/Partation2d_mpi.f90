! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
!OpenCFD ver 2, CopyRight by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine partition2d_mpi
    use flow_para
    use mesh_data
    implicit none

    integer :: k, ka, npx1, npx2, npy1, npy2, my_mod1
    integer :: IERR

    call compute_partition_2d
    npx = mod(my_id, npx0) ! Block index in i-direction
    npy = my_id/npx0       ! Block index in j-direction
    !> Create communicators for i- and j-directions
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, npy, npx, MPI_COMM_X, IERR)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, npx, npy, MPI_COMM_Y, IERR)

    nx = nx_global/npx0
    ny = ny_global/npy0
    if (npx < mod(nx_global, npx0)) nx = nx+1
    if (npy < mod(ny_global, npy0)) ny = ny+1

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

    npx1 = my_mod1(npx-1, npx0)
    npx2 = my_mod1(npx+1, npx0)
    ID_XM1 = npy*npx0+npx1    ! -1 proc in x-direction
    ID_XP1 = npy*npx0+npx2    ! +1 proc in x-direction
    if (Para%Iperiodic_X == 0 .and. npx == 0) ID_XM1 = MPI_PROC_NULL     ! if not periodic, 0 node do not send mesg to npx0-1 node
    if (Para%Iperiodic_X == 0 .and. npx == npx0-1) ID_XP1 = MPI_PROC_NULL

    npy1 = my_mod1(npy-1, npy0)
    npy2 = my_mod1(npy+1, npy0)
    ID_YM1 = npy1*npx0+npx
    ID_YP1 = npy2*npx0+npx
    if (Para%Iperiodic_Y == 0 .and. npy == 0) ID_YM1 = MPI_PROC_NULL     ! if not periodic, 0 node do not send mesg to npy0-1 node
    if (Para%Iperiodic_Y == 0 .and. npy == npy0-1) ID_YP1 = MPI_PROC_NULL

end subroutine partition2d_mpi

subroutine compute_partition_2d
    use flow_para
    use mesh_data
    implicit none

    integer :: cand_npx0
    integer :: load_x, load_y, load_max
    integer :: best_npx0, best_npy0, best_load_max
    real(dp) :: ratio_xy, cand_ratio, shape_err, best_shape_err

    real(dp) :: ratio_diff
    integer :: limit, k, IERR
    integer :: test_npx0(2), test_npy0(2)
    integer :: min_block, max_npx0, max_npy0
    
    min_block = max(LAP, 5)

    ratio_xy = real(nx_global, dp)/real(max(1, ny_global), dp)
    best_npx0 = 0
    best_npy0 = 0
    best_load_max = huge(best_load_max)
    best_shape_err = huge(best_shape_err)

    limit = int(sqrt(real(np_size)))

    do cand_npx0 = 1, limit
        if (mod(np_size, cand_npx0) /= 0) cycle

        test_npx0(1) = cand_npx0
        test_npy0(1) = np_size/cand_npx0

        test_npx0(2) = np_size/cand_npx0
        test_npy0(2) = cand_npx0

        do k = 1, 2
            if (k == 2 .and. test_npx0(1) == test_npx0(2)) cycle

            !> reject decompositions whose smallest block is too thin
            if (nx_global/test_npx0(k) < min_block) cycle
            if (ny_global/test_npy0(k) < min_block) cycle

            load_x = (nx_global+test_npx0(k)-1)/test_npx0(k)
            load_y = (ny_global+test_npy0(k)-1)/test_npy0(k)
            load_max = load_x*load_y

            cand_ratio = real(test_npx0(k), dp)/real(test_npy0(k), dp)

            ratio_diff = cand_ratio/ratio_xy
            shape_err = max(ratio_diff, 1.0_dp/ratio_diff)

            if (load_max < best_load_max) then
                best_load_max = load_max
                best_shape_err = shape_err
                best_npx0 = test_npx0(k)
                best_npy0 = test_npy0(k)
            else if (load_max == best_load_max .and. shape_err < best_shape_err) then
                best_shape_err = shape_err
                best_npx0 = test_npx0(k)
                best_npy0 = test_npy0(k)
            end if
        end do
    end do

    if (best_npx0 == 0) then
        max_npx0 = nx_global/min_block
        max_npy0 = ny_global/min_block
        if (my_id == 0) then
            print *, "*** Error: cannot decompose ", nx_global, "x", ny_global, &
                " over ", np_size, " processes ***"
            print *, "Each local block must be at least ", min_block, &
                " points wide in every direction (LAP = ", LAP, ")."
            print *, "Max processes in i-direction: ", max_npx0, &
                ", in j-direction: ", max_npy0
            print *, "Use at most ", max_npx0*max_npy0, " processes for this grid."
        end if
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end if

    npx0 = best_npx0
    npy0 = best_npy0

    if (my_id == 0) then
        print *, "MPI decomposition: npx0 = ", npx0, ", npy0 = ", npy0
    end if

end subroutine compute_partition_2d

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

subroutine mpi_exchange_xy2d(f)
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)

    call mpi_exchange_faces2d(f, .true., .true.)
end subroutine mpi_exchange_xy2d

subroutine mpi_exchange_x2d(f)
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)

    call mpi_exchange_faces2d(f, .true., .false.)
end subroutine mpi_exchange_x2d

subroutine mpi_exchange_y2d(f)
    use mesh_data
    implicit none
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)

    call mpi_exchange_faces2d(f, .false., .true.)
end subroutine mpi_exchange_y2d

subroutine mpi_exchange_faces2d(f, exchange_x, exchange_y)
    use mesh_data
    implicit none

    logical, intent(IN) :: exchange_x, exchange_y
    real(dp), intent(INOUT) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)
    real(dp), asynchronous :: send_x(LAP*merge(ny, 0, exchange_x .and. ID_XM1 /= my_id), 2)
    real(dp), asynchronous :: recv_x(LAP*merge(ny, 0, exchange_x .and. ID_XM1 /= my_id), 2)
    real(dp), asynchronous :: send_y(LAP*merge(nx, 0, exchange_y .and. ID_YM1 /= my_id), 2)
    real(dp), asynchronous :: recv_y(LAP*merge(nx, 0, exchange_y .and. ID_YM1 /= my_id), 2)
    integer, parameter :: TAG_XM = 9000, TAG_XP = 8000, TAG_YM = 7000, TAG_YP = 6000
    integer :: requests(8), statuses(MPI_STATUS_SIZE, 8), ierr
    integer :: i, j, k1
    logical :: mpi_x, mpi_y

    mpi_x = exchange_x .and. ID_XM1 /= my_id
    mpi_y = exchange_y .and. ID_YM1 /= my_id

    requests = MPI_REQUEST_NULL
    if (mpi_x) then
        do j = 1, ny
        do i = 1, LAP
            k1 = (j-1)*LAP+i
            send_x(k1, 1) = f(i, j)
            send_x(k1, 2) = f(nx-LAP+i, j)
        end do
        end do
    end if
    if (mpi_y) then
        do j = 1, LAP
        do i = 1, nx
            k1 = (j-1)*nx+i
            send_y(k1, 1) = f(i, j)
            send_y(k1, 2) = f(i, ny-LAP+j)
        end do
        end do
    end if

    if (mpi_x) then
        call MPI_IRECV(recv_x(1, 1), LAP*ny, OCFD_DATA_TYPE, ID_XM1, TAG_XP, &
                       MPI_COMM_WORLD, requests(1), ierr)
        call MPI_IRECV(recv_x(1, 2), LAP*ny, OCFD_DATA_TYPE, ID_XP1, TAG_XM, &
                       MPI_COMM_WORLD, requests(2), ierr)
    end if
    if (mpi_y) then
        call MPI_IRECV(recv_y(1, 1), LAP*nx, OCFD_DATA_TYPE, ID_YM1, TAG_YP, &
                       MPI_COMM_WORLD, requests(3), ierr)
        call MPI_IRECV(recv_y(1, 2), LAP*nx, OCFD_DATA_TYPE, ID_YP1, TAG_YM, &
                       MPI_COMM_WORLD, requests(4), ierr)
    end if
    if (mpi_x) then
        call MPI_ISEND(send_x(1, 1), LAP*ny, OCFD_DATA_TYPE, ID_XM1, TAG_XM, &
                       MPI_COMM_WORLD, requests(5), ierr)
        call MPI_ISEND(send_x(1, 2), LAP*ny, OCFD_DATA_TYPE, ID_XP1, TAG_XP, &
                       MPI_COMM_WORLD, requests(6), ierr)
    end if
    if (mpi_y) then
        call MPI_ISEND(send_y(1, 1), LAP*nx, OCFD_DATA_TYPE, ID_YM1, TAG_YM, &
                       MPI_COMM_WORLD, requests(7), ierr)
        call MPI_ISEND(send_y(1, 2), LAP*nx, OCFD_DATA_TYPE, ID_YP1, TAG_YP, &
                       MPI_COMM_WORLD, requests(8), ierr)
    end if
    call MPI_WAITALL(8, requests, statuses, ierr)

    if (mpi_x) then
        do j = 1, ny
        do i = 1, LAP
            k1 = (j-1)*LAP+i
            if (ID_XM1 /= MPI_PROC_NULL) f(-LAP+i, j) = recv_x(k1, 1)
            if (ID_XP1 /= MPI_PROC_NULL) f(nx+i, j) = recv_x(k1, 2)
        end do
        end do
    end if
    if (mpi_y) then
        do j = 1, LAP
        do i = 1, nx
            k1 = (j-1)*nx+i
            if (ID_YM1 /= MPI_PROC_NULL) f(i, -LAP+j) = recv_y(k1, 1)
            if (ID_YP1 /= MPI_PROC_NULL) f(i, ny+j) = recv_y(k1, 2)
        end do
        end do
    end if
    if (exchange_x .and. .not. mpi_x) then
        f(1-LAP:0, 1:ny) = f(nx-LAP+1:nx, 1:ny)
        f(nx+1:nx+LAP, 1:ny) = f(1:LAP, 1:ny)
    end if
    if (exchange_y .and. .not. mpi_y) then
        f(1:nx, 1-LAP:0) = f(1:nx, ny-LAP+1:ny)
        f(1:nx, ny+1:ny+LAP) = f(1:nx, 1:LAP)
    end if
end subroutine mpi_exchange_faces2d

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
