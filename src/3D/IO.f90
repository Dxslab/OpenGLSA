! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
!Copyright by Li Xinliang, lixl@lnm.imech.ac.cn
! Copyright by Li Xinliang
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module IO
    use flow_para
    use mesh_data
    implicit none

    integer, parameter :: NCOL_INT_PPLIM = 9

contains

    subroutine write_flow
        use flow_data
        use fluc_data
        implicit none
        character(LEN=100) :: filename

        if (my_id == 0) then
            if (Para%RunMode == 0) then
                write (filename, "('flow3d-',I8.8,'.dat')") Istep
                print *, "Writing flow data to file: ", filename
            else if (Para%RunMode == 2) then
                write (filename, "('fluc3d-',I8.8,'.dat')") Istep
                print *, "Writing fluc data to file: ", filename
            end if
            open (55, FILE=filename, FORM='UNFORMATTED')
            write (55) Istep, tt
        end if
        if (Para%RunMode == 0) then
            call write_data(55, d)
            call write_data(55, u)
            call write_data(55, v)
            call write_data(55, w)
            call write_data(55, T)
        else if (Para%RunMode == 2) then
            call write_data(55, d_fluc)
            call write_data(55, u_fluc)
            call write_data(55, v_fluc)
            call write_data(55, w_fluc)
            call write_data(55, T_fluc)
        end if

        if (my_id == 0) then
            close (55)
        end if

        if (Para%RunMode == 0 .and. Para%IF_Character == 2 .and. Para%pp_limiter /= 0) then
            call write_pp_limiter_map
        end if
    end subroutine write_flow

    subroutine write_pp_limiter_map
        use flow_data, only: pp_mark_rho_x, pp_mark_p_x, pp_mark_rho_y, pp_mark_p_y, pp_mark_rho_z, pp_mark_p_z
        implicit none
        integer :: Status(MPI_status_Size)
        character(LEN=100) :: filename
        integer :: i, j, k, n, np, n_loc, n_rec, IERR
        integer, allocatable :: ibuf(:), ibuf_recv(:)
        real(dp), allocatable :: rbuf(:), rbuf_recv(:)

        if (.not. allocated(pp_mark_rho_x)) return

        n_loc = 0
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            if (pp_mark_rho_x(i, j, k)+pp_mark_p_x(i, j, k)+ &
                pp_mark_rho_y(i, j, k)+pp_mark_p_y(i, j, k)+ &
                pp_mark_rho_z(i, j, k)+pp_mark_p_z(i, j, k) > 0) n_loc = n_loc+1
        end do
        end do
        end do

        allocate (ibuf(max(1, NCOL_INT_PPLIM*n_loc)), rbuf(max(1, 3*n_loc)))

        n = 0
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            if (pp_mark_rho_x(i, j, k)+pp_mark_p_x(i, j, k)+ &
                pp_mark_rho_y(i, j, k)+pp_mark_p_y(i, j, k)+ &
                pp_mark_rho_z(i, j, k)+pp_mark_p_z(i, j, k) > 0) then
                ibuf(NCOL_INT_PPLIM*n+1) = i_offset(npx)+i-1
                ibuf(NCOL_INT_PPLIM*n+2) = j_offset(npy)+j-1
                ibuf(NCOL_INT_PPLIM*n+3) = k_offset(npz)+k-1
                ibuf(NCOL_INT_PPLIM*n+4) = pp_mark_rho_x(i, j, k)
                ibuf(NCOL_INT_PPLIM*n+5) = pp_mark_p_x(i, j, k)
                ibuf(NCOL_INT_PPLIM*n+6) = pp_mark_rho_y(i, j, k)
                ibuf(NCOL_INT_PPLIM*n+7) = pp_mark_p_y(i, j, k)
                ibuf(NCOL_INT_PPLIM*n+8) = pp_mark_rho_z(i, j, k)
                ibuf(NCOL_INT_PPLIM*n+9) = pp_mark_p_z(i, j, k)
                rbuf(3*n+1) = Axx(i, j, k)
                rbuf(3*n+2) = Ayy(i, j, k)
                rbuf(3*n+3) = Azz(i, j, k)
                n = n+1
            end if
        end do
        end do
        end do

        if (my_id == 0) then
            filename = 'pplim.txt'
            print *, "Writing PP limiter map to file: ", trim(filename)
            open (57, FILE=filename, STATUS='REPLACE')
            write (57, '(A,I8,A,ES20.12)') '# PP limiter activations   Istep = ', Istep, '   tt = ', tt
            write (57, '(A)') '# i_glob  j_glob  k_glob'// &
                '                x                y                z'// &
                '  n_rho_x    n_p_x  n_rho_y    n_p_y  n_rho_z    n_p_z'
            call write_pp_records(57, n_loc, ibuf, rbuf)
            do np = 1, np_size-1
                call MPI_RECV(n_rec, 1, MPI_INTEGER, np, 9100, MPI_COMM_WORLD, Status, IERR)
                if (n_rec > 0) then
                    allocate (ibuf_recv(NCOL_INT_PPLIM*n_rec), rbuf_recv(3*n_rec))
                    call MPI_RECV(ibuf_recv, NCOL_INT_PPLIM*n_rec, MPI_INTEGER, np, 9101, MPI_COMM_WORLD, Status, IERR)
                    call MPI_RECV(rbuf_recv, 3*n_rec, MPI_DOUBLE_PRECISION, np, 9102, MPI_COMM_WORLD, Status, IERR)
                    call write_pp_records(57, n_rec, ibuf_recv, rbuf_recv)
                    deallocate (ibuf_recv, rbuf_recv)
                end if
            end do
            close (57)
        else
            call MPI_SEND(n_loc, 1, MPI_INTEGER, 0, 9100, MPI_COMM_WORLD, IERR)
            if (n_loc > 0) then
                call MPI_SEND(ibuf, NCOL_INT_PPLIM*n_loc, MPI_INTEGER, 0, 9101, MPI_COMM_WORLD, IERR)
                call MPI_SEND(rbuf, 3*n_loc, MPI_DOUBLE_PRECISION, 0, 9102, MPI_COMM_WORLD, IERR)
            end if
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        deallocate (ibuf, rbuf)
    end subroutine write_pp_limiter_map

    subroutine write_pp_records(file_no, n, ibuf, rbuf)
        implicit none
        integer, intent(in) :: file_no, n
        integer, intent(in) :: ibuf(*)
        real(dp), intent(in) :: rbuf(*)
        integer :: k

        do k = 0, n-1
            write (file_no, '(3I8,3ES17.9,6I9)') &
                ibuf(NCOL_INT_PPLIM*k+1:NCOL_INT_PPLIM*k+3), rbuf(3*k+1:3*k+3), &
                ibuf(NCOL_INT_PPLIM*k+4:NCOL_INT_PPLIM*k+9)
        end do
    end subroutine write_pp_records

    subroutine write_data(file_no, U)
        integer :: file_no, i, j, k
        real(dp) :: U(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp), allocatable :: U1(:, :, :)
        allocate (U1(nx, ny, nz))

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            U1(i, j, k) = U(i, j, k)
        end do
        end do
        end do
        call write_mpi(file_no, U1)
        deallocate (U1)
    end subroutine write_data

    subroutine read_data(file_no, U)
        integer :: file_no, i, j, k
        real(dp) :: U(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp), allocatable :: U1(:, :, :)
        allocate (U1(nx, ny, nz))

        call read_mpi(file_no, U1)
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            U(i, j, k) = U1(i, j, k)
        end do
        end do
        end do
        deallocate (U1)
    end subroutine read_data

    subroutine write_mpi(file_no, U)
        integer :: Status(MPI_status_Size)
        integer :: file_no, i, kk, nk, i0, j0, ia, i1, j1, k1, npx1, npy1, IERR
        real(dp) :: U(nx, ny, nz)
        real(dp), allocatable :: buff2d(:, :), buff1(:, :), buff2(:), buff_send(:)
        integer, allocatable :: recvcounts1(:), displs1(:), recvcounts2(:), displs2(:)

        allocate (buff2d(nx_global, ny_global), buff1(nx_global, ny))
        allocate (buff2(nx_global*ny), buff_send(nx*ny))
        allocate (recvcounts1(npy0), displs1(npy0), recvcounts2(npx0), displs2(npx0))

        if (my_id == 0) print *, 'write 3d data ...'

        do npy1 = 0, npy0-1
            recvcounts1(npy1+1) = nx_global*j_nn(npy1)
            if (npy1 == 0) then
                displs1(npy1+1) = 0
            else
                displs1(npy1+1) = displs1(npy1)+recvcounts1(npy1)
            end if
        end do

        do npx1 = 0, npx0-1
            recvcounts2(npx1+1) = i_nn(npx1)*ny
            if (npx1 == 0) then
                displs2(npx1+1) = 0
            else
                displs2(npx1+1) = displs2(npx1)+recvcounts2(npx1)
            end if
        end do

        do kk = 1, nz_global
            call get_k_node(kk, nk, k1)

            if (npz == nk) then
                do j1 = 1, ny
                do i1 = 1, nx
                    buff_send(i1+(j1-1)*nx) = U(i1, j1, k1)
                end do
                end do

                call MPI_GATHERV(buff_send, nx*ny, OCFD_DATA_TYPE, &
                                 buff2, recvcounts2, displs2, OCFD_DATA_TYPE, &
                                 0, MPI_COMM_X, IERR)

                if (npx == 0) then
                    ia = 0
                    do npx1 = 0, npx0-1
                        do i = 1, i_nn(npx1)*ny
                            i0 = i_offset(npx1)+mod(i-1, i_nn(npx1))
                            j0 = (i-1)/i_nn(npx1)+1
                            buff1(i0, j0) = buff2(ia+i)
                        end do
                        ia = ia+i_nn(npx1)*ny
                    end do

                    call MPI_GATHERV(buff1, nx_global*ny, OCFD_DATA_TYPE, &
                                     buff2d, recvcounts1, displs1, OCFD_DATA_TYPE, &
                                     0, MPI_COMM_Y, IERR)
                end if
            end if

            if (nk /= 0) then
                if (npx == 0 .and. npy == 0 .and. npz == nk) then
                    call MPI_SEND(buff2d, nx_global*ny_global, OCFD_DATA_TYPE, &
                                  0, 6666, MPI_COMM_WORLD, IERR)
                end if
                if (my_id == 0) then
                    call MPI_RECV(buff2d, nx_global*ny_global, OCFD_DATA_TYPE, &
                                  nk*(npx0*npy0), 6666, MPI_COMM_WORLD, Status, IERR)
                end if
            end if

            if (my_id == 0) then
                write (file_no) buff2d
            end if
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        deallocate (buff2d, buff1, buff2, buff_send)
        deallocate (recvcounts1, displs1, recvcounts2, displs2)
    end subroutine write_mpi

    subroutine read_mpi(file_no, U)
        integer :: Status(MPI_status_Size)
        integer :: file_no, i, kk, nk, i0, j0, ia, i1, j1, k1, npx1, npy1, IERR
        real(dp) :: U(nx, ny, nz)
        real(dp), allocatable :: buff2d(:, :), buff1(:, :), buff2(:), buff_recv(:)
        integer, allocatable :: sendcounts1(:), displs1(:), sendcounts2(:), displs2(:)

        allocate (buff2d(nx_global, ny_global), buff1(nx_global, ny))
        allocate (buff2(nx_global*ny), buff_recv(nx*ny))
        allocate (sendcounts1(npy0), displs1(npy0), sendcounts2(npx0), displs2(npx0))

        if (my_id == 0) print *, 'read 3d data ...'

        do npy1 = 0, npy0-1
            sendcounts1(npy1+1) = nx_global*j_nn(npy1)
            if (npy1 == 0) then
                displs1(npy1+1) = 0
            else
                displs1(npy1+1) = displs1(npy1)+sendcounts1(npy1)
            end if
        end do

        do npx1 = 0, npx0-1
            sendcounts2(npx1+1) = i_nn(npx1)*ny
            if (npx1 == 0) then
                displs2(npx1+1) = 0
            else
                displs2(npx1+1) = displs2(npx1)+sendcounts2(npx1)
            end if
        end do

        do kk = 1, nz_global
            call get_k_node(kk, nk, k1)

            if (my_id == 0) then
                read (file_no) buff2d
            end if

            if (nk /= 0) then
                if (my_id == 0) then
                    call MPI_SEND(buff2d, nx_global*ny_global, OCFD_DATA_TYPE, &
                                  nk*(npx0*npy0), 6666, MPI_COMM_WORLD, IERR)
                end if
                if (npx == 0 .and. npy == 0 .and. npz == nk) then
                    call MPI_RECV(buff2d, nx_global*ny_global, OCFD_DATA_TYPE, &
                                  0, 6666, MPI_COMM_WORLD, Status, IERR)
                end if
            end if

            if (npz == nk) then
                if (npx == 0) then
                    call MPI_SCATTERV(buff2d, sendcounts1, displs1, OCFD_DATA_TYPE, &
                                      buff1, nx_global*ny, OCFD_DATA_TYPE, &
                                      0, MPI_COMM_Y, IERR)

                    ia = 0
                    do npx1 = 0, npx0-1
                        do i = 1, i_nn(npx1)*ny
                            i0 = i_offset(npx1)+mod(i-1, i_nn(npx1))
                            j0 = (i-1)/i_nn(npx1)+1
                            buff2(ia+i) = buff1(i0, j0)
                        end do
                        ia = ia+i_nn(npx1)*ny
                    end do
                end if

                call MPI_SCATTERV(buff2, sendcounts2, displs2, OCFD_DATA_TYPE, &
                                  buff_recv, nx*ny, OCFD_DATA_TYPE, &
                                  0, MPI_COMM_X, IERR)

                do j1 = 1, ny
                do i1 = 1, nx
                    U(i1, j1, k1) = buff_recv(i1+(j1-1)*nx)
                end do
                end do
            end if
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        deallocate (buff2d, buff1, buff2, buff_recv)
        deallocate (sendcounts1, displs1, sendcounts2, displs2)
    end subroutine read_mpi

end module IO
