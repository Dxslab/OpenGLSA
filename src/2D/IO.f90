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

    integer, parameter :: NCOL_INT_PPLIM = 6

contains

    subroutine write_flow
        use flow_data
        use fluc_data
        implicit none
        character(LEN=100) :: filename

        if (my_id == 0) then
            if (Para%RunMode == 0) then
                write (filename, "('flow-',I8.8,'.dat')") Istep
                print *, "Writing flow data to file: ", filename
            else if (Para%RunMode == 2) then
                write (filename, "('fluc-',I8.8,'.dat')") Istep
                print *, "Writing fluc data to file: ", filename
            end if
            open (55, FILE=filename, FORM='UNFORMATTED')
            write (55) Istep, tt
        end if
        if (Para%RunMode == 0) then
            call write_data(55, d)
            call write_data(55, u)
            call write_data(55, v)
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

    subroutine write_resid_field
        use flow_data, only: f, fn
        implicit none
        integer :: i, j, m
        real(dp), allocatable :: R(:, :), xg(:, :), yg(:, :), Rg(:, :, :)

        allocate (R(nx, ny), xg(nx_global, ny_global), yg(nx_global, ny_global), &
                  Rg(nx_global, ny_global, 4))

        R = Axx(1:nx, 1:ny)
        call gather_mpi(R, xg)
        R = Ayy(1:nx, 1:ny)
        call gather_mpi(R, yg)
        do m = 1, 4
            do j = 1, ny
            do i = 1, nx
                R(i, j) = abs(f(i, j, m)-fn(i, j, m))/Para%dt
            end do
            end do
            call gather_mpi(R, Rg(:, :, m))
        end do

        if (my_id == 0) then
            print *, "Writing residual field to file: resid.dat"
            open (56, FILE='resid.dat')
            write (56, '(A)') 'TITLE = "Residual"'
            write (56, '(A)') 'VARIABLES = "X" "Y" "R_mean" "R_rho" "R_rhou" "R_rhov" "R_rhoE"'
            write (56, '(A,I0,A,I0,A,I0,A,ES20.12,A)') 'ZONE T="Istep=', Istep, '", I=', nx_global, &
                ', J=', ny_global, ', SOLUTIONTIME=', tt, ', DATAPACKING=POINT'
            do j = 1, ny_global
            do i = 1, nx_global
                write (56, '(7ES17.8E3)') xg(i, j), yg(i, j), sum(Rg(i, j, :))/4.0d0, Rg(i, j, 1:4)
            end do
            end do
            close (56)
        end if

        deallocate (R, xg, yg, Rg)
    end subroutine write_resid_field

    subroutine write_pp_limiter_map
        use flow_data, only: pp_mark_rho_x, pp_mark_p_x, pp_mark_rho_y, pp_mark_p_y
        implicit none
        integer :: Status(MPI_status_Size)
        character(LEN=100) :: filename
        integer :: i, j, k, np, n_loc, n_rec, IERR
        integer, allocatable :: ibuf(:), ibuf_recv(:)
        real(dp), allocatable :: rbuf(:), rbuf_recv(:)

        if (.not. allocated(pp_mark_rho_x)) return

        n_loc = 0
        do j = 1, ny
        do i = 1, nx
            if (pp_mark_rho_x(i, j)+pp_mark_p_x(i, j)+ &
                pp_mark_rho_y(i, j)+pp_mark_p_y(i, j) > 0) n_loc = n_loc+1
        end do
        end do

        allocate (ibuf(max(1, NCOL_INT_PPLIM*n_loc)), rbuf(max(1, 2*n_loc)))

        k = 0
        do j = 1, ny
        do i = 1, nx
            if (pp_mark_rho_x(i, j)+pp_mark_p_x(i, j)+ &
                pp_mark_rho_y(i, j)+pp_mark_p_y(i, j) > 0) then
                ibuf(NCOL_INT_PPLIM*k+1) = i_offset(npx)+i-1
                ibuf(NCOL_INT_PPLIM*k+2) = j_offset(npy)+j-1
                ibuf(NCOL_INT_PPLIM*k+3) = pp_mark_rho_x(i, j)
                ibuf(NCOL_INT_PPLIM*k+4) = pp_mark_p_x(i, j)
                ibuf(NCOL_INT_PPLIM*k+5) = pp_mark_rho_y(i, j)
                ibuf(NCOL_INT_PPLIM*k+6) = pp_mark_p_y(i, j)
                rbuf(2*k+1) = Axx(i, j)
                rbuf(2*k+2) = Ayy(i, j)
                k = k+1
            end if
        end do
        end do

        if (my_id == 0) then
            filename = 'pplim.txt'
            print *, "Writing PP limiter map to file: ", trim(filename)
            open (57, FILE=filename, STATUS='REPLACE')
            write (57, '(A,I8,A,ES20.12)') '# PP limiter activations   Istep = ', Istep, '   tt = ', tt
            write (57, '(A)') '# i_glob  j_glob'// &
                '                x                y  n_rho_x    n_p_x  n_rho_y    n_p_y'
            call write_pp_records(57, n_loc, ibuf, rbuf)
            do np = 1, np_size-1
                call MPI_RECV(n_rec, 1, MPI_INTEGER, np, 9100, MPI_COMM_WORLD, Status, IERR)
                if (n_rec > 0) then
                    allocate (ibuf_recv(NCOL_INT_PPLIM*n_rec), rbuf_recv(2*n_rec))
                    call MPI_RECV(ibuf_recv, NCOL_INT_PPLIM*n_rec, MPI_INTEGER, np, 9101, MPI_COMM_WORLD, Status, IERR)
                    call MPI_RECV(rbuf_recv, 2*n_rec, MPI_DOUBLE_PRECISION, np, 9102, MPI_COMM_WORLD, Status, IERR)
                    call write_pp_records(57, n_rec, ibuf_recv, rbuf_recv)
                    deallocate (ibuf_recv, rbuf_recv)
                end if
            end do
            close (57)
        else
            call MPI_SEND(n_loc, 1, MPI_INTEGER, 0, 9100, MPI_COMM_WORLD, IERR)
            if (n_loc > 0) then
                call MPI_SEND(ibuf, NCOL_INT_PPLIM*n_loc, MPI_INTEGER, 0, 9101, MPI_COMM_WORLD, IERR)
                call MPI_SEND(rbuf, 2*n_loc, MPI_DOUBLE_PRECISION, 0, 9102, MPI_COMM_WORLD, IERR)
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
            write (file_no, '(2I8,2ES17.9,4I9)') &
                ibuf(NCOL_INT_PPLIM*k+1), ibuf(NCOL_INT_PPLIM*k+2), rbuf(2*k+1), rbuf(2*k+2), &
                ibuf(NCOL_INT_PPLIM*k+3), ibuf(NCOL_INT_PPLIM*k+4), ibuf(NCOL_INT_PPLIM*k+5), ibuf(NCOL_INT_PPLIM*k+6)
        end do
    end subroutine write_pp_records

    subroutine write_data(file_no, U)
        integer :: file_no, i, j
        real(dp) :: U(1-LAP:nx+LAP, 1-LAP:ny+LAP)
        real(dp), allocatable :: U1(:, :)
        allocate (U1(nx, ny))

        do j = 1, ny
        do i = 1, nx
            U1(i, j) = U(i, j)
        end do
        end do
        call write_mpi(file_no, U1)
        deallocate (U1)
    end subroutine write_data

    subroutine read_data(file_no, U)
        integer :: file_no, i, j
        real(dp) :: U(1-LAP:nx+LAP, 1-LAP:ny+LAP)
        real(dp), allocatable :: U1(:, :)
        allocate (U1(nx, ny))

        call read_mpi(file_no, U1)
        do j = 1, ny
        do i = 1, nx
            U(i, j) = U1(i, j)
        end do
        end do
        deallocate (U1)
    end subroutine read_data

    subroutine write_mpi(file_no, U)
        integer :: file_no
        real(dp) :: U(nx, ny)
        real(dp), allocatable :: buff2d(:, :)

        allocate (buff2d(nx_global, ny_global))
        call gather_mpi(U, buff2d)
        if (my_id == 0) write (file_no) buff2d
        deallocate (buff2d)
    end subroutine write_mpi

    subroutine gather_mpi(U, buff2d)
        integer :: Status(MPI_status_Size)
        integer :: i, j, ki, kj, np_recv, i1, j1, ia, IERR
        real(dp) :: U(nx, ny), buff2d(nx_global, ny_global)
        real(dp), allocatable :: buff1(:), buff2(:, :)

        allocate (buff1(nx*ny), buff2(nx, ny))

        if (my_id == 0) then
            do kj = 0, npy0-1
            do ki = 0, npx0-1

                np_recv = kj*npx0+ki
                if (np_recv == 0) then
                    do i1 = 1, nx
                    do j1 = 1, ny
                        buff2d(i1, j1) = U(i1, j1)
                    end do
                    end do
                else
                    call MPI_RECV(buff1, i_nn(ki)*j_nn(kj), MPI_DOUBLE_PRECISION, np_recv, 9000, MPI_COMM_WORLD, status, IERR)
                    do j1 = 1, j_nn(kj)
                    do i1 = 1, i_nn(ki)
                        i = i_offset(ki)+i1-1
                        j = j_offset(kj)+j1-1
                        ia = (j1-1)*i_nn(ki)+i1
                        buff2d(i, j) = buff1(ia)
                    end do
                    end do
                end if
            end do
            end do
        else
            do i1 = 1, nx
            do j1 = 1, ny
                buff2(i1, j1) = U(i1, j1)
            end do
            end do
            call MPI_SEND(buff2(1, 1), nx*ny, MPI_DOUBLE_PRECISION, 0, 9000, MPI_COMM_WORLD, IERR)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        deallocate (buff1, buff2)
    end subroutine gather_mpi

    subroutine read_mpi(file_no, U)
        integer :: Status(MPI_status_Size)
        integer :: file_no, i, j, ki, kj, np_send, i1, j1, ia, IERR
        real(dp) :: U(nx, ny)
        real(dp), allocatable :: buff2d(:, :), buff1(:), buff2(:, :)

        allocate (buff2d(nx_global, ny_global), buff1(nx*ny), buff2(nx, ny))

        if (my_id == 0) then
            read (file_no) buff2d

            do kj = 0, npy0-1
            do ki = 0, npx0-1
                np_send = kj*npx0+ki
                if (np_send == 0) then
                    do j1 = 1, ny
                    do i1 = 1, nx
                        U(i1, j1) = buff2d(i1, j1)
                    end do
                    end do
                else
                    do j1 = 1, j_nn(kj)
                    do i1 = 1, i_nn(ki)
                        i = i_offset(ki)+i1-1
                        j = j_offset(kj)+j1-1
                        ia = (j1-1)*i_nn(ki)+i1
                        buff1(ia) = buff2d(i, j)
                    end do
                    end do
                    call MPI_SEND(buff1, i_nn(ki)*j_nn(kj), MPI_DOUBLE_PRECISION, np_send, 9000, MPI_COMM_WORLD, IERR)
                end if
            end do
            end do
        else
            call MPI_RECV(buff2, nx*ny, MPI_DOUBLE_PRECISION, 0, 9000, MPI_COMM_WORLD, Status, IERR)
            do i1 = 1, nx
            do j1 = 1, ny
                U(i1, j1) = buff2(i1, j1)
            end do
            end do
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        deallocate (buff2d, buff1, buff2)
    end subroutine read_mpi

end module IO
