! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
!Copyright by Li Xinliang, lixl@lnm.imech.ac.cn
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine init_flow
    use flow_data
    use fluc_data
    use mesh_data
    use IO
    use Flux_Viscous, only: compute_mu
    implicit none

    integer :: i, j, IERR
    real(dp) :: r1, r2, r3, r4, r5, norm, norm_g

    select case (Para%RunMode)
    case (0)
        if (Para%IF_init == 0) then
            if (my_id == 0) then
                print *, "Initializing uniform flow"
            end if
            Istep = 0; tt = 0.0d0
            do j = 1, ny
            do i = 1, nx
                d(i, j) = 1.0d0
                u(i, j) = 1.0d0
                v(i, j) = 0.0d0
                T(i, j) = 1.0d0
            end do
            end do
        else if (Para%IF_init == 2) then
            if (my_id == 0) then
                print *, "Initializing quiescent flow"
            end if
            Istep = 0; tt = 0.0d0
            do j = 1, ny
            do i = 1, nx
                d(i, j) = 1.0d0
                u(i, j) = 0.0d0
                v(i, j) = 0.0d0
                T(i, j) = 1.0d0
            end do
            end do
        else if (Para%IF_init == 1) then
            if (my_id == 0) then
                print *, "Reading flow field from file"
                open (55, FILE='flow.dat', FORM='UNFORMATTED')
                read (55) Istep, tt
            end if
            !
            call MPI_BCAST(Istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
            call MPI_BCAST(tt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
            !
            call read_data(55, d)
            call read_data(55, u)
            call read_data(55, v)
            call read_data(55, T)
            !
            if (my_id == 0) then
                print *, "Read flow field ok!"
                close (55)
            end if
        else
            if (my_id == 0) then
                print *, "Initializing zero flow"
            end if
            Istep = 0; tt = 0.0d0
            do j = 1, ny
            do i = 1, nx
                d(i, j) = 1.0d0
                u(i, j) = 0.0d0
                v(i, j) = 0.0d0
                T(i, j) = 1.0d0
            end do
            end do
        end if
        !
        call apply_bc
        do j = 1, ny
        do i = 1, nx
            f(i, j, 1) = d(i, j)
            f(i, j, 2) = d(i, j)*u(i, j)
            f(i, j, 3) = d(i, j)*v(i, j)
            f(i, j, 4) = Para%Cv*d(i, j)*T(i, j)+d(i, j)*(u(i, j)**2+v(i, j)**2)/2.0d0
        end do
        end do
        !
        do j = 1, ny
        do i = 1, nx
            c(i, j) = sqrt(T(i, j))/Para%Ma
            p(i, j) = d(i, j)*c(i, j)**2/Para%gamma
            rhoe(i, j) = p(i, j)/(Para%gamma-1.0d0)+0.5d0*d(i, j)*(u(i, j)**2+v(i, j)**2)
        end do
        end do
        !
        if (Para%IF_Viscous == 1) call compute_mu
        !
        call write_flow
        !
    case (1)
        !
        call init_baseflow
        !
    case (2)
        call init_baseflow
        if (Para%IF_init == 0) then
            if (my_id == 0) then
                print *, "Initializing random fluc field"
            end if
            Istep = 0; tt = 0.0d0
            call random_seed()
            do j = 1, ny
            do i = 1, nx
                call random_number(r1)
                call random_number(r2)
                call random_number(r3)
                call random_number(r4)
                call random_number(r5)
                d_fluc(i, j) = r1
                u_fluc(i, j) = r2
                v_fluc(i, j) = r3
                T_fluc(i, j) = r4
                if (Para%wavenumber < 1.e-10) then
                    w_fluc(i, j) = 0.0d0
                else
                    w_fluc(i, j) = r5
                end if
            end do
            end do
            !
            norm = 0.0d0
            do j = 1, ny
            do i = 1, nx
                norm = norm+d_fluc(i, j)**2+u_fluc(i, j)**2+v_fluc(i, j)**2+w_fluc(i, j)**2+T_fluc(i, j)**2
            end do
            end do
            !
            call MPI_ALLREDUCE(norm, norm_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
            norm_g = sqrt(norm_g)
            !
            do j = 1, ny
            do i = 1, nx
                d_fluc(i, j) = d_fluc(i, j)/norm_g
                u_fluc(i, j) = u_fluc(i, j)/norm_g
                v_fluc(i, j) = v_fluc(i, j)/norm_g
                T_fluc(i, j) = T_fluc(i, j)/norm_g
                if (Para%wavenumber < 1.e-10) then
                    w_fluc(i, j) = 0.0d0
                else
                    w_fluc(i, j) = w_fluc(i, j)/norm_g
                end if
            end do
            end do
            !
        else if (Para%IF_init == 1) then
            if (my_id == 0) then
                print *, "Reading fluc field from file"
                open (55, FILE='fluc.dat', FORM='UNFORMATTED')
                read (55) Istep, tt
            end if
            !
            call MPI_BCAST(Istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
            call MPI_BCAST(tt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
            !
            call read_data(55, d_fluc)
            call read_data(55, u_fluc)
            call read_data(55, v_fluc)
            call read_data(55, w_fluc)
            call read_data(55, T_fluc)
            !
            if (my_id == 0) then
                print *, "Read fluc field ok!"
                close (55)
            end if
        else
            if (my_id == 0) then
                print *, "*** Error: please choose a valid initialization method! ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if
        end if
        !
        do j = 1, ny
        do i = 1, nx
            f(i, j, 1) = d_fluc(i, j)
            f(i, j, 2) = d(i, j)*u_fluc(i, j)+d_fluc(i, j)*u(i, j)
            f(i, j, 3) = d(i, j)*v_fluc(i, j)+d_fluc(i, j)*v(i, j)
            f(i, j, 4) = d(i, j)*w_fluc(i, j)
            f(i, j, 5) = d(i, j)*(u(i, j)*u_fluc(i, j)+v(i, j)*v_fluc(i, j)+Para%Cv*T_fluc(i, j))+ &
                         d_fluc(i, j)*((u(i, j)**2+v(i, j)**2)*0.5d0+Para%Cv*T(i, j))
        end do
        end do
        !
        call write_flow
        !
    case DEFAULT
        if (my_id == 0) print *, "*** Error: RunMode not supported! ***"
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select

end subroutine init_flow

subroutine init_baseflow
    use flow_data
    use mesh_data
    use Centrals
    use IO
    use Flux_Viscous, only: compute_mu, diff_mu
    implicit none

    integer :: i, j

    if (my_id == 0) then
        print *, "Initializing base flow"
        open (55, FILE='flow.dat', FORM='UNFORMATTED')
        read (55) Istep, tt
    end if

    call read_data(55, d)
    call read_data(55, u)
    call read_data(55, v)
    call read_data(55, T)

    if (my_id == 0) then
        print *, "Read base flow field ok!"
        close (55)
    end if
    Istep = 0; tt = 0.0d0

    call apply_bc

    call mpi_exchange_xy2d(d)
    call mpi_exchange_xy2d(u)
    call mpi_exchange_xy2d(v)
    call mpi_exchange_xy2d(T)

    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j) = sqrt(T(i, j))/Para%Ma
        p(i, j) = d(i, j)*c(i, j)**2/Para%gamma
        rhoe(i, j) = p(i, j)/(Para%gamma-1.0d0)+0.5d0*d(i, j)*(u(i, j)**2+v(i, j)**2)
    end do
    end do

    call mpi_exchange_xy2d(c)
    call mpi_exchange_xy2d(p)
    call mpi_exchange_xy2d(rhoe)

    if (Para%IF_Viscous == 1) then
        !
        call OCFD2d_dx0(u, dudk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(v, dvdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(T, dTdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
        !
        call OCFD2d_dy0(u, dudi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(v, dvdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(T, dTdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
        !
        call mpi_exchange_y2d(dudk)
        call mpi_exchange_y2d(dvdk)
        call mpi_exchange_y2d(dTdk)
        !
        call OCFD2d_dy0(dudk, d2udkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(dvdk, d2vdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(dTdk, d2Tdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
        !
        call OCFD2d_dx2(u, d2udk2, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx2(v, d2vdk2, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx2(T, d2Tdk2, nx, ny, Scheme%Scheme_Vis, BcI)
        !
        call OCFD2d_dy2(u, d2udi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy2(v, d2vdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy2(T, d2Tdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        !
        call compute_mu
        dmudk = diff_mu(dTdk, T)
        dmudi = diff_mu(dTdi, T)
        !
    end if

end subroutine init_baseflow
