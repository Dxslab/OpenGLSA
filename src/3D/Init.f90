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

    integer :: i, j, k, IERR
    real(dp) :: r1, r2, r3, r4, r5, norm, norm_g

    select case (Para%RunMode)
    case (0)
        if (Para%IF_init == 0) then
            if (my_id == 0) then
                print *, "Initializing uniform flow"
            end if
            Istep = 0; tt = 0.0d0
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                d(i, j, k) = 1.0d0
                u(i, j, k) = 1.0d0
                v(i, j, k) = 0.0d0
                w(i, j, k) = 0.0d0
                T(i, j, k) = 1.0d0
            end do
            end do
            end do
        else if (Para%IF_init == 2) then
            if (my_id == 0) then
                print *, "Initializing quiescent flow"
            end if
            Istep = 0; tt = 0.0d0
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                d(i, j, k) = 1.0d0
                u(i, j, k) = 0.0d0
                v(i, j, k) = 0.0d0
                w(i, j, k) = 0.0d0
                T(i, j, k) = 1.0d0
            end do
            end do
            end do
        else if (Para%IF_init == 1) then
            if (my_id == 0) then
                print *, "Reading flow field from file"
                open (55, FILE='flow3d.dat', FORM='UNFORMATTED')
                read (55) Istep, tt
            end if
            !
            call MPI_BCAST(Istep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
            call MPI_BCAST(tt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
            !
            call read_data(55, d)
            call read_data(55, u)
            call read_data(55, v)
            call read_data(55, w)
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
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                d(i, j, k) = 1.0d0
                u(i, j, k) = 0.0d0
                v(i, j, k) = 0.0d0
                w(i, j, k) = 0.0d0
                T(i, j, k) = 1.0d0
            end do
            end do
            end do
        end if
        !
        call apply_bc
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            f(i, j, k, 1) = d(i, j, k)
            f(i, j, k, 2) = d(i, j, k)*u(i, j, k)
            f(i, j, k, 3) = d(i, j, k)*v(i, j, k)
            f(i, j, k, 4) = d(i, j, k)*w(i, j, k)
            f(i, j, k, 5) = Para%Cv*d(i, j, k)*T(i, j, k)+d(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)*0.5d0
        end do
        end do
        end do
        !
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            c(i, j, k) = sqrt(T(i, j, k))/Para%Ma
            p(i, j, k) = d(i, j, k)*c(i, j, k)**2/Para%gamma
            rhoe(i, j, k) = p(i, j, k)/(Para%gamma-1.0d0)+ &
                            0.5d0*d(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)
        end do
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
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                call random_number(r1)
                call random_number(r2)
                call random_number(r3)
                call random_number(r4)
                call random_number(r5)
                d_fluc(i, j, k) = r1
                u_fluc(i, j, k) = r2
                v_fluc(i, j, k) = r3
                T_fluc(i, j, k) = r4
                w_fluc(i, j, k) = r5
            end do
            end do
            end do
            !
            norm = 0.0d0
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                norm = norm+d_fluc(i, j, k)**2+u_fluc(i, j, k)**2+v_fluc(i, j, k)**2+w_fluc(i, j, k)**2+T_fluc(i, j, k)**2
            end do
            end do
            end do
            !
            call MPI_ALLREDUCE(norm, norm_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
            norm_g = sqrt(norm_g)
            !
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                d_fluc(i, j, k) = d_fluc(i, j, k)/norm_g
                u_fluc(i, j, k) = u_fluc(i, j, k)/norm_g
                v_fluc(i, j, k) = v_fluc(i, j, k)/norm_g
                T_fluc(i, j, k) = T_fluc(i, j, k)/norm_g
                w_fluc(i, j, k) = w_fluc(i, j, k)/norm_g
            end do
            end do
            end do
            !
        else if (Para%IF_init == 1) then
            if (my_id == 0) then
                print *, "Reading fluc field from file"
                open (55, FILE='fluc3d.dat', FORM='UNFORMATTED')
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
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            f(i, j, k, 1) = d_fluc(i, j, k)
            f(i, j, k, 2) = d_fluc(i, j, k)*u(i, j, k)+d(i, j, k)*u_fluc(i, j, k)
            f(i, j, k, 3) = d_fluc(i, j, k)*v(i, j, k)+d(i, j, k)*v_fluc(i, j, k)
            f(i, j, k, 4) = d_fluc(i, j, k)*w(i, j, k)+d(i, j, k)*w_fluc(i, j, k)
            f(i, j, k, 5) = d_fluc(i, j, k)*(Para%Cv*T(i, j, k)+(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)*0.5d0)+ &
                            d(i, j, k)*(Para%Cv*T_fluc(i, j, k)+ &
                                        (u(i, j, k)*u_fluc(i, j, k)+ &
                                         v(i, j, k)*v_fluc(i, j, k)+ &
                                         w(i, j, k)*w_fluc(i, j, k)))
        end do
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

    integer :: i, j, k

    if (my_id == 0) then
        print *, "Initializing base flow"
        open (55, FILE='flow3d.dat', FORM='UNFORMATTED')
        read (55) Istep, tt
    end if

    call read_data(55, d)
    call read_data(55, u)
    call read_data(55, v)
    call read_data(55, w)
    call read_data(55, T)

    if (my_id == 0) then
        print *, "Read base flow field ok!"
        close (55)
    end if
    Istep = 0; tt = 0.0d0

    call apply_bc

    call mpi_exchange_xyz3d(d)
    call mpi_exchange_xyz3d(u)
    call mpi_exchange_xyz3d(v)
    call mpi_exchange_xyz3d(w)
    call mpi_exchange_xyz3d(T)

    do k = 1-LAP, nz+LAP
    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j, k) = sqrt(T(i, j, k))/Para%Ma
        p(i, j, k) = d(i, j, k)*c(i, j, k)**2/Para%gamma
        rhoe(i, j, k) = p(i, j, k)/(Para%gamma-1.0d0)+0.5d0*d(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)
    end do
    end do
    end do

    call mpi_exchange_xyz3d(c)
    call mpi_exchange_xyz3d(p)
    call mpi_exchange_xyz3d(rhoe)

    if (Para%IF_Viscous == 1) then
        !
        call compute_mu
        !
        call OCFD_dx0(u, dudk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(v, dvdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(w, dwdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(T, dTdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        !
        call OCFD_dy0(u, dudi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(v, dvdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(w, dwdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(T, dTdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        !
        call OCFD_dz0(u, duds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(v, dvds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(w, dwds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(T, dTds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        !
        call mpi_exchange_yz3d(dudk)
        call mpi_exchange_yz3d(dvdk)
        call mpi_exchange_yz3d(dwdk)
        call mpi_exchange_yz3d(dTdk)
        call mpi_exchange_z3d(dudi)
        call mpi_exchange_z3d(dvdi)
        call mpi_exchange_z3d(dwdi)
        call mpi_exchange_z3d(dTdi)
        !
        call OCFD_dy0(dudk, d2udkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dvdk, d2vdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dwdk, d2wdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dTdk, d2Tdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dz0(dudk, d2udkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dvdk, d2vdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dwdk, d2wdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dTdk, d2Tdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dudi, d2udids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dvdi, d2vdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dwdi, d2wdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dTdi, d2Tdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        !
        call OCFD_dx2(u, d2udk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(v, d2vdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(w, d2wdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(T, d2Tdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dy2(u, d2udi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(v, d2vdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(w, d2wdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(T, d2Tdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dz2(u, d2uds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(v, d2vds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(w, d2wds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(T, d2Tds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        !
        dmudk = diff_mu(dTdk, T)
        dmudi = diff_mu(dTdi, T)
        dmuds = diff_mu(dTds, T)
        !
    end if

end subroutine init_baseflow
