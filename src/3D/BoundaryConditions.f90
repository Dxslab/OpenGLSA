! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module LDC_Wall
    use mesh_data
    implicit none

    real(dp), parameter :: c_os4 = 12.0d0/25.0d0

contains

    subroutine wall_index(dir, n0, a, b, i0, j0, k0)
        integer, intent(IN) :: dir, n0, a, b
        integer, intent(OUT) :: i0, j0, k0

        select case (dir)
        case (1)
            i0 = n0; j0 = a; k0 = b
        case (2)
            i0 = a; j0 = n0; k0 = b
        case default
            i0 = a; j0 = b; k0 = n0
        end select
    end subroutine wall_index

    subroutine shift_index(dir, m, i0, j0, k0, ii, jj, kk)
        integer, intent(IN) :: dir, m, i0, j0, k0
        integer, intent(OUT) :: ii, jj, kk

        ii = i0; jj = j0; kk = k0
        select case (dir)
        case (1)
            ii = i0+m
        case (2)
            jj = j0+m
        case default
            kk = k0+m
        end select
    end subroutine shift_index

    real(dp) function metric_diag(dir, i, j, k)
        integer, intent(IN) :: dir, i, j, k

        select case (dir)
        case (1)
            metric_diag = Akx(i, j, k)
        case (2)
            metric_diag = Aiy(i, j, k)
        case default
            metric_diag = Asz(i, j, k)
        end select
    end function metric_diag

    real(dp) function wall_pressure(pp, pn, inward)
        real(dp), intent(IN) :: pp(4), pn
        integer, intent(IN) :: inward

        wall_pressure = (48.0d0*pp(1)-36.0d0*pp(2)+16.0d0*pp(3)-3.0d0*pp(4))/25.0d0-inward*c_os4*pn
    end function wall_pressure

    real(dp) function normal_stress_div(n, mu, gmu, g, h)
        integer, intent(IN) :: n
        real(dp), intent(IN) :: mu, gmu(3), g(3, 3), h(3, 3, 3)
        integer :: b
        real(dp) :: lapp, ddiv, divu, sgm

        lapp = 0.0d0; ddiv = 0.0d0; divu = 0.0d0; sgm = 0.0d0
        do b = 1, 3
            lapp = lapp+h(n, b, b)
            ddiv = ddiv+h(b, b, n)
            divu = divu+g(b, b)
            sgm = sgm+gmu(b)*(g(n, b)+g(b, n))
        end do
        normal_stress_div = mu*(lapp+ddiv/3.0d0)+sgm-2.0d0/3.0d0*gmu(n)*divu
    end function normal_stress_div

    subroutine physical_gradients(A, dA, d1, d2, dkdi, dkds, dids, g, h)
        real(dp), intent(IN) :: A(3), dA(3), d1(3, 3), d2(3, 3), dkdi(3), dkds(3), dids(3)
        real(dp), intent(OUT) :: g(3, 3), h(3, 3, 3)
        integer :: ia, ib

        do ia = 1, 3
            do ib = 1, 3
                g(ia, ib) = A(ib)*d1(ia, ib)
                h(ia, ib, ib) = A(ib)**2*d2(ia, ib)+A(ib)*dA(ib)*d1(ia, ib)
            end do
            h(ia, 1, 2) = A(1)*A(2)*dkdi(ia); h(ia, 2, 1) = h(ia, 1, 2)
            h(ia, 1, 3) = A(1)*A(3)*dkds(ia); h(ia, 3, 1) = h(ia, 1, 3)
            h(ia, 2, 3) = A(2)*A(3)*dids(ia); h(ia, 3, 2) = h(ia, 2, 3)
        end do
    end subroutine physical_gradients

    subroutine base_gradients(i, j, k, g, h, gmu, mu0)
        use flow_data
        use Flux_Viscous, only: diff_mu
        integer, intent(IN) :: i, j, k
        real(dp), intent(OUT) :: g(3, 3), h(3, 3, 3), gmu(3), mu0
        real(dp) :: A(3), dA(3), d1(3, 3), d2(3, 3), dkdi(3), dkds(3), dids(3), Tsb

        A = (/Akx(i, j, k), Aiy(i, j, k), Asz(i, j, k)/)
        dA = (/Akxk(i, j, k), Aiyi(i, j, k), Aszs(i, j, k)/)
        d1(1, :) = (/dudk(i, j, k), dudi(i, j, k), duds(i, j, k)/)
        d1(2, :) = (/dvdk(i, j, k), dvdi(i, j, k), dvds(i, j, k)/)
        d1(3, :) = (/dwdk(i, j, k), dwdi(i, j, k), dwds(i, j, k)/)
        d2(1, :) = (/d2udk2(i, j, k), d2udi2(i, j, k), d2uds2(i, j, k)/)
        d2(2, :) = (/d2vdk2(i, j, k), d2vdi2(i, j, k), d2vds2(i, j, k)/)
        d2(3, :) = (/d2wdk2(i, j, k), d2wdi2(i, j, k), d2wds2(i, j, k)/)
        dkdi = (/d2udkdi(i, j, k), d2vdkdi(i, j, k), d2wdkdi(i, j, k)/)
        dkds = (/d2udkds(i, j, k), d2vdkds(i, j, k), d2wdkds(i, j, k)/)
        dids = (/d2udids(i, j, k), d2vdids(i, j, k), d2wdids(i, j, k)/)
        call physical_gradients(A, dA, d1, d2, dkdi, dkds, dids, g, h)

        Tsb = 110.4d0/Para%T_inf
        mu0 = (1.0d0+Tsb)*sqrt(T(i, j, k)**3)/(Para%Re*(Tsb+T(i, j, k)))
        gmu = A*diff_mu((/dTdk(i, j, k), dTdi(i, j, k), dTds(i, j, k)/), T(i, j, k))
    end subroutine base_gradients

    subroutine fluc_gradients(i, j, k, g, h, gmu, mu0)
        use flow_data
        use fluc_data
        use Flux_Viscous_d, only: diff_mu_d
        integer, intent(IN) :: i, j, k
        real(dp), intent(OUT) :: g(3, 3), h(3, 3, 3), gmu(3), mu0
        real(dp) :: A(3), dA(3), d1(3, 3), d2(3, 3), dkdi(3), dkds(3), dids(3), Tsb

        A = (/Akx(i, j, k), Aiy(i, j, k), Asz(i, j, k)/)
        dA = (/Akxk(i, j, k), Aiyi(i, j, k), Aszs(i, j, k)/)
        d1(1, :) = (/dufdk(i, j, k), dufdi(i, j, k), dufds(i, j, k)/)
        d1(2, :) = (/dvfdk(i, j, k), dvfdi(i, j, k), dvfds(i, j, k)/)
        d1(3, :) = (/dwfdk(i, j, k), dwfdi(i, j, k), dwfds(i, j, k)/)
        d2(1, :) = (/d2ufdk2(i, j, k), d2ufdi2(i, j, k), d2ufds2(i, j, k)/)
        d2(2, :) = (/d2vfdk2(i, j, k), d2vfdi2(i, j, k), d2vfds2(i, j, k)/)
        d2(3, :) = (/d2wfdk2(i, j, k), d2wfdi2(i, j, k), d2wfds2(i, j, k)/)
        dkdi = (/d2ufdkdi(i, j, k), d2vfdkdi(i, j, k), d2wfdkdi(i, j, k)/)
        dkds = (/d2ufdkds(i, j, k), d2vfdkds(i, j, k), d2wfdkds(i, j, k)/)
        dids = (/d2ufdids(i, j, k), d2vfdids(i, j, k), d2wfdids(i, j, k)/)
        call physical_gradients(A, dA, d1, d2, dkdi, dkds, dids, g, h)

        ! mu' = dmu/dT*T' (zero on isothermal walls, kept for generality)
        Tsb = 110.4d0/Para%T_inf
        mu0 = (1.0d0+Tsb)/(2.0d0*Para%Re)*sqrt(T(i, j, k))*(T(i, j, k)+3.0d0*Tsb)*T_fluc(i, j, k)/((T(i, j, k)+Tsb)**2)
        gmu = A*diff_mu_d((/dTfdk(i, j, k), dTfdi(i, j, k), dTfds(i, j, k)/), T_fluc(i, j, k), &
                          (/dTdk(i, j, k), dTdi(i, j, k), dTds(i, j, k)/), T(i, j, k))
    end subroutine fluc_gradients

end module LDC_Wall

subroutine ocfd_bc_Lid_Driven_Cavity
    use mesh_data
    use flow_data
    use Flux_Viscous, only: compute_derivatives
    use LDC_Wall
    implicit none

    integer :: i, j, k
    real(dp) :: Rgas
    real(dp), allocatable :: psum(:, :, :)
    integer, allocatable :: pcnt(:, :, :)

    Rgas = 1.0d0/(Para%Ma**2*Para%gamma)

    if (npy == npy0-1) then
        do k = 1, nz
        do i = 1, nx
            u(i, ny, k) = (1.0d0-(2.0d0*Axx(i, ny, k)-1.0d0)**18)**2*(1.0d0-(2.0d0*Azz(i, ny, k)-1.0d0)**18)**2
            v(i, ny, k) = 0.0d0
            w(i, ny, k) = 0.0d0
            T(i, ny, k) = 1.0d0
        end do
        end do
    end if
    if (npy == 0) call set_wall_uvwT(2, 1)
    if (npx == 0) call set_wall_uvwT(1, 1)
    if (npx == npx0-1) call set_wall_uvwT(1, nx)
    if (npz == 0) call set_wall_uvwT(3, 1)
    if (npz == npz0-1) call set_wall_uvwT(3, nz)

    call mpi_exchange_xyz3d(u)
    call mpi_exchange_xyz3d(v)
    call mpi_exchange_xyz3d(w)
    call mpi_exchange_xyz3d(T)

    call compute_derivatives

    allocate (psum(nx, ny, nz), pcnt(nx, ny, nz))
    psum = 0.0d0
    pcnt = 0
    if (npx == 0) call add_wall_pressure(1, 1, 1)
    if (npx == npx0-1) call add_wall_pressure(1, nx, -1)
    if (npy == 0) call add_wall_pressure(2, 1, 1)
    if (npy == npy0-1) call add_wall_pressure(2, ny, -1)
    if (npz == 0) call add_wall_pressure(3, 1, 1)
    if (npz == npz0-1) call add_wall_pressure(3, nz, -1)

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
        if (pcnt(i, j, k) > 0) d(i, j, k) = psum(i, j, k)/pcnt(i, j, k)/(Rgas*T(i, j, k))
    end do
    end do
    end do
    deallocate (psum, pcnt)

contains

    subroutine set_wall_uvwT(dir, n0)
        integer, intent(IN) :: dir, n0
        integer :: a, b, i0, j0, k0

        do b = 1, merge(ny, nz, dir == 3)
        do a = 1, merge(nx, ny, dir /= 1)
            call wall_index(dir, n0, a, b, i0, j0, k0)
            u(i0, j0, k0) = 0.0d0
            v(i0, j0, k0) = 0.0d0
            w(i0, j0, k0) = 0.0d0
            T(i0, j0, k0) = 1.0d0
        end do
        end do
    end subroutine set_wall_uvwT

    subroutine add_wall_pressure(dir, n0, inward)
        integer, intent(IN) :: dir, n0, inward
        integer :: a, b, m, i0, j0, k0, ii, jj, kk
        real(dp) :: pp(4), pn, g(3, 3), h(3, 3, 3), gmu(3), muw

        do b = 1, merge(ny, nz, dir == 3)
        do a = 1, merge(nx, ny, dir /= 1)
            call wall_index(dir, n0, a, b, i0, j0, k0)
            do m = 1, 4
                call shift_index(dir, m*inward, i0, j0, k0, ii, jj, kk)
                pp(m) = d(ii, jj, kk)*T(ii, jj, kk)*Rgas
            end do
            call base_gradients(i0, j0, k0, g, h, gmu, muw)
            pn = normal_stress_div(dir, muw, gmu, g, h)/metric_diag(dir, i0, j0, k0)
            psum(i0, j0, k0) = psum(i0, j0, k0)+wall_pressure(pp, pn, inward)
            pcnt(i0, j0, k0) = pcnt(i0, j0, k0)+1
        end do
        end do
    end subroutine add_wall_pressure

end subroutine ocfd_bc_Lid_Driven_Cavity

subroutine ocfd_bc_Lid_Driven_Cavity_d
    ! Linearized counterpart of ocfd_bc_Lid_Driven_Cavity
    use mesh_data
    use flow_data
    use fluc_data
    use Flux_Viscous_d, only: compute_derivatives_d
    use LDC_Wall
    implicit none

    integer :: i, j, k
    real(dp) :: Rgas
    real(dp), allocatable :: psum(:, :, :)
    integer, allocatable :: pcnt(:, :, :)

    Rgas = 1.0d0/(Para%Ma**2*Para%gamma)

    if (npx == 0) call set_wall_fluc(1, 1)
    if (npx == npx0-1) call set_wall_fluc(1, nx)
    if (npy == 0) call set_wall_fluc(2, 1)
    if (npy == npy0-1) call set_wall_fluc(2, ny)
    if (npz == 0) call set_wall_fluc(3, 1)
    if (npz == npz0-1) call set_wall_fluc(3, nz)

    call mpi_exchange_xyz3d(u_fluc)
    call mpi_exchange_xyz3d(v_fluc)
    call mpi_exchange_xyz3d(w_fluc)
    call mpi_exchange_xyz3d(T_fluc)

    call compute_derivatives_d

    allocate (psum(nx, ny, nz), pcnt(nx, ny, nz))
    psum = 0.0d0
    pcnt = 0
    if (npx == 0) call add_wall_pressure_d(1, 1, 1)
    if (npx == npx0-1) call add_wall_pressure_d(1, nx, -1)
    if (npy == 0) call add_wall_pressure_d(2, 1, 1)
    if (npy == npy0-1) call add_wall_pressure_d(2, ny, -1)
    if (npz == 0) call add_wall_pressure_d(3, 1, 1)
    if (npz == npz0-1) call add_wall_pressure_d(3, nz, -1)

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
        if (pcnt(i, j, k) > 0) &
            d_fluc(i, j, k) = (psum(i, j, k)/pcnt(i, j, k)/Rgas-d(i, j, k)*T_fluc(i, j, k))/T(i, j, k)
    end do
    end do
    end do
    deallocate (psum, pcnt)

contains

    subroutine set_wall_fluc(dir, n0)
        integer, intent(IN) :: dir, n0
        integer :: a, b, i0, j0, k0

        do b = 1, merge(ny, nz, dir == 3)
        do a = 1, merge(nx, ny, dir /= 1)
            call wall_index(dir, n0, a, b, i0, j0, k0)
            u_fluc(i0, j0, k0) = 0.0d0
            v_fluc(i0, j0, k0) = 0.0d0
            w_fluc(i0, j0, k0) = 0.0d0
            T_fluc(i0, j0, k0) = 0.0d0
        end do
        end do
    end subroutine set_wall_fluc

    subroutine add_wall_pressure_d(dir, n0, inward)
        integer, intent(IN) :: dir, n0, inward
        integer :: a, b, m, i0, j0, k0, ii, jj, kk
        real(dp) :: pp(4), pn, muw, muwf
        real(dp) :: g(3, 3), h(3, 3, 3), gmu(3), gf(3, 3), hf(3, 3, 3), gmuf(3)

        do b = 1, merge(ny, nz, dir == 3)
        do a = 1, merge(nx, ny, dir /= 1)
            call wall_index(dir, n0, a, b, i0, j0, k0)
            do m = 1, 4
                call shift_index(dir, m*inward, i0, j0, k0, ii, jj, kk)
                pp(m) = (d_fluc(ii, jj, kk)*T(ii, jj, kk)+d(ii, jj, kk)*T_fluc(ii, jj, kk))*Rgas
            end do
            call base_gradients(i0, j0, k0, g, h, gmu, muw)
            call fluc_gradients(i0, j0, k0, gf, hf, gmuf, muwf)
            ! the stress divergence is bilinear in (mu, grad mu) and the velocity field
            pn = (normal_stress_div(dir, muwf, gmuf, g, h)+normal_stress_div(dir, muw, gmu, gf, hf)) &
                 /metric_diag(dir, i0, j0, k0)
            psum(i0, j0, k0) = psum(i0, j0, k0)+wall_pressure(pp, pn, inward)
            pcnt(i0, j0, k0) = pcnt(i0, j0, k0)+1
        end do
        end do
    end subroutine add_wall_pressure_d

end subroutine ocfd_bc_Lid_Driven_Cavity_d

subroutine ocfd_bc_Flat_Plate
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, k

    real(dp) :: x0, Tw, p1, p2, pb

    x0 = Para%BC_para(1)
    Tw = Para%BC_para(2)

    if (npx == 0) then
        do k = 1, nz
        do j = 1, ny
            d(1, j, k) = 1.0d0
            u(1, j, k) = 1.0d0
            v(1, j, k) = 0.0d0
            w(1, j, k) = 0.0d0
            T(1, j, k) = 1.0d0
        end do
        end do
    end if

    if (npy == npy0-1) then
        do k = 1, nz
        do i = 1, nx
            d(i, ny, k) = 1.0d0
            u(i, ny, k) = 1.0d0
            v(i, ny, k) = 0.0d0
            w(i, ny, k) = 0.0d0
            T(i, ny, k) = 1.0d0
        end do
        end do
    end if

    if (npy == 0) then
        do k = 1, nz
        do i = 1, nx
            if (Axx(i, 1, k) < x0) then
                d(i, 1, k) = (4.0d0*d(i, 2, k)-d(i, 3, k))/3.0d0
                u(i, 1, k) = (4.0d0*u(i, 2, k)-u(i, 3, k))/3.0d0
                v(i, 1, k) = 0.0d0
                w(i, 1, k) = 0.0d0
                T(i, 1, k) = (4.0d0*T(i, 2, k)-T(i, 3, k))/3.0d0
            else
                u(i, 1, k) = 0.0d0
                v(i, 1, k) = 0.0d0
                w(i, 1, k) = 0.0d0
                if (Tw <= 0) then
                    T(i, 1, k) = (4.0d0*T(i, 2, k)-T(i, 3, k))/3.0d0
                    p1 = d(i, 2, k)*T(i, 2, k)
                    p2 = d(i, 3, k)*T(i, 3, k)
                    pb = (4.0d0*p1-p2)/3.0d0
                    d(i, 1, k) = pb/T(i, 1, k)
                else
                    T(i, 1, k) = Tw
                    p1 = d(i, 2, k)*T(i, 2, k)
                    p2 = d(i, 3, k)*T(i, 3, k)
                    pb = (4.0d0*p1-p2)/3.0d0
                    d(i, 1, k) = pb/Tw
                end if
            end if
        end do
        end do
    end if

    if (npx == npx0-1) then
        do k = 1, nz
        do j = 1, ny
            d(nx, j, k) = 2.0d0*d(nx-1, j, k)-d(nx-2, j, k)
            u(nx, j, k) = 2.0d0*u(nx-1, j, k)-u(nx-2, j, k)
            v(nx, j, k) = 2.0d0*v(nx-1, j, k)-v(nx-2, j, k)
            w(nx, j, k) = 2.0d0*w(nx-1, j, k)-w(nx-2, j, k)
            T(nx, j, k) = 2.0d0*T(nx-1, j, k)-T(nx-2, j, k)
        end do
        end do
    end if

end subroutine ocfd_bc_Flat_Plate

subroutine ocfd_bc_Flat_Plate_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    integer :: i, j, k

    real(dp) :: x0, Tw, p1_fluc, p2_fluc, pb_fluc

    x0 = Para%BC_para(1)
    Tw = Para%BC_para(2)

    if (npx == 0) then
        do k = 1, nz
        do j = 1, ny
            d_fluc(1, j, k) = 0.0d0
            u_fluc(1, j, k) = 0.0d0
            v_fluc(1, j, k) = 0.0d0
            w_fluc(1, j, k) = 0.0d0
            T_fluc(1, j, k) = 0.0d0
        end do
        end do
    end if

    if (npy == npy0-1) then
        do k = 1, nz
        do i = 1, nx
            d_fluc(i, ny, k) = 0.0d0
            u_fluc(i, ny, k) = 0.0d0
            v_fluc(i, ny, k) = 0.0d0
            w_fluc(i, ny, k) = 0.0d0
            T_fluc(i, ny, k) = 0.0d0
        end do
        end do
    end if

    if (npy == 0) then
        do k = 1, nz
        do i = 1, nx
            if (Axx(i, 1, k) < x0) then
                d_fluc(i, 1, k) = (4.0d0*d_fluc(i, 2, k)-d_fluc(i, 3, k))/3.0d0
                u_fluc(i, 1, k) = (4.0d0*u_fluc(i, 2, k)-u_fluc(i, 3, k))/3.0d0
                v_fluc(i, 1, k) = 0.0d0
                w_fluc(i, 1, k) = 0.0d0
                T_fluc(i, 1, k) = (4.0d0*T_fluc(i, 2, k)-T_fluc(i, 3, k))/3.0d0
            else
                u_fluc(i, 1, k) = 0.0d0
                v_fluc(i, 1, k) = 0.0d0
                w_fluc(i, 1, k) = 0.0d0
                if (Tw <= 0) then
                    T_fluc(i, 1, k) = (4.0d0*T_fluc(i, 2, k)-T_fluc(i, 3, k))/3.0d0
                    p1_fluc = d_fluc(i, 2, k)*T(i, 2, k)+d(i, 2, k)*T_fluc(i, 2, k)
                    p2_fluc = d_fluc(i, 3, k)*T(i, 3, k)+d(i, 3, k)*T_fluc(i, 3, k)
                    pb_fluc = (4.0d0*p1_fluc-p2_fluc)/3.0d0
                    d_fluc(i, 1, k) = (pb_fluc-T_fluc(i, 1, k)*d(i, 1, k))/T(i, 1, k)
                else
                    T_fluc(i, 1, k) = 0.0d0
                    p1_fluc = d_fluc(i, 2, k)*T(i, 2, k)+d(i, 2, k)*T_fluc(i, 2, k)
                    p2_fluc = d_fluc(i, 3, k)*T(i, 3, k)+d(i, 3, k)*T_fluc(i, 3, k)
                    pb_fluc = (4.0d0*p1_fluc-p2_fluc)/3.0d0
                    d_fluc(i, 1, k) = (pb_fluc-T_fluc(i, 1, k)*d(i, 1, k))/Tw
                end if
            end if
        end do
        end do
    end if

    if (npx == npx0-1) then
        do k = 1, nz
        do j = 1, ny
            d_fluc(nx, j, k) = 2.0d0*d_fluc(nx-1, j, k)-d_fluc(nx-2, j, k)
            u_fluc(nx, j, k) = 2.0d0*u_fluc(nx-1, j, k)-u_fluc(nx-2, j, k)
            v_fluc(nx, j, k) = 2.0d0*v_fluc(nx-1, j, k)-v_fluc(nx-2, j, k)
            w_fluc(nx, j, k) = 2.0d0*w_fluc(nx-1, j, k)-w_fluc(nx-2, j, k)
            T_fluc(nx, j, k) = 2.0d0*T_fluc(nx-1, j, k)-T_fluc(nx-2, j, k)
        end do
        end do
    end if

end subroutine ocfd_bc_Flat_Plate_d

subroutine ocfd_bc_User_Def
    use mesh_data
    use flow_data
    implicit none

end subroutine ocfd_bc_User_Def

subroutine ocfd_bc_User_Def_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

end subroutine ocfd_bc_User_Def_d
