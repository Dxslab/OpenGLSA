! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid_character_pp
    use flow_data
    use mesh_data
    implicit none

    integer :: i, j, k, m, ierr
    real(dp) :: eps_rho_local, eps_rho_global
    integer :: Ka1, Kb1, Ka2, Kb2
    integer :: Bc_type(2)
    integer :: mark_rho_x(0:nx+1), mark_p_x(0:nx+1)
    integer :: mark_rho_y(0:ny+1), mark_p_y(0:ny+1)
    integer :: mark_rho_z(0:nz+1), mark_p_z(0:nz+1)
    logical :: track_limiter

    Ka1 = Scheme%Ka1
    Kb1 = Scheme%Kb1
    Ka2 = Scheme%Ka2
    Kb2 = Scheme%Kb2

    ! Initialize limiter statistics
    limiter_count_rho = 0
    limiter_count_pressure = 0
    limiter_count_total = 0

    ! The activation map is a snapshot of this sub-stage only, so it is reset here too
    track_limiter = allocated(pp_mark_rho_x)
    if (track_limiter) then
        pp_mark_rho_x = 0
        pp_mark_p_x = 0
        pp_mark_rho_y = 0
        pp_mark_p_y = 0
        pp_mark_rho_z = 0
        pp_mark_p_z = 0
    end if

    do k = 1-LAP, nz+LAP
    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j, k) = sqrt(T(i, j, k))/Para%Ma
        Hc(i, j, k) = c(i, j, k)**2/(Para%gamma-1.d0)+0.5d0*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)
    end do
    end do
    end do

    eps_rho_local = 1.0d-13
    eps_rho_global = 1.0d-13
    if (Para%pp_limiter /= 0) then
        do k = 1, nz
        do j = 1, ny
            call split_density_floor_pp(nx, LAP, Ka1, Kb2, &
                                        Akx1(1-LAP:nx+LAP, j, k), Aky1(1-LAP:nx+LAP, j, k), Akz1(1-LAP:nx+LAP, j, k), &
                                        c(1-LAP:nx+LAP, j, k), d(1-LAP:nx+LAP, j, k), &
                                        u(1-LAP:nx+LAP, j, k), v(1-LAP:nx+LAP, j, k), w(1-LAP:nx+LAP, j, k), &
                                        BcI(j, k, 1:2), eps_rho_local)
        end do
        end do
        do k = 1, nz
        do i = 1, nx
            call split_density_floor_pp(ny, LAP, Ka1, Kb2, &
                                        Aix1(i, 1-LAP:ny+LAP, k), Aiy1(i, 1-LAP:ny+LAP, k), Aiz1(i, 1-LAP:ny+LAP, k), &
                                        c(i, 1-LAP:ny+LAP, k), d(i, 1-LAP:ny+LAP, k), &
                                        u(i, 1-LAP:ny+LAP, k), v(i, 1-LAP:ny+LAP, k), w(i, 1-LAP:ny+LAP, k), &
                                        BcJ(i, k, 1:2), eps_rho_local)
        end do
        end do
        do j = 1, ny
        do i = 1, nx
            call split_density_floor_pp(nz, LAP, Ka1, Kb2, &
                                        Asx1(i, j, 1-LAP:nz+LAP), Asy1(i, j, 1-LAP:nz+LAP), Asz1(i, j, 1-LAP:nz+LAP), &
                                        c(i, j, 1-LAP:nz+LAP), d(i, j, 1-LAP:nz+LAP), &
                                        u(i, j, 1-LAP:nz+LAP), v(i, j, 1-LAP:nz+LAP), w(i, j, 1-LAP:nz+LAP), &
                                        BcK(i, j, 1:2), eps_rho_local)
        end do
        end do
        call MPI_ALLREDUCE(eps_rho_local, eps_rho_global, 1, OCFD_DATA_TYPE, MPI_MIN, MPI_COMM_WORLD, ierr)
        if (eps_rho_global <= 0.0d0) then
            if (my_id == 0) print *, 'PP limiter: nonpositive split base density.'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            stop 1
        end if
    end if

    !> x direction
    do k = 1, nz
    do j = 1, ny
        do i = 1-LAP, nx+LAP
            Hc1x(i) = Hc(i, j, k)
            c1x(i) = c(i, j, k)
            d1x(i) = d(i, j, k)
            u1x(i) = u(i, j, k)
            v1x(i) = v(i, j, k)
            w1x(i) = w(i, j, k)
            A1x(i) = Akx1(i, j, k)
            A2x(i) = Aky1(i, j, k)
            A3x(i) = Akz1(i, j, k)
        end do
        !
        Bc_type = BcI(j, k, 1:2)
        call flux_characteristic_pp(nx, LAP, Ka1, Kb1, Ka2, Kb2, &
                                    hhx(0:nx, 1:5), A1x(1-LAP:nx+LAP), A2x(1-LAP:nx+LAP), A3x(1-LAP:nx+LAP), &
                                    Hc1x(1-LAP:nx+LAP), c1x(1-LAP:nx+LAP), d1x(1-LAP:nx+LAP), &
                                    u1x(1-LAP:nx+LAP), v1x(1-LAP:nx+LAP), w1x(1-LAP:nx+LAP), &
                                    Bc_type, eps_rho_global, mark_rho_x, mark_p_x)
        !
        if (track_limiter) then
            pp_mark_rho_x(1:nx, j, k) = mark_rho_x(1:nx)
            pp_mark_p_x(1:nx, j, k) = mark_p_x(1:nx)
        end if
        !
        do m = 1, 5
            do i = 1, nx
                du(i, j, k, m) = hhx(i, m)-hhx(i-1, m)
            end do
        end do
    end do
    end do
    !> y direction
    do k = 1, nz
    do i = 1, nx
        do j = 1-LAP, ny+LAP
            Hc1y(j) = Hc(i, j, k)
            c1y(j) = c(i, j, k)
            d1y(j) = d(i, j, k)
            u1y(j) = u(i, j, k)
            v1y(j) = v(i, j, k)
            w1y(j) = w(i, j, k)
            A1y(j) = Aix1(i, j, k)
            A2y(j) = Aiy1(i, j, k)
            A3y(j) = Aiz1(i, j, k)
        end do
        !
        Bc_type = BcJ(i, k, 1:2)
        call flux_characteristic_pp(ny, LAP, Ka1, Kb1, Ka2, Kb2, &
                                    hhy(0:ny, 1:5), A1y(1-LAP:ny+LAP), A2y(1-LAP:ny+LAP), A3y(1-LAP:ny+LAP), &
                                    Hc1y(1-LAP:ny+LAP), c1y(1-LAP:ny+LAP), d1y(1-LAP:ny+LAP), &
                                    u1y(1-LAP:ny+LAP), v1y(1-LAP:ny+LAP), w1y(1-LAP:ny+LAP), &
                                    Bc_type, eps_rho_global, mark_rho_y, mark_p_y)
        !
        if (track_limiter) then
            pp_mark_rho_y(i, 1:ny, k) = mark_rho_y(1:ny)
            pp_mark_p_y(i, 1:ny, k) = mark_p_y(1:ny)
        end if
        !
        do m = 1, 5
            do j = 1, ny
                du(i, j, k, m) = du(i, j, k, m)+hhy(j, m)-hhy(j-1, m)
            end do
        end do
    end do
    end do
    !> z direction
    do j = 1, ny
    do i = 1, nx
        do k = 1-LAP, nz+LAP
            Hc1z(k) = Hc(i, j, k)
            c1z(k) = c(i, j, k)
            d1z(k) = d(i, j, k)
            u1z(k) = u(i, j, k)
            v1z(k) = v(i, j, k)
            w1z(k) = w(i, j, k)
            A1z(k) = Asx1(i, j, k)
            A2z(k) = Asy1(i, j, k)
            A3z(k) = Asz1(i, j, k)
        end do
        !
        Bc_type = BcK(i, j, 1:2)
        call flux_characteristic_pp(nz, LAP, Ka1, Kb1, Ka2, Kb2, &
                                    hhz(0:nz, 1:5), A1z(1-LAP:nz+LAP), A2z(1-LAP:nz+LAP), A3z(1-LAP:nz+LAP), &
                                    Hc1z(1-LAP:nz+LAP), c1z(1-LAP:nz+LAP), d1z(1-LAP:nz+LAP), &
                                    u1z(1-LAP:nz+LAP), v1z(1-LAP:nz+LAP), w1z(1-LAP:nz+LAP), &
                                    Bc_type, eps_rho_global, mark_rho_z, mark_p_z)
        !
        if (track_limiter) then
            pp_mark_rho_z(i, j, 1:nz) = mark_rho_z(1:nz)
            pp_mark_p_z(i, j, 1:nz) = mark_p_z(1:nz)
        end if
        !
        do m = 1, 5
            do k = 1, nz
                du(i, j, k, m) = -(du(i, j, k, m)+hhz(k, m)-hhz(k-1, m))*Ajac(i, j, k)
            end do
        end do
    end do
    end do

end subroutine flux_inviscid_character_pp

subroutine flux_characteristic_pp(nn, LAP1, Ka1, Kb1, Ka2, Kb2, hh, AA1, AA2, AA3, Hc, cc, dd, uu, vv, ww, &
                                  Bc_type, eps_rho, mark_rho, mark_p)
    use flow_para
    implicit none

    integer :: Bc_type(2), ib, ie
    integer :: nn, LAP1, Ka1, Kb1, Ka2, Kb2, i, m, mm, ms
    !> Limiter hits, attributed to the cell whose base state was limited
    integer, intent(out) :: mark_rho(0:nn+1), mark_p(0:nn+1)
    logical :: rho_limited, pressure_limited, base_invalid
    real(dp), dimension(1-LAP1:nn+LAP1):: Hc, cc, dd, uu, vv, ww, AA1, AA2, AA3
    real(dp), intent(in) :: eps_rho
    real(dp):: hh(0:nn, 5)
    real(dp):: fpc(Ka1:Kb1, 5)
    real(dp):: fmc(Ka2:Kb2, 5)
    real(dp):: h1(5), h2(5), hp(5), hm(5), qbase(5)
    real(dp):: S(5, 5), S1(5, 5)
    real(dp):: d1, u1, v1, w1, c1, V2, un, ul, um, dtmp, Hc1
    real(dp):: a1, a2, a3, ss, n1, n2, n3, l1, l2, l3, m1, m2, m3, KK1, KK, X1, H
    real(dp):: rho, rhou, rhov, rhow, rhoe, u2, pp, lambda(5), lambda_max
    real(dp):: sc, eps2, wc, qc
    real(dp):: Ustate(5), Fstate(5)
    real(dp), parameter :: delta_ef = 0.01d0

    mark_rho = 0
    mark_p = 0

    if (Bc_type(1) == 1) then
        ib = 1
    else
        ib = 0
    end if

    if (Bc_type(2) == 1) then
        ie = nn-1
    else
        ie = nn
    end if

    do i = ib, ie
        ! Roe Average
        d1 = sqrt(dd(i)*dd(i+1))
        dtmp = 1.0d0/(sqrt(dd(i))+sqrt(dd(i+1)))
        u1 = (sqrt(dd(i))*uu(i)+sqrt(dd(i+1))*uu(i+1))*dtmp
        v1 = (sqrt(dd(i))*vv(i)+sqrt(dd(i+1))*vv(i+1))*dtmp
        w1 = (sqrt(dd(i))*ww(i)+sqrt(dd(i+1))*ww(i+1))*dtmp
        Hc1 = (sqrt(dd(i))*Hc(i)+sqrt(dd(i+1))*Hc(i+1))*dtmp
        c1 = sqrt((Para%gamma-1.0d0)*(Hc1-0.5d0*(u1*u1+v1*v1+w1*w1)))

        a1 = (AA1(i)+AA1(i+1))*0.5d0
        a2 = (AA2(i)+AA2(i+1))*0.5d0
        a3 = (AA3(i)+AA3(i+1))*0.5d0

        ss = sqrt(a1*a1+a2*a2+a3*a3)
        n1 = a1/ss; n2 = a2/ss; n3 = a3/ss

        if (abs(n3) <= abs(n2)) then
            ss = sqrt(n1*n1+n2*n2)
            l1 = -n2/ss; l2 = n1/ss; l3 = 0.d0
        else
            ss = sqrt(n1*n1+n3*n3)
            l1 = -n3/ss; l2 = 0.d0; l3 = n1/ss
        end if

        m1 = n2*l3-n3*l2
        m2 = n3*l1-n1*l3
        m3 = n1*l2-n2*l1
        un = u1*n1+v1*n2+w1*n3
        ul = u1*l1+v1*l2+w1*l3
        um = u1*m1+v1*m2+w1*m3
        V2 = (u1*u1+v1*v1+w1*w1)*0.5d0
        KK = (para%gamma-1.d0)/(c1*c1)
        KK1 = KK*0.5d0
        X1 = 1.d0/(2.d0*c1)
        H = V2+1.d0/KK

        S(1, 1) = 1.d0-KK*V2; S(1, 2) = KK*u1; S(1, 3) = KK*v1; S(1, 4) = KK*w1; S(1, 5) = -KK
        S(2, 1) = -ul; S(2, 2) = l1; S(2, 3) = l2; S(2, 4) = l3; S(2, 5) = 0.d0
        S(3, 1) = -um; S(3, 2) = m1; S(3, 3) = m2; S(3, 4) = m3; S(3, 5) = 0.d0
        S(4, 1) = KK1*V2+X1*un; S(4, 2) = -X1*n1-KK1*u1; S(4, 3) = -X1*n2-KK1*v1; S(4, 4) = -X1*n3-KK1*w1; S(4, 5) = KK1
        S(5, 1) = KK1*V2-X1*un; S(5, 2) = X1*n1-KK1*u1; S(5, 3) = X1*n2-KK1*v1; S(5, 4) = X1*n3-KK1*w1; S(5, 5) = KK1

        S1(1, 1) = 1.d0; S1(1, 2) = 0.d0; S1(1, 3) = 0.d0; S1(1, 4) = 1.d0; S1(1, 5) = 1.d0
        S1(2, 1) = u1; S1(2, 2) = l1; S1(2, 3) = m1; S1(2, 4) = u1-c1*n1; S1(2, 5) = u1+c1*n1
        S1(3, 1) = v1; S1(3, 2) = l2; S1(3, 3) = m2; S1(3, 4) = v1-c1*n2; S1(3, 5) = v1+c1*n2
        S1(4, 1) = w1; S1(4, 2) = l3; S1(4, 3) = m3; S1(4, 4) = w1-c1*n3; S1(4, 5) = w1+c1*n3
        S1(5, 1) = V2; S1(5, 2) = ul; S1(5, 3) = um; S1(5, 4) = H-c1*un; S1(5, 5) = H+c1*un

        ! Characteristic fields: 1-3 convect with u_n, 4 is u_n-c, 5 is u_n+c
        lambda = 0.0d0
        lambda_max = 1.0d-14
        do mm = Ka1, Kb2
            m = i+mm
            ms = m
            if (m < 1 .and. Bc_type(1) == 1) ms = 1
            if (m > nn .and. Bc_type(2) == 1) ms = nn

            u2 = uu(ms)*AA1(ms)+vv(ms)*AA2(ms)+ww(ms)*AA3(ms)
            sc = cc(ms)*sqrt(AA1(ms)*AA1(ms)+AA2(ms)*AA2(ms)+AA3(ms)*AA3(ms))
            if (Para%pp_limiter /= 0) then
                lambda_max = max(lambda_max, abs(u2)+sc)
            else
                eps2 = (max(1.0d-12, delta_ef*(abs(u2)+sc)))**2
                lambda(1) = max(lambda(1), sqrt(u2*u2+eps2))
                lambda(4) = max(lambda(4), sqrt((u2-sc)**2+eps2))
                lambda(5) = max(lambda(5), sqrt((u2+sc)**2+eps2))
            end if
        end do
        lambda(2) = lambda(1)
        lambda(3) = lambda(1)

        do mm = Ka1, Kb2
            m = i+mm
            ms = m
            if (m < 1 .and. Bc_type(1) == 1) ms = 1
            if (m > nn .and. Bc_type(2) == 1) ms = nn

            rho = dd(ms)
            rhou = dd(ms)*uu(ms)
            rhov = dd(ms)*vv(ms)
            rhow = dd(ms)*ww(ms)
            rhoe = dd(ms)*cc(ms)*cc(ms)/(Para%gamma*(Para%gamma-1.0d0))+ &
                   0.5d0*dd(ms)*(uu(ms)*uu(ms)+vv(ms)*vv(ms)+ww(ms)*ww(ms))
            u2 = uu(ms)*AA1(ms)+vv(ms)*AA2(ms)+ww(ms)*AA3(ms)
            pp = dd(ms)*cc(ms)*cc(ms)/Para%gamma

            Ustate = (/rho, rhou, rhov, rhow, rhoe/)
            Fstate = (/rho*u2, rhou*u2+AA1(ms)*pp, rhov*u2+AA2(ms)*pp, rhow*u2+AA3(ms)*pp, rho*u2*Hc(ms)/)

            if (Para%pp_limiter /= 0) then
                hp = 0.5d0*(Ustate+Fstate/lambda_max)
                hm = 0.5d0*(Ustate-Fstate/lambda_max)

                if (mm >= Ka1 .and. mm <= Kb1) then
                    do m = 1, 5
                        fpc(mm, m) = S(m, 1)*hp(1)+S(m, 2)*hp(2)+S(m, 3)*hp(3)+S(m, 4)*hp(4)+S(m, 5)*hp(5)
                    end do
                end if

                if (mm >= Ka2 .and. mm <= Kb2) then
                    do m = 1, 5
                        fmc(mm, m) = S(m, 1)*hm(1)+S(m, 2)*hm(2)+S(m, 3)*hm(3)+S(m, 4)*hm(4)+S(m, 5)*hm(5)
                    end do
                end if
            else
                do m = 1, 5
                    wc = S(m, 1)*Fstate(1)+S(m, 2)*Fstate(2)+S(m, 3)*Fstate(3)+S(m, 4)*Fstate(4)+S(m, 5)*Fstate(5)
                    qc = S(m, 1)*Ustate(1)+S(m, 2)*Ustate(2)+S(m, 3)*Ustate(3)+S(m, 4)*Ustate(4)+S(m, 5)*Ustate(5)
                    if (mm >= Ka1 .and. mm <= Kb1) fpc(mm, m) = 0.5d0*(wc+lambda(m)*qc)
                    if (mm >= Ka2 .and. mm <= Kb2) fmc(mm, m) = 0.5d0*(wc-lambda(m)*qc)
                end do
            end if
        end do

        do m = 1, 5
            call flux_x1(Ka1, Kb1, fpc(:, m), h1(m), Scheme%Scheme_Invis, Bc_type, i, nn)
            call flux_x2(Ka2, Kb2, fmc(:, m), h2(m), Scheme%Scheme_Invis, Bc_type, i, nn)
        end do

        do m = 1, 5
            hp(m) = S1(m, 1)*h1(1)+S1(m, 2)*h1(2)+S1(m, 3)*h1(3)+S1(m, 4)*h1(4)+S1(m, 5)*h1(5)
            hm(m) = S1(m, 1)*h2(1)+S1(m, 2)*h2(2)+S1(m, 3)*h2(3)+S1(m, 4)*h2(4)+S1(m, 5)*h2(5)
        end do

        if (Para%pp_limiter /= 0) then
            call split_state_pp(i, LAP1, 1, lambda_max, AA1, AA2, AA3, Hc, cc, dd, uu, vv, ww, qbase)
            call limit_split_endpoint_pp(hp, qbase, 1.0d0/12.0d0, eps_rho, rho_limited, pressure_limited, base_invalid)
            if (rho_limited) mark_rho(i) = mark_rho(i)+1
            if (pressure_limited .or. base_invalid) mark_p(i) = mark_p(i)+1

            call split_state_pp(i+1, LAP1, -1, lambda_max, AA1, AA2, AA3, Hc, cc, dd, uu, vv, ww, qbase)
            call limit_split_endpoint_pp(hm, qbase, 1.0d0/12.0d0, eps_rho, rho_limited, pressure_limited, base_invalid)
            if (rho_limited) mark_rho(i+1) = mark_rho(i+1)+1
            if (pressure_limited .or. base_invalid) mark_p(i+1) = mark_p(i+1)+1

            do m = 1, 5
                hh(i, m) = lambda_max*(hp(m)-hm(m))
            end do
        else
            do m = 1, 5
                hh(i, m) = hp(m)+hm(m)
            end do
        end if
    end do

    if (Bc_type(1) == 1) then
        do m = 1, 5
            hh(0, m) = hh(1, m)
        end do
    end if

    if (Bc_type(2) == 1) then
        do m = 1, 5
            hh(nn, m) = hh(nn-1, m)
        end do
    end if

end subroutine flux_characteristic_pp

subroutine split_density_floor_pp(nn, LAP1, Ka1, Kb2, A1, A2, A3, cc, dd, uu, vv, ww, Bc_type, floor_rho)
    use flow_para
    implicit none
    integer, intent(in) :: nn, LAP1, Ka1, Kb2, Bc_type(2)
    real(dp), intent(in) :: A1(1-LAP1:nn+LAP1), A2(1-LAP1:nn+LAP1), A3(1-LAP1:nn+LAP1)
    real(dp), intent(in) :: cc(1-LAP1:nn+LAP1), dd(1-LAP1:nn+LAP1)
    real(dp), intent(in) :: uu(1-LAP1:nn+LAP1), vv(1-LAP1:nn+LAP1), ww(1-LAP1:nn+LAP1)
    real(dp), intent(inout) :: floor_rho
    real(dp) :: speed, u_normal, sound_normal, rho_plus, rho_minus
    integer :: i, k, ms, ib, ie

    ib = 0
    ie = nn
    if (Bc_type(1) == 1) ib = 1
    if (Bc_type(2) == 1) ie = nn-1
    do i = ib, ie
        speed = 1.0d-14
        do k = Ka1, Kb2
            ms = i+k
            if (ms < 1 .and. Bc_type(1) == 1) ms = 1
            if (ms > nn .and. Bc_type(2) == 1) ms = nn
            u_normal = uu(ms)*A1(ms)+vv(ms)*A2(ms)+ww(ms)*A3(ms)
            sound_normal = cc(ms)*sqrt(A1(ms)**2+A2(ms)**2+A3(ms)**2)
            speed = max(speed, abs(u_normal)+sound_normal)
        end do
        u_normal = uu(i)*A1(i)+vv(i)*A2(i)+ww(i)*A3(i)
        rho_plus = 0.5d0*dd(i)*(1.0d0+u_normal/speed)
        u_normal = uu(i+1)*A1(i+1)+vv(i+1)*A2(i+1)+ww(i+1)*A3(i+1)
        rho_minus = 0.5d0*dd(i+1)*(1.0d0-u_normal/speed)
        floor_rho = min(floor_rho, rho_plus, rho_minus)
    end do
end subroutine split_density_floor_pp

subroutine split_state_pp(cell_index, LAP1, split_sign, lambda, A1, A2, A3, Hc, cc, dd, uu, vv, ww, state)
    use flow_para
    implicit none
    integer, intent(in) :: cell_index, LAP1, split_sign
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: A1(*), A2(*), A3(*), Hc(*), cc(*), dd(*), uu(*), vv(*), ww(*)
    real(dp), intent(out) :: state(5)
    real(dp) :: rho, u2, pp, energy
    integer :: ic

    ic = cell_index+LAP1
    rho = dd(ic)
    u2 = uu(ic)*A1(ic)+vv(ic)*A2(ic)+ww(ic)*A3(ic)
    pp = dd(ic)*cc(ic)*cc(ic)/Para%gamma
    energy = dd(ic)*Hc(ic)-pp    ! rho*E = rho*H - p
    state(1) = 0.5d0*(rho+split_sign*rho*u2/lambda)
    state(2) = 0.5d0*(rho*uu(ic)+split_sign*(rho*uu(ic)*u2+A1(ic)*pp)/lambda)
    state(3) = 0.5d0*(rho*vv(ic)+split_sign*(rho*vv(ic)*u2+A2(ic)*pp)/lambda)
    state(4) = 0.5d0*(rho*ww(ic)+split_sign*(rho*ww(ic)*u2+A3(ic)*pp)/lambda)
    state(5) = 0.5d0*(energy+split_sign*(rho*u2*Hc(ic))/lambda)
end subroutine split_state_pp

subroutine limit_split_endpoint_pp(endpoint, base, endpoint_weight, eps_rho, rho_limited, pressure_limited, base_invalid)
    use flow_para
    use flow_data, only: limiter_count_rho, limiter_count_pressure, limiter_count_total
    implicit none
    real(dp), intent(inout) :: endpoint(5)
    real(dp), intent(in) :: base(5), endpoint_weight, eps_rho
    logical, intent(out) :: rho_limited, pressure_limited, base_invalid
    real(dp) :: qstar(5), modified(5), rho_min, theta_rho, theta_p, tm
    integer :: ierr
    real(dp) :: pressure_split_pp
    integer :: branch_index

    rho_limited = .false.
    pressure_limited = .false.
    base_invalid = .false.

    if (pressure_split_pp(base) <= 0.0d0) then
        base_invalid = .true.
        print *, 'PP limiter: inadmissible split base on rank ', my_id, base
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        stop 1
    end if
    qstar = (base-endpoint_weight*endpoint)/(1.0d0-endpoint_weight)
    rho_min = min(endpoint(1), qstar(1))
    if (rho_min < eps_rho) then
        if (abs(base(1)-rho_min) > 1.0d-30) then
            theta_rho = min((base(1)-eps_rho)/(base(1)-rho_min), 1.0d0)
        else
            theta_rho = 0.0d0
        end if
        endpoint(1) = base(1)+theta_rho*(endpoint(1)-base(1))
        qstar = (base-endpoint_weight*endpoint)/(1.0d0-endpoint_weight)
        rho_limited = .true.
    end if

    theta_p = 1.0d0
    do branch_index = 1, 2
        if (branch_index == 1) then
            modified = endpoint
        else
            modified = qstar
        end if
        if (pressure_split_pp(modified) < 0.0d0) then
            call pressure_root_pp(base, modified, tm)
            theta_p = min(theta_p, tm)
            pressure_limited = .true.
        end if
    end do
    endpoint = base+theta_p*(endpoint-base)

    ! Update statistics
    if (rho_limited) limiter_count_rho = limiter_count_rho+1
    if (pressure_limited) limiter_count_pressure = limiter_count_pressure+1
    if (rho_limited .or. pressure_limited) limiter_count_total = limiter_count_total+1
end subroutine limit_split_endpoint_pp

subroutine pressure_root_pp(base, target, theta)
    use flow_para
    implicit none
    real(dp), intent(in) :: base(5), target(5)
    real(dp), intent(out) :: theta
    real(dp) :: left, right, mid, p_mid
    real(dp) :: dq(5), a, b, c, sqrt_disc
    integer :: iter
    real(dp) :: pressure_split_pp

    dq = target-base
    a = dq(1)*dq(5)-0.5d0*(dq(2)*dq(2)+dq(3)*dq(3)+dq(4)*dq(4))
    b = base(1)*dq(5)+base(5)*dq(1)-(base(2)*dq(2)+base(3)*dq(3)+base(4)*dq(4))
    c = base(1)*base(5)-0.5d0*(base(2)*base(2)+base(3)*base(3)+base(4)*base(4))
    sqrt_disc = sqrt(max(b*b-4.0d0*a*c, 0.0d0))
    theta = -1.0d0
    if (b <= 0.0d0) then
        if (sqrt_disc-b > 0.0d0) theta = 2.0d0*c/(sqrt_disc-b)
    else if (a < 0.0d0) then
        theta = (b+sqrt_disc)/(-2.0d0*a)
    end if
    if (theta >= 0.0d0) then
        theta = min(theta, 1.0d0)
        ! Stay just inside the admissible side when rebuilding the auxiliary state.
        theta = theta*(1.0d0-16.0d0*epsilon(1.0d0))
        if (pressure_split_pp(base+theta*dq) > 0.0d0) return
    end if

    left = 0.0d0
    right = 1.0d0
    do iter = 1, 60
        mid = 0.5d0*(left+right)
        p_mid = pressure_split_pp((1.0d0-mid)*base+mid*target)
        if (p_mid > 0.0d0) then
            left = mid
        else
            right = mid
        end if
    end do
    theta = left*(1.0d0-16.0d0*epsilon(1.0d0))
end subroutine pressure_root_pp

real(dp) function pressure_split_pp(state)
    use flow_para
    implicit none
    real(dp), intent(in) :: state(5)
    if (state(1) <= 0.0d0) then
        pressure_split_pp = -huge(1.0d0)
    else
        pressure_split_pp = (Para%gamma-1.0d0)*(state(5)-0.5d0*(state(2)*state(2)+state(3)*state(3)+state(4)*state(4))/state(1))
    end if
end function pressure_split_pp
