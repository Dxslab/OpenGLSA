! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid2d_character_pp
    use flow_data
    use mesh_data
    implicit none

    integer :: i, j, m, ierr
    real(dp) :: eps_rho_local, eps_rho_global
    integer :: Ka1, Kb1, Ka2, Kb2
    integer :: Bc_type(2)
    integer :: mark_rho_x(0:nx+1), mark_p_x(0:nx+1)
    integer :: mark_rho_y(0:ny+1), mark_p_y(0:ny+1)
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
    end if

    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j) = sqrt(T(i, j))/Para%Ma
        Hc(i, j) = c(i, j)**2/(Para%gamma-1.d0)+0.5d0*(u(i, j)**2+v(i, j)**2)
    end do
    end do

    eps_rho_local = 1.0d-13
    eps_rho_global = 1.0d-13
    if (Para%pp_limiter /= 0) then
        do j = 1, ny
            call split_density_floor_pp(nx, LAP, Ka1, Kb2, &
                                        Akx1(1-LAP:nx+LAP, j), Aky1(1-LAP:nx+LAP, j), &
                                        c(1-LAP:nx+LAP, j), d(1-LAP:nx+LAP, j), &
                                        u(1-LAP:nx+LAP, j), v(1-LAP:nx+LAP, j), BcI(j, 1:2), eps_rho_local)
        end do
        do i = 1, nx
            call split_density_floor_pp(ny, LAP, Ka1, Kb2, &
                                        Aix1(i, 1-LAP:ny+LAP), Aiy1(i, 1-LAP:ny+LAP), &
                                        c(i, 1-LAP:ny+LAP), d(i, 1-LAP:ny+LAP), &
                                        u(i, 1-LAP:ny+LAP), v(i, 1-LAP:ny+LAP), BcJ(i, 1:2), eps_rho_local)
        end do
        call MPI_ALLREDUCE(eps_rho_local, eps_rho_global, 1, OCFD_DATA_TYPE, MPI_MIN, MPI_COMM_WORLD, ierr)
        if (eps_rho_global <= 0.0d0) then
            if (my_id == 0) print *, 'PP limiter: nonpositive split base density.'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            stop 1
        end if
    end if

    !> x direction
    do j = 1, ny
        do i = 1-LAP, nx+LAP
            Hc1x(i) = Hc(i, j)
            c1x(i) = c(i, j)
            d1x(i) = d(i, j)
            u1x(i) = u(i, j)
            v1x(i) = v(i, j)
            A1x(i) = Akx1(i, j)
            A2x(i) = Aky1(i, j)
        end do
        !
        Bc_type = BcI(j, 1:2)
        call flux_characteristic_pp(nx, LAP, Ka1, Kb1, Ka2, Kb2, &
                                    hhx(0:nx, 1:4), A1x(1-LAP:nx+LAP), A2x(1-LAP:nx+LAP), &
                                    Hc1x(1-LAP:nx+LAP), c1x(1-LAP:nx+LAP), &
                                    d1x(1-LAP:nx+LAP), u1x(1-LAP:nx+LAP), v1x(1-LAP:nx+LAP), &
                                    Bc_type, eps_rho_global, mark_rho_x, mark_p_x)
        !
        if (track_limiter) then
            pp_mark_rho_x(1:nx, j) = mark_rho_x(1:nx)
            pp_mark_p_x(1:nx, j) = mark_p_x(1:nx)
        end if
        !
        do m = 1, 4
            do i = 1, nx
                du(i, j, m) = hhx(i, m)-hhx(i-1, m)
            end do
        end do
    end do
    !> y direction
    do i = 1, nx
        do j = 1-LAP, ny+LAP
            Hc1y(j) = Hc(i, j)
            c1y(j) = c(i, j)
            d1y(j) = d(i, j)
            u1y(j) = u(i, j)
            v1y(j) = v(i, j)
            A1y(j) = Aix1(i, j)
            A2y(j) = Aiy1(i, j)
        end do
        !
        Bc_type = BcJ(i, 1:2)
        call flux_characteristic_pp(ny, LAP, Ka1, Kb1, Ka2, Kb2, &
                                    hhy(0:ny, 1:4), A1y(1-LAP:ny+LAP), A2y(1-LAP:ny+LAP), &
                                    Hc1y(1-LAP:ny+LAP), c1y(1-LAP:ny+LAP), &
                                    d1y(1-LAP:ny+LAP), u1y(1-LAP:ny+LAP), v1y(1-LAP:ny+LAP), &
                                    Bc_type, eps_rho_global, mark_rho_y, mark_p_y)
        !
        if (track_limiter) then
            pp_mark_rho_y(i, 1:ny) = mark_rho_y(1:ny)
            pp_mark_p_y(i, 1:ny) = mark_p_y(1:ny)
        end if
        !
        do m = 1, 4
            do j = 1, ny
                du(i, j, m) = -(du(i, j, m)+hhy(j, m)-hhy(j-1, m))*Ajac(i, j)
            end do
        end do
    end do

    if (Para%IF_Cylindrical_Coordinate == 1) then
        do j = 1, ny
        do i = 1, nx
            if (abs(Ayy(i, j)) > 1.0d-10) then
                Source_Invis(i, j, 1) = d(i, j)*v(i, j)/Ayy(i, j)
                Source_Invis(i, j, 2) = d(i, j)*v(i, j)*u(i, j)/Ayy(i, j)
                Source_Invis(i, j, 3) = d(i, j)*v(i, j)*v(i, j)/Ayy(i, j)
                Source_Invis(i, j, 4) = d(i, j)*v(i, j)*Hc(i, j)/Ayy(i, j)
            else
                Source_Invis(i, j, 1) = 0.0d0
                Source_Invis(i, j, 2) = 0.0d0
                Source_Invis(i, j, 3) = 0.0d0
                Source_Invis(i, j, 4) = 0.0d0
            end if
        end do
        end do

        if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
            do j = 1, ny
                Source_Invis(1, j, 1) = (56.0d0*Source_Invis(2, j, 1)-28.0d0*Source_Invis(3, j, 1)+ &
                                         8.0d0*Source_Invis(4, j, 1)-Source_Invis(5, j, 1))/35.0d0
                Source_Invis(1, j, 2) = (56.0d0*Source_Invis(2, j, 2)-28.0d0*Source_Invis(3, j, 2)+ &
                                         8.0d0*Source_Invis(4, j, 2)-Source_Invis(5, j, 2))/35.0d0
                Source_Invis(1, j, 3) = 0.0d0
                Source_Invis(1, j, 4) = (56.0d0*Source_Invis(2, j, 4)-28.0d0*Source_Invis(3, j, 4)+ &
                                         8.0d0*Source_Invis(4, j, 4)-Source_Invis(5, j, 4))/35.0d0
            end do
        end if

        if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
            do j = 1, ny
                Source_Invis(nx, j, 1) = (56.0d0*Source_Invis(nx-1, j, 1)-28.0d0*Source_Invis(nx-2, j, 1)+ &
                                          8.0d0*Source_Invis(nx-3, j, 1)-Source_Invis(nx-4, j, 1))/35.0d0
                Source_Invis(nx, j, 2) = (56.0d0*Source_Invis(nx-1, j, 2)-28.0d0*Source_Invis(nx-2, j, 2)+ &
                                          8.0d0*Source_Invis(nx-3, j, 2)-Source_Invis(nx-4, j, 2))/35.0d0
                Source_Invis(nx, j, 3) = 0.0d0
                Source_Invis(nx, j, 4) = (56.0d0*Source_Invis(nx-1, j, 4)-28.0d0*Source_Invis(nx-2, j, 4)+ &
                                          8.0d0*Source_Invis(nx-3, j, 4)-Source_Invis(nx-4, j, 4))/35.0d0
            end do
        end if

        if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
            do i = 1, nx
                Source_Invis(i, 1, 1) = (56.0d0*Source_Invis(i, 2, 1)-28.0d0*Source_Invis(i, 3, 1)+ &
                                         8.0d0*Source_Invis(i, 4, 1)-Source_Invis(i, 5, 1))/35.0d0
                Source_Invis(i, 1, 2) = (56.0d0*Source_Invis(i, 2, 2)-28.0d0*Source_Invis(i, 3, 2)+ &
                                         8.0d0*Source_Invis(i, 4, 2)-Source_Invis(i, 5, 2))/35.0d0
                Source_Invis(i, 1, 3) = 0.0d0
                Source_Invis(i, 1, 4) = (56.0d0*Source_Invis(i, 2, 4)-28.0d0*Source_Invis(i, 3, 4)+ &
                                         8.0d0*Source_Invis(i, 4, 4)-Source_Invis(i, 5, 4))/35.0d0
            end do
        end if

        if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
            do i = 1, nx
                Source_Invis(i, ny, 1) = (56.0d0*Source_Invis(i, ny-1, 1)-28.0d0*Source_Invis(i, ny-2, 1)+ &
                                          8.0d0*Source_Invis(i, ny-3, 1)-Source_Invis(i, ny-4, 1))/35.0d0
                Source_Invis(i, ny, 2) = (56.0d0*Source_Invis(i, ny-1, 2)-28.0d0*Source_Invis(i, ny-2, 2)+ &
                                          8.0d0*Source_Invis(i, ny-3, 2)-Source_Invis(i, ny-4, 2))/35.0d0
                Source_Invis(i, ny, 3) = 0.0d0
                Source_Invis(i, ny, 4) = (56.0d0*Source_Invis(i, ny-1, 4)-28.0d0*Source_Invis(i, ny-2, 4)+ &
                                          8.0d0*Source_Invis(i, ny-3, 4)-Source_Invis(i, ny-4, 4))/35.0d0
            end do
        end if

        do j = 1, ny
        do i = 1, nx
            du(i, j, 1) = du(i, j, 1)-Source_Invis(i, j, 1)
            du(i, j, 2) = du(i, j, 2)-Source_Invis(i, j, 2)
            du(i, j, 3) = du(i, j, 3)-Source_Invis(i, j, 3)
            du(i, j, 4) = du(i, j, 4)-Source_Invis(i, j, 4)
        end do
        end do

    end if

end subroutine flux_inviscid2d_character_pp

subroutine flux_characteristic_pp(nn, LAP1, Ka1, Kb1, Ka2, Kb2, hh, A1, A2, Hc, cc, dd, uu, vv, Bc_type, eps_rho, mark_rho, mark_p)
    use flow_para
    implicit none

    integer :: Bc_type(2), ib, ie
    integer :: nn, LAP1, Ka1, Kb1, Ka2, Kb2, i, m, m1, ms
    !> Limiter hits, attributed to the cell whose base state was limited
    integer, intent(out) :: mark_rho(0:nn+1), mark_p(0:nn+1)
    logical :: rho_limited, pressure_limited, base_invalid
    real(dp), dimension(1-LAP1:nn+LAP1):: Hc, cc, dd, uu, vv, A1, A2
    real(dp), intent(in) :: eps_rho
    real(dp):: hh(0:nn, 4)
    real(dp):: fpc(Ka1:Kb1, 4)
    real(dp):: fmc(Ka2:Kb2, 4)
    real(dp):: h1(4), h2(4), hp(4), hm(4), qbase(4)
    real(dp):: S(4, 4), S1(4, 4)
    real(dp):: dtmp
    real(dp):: d1, u1, v1, c1, Hc1, alfa, beta, ss
    real(dp):: rho, rhou, rhov, rhoe, u2, pp, lambda(4), lambda_max
    real(dp):: sc, eps2, wc, qc
    real(dp):: Ustate(4), Fstate(4)
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
        Hc1 = (sqrt(dd(i))*Hc(i)+sqrt(dd(i+1))*Hc(i+1))*dtmp
        c1 = sqrt((Para%gamma-1.0d0)*(Hc1-0.5d0*(u1*u1+v1*v1)))

        alfa = (A1(i)+A1(i+1))*0.5d0
        beta = (A2(i)+A2(i+1))*0.5d0

        ss = sqrt(alfa*alfa+beta*beta)
        alfa = alfa/ss
        beta = beta/ss

        S(1, 1) = 1.0d0-(Para%gamma-1.0d0)*(u1*u1+v1*v1)/(2.0d0*c1*c1)
        S(1, 2) = (Para%gamma-1.0d0)*u1/(c1*c1)
        S(1, 3) = (Para%gamma-1.0d0)*v1/(c1*c1)
        S(1, 4) = -(Para%gamma-1.0d0)/(c1*c1)

        S(2, 1) = (alfa*v1-beta*u1)/d1
        S(2, 2) = beta/d1
        S(2, 3) = -alfa/d1
        S(2, 4) = 0.0d0

        S(3, 1) = ((Para%gamma-1.0d0)*(u1*u1+v1*v1)/2.0d0-c1*(alfa*u1+beta*v1))/(sqrt(2.0d0)*d1*c1)
        S(3, 2) = (alfa*c1-u1*(Para%gamma-1.0d0))/(sqrt(2.0d0)*d1*c1)
        S(3, 3) = (beta*c1-v1*(Para%gamma-1.0d0))/(sqrt(2.0d0)*d1*c1)
        S(3, 4) = (Para%gamma-1.0d0)/(sqrt(2.0d0)*d1*c1)

        S(4, 1) = ((Para%gamma-1.0d0)*(u1*u1+v1*v1)/2.0d0+c1*(alfa*u1+beta*v1))/(sqrt(2.0d0)*d1*c1)
        S(4, 2) = (-alfa*c1-u1*(Para%gamma-1.0d0))/(sqrt(2.0d0)*d1*c1)
        S(4, 3) = (-beta*c1-v1*(Para%gamma-1.0d0))/(sqrt(2.0d0)*d1*c1)
        S(4, 4) = (Para%gamma-1.0d0)/(sqrt(2.0d0)*d1*c1)

        S1(1, 1) = 1.0d0
        S1(1, 2) = 0.0d0
        S1(1, 3) = d1/(sqrt(2.0d0)*c1)
        S1(1, 4) = d1/(sqrt(2.0d0)*c1)

        S1(2, 1) = u1
        S1(2, 2) = beta*d1
        S1(2, 3) = (alfa*c1+u1)*d1/(sqrt(2.0d0)*c1)
        S1(2, 4) = (-alfa*c1+u1)*d1/(sqrt(2.0d0)*c1)

        S1(3, 1) = v1
        S1(3, 2) = -alfa*d1
        S1(3, 3) = (beta*c1+v1)*d1/(sqrt(2.0d0)*c1)
        S1(3, 4) = (-beta*c1+v1)*d1/(sqrt(2.0d0)*c1)

        S1(4, 1) = (u1*u1+v1*v1)/2.0d0
        S1(4, 2) = (beta*u1-alfa*v1)*d1
        S1(4, 3) = ((c1*c1/(Para%gamma-1.0d0)+(u1*u1+v1*v1)/2.0d0)+(alfa*u1+beta*v1)*c1)*d1/(sqrt(2.0d0)*c1)
        S1(4, 4) = ((c1*c1/(Para%gamma-1.0d0)+(u1*u1+v1*v1)/2.0d0)-(alfa*u1+beta*v1)*c1)*d1/(sqrt(2.0d0)*c1)

        lambda = 0.0d0
        lambda_max = 1.0d-14
        do m1 = Ka1, Kb2
            m = i+m1
            ms = m
            if (m < 1 .and. Bc_type(1) == 1) ms = 1
            if (m > nn .and. Bc_type(2) == 1) ms = nn

            u2 = uu(ms)*A1(ms)+vv(ms)*A2(ms)
            sc = cc(ms)*sqrt(A1(ms)*A1(ms)+A2(ms)*A2(ms))
            if (Para%pp_limiter /= 0) then
                lambda_max = max(lambda_max, abs(u2)+sc)
            else
                eps2 = (max(1.0d-12, delta_ef*(abs(u2)+sc)))**2
                lambda(1) = max(lambda(1), sqrt(u2*u2+eps2))
                lambda(3) = max(lambda(3), sqrt((u2+sc)**2+eps2))
                lambda(4) = max(lambda(4), sqrt((u2-sc)**2+eps2))
            end if
        end do
        lambda(2) = lambda(1)

        do m1 = Ka1, Kb2
            m = i+m1
            ms = m
            if (m < 1 .and. Bc_type(1) == 1) ms = 1
            if (m > nn .and. Bc_type(2) == 1) ms = nn

            rho = dd(ms)
            rhou = dd(ms)*uu(ms)
            rhov = dd(ms)*vv(ms)
            rhoe = dd(ms)*cc(ms)*cc(ms)/(Para%gamma*(Para%gamma-1.0d0))+0.5d0*dd(ms)*(uu(ms)*uu(ms)+vv(ms)*vv(ms))
            u2 = uu(ms)*A1(ms)+vv(ms)*A2(ms)
            pp = dd(ms)*cc(ms)*cc(ms)/Para%gamma

            Ustate = (/rho, rhou, rhov, rhoe/)
            Fstate = (/rho*u2, rhou*u2+A1(ms)*pp, rhov*u2+A2(ms)*pp, rho*u2*Hc(ms)/)

            if (Para%pp_limiter /= 0) then
                hp = 0.5d0*(Ustate+Fstate/lambda_max)
                hm = 0.5d0*(Ustate-Fstate/lambda_max)

                if (m1 >= Ka1 .and. m1 <= Kb1) then
                    do m = 1, 4
                        fpc(m1, m) = S(m, 1)*hp(1)+S(m, 2)*hp(2)+S(m, 3)*hp(3)+S(m, 4)*hp(4)
                    end do
                end if

                if (m1 >= Ka2 .and. m1 <= Kb2) then
                    do m = 1, 4
                        fmc(m1, m) = S(m, 1)*hm(1)+S(m, 2)*hm(2)+S(m, 3)*hm(3)+S(m, 4)*hm(4)
                    end do
                end if
            else
                do m = 1, 4
                    wc = S(m, 1)*Fstate(1)+S(m, 2)*Fstate(2)+S(m, 3)*Fstate(3)+S(m, 4)*Fstate(4)
                    qc = S(m, 1)*Ustate(1)+S(m, 2)*Ustate(2)+S(m, 3)*Ustate(3)+S(m, 4)*Ustate(4)
                    if (m1 >= Ka1 .and. m1 <= Kb1) fpc(m1, m) = 0.5d0*(wc+lambda(m)*qc)
                    if (m1 >= Ka2 .and. m1 <= Kb2) fmc(m1, m) = 0.5d0*(wc-lambda(m)*qc)
                end do
            end if
        end do

        do m = 1, 4
            call flux_x1(Ka1, Kb1, fpc(:, m), h1(m), Scheme%Scheme_Invis, Bc_type, i, nn)
            call flux_x2(Ka2, Kb2, fmc(:, m), h2(m), Scheme%Scheme_Invis, Bc_type, i, nn)
        end do

        do m = 1, 4
            hp(m) = S1(m, 1)*h1(1)+S1(m, 2)*h1(2)+S1(m, 3)*h1(3)+S1(m, 4)*h1(4)
            hm(m) = S1(m, 1)*h2(1)+S1(m, 2)*h2(2)+S1(m, 3)*h2(3)+S1(m, 4)*h2(4)
        end do

        if (Para%pp_limiter /= 0) then
            call split_state_pp(i, LAP1, 1, lambda_max, A1, A2, Hc, cc, dd, uu, vv, qbase)
            call limit_split_endpoint_pp(hp, qbase, 1.0d0/12.0d0, eps_rho, rho_limited, pressure_limited, base_invalid)
            if (rho_limited) mark_rho(i) = mark_rho(i)+1
            if (pressure_limited .or. base_invalid) mark_p(i) = mark_p(i)+1

            call split_state_pp(i+1, LAP1, -1, lambda_max, A1, A2, Hc, cc, dd, uu, vv, qbase)
            call limit_split_endpoint_pp(hm, qbase, 1.0d0/12.0d0, eps_rho, rho_limited, pressure_limited, base_invalid)
            if (rho_limited) mark_rho(i+1) = mark_rho(i+1)+1
            if (pressure_limited .or. base_invalid) mark_p(i+1) = mark_p(i+1)+1

            do m = 1, 4
                hh(i, m) = lambda_max*(hp(m)-hm(m))
            end do
        else
            do m = 1, 4
                hh(i, m) = hp(m)+hm(m)
            end do
        end if
    end do

    if (Bc_type(1) == 1) then
        do m = 1, 4
            hh(0, m) = hh(1, m)
        end do
    end if

    if (Bc_type(2) == 1) then
        do m = 1, 4
            hh(nn, m) = hh(nn-1, m)
        end do
    end if

end subroutine flux_characteristic_pp

subroutine split_density_floor_pp(nn, LAP1, Ka1, Kb2, A1, A2, cc, dd, uu, vv, Bc_type, floor_rho)
    use flow_para
    implicit none
    integer, intent(in) :: nn, LAP1, Ka1, Kb2, Bc_type(2)
    real(dp), intent(in) :: A1(1-LAP1:nn+LAP1), A2(1-LAP1:nn+LAP1)
    real(dp), intent(in) :: cc(1-LAP1:nn+LAP1), dd(1-LAP1:nn+LAP1)
    real(dp), intent(in) :: uu(1-LAP1:nn+LAP1), vv(1-LAP1:nn+LAP1)
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
            u_normal = uu(ms)*A1(ms)+vv(ms)*A2(ms)
            sound_normal = cc(ms)*sqrt(A1(ms)**2+A2(ms)**2)
            speed = max(speed, abs(u_normal)+sound_normal)
        end do
        u_normal = uu(i)*A1(i)+vv(i)*A2(i)
        rho_plus = 0.5d0*dd(i)*(1.0d0+u_normal/speed)
        u_normal = uu(i+1)*A1(i+1)+vv(i+1)*A2(i+1)
        rho_minus = 0.5d0*dd(i+1)*(1.0d0-u_normal/speed)
        floor_rho = min(floor_rho, rho_plus, rho_minus)
    end do
end subroutine split_density_floor_pp

subroutine split_state_pp(cell_index, LAP1, split_sign, lambda, A1, A2, Hc, cc, dd, uu, vv, state)
    use flow_para
    implicit none
    integer, intent(in) :: cell_index, LAP1, split_sign
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: A1(*), A2(*), Hc(*), cc(*), dd(*), uu(*), vv(*)
    real(dp), intent(out) :: state(4)
    real(dp) :: rho, u2, pp, energy

    rho = dd(cell_index+LAP1)
    u2 = uu(cell_index+LAP1)*A1(cell_index+LAP1)+vv(cell_index+LAP1)*A2(cell_index+LAP1)
    pp = dd(cell_index+LAP1)*cc(cell_index+LAP1)*cc(cell_index+LAP1)/Para%gamma
    energy = dd(cell_index+LAP1)*Hc(cell_index+LAP1)-pp    ! rho*E = rho*H - p
    state(1) = 0.5d0*(rho+split_sign*rho*u2/lambda)
    state(2) = 0.5d0*(rho*uu(cell_index+LAP1)+split_sign*(rho*uu(cell_index+LAP1)*u2+A1(cell_index+LAP1)*pp)/lambda)
    state(3) = 0.5d0*(rho*vv(cell_index+LAP1)+split_sign*(rho*vv(cell_index+LAP1)*u2+A2(cell_index+LAP1)*pp)/lambda)
    state(4) = 0.5d0*(energy+split_sign*(rho*u2*Hc(cell_index+LAP1))/lambda)
end subroutine split_state_pp

subroutine limit_split_endpoint_pp(endpoint, base, endpoint_weight, eps_rho, rho_limited, pressure_limited, base_invalid)
    use flow_para
    use flow_data, only: limiter_count_rho, limiter_count_pressure, limiter_count_total
    implicit none
    real(dp), intent(inout) :: endpoint(4)
    real(dp), intent(in) :: base(4), endpoint_weight, eps_rho
    logical, intent(out) :: rho_limited, pressure_limited, base_invalid
    real(dp) :: qstar(4), modified(4), rho_min, theta_rho, theta_p, tm
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
    real(dp), intent(in) :: base(4), target(4)
    real(dp), intent(out) :: theta
    real(dp) :: left, right, mid, p_mid
    real(dp) :: dq(4), a, b, c, sqrt_disc
    integer :: iter
    real(dp) :: pressure_split_pp

    dq = target-base
    a = dq(1)*dq(4)-0.5d0*(dq(2)*dq(2)+dq(3)*dq(3))
    b = base(1)*dq(4)+base(4)*dq(1)-(base(2)*dq(2)+base(3)*dq(3))
    c = base(1)*base(4)-0.5d0*(base(2)*base(2)+base(3)*base(3))
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
    real(dp), intent(in) :: state(4)
    if (state(1) <= 0.0d0) then
        pressure_split_pp = -huge(1.0d0)
    else
        pressure_split_pp = (Para%gamma-1.0d0)*(state(4)-0.5d0*(state(2)*state(2)+state(3)*state(3))/state(1))
    end if
end function pressure_split_pp
