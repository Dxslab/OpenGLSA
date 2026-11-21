! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid_character
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, k, m
    integer :: Ka1, Kb1, Ka2, Kb2
    integer :: Bc_type(2)

    Ka1 = Scheme%Ka1
    Kb1 = Scheme%Kb1
    Ka2 = Scheme%Ka2
    Kb2 = Scheme%Kb2

    do k = 1-LAP, nz+LAP
    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j, k) = sqrt(T(i, j, k))/Para%Ma
        Hc(i, j, k) = c(i, j, k)**2/(Para%gamma-1.d0)+0.5d0*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)
    end do
    end do
    end do

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
        call split_jac(nx, LAP, c1x, d1x, u1x, v1x, w1x, A1x, A2x, A3x, fpx, fmx)
        Bc_type = BcI(j, k, 1:2)
        call flux_characteristic(nx, LAP, Ka1, Kb1, Ka2, Kb2, &
                                 hhx(0:nx, 1:5), fpx(1-LAP:nx+LAP, 1:5), fmx(1-LAP:nx+LAP, 1:5), &
                                 A1x(1-LAP:nx+LAP), A2x(1-LAP:nx+LAP), A3x(1-LAP:nx+LAP), &
                                 Hc1x(1-LAP:nx+LAP), d1x(1-LAP:nx+LAP), u1x(1-LAP:nx+LAP), v1x(1-LAP:nx+LAP), w1x(1-LAP:nx+LAP), &
                                 Bc_type)
        !
        do m = 1, 5
            do i = 1, nx
                du(i, j, k, m) = hhx(i, m)-hhx(i-1, m)
            end do
        end do
    end do
    end do

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
        call split_jac(ny, LAP, c1y, d1y, u1y, v1y, w1y, A1y, A2y, A3y, fpy, fmy)
        Bc_type = BcJ(i, k, 1:2)
        call flux_characteristic(ny, LAP, Ka1, Kb1, Ka2, Kb2, &
                                 hhy(0:ny, 1:5), fpy(1-LAP:ny+LAP, 1:5), fmy(1-LAP:ny+LAP, 1:5), &
                                 A1y(1-LAP:ny+LAP), A2y(1-LAP:ny+LAP), A3y(1-LAP:ny+LAP), &
                                 Hc1y(1-LAP:ny+LAP), d1y(1-LAP:ny+LAP), u1y(1-LAP:ny+LAP), v1y(1-LAP:ny+LAP), w1y(1-LAP:ny+LAP), &
                                 Bc_type)
        !
        do m = 1, 5
            do j = 1, ny
                du(i, j, k, m) = du(i, j, k, m)+hhy(j, m)-hhy(j-1, m)
            end do
        end do
    end do
    end do

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
        call split_jac(nz, LAP, c1z, d1z, u1z, v1z, w1z, A1z, A2z, A3z, fpz, fmz)
        Bc_type = BcK(i, j, 1:2)
        call flux_characteristic(nz, LAP, Ka1, Kb1, Ka2, Kb2, &
                                 hhz(0:nz, 1:5), fpz(1-LAP:nz+LAP, 1:5), fmz(1-LAP:nz+LAP, 1:5), &
                                 A1z(1-LAP:nz+LAP), A2z(1-LAP:nz+LAP), A3z(1-LAP:nz+LAP), &
                                 Hc1z(1-LAP:nz+LAP), d1z(1-LAP:nz+LAP), u1z(1-LAP:nz+LAP), v1z(1-LAP:nz+LAP), w1z(1-LAP:nz+LAP), &
                                 Bc_type)
        !
        do m = 1, 5
            do k = 1, nz
                du(i, j, k, m) = -(du(i, j, k, m)+hhz(k, m)-hhz(k-1, m))*Ajac(i, j, k)
            end do
        end do
    end do
    end do

end subroutine flux_inviscid_character

subroutine flux_characteristic(nn, LAP1, Ka1, Kb1, Ka2, Kb2, hh, fp, fm, AA1, AA2, AA3, Hc, dd, uu, vv, ww, Bc_type)
    use flow_para
    implicit none

    integer :: Bc_type(2), nn, LAP1, Ka1, Kb1, Ka2, Kb2, i, m, mm, ib, ie
    real(dp), dimension(1-LAP1:nn+LAP1):: dd, uu, vv, ww, AA1, AA2, AA3, Hc
    real(dp), dimension(1-LAP1:nn+LAP1, 5):: fp, fm
    real(dp):: hh(0:nn, 5)
    real(dp):: fpc(Ka1:Kb1, 5), fmc(Ka2:Kb2, 5), h1(5), h2(5), h0(5)
    real(dp):: S(5, 5), S1(5, 5)
    real(dp):: d1, u1, v1, w1, c1, V2, un, ul, um, dtmp, Hc1
    real(dp):: a1, a2, a3, ss, n1, n2, n3, l1, l2, l3, m1, m2, m3, KK1, KK, X1, H

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

        do mm = Ka1, Kb1
            m = i+mm
            fpc(mm, 1) = S(1, 1)*fp(m, 1)+S(1, 2)*fp(m, 2)+S(1, 3)*fp(m, 3)+S(1, 4)*fp(m, 4)+S(1, 5)*fp(m, 5)
            fpc(mm, 2) = S(2, 1)*fp(m, 1)+S(2, 2)*fp(m, 2)+S(2, 3)*fp(m, 3)+S(2, 4)*fp(m, 4)+S(2, 5)*fp(m, 5)
            fpc(mm, 3) = S(3, 1)*fp(m, 1)+S(3, 2)*fp(m, 2)+S(3, 3)*fp(m, 3)+S(3, 4)*fp(m, 4)+S(3, 5)*fp(m, 5)
            fpc(mm, 4) = S(4, 1)*fp(m, 1)+S(4, 2)*fp(m, 2)+S(4, 3)*fp(m, 3)+S(4, 4)*fp(m, 4)+S(4, 5)*fp(m, 5)
            fpc(mm, 5) = S(5, 1)*fp(m, 1)+S(5, 2)*fp(m, 2)+S(5, 3)*fp(m, 3)+S(5, 4)*fp(m, 4)+S(5, 5)*fp(m, 5)
        end do

        do mm = Ka2, Kb2
            m = i+mm
            fmc(mm, 1) = S(1, 1)*fm(m, 1)+S(1, 2)*fm(m, 2)+S(1, 3)*fm(m, 3)+S(1, 4)*fm(m, 4)+S(1, 5)*fm(m, 5)
            fmc(mm, 2) = S(2, 1)*fm(m, 1)+S(2, 2)*fm(m, 2)+S(2, 3)*fm(m, 3)+S(2, 4)*fm(m, 4)+S(2, 5)*fm(m, 5)
            fmc(mm, 3) = S(3, 1)*fm(m, 1)+S(3, 2)*fm(m, 2)+S(3, 3)*fm(m, 3)+S(3, 4)*fm(m, 4)+S(3, 5)*fm(m, 5)
            fmc(mm, 4) = S(4, 1)*fm(m, 1)+S(4, 2)*fm(m, 2)+S(4, 3)*fm(m, 3)+S(4, 4)*fm(m, 4)+S(4, 5)*fm(m, 5)
            fmc(mm, 5) = S(5, 1)*fm(m, 1)+S(5, 2)*fm(m, 2)+S(5, 3)*fm(m, 3)+S(5, 4)*fm(m, 4)+S(5, 5)*fm(m, 5)
        end do

        do m = 1, 5
            call flux_x1(Ka1, Kb1, fpc(:, m), h1(m), Scheme%Scheme_Invis, Bc_type, i, nn)
            call flux_x2(Ka2, Kb2, fmc(:, m), h2(m), Scheme%Scheme_Invis, Bc_type, i, nn)
        end do

        h0 = h1+h2

        hh(i, 1) = S1(1, 1)*h0(1)+S1(1, 2)*h0(2)+S1(1, 3)*h0(3)+S1(1, 4)*h0(4)+S1(1, 5)*h0(5)
        hh(i, 2) = S1(2, 1)*h0(1)+S1(2, 2)*h0(2)+S1(2, 3)*h0(3)+S1(2, 4)*h0(4)+S1(2, 5)*h0(5)
        hh(i, 3) = S1(3, 1)*h0(1)+S1(3, 2)*h0(2)+S1(3, 3)*h0(3)+S1(3, 4)*h0(4)+S1(3, 5)*h0(5)
        hh(i, 4) = S1(4, 1)*h0(1)+S1(4, 2)*h0(2)+S1(4, 3)*h0(3)+S1(4, 4)*h0(4)+S1(4, 5)*h0(5)
        hh(i, 5) = S1(5, 1)*h0(1)+S1(5, 2)*h0(2)+S1(5, 3)*h0(3)+S1(5, 4)*h0(4)+S1(5, 5)*h0(5)
    end do

    if (Bc_type(1) == 1) then
        do m = 1, 5
            hh(0, m) = (2.d0*fp(1, m)+5.d0*fp(2, m)-fp(3, m)+2.d0*fm(3, m)-7.d0*fm(2, m)+11.d0*fm(1, m))/6.d0
        end do
    end if

    if (Bc_type(2) == 1) then
        do m = 1, 5
            hh(nn, m) = (2.d0*fp(nn-2, m)-7.d0*fp(nn-1, m)+11.d0*fp(nn, m)+2.d0*fm(nn, m)+5.d0*fm(nn-1, m)-fm(nn-2, m))/6.d0
        end do
    end if

end subroutine flux_characteristic

subroutine flux_x1(Ka, Kb, ff, hh, Num_Scheme, Bc_type, i, nn)
    use Constants
    use WENO_Reconstruction
    implicit none

    integer :: Bc_type(2), IERR
    integer :: Num_Scheme, Ka, Kb, i, nn
    real(dp) :: ff(Ka:Kb), hh

    if ((Bc_type(1) == 0 .or. i > -Ka) .and. (Bc_type(2) == 0 .or. i <= nn-Kb)) then
        !
        select case (Num_Scheme)
        case (OCFD_Scheme_WENOZ5)
            call hh_WENO5ZP(Ka, Kb, ff, hh)
        case (OCFD_Scheme_WENOZQ5)
            call hh_WENOZQ5P(Ka, Kb, ff, hh)
        case (OCFD_Scheme_WENOMR5)
            call hh_WENOMR5P(Ka, Kb, ff, hh)
        case (OCFD_Scheme_LINEAR5)
            call hh_LINEAR5P(Ka, Kb, ff, hh)
        case DEFAULT
            print *, '*** Error: Unknown scheme number in flux_x1 ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        end select
        !
    else if (i <= -Ka .and. Bc_type(1) == 1) then
        !
        if (i == 1) then
            hh = (2.0d0*ff(0)+5.0d0*ff(1)-ff(2))/6.0d0
        else if (i == 2) then
            if (Num_Scheme == OCFD_Scheme_LINEAR5) then
                call hh_LINEAR5P_boundary(Ka, Kb, ff, hh, del_left)
            else
                call hh_WENOMR3P(Ka, Kb, ff, hh)
            end if
        else
            call hh_WENOMR5P(Ka, Kb, ff, hh)
        end if
        !
    else if (i > nn-Kb .and. Bc_type(2) == 1) then
        !
        if (i == nn-1) then
            if (Num_Scheme == OCFD_Scheme_LINEAR5) then
                call hh_LINEAR5P_boundary(Ka, Kb, ff, hh, del_right)
            else
                call hh_WENOMR3P(Ka, Kb, ff, hh)
            end if
        else
            call hh_WENOMR5P(Ka, Kb, ff, hh)
        end if
        !
    end if

end subroutine flux_x1

subroutine flux_x2(Ka, Kb, ff, hh, Num_Scheme, Bc_type, i, nn)
    use Constants
    use WENO_Reconstruction
    implicit none

    integer :: Bc_type(2), IERR
    integer :: Num_Scheme, Ka, Kb, i, nn
    real(dp) ::  ff(Ka:Kb), hh

    if ((Bc_type(1) == 0 .or. i > -Ka) .and. (Bc_type(2) == 0 .or. i <= nn-Kb)) then
        !
        select case (Num_Scheme)
        case (OCFD_Scheme_WENOZ5)
            call hh_WENO5ZM(Ka, Kb, ff, hh)
        case (OCFD_Scheme_WENOZQ5)
            call hh_WENOZQ5M(Ka, Kb, ff, hh)
        case (OCFD_Scheme_WENOMR5)
            call hh_WENOMR5M(Ka, Kb, ff, hh)
        case (OCFD_Scheme_LINEAR5)
            call hh_LINEAR5M(Ka, Kb, ff, hh)
        case DEFAULT
            print *, '*** Error: Unknown scheme number in flux_x2 ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        end select
        !
    else if (i <= -Ka .and. Bc_type(1) == 1) then
        !
        if (i == 1) then
            if (Num_Scheme == OCFD_Scheme_LINEAR5) then
                call hh_LINEAR5M_boundary(Ka, Kb, ff, hh, del_left)
            else
                call hh_WENOMR3M(Ka, Kb, ff, hh)
            end if
        else
            call hh_WENOMR5M(Ka, Kb, ff, hh)
        end if
        !
    else if (i > nn-Kb .and. Bc_type(2) == 1) then
        !
        if (i == nn-1) then
            hh = (2.0d0*ff(1)+5.0d0*ff(0)-ff(-1))/6.0d0
        else if (i == nn-2) then
            if (Num_Scheme == OCFD_Scheme_LINEAR5) then
                call hh_LINEAR5M_boundary(Ka, Kb, ff, hh, del_right)
            else
                call hh_WENOMR3M(Ka, Kb, ff, hh)
            end if
        else
            call hh_WENOMR5M(Ka, Kb, ff, hh)
        end if
        !
    end if

end subroutine flux_x2
