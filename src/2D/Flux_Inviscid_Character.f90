! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid2d_character
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, m
    integer :: Ka1, Kb1, Ka2, Kb2
    integer :: Bc_type(2)

    Ka1 = Scheme%Ka1
    Kb1 = Scheme%Kb1
    Ka2 = Scheme%Ka2
    Kb2 = Scheme%Kb2

    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j) = sqrt(T(i, j))/Para%Ma
        Hc(i, j) = c(i, j)**2/(Para%gamma-1.d0)+0.5d0*(u(i, j)**2+v(i, j)**2)
    end do
    end do

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
        call split2d_jac(nx, LAP, c1x, d1x, u1x, v1x, A1x, A2x, fpx, fmx)
        Bc_type = BcI(j, 1:2)
        call flux_characteristic(nx, LAP, Ka1, Kb1, Ka2, Kb2, &
                                 hhx(0:nx, 1:4), fpx(1-LAP:nx+LAP, 1:4), fmx(1-LAP:nx+LAP, 1:4), &
                                 A1x(1-LAP:nx+LAP), A2x(1-LAP:nx+LAP), &
                                 Hc1x(1-LAP:nx+LAP), d1x(1-LAP:nx+LAP), u1x(1-LAP:nx+LAP), v1x(1-LAP:nx+LAP), &
                                 Bc_type)
        !
        do m = 1, 4
            do i = 1, nx
                du(i, j, m) = hhx(i, m)-hhx(i-1, m)
            end do
        end do
    end do

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
        call split2d_jac(ny, LAP, c1y, d1y, u1y, v1y, A1y, A2y, fpy, fmy)
        Bc_type = BcJ(i, 1:2)
        call flux_characteristic(ny, LAP, Ka1, Kb1, Ka2, Kb2, &
                                 hhy(0:ny, 1:4), fpy(1-LAP:ny+LAP, 1:4), fmy(1-LAP:ny+LAP, 1:4), &
                                 A1y(1-LAP:ny+LAP), A2y(1-LAP:ny+LAP), &
                                 Hc1y(1-LAP:ny+LAP), d1y(1-LAP:ny+LAP), u1y(1-LAP:ny+LAP), v1y(1-LAP:ny+LAP), &
                                 Bc_type)
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

end subroutine flux_inviscid2d_character

subroutine flux_characteristic(nn, LAP1, Ka1, Kb1, Ka2, Kb2, hh, fp, fm, A1, A2, Hc, dd, uu, vv, Bc_type)
    use flow_para
    implicit none

    integer :: Bc_type(2), ib, ie
    integer :: nn, LAP1, Ka1, Kb1, Ka2, Kb2, i, m, m1
    real(dp), dimension(1-LAP1:nn+LAP1):: Hc, dd, uu, vv, A1, A2
    real(dp), dimension(1-LAP1:nn+LAP1, 4):: fp, fm
    real(dp):: hh(0:nn, 4)
    real(dp):: fpc(Ka1:Kb1, 4)
    real(dp):: fmc(Ka2:Kb2, 4)
    real(dp):: h1(4), h2(4), h0(4)
    real(dp):: S(4, 4), S1(4, 4)
    real(dp):: dtmp
    real(dp):: d1, u1, v1, c1, Hc1, alfa, beta, ss

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

        do m1 = Ka1, Kb1
            m = i+m1
            fpc(m1, 1) = S(1, 1)*fp(m, 1)+S(1, 2)*fp(m, 2)+S(1, 3)*fp(m, 3)+S(1, 4)*fp(m, 4)
            fpc(m1, 2) = S(2, 1)*fp(m, 1)+S(2, 2)*fp(m, 2)+S(2, 3)*fp(m, 3)+S(2, 4)*fp(m, 4)
            fpc(m1, 3) = S(3, 1)*fp(m, 1)+S(3, 2)*fp(m, 2)+S(3, 3)*fp(m, 3)+S(3, 4)*fp(m, 4)
            fpc(m1, 4) = S(4, 1)*fp(m, 1)+S(4, 2)*fp(m, 2)+S(4, 3)*fp(m, 3)+S(4, 4)*fp(m, 4)
        end do

        do m1 = Ka2, Kb2
            m = i+m1
            fmc(m1, 1) = S(1, 1)*fm(m, 1)+S(1, 2)*fm(m, 2)+S(1, 3)*fm(m, 3)+S(1, 4)*fm(m, 4)
            fmc(m1, 2) = S(2, 1)*fm(m, 1)+S(2, 2)*fm(m, 2)+S(2, 3)*fm(m, 3)+S(2, 4)*fm(m, 4)
            fmc(m1, 3) = S(3, 1)*fm(m, 1)+S(3, 2)*fm(m, 2)+S(3, 3)*fm(m, 3)+S(3, 4)*fm(m, 4)
            fmc(m1, 4) = S(4, 1)*fm(m, 1)+S(4, 2)*fm(m, 2)+S(4, 3)*fm(m, 3)+S(4, 4)*fm(m, 4)
        end do

        do m = 1, 4
            call flux_x1(Ka1, Kb1, fpc(:, m), h1(m), Scheme%Scheme_Invis, Bc_type, i, nn)
            call flux_x2(Ka2, Kb2, fmc(:, m), h2(m), Scheme%Scheme_Invis, Bc_type, i, nn)
        end do

        h0 = h1+h2
        hh(i, 1) = S1(1, 1)*h0(1)+S1(1, 2)*h0(2)+S1(1, 3)*h0(3)+S1(1, 4)*h0(4)
        hh(i, 2) = S1(2, 1)*h0(1)+S1(2, 2)*h0(2)+S1(2, 3)*h0(3)+S1(2, 4)*h0(4)
        hh(i, 3) = S1(3, 1)*h0(1)+S1(3, 2)*h0(2)+S1(3, 3)*h0(3)+S1(3, 4)*h0(4)
        hh(i, 4) = S1(4, 1)*h0(1)+S1(4, 2)*h0(2)+S1(4, 3)*h0(3)+S1(4, 4)*h0(4)
    end do

    if (Bc_type(1) == 1) then
        do m = 1, 4
            hh(0, m) = (2.d0*fp(1, m)+5.d0*fp(2, m)-fp(3, m)+ &
                        2.d0*fm(3, m)-7.d0*fm(2, m)+11.d0*fm(1, m))/6.d0
        end do
    end if

    if (Bc_type(2) == 1) then
        do m = 1, 4
            hh(nn, m) = (2.d0*fp(nn-2, m)-7.d0*fp(nn-1, m)+11.d0*fp(nn, m)+ &
                         2.d0*fm(nn, m)+5.d0*fm(nn-1, m)-fm(nn-2, m))/6.d0
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
                ! call hh_WENO5P_boundary(Ka, Kb, ff, hh, del_left)
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
                ! call hh_WENO5P_boundary(Ka, Kb, ff, hh, del_right)
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
                ! call hh_WENO5M_boundary(Ka, Kb, ff, hh, del_left)
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
                ! call hh_WENO5M_boundary(Ka, Kb, ff, hh, del_right)
                call hh_WENOMR3M(Ka, Kb, ff, hh)
            end if
        else
            call hh_WENOMR5M(Ka, Kb, ff, hh)
        end if
        !
    end if

end subroutine flux_x2
