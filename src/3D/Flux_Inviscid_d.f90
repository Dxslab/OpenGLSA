! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid_d
    use flow_data
    use fluc_data
    use mesh_data
    implicit none

    integer :: i, j, k, m
    integer :: Ka1, Kb1, Ka2, Kb2, Bc_type(2)

    Ka1 = Scheme%Ka1
    Kb1 = Scheme%Kb1
    Ka2 = Scheme%Ka2
    Kb2 = Scheme%Kb2

    do k = 1-LAP, nz+LAP
    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        p_fluc(i, j, k) = (T_fluc(i, j, k)/T(i, j, k)+d_fluc(i, j, k)/d(i, j, k))*p(i, j, k)
        rhoe_fluc(i, j, k) = p_fluc(i, j, k)/(Para%gamma-1.0d0)+ &
                             0.5d0*d_fluc(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)+ &
                             d(i, j, k)*(u(i, j, k)*u_fluc(i, j, k)+ &
                                         v(i, j, k)*v_fluc(i, j, k)+ &
                                         w(i, j, k)*w_fluc(i, j, k))
    end do
    end do
    end do

    do k = 1, nz
    do j = 1, ny
        do i = 1-LAP, nx+LAP
            c1x(i) = c(i, j, k)
            d1x(i) = d(i, j, k)
            u1x(i) = u(i, j, k)
            v1x(i) = v(i, j, k)
            w1x(i) = w(i, j, k)
            p2x(i) = p_fluc(i, j, k)
            d2x(i) = d_fluc(i, j, k)
            u2x(i) = u_fluc(i, j, k)
            v2x(i) = v_fluc(i, j, k)
            w2x(i) = w_fluc(i, j, k)
            A1x(i) = Akx1(i, j, k)
            A2x(i) = Aky1(i, j, k)
            A3x(i) = Akz1(i, j, k)
        end do

        if (Para%IF_Character == 2) then
            call split_jac_llf_d(nx, LAP, c1x, d1x, u1x, v1x, w1x, &
                                 p2x, d2x, u2x, v2x, w2x, A1x, A2x, A3x, fpx, fmx)
        else
            call split_jac_d(nx, LAP, c1x, d1x, u1x, v1x, w1x, &
                             p2x, d2x, u2x, v2x, w2x, A1x, A2x, A3x, fpx, fmx)
        end if

        Bc_type = BcI(j, k, :)
        do m = 1, 5
            call OCFD_flux1(fpx(1-LAP, m), hhx1, nx, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD_flux2(fmx(1-LAP, m), hhx2, nx, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)
            !
            do i = 1, nx
                du(i, j, k, m) = hhx1(i)-hhx1(i-1)+hhx2(i)-hhx2(i-1)
            end do
        end do
    end do
    end do

    do k = 1, nz
    do i = 1, nx
        do j = 1-LAP, ny+LAP
            c1y(j) = c(i, j, k)
            d1y(j) = d(i, j, k)
            u1y(j) = u(i, j, k)
            v1y(j) = v(i, j, k)
            w1y(j) = w(i, j, k)
            p2y(j) = p_fluc(i, j, k)
            d2y(j) = d_fluc(i, j, k)
            u2y(j) = u_fluc(i, j, k)
            v2y(j) = v_fluc(i, j, k)
            w2y(j) = w_fluc(i, j, k)
            A1y(j) = Aix1(i, j, k)
            A2y(j) = Aiy1(i, j, k)
            A3y(j) = Aiz1(i, j, k)
        end do

        if (Para%IF_Character == 2) then
            call split_jac_llf_d(ny, LAP, c1y, d1y, u1y, v1y, w1y, &
                                 p2y, d2y, u2y, v2y, w2y, A1y, A2y, A3y, fpy, fmy)
        else
            call split_jac_d(ny, LAP, c1y, d1y, u1y, v1y, w1y, &
                             p2y, d2y, u2y, v2y, w2y, A1y, A2y, A3y, fpy, fmy)
        end if

        Bc_type = BcJ(i, k, :)
        do m = 1, 5
            call OCFD_flux1(fpy(1-LAP, m), hhy1, ny, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD_flux2(fmy(1-LAP, m), hhy2, ny, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)
            !
            do j = 1, ny
                du(i, j, k, m) = du(i, j, k, m)+hhy1(j)-hhy1(j-1)+hhy2(j)-hhy2(j-1)
            end do
        end do
    end do
    end do

    do j = 1, ny
    do i = 1, nx
        do k = 1-LAP, nz+LAP
            c1z(k) = c(i, j, k)
            d1z(k) = d(i, j, k)
            u1z(k) = u(i, j, k)
            v1z(k) = v(i, j, k)
            w1z(k) = w(i, j, k)
            p2z(k) = p_fluc(i, j, k)
            d2z(k) = d_fluc(i, j, k)
            u2z(k) = u_fluc(i, j, k)
            v2z(k) = v_fluc(i, j, k)
            w2z(k) = w_fluc(i, j, k)
            A1z(k) = Asx1(i, j, k)
            A2z(k) = Asy1(i, j, k)
            A3z(k) = Asz1(i, j, k)
        end do

        if (Para%IF_Character == 2) then
            call split_jac_llf_d(nz, LAP, c1z, d1z, u1z, v1z, w1z, &
                                 p2z, d2z, u2z, v2z, w2z, A1z, A2z, A3z, fpz, fmz)
        else
            call split_jac_d(nz, LAP, c1z, d1z, u1z, v1z, w1z, &
                             p2z, d2z, u2z, v2z, w2z, A1z, A2z, A3z, fpz, fmz)
        end if

        Bc_type = BcK(i, j, :)
        do m = 1, 5
            call OCFD_flux1(fpz(1-LAP, m), hhz1, nz, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD_flux2(fmz(1-LAP, m), hhz2, nz, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)
            !
            do k = 1, nz
                du(i, j, k, m) = -(du(i, j, k, m)+hhz1(k)-hhz1(k-1)+hhz2(k)-hhz2(k-1))*Ajac(i, j, k)
            end do
        end do
    end do
    end do

end subroutine flux_inviscid_d
