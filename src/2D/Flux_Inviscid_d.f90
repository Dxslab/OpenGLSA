! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid2d_d
    use flow_data
    use fluc_data
    use mesh_data
    implicit none

    integer :: i, j, m
    integer :: Ka1, Kb1, Ka2, Kb2, Bc_type(2)

    integer :: m_int
    real(dp) :: inv_y

    Ka1 = Scheme%Ka1
    Kb1 = Scheme%Kb1
    Ka2 = Scheme%Ka2
    Kb2 = Scheme%Kb2

    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        p_fluc(i, j) = (T_fluc(i, j)/T(i, j)+d_fluc(i, j)/d(i, j))*p(i, j)
        rhoe_fluc(i, j) = p_fluc(i, j)/(Para%gamma-1.0d0)+0.5d0*d_fluc(i, j)*(u(i, j)**2+v(i, j)**2)+ &
                          d(i, j)*(u(i, j)*u_fluc(i, j)+v(i, j)*v_fluc(i, j))
    end do
    end do

    do j = 1, ny
        do i = 1-LAP, nx+LAP
            c1x(i) = c(i, j)
            d1x(i) = d(i, j)
            u1x(i) = u(i, j)
            v1x(i) = v(i, j)
            p2x(i) = p_fluc(i, j)
            d2x(i) = d_fluc(i, j)
            u2x(i) = u_fluc(i, j)
            v2x(i) = v_fluc(i, j)
            w2x(i) = w_fluc(i, j)
            A1x(i) = Akx1(i, j)
            A2x(i) = Aky1(i, j)
        end do

        Bc_type = BcI(j, :)
        if (Para%IF_Character == 2) then
            call split2d_jac_llf_d(nx, LAP, c1x, d1x, u1x, v1x, &
                                   p2x, d2x, u2x, v2x, w2x, A1x, A2x, fpx, fmx)
        else
            call split2d_jac_d(nx, LAP, c1x, d1x, u1x, v1x, &
                               p2x, d2x, u2x, v2x, w2x, A1x, A2x, fpx, fmx)
        end if

        do m = 1, 5
            call OCFD2d_flux1(fpx(1-LAP, m), hhx1, nx, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD2d_flux2(fmx(1-LAP, m), hhx2, nx, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)
            !
            do i = 1, nx
                du(i, j, m) = hhx1(i)-hhx1(i-1)+hhx2(i)-hhx2(i-1)
            end do
        end do
    end do

    do i = 1, nx
        do j = 1-LAP, ny+LAP
            c1y(j) = c(i, j)
            d1y(j) = d(i, j)
            u1y(j) = u(i, j)
            v1y(j) = v(i, j)
            p2y(j) = p_fluc(i, j)
            d2y(j) = d_fluc(i, j)
            u2y(j) = u_fluc(i, j)
            v2y(j) = v_fluc(i, j)
            w2y(j) = w_fluc(i, j)
            A1y(j) = Aix1(i, j)
            A2y(j) = Aiy1(i, j)
        end do

        Bc_type = BcJ(i, :)
        if (Para%IF_Character == 2) then
            call split2d_jac_llf_d(ny, LAP, c1y, d1y, u1y, v1y, &
                                   p2y, d2y, u2y, v2y, w2y, A1y, A2y, fpy, fmy)
        else
            call split2d_jac_d(ny, LAP, c1y, d1y, u1y, v1y, &
                               p2y, d2y, u2y, v2y, w2y, A1y, A2y, fpy, fmy)
        end if

        do m = 1, 5
            call OCFD2d_flux1(fpy(1-LAP, m), hhy1, ny, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD2d_flux2(fmy(1-LAP, m), hhy2, ny, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)
            !
            do j = 1, ny
                du(i, j, m) = -(du(i, j, m)+hhy1(j)-hhy1(j-1)+hhy2(j)-hhy2(j-1))*Ajac(i, j)
            end do
        end do

    end do

    m_int = abs(nint(Para%Wavenumber))
    if (Para%IF_Cylindrical_Coordinate == 1) then
        do j = 1, ny
        do i = 1, nx
            if (abs(Ayy(i, j)) > 1.0d-10) then
                inv_y = 1.0d0/Ayy(i, j)
                Source_Invis(i, j, 1) = Para%Wavenumber*d(i, j)*w_fluc(i, j)*inv_y &
                                        +(d(i, j)*v_fluc(i, j)+d_fluc(i, j)*v(i, j))*inv_y

                Source_Invis(i, j, 2) = Para%Wavenumber*d(i, j)*u(i, j)*w_fluc(i, j)*inv_y &
                                        +(d(i, j)*u_fluc(i, j)*v(i, j)+ &
                                          d(i, j)*u(i, j)*v_fluc(i, j)+ &
                                          d_fluc(i, j)*u(i, j)*v(i, j))*inv_y

                Source_Invis(i, j, 3) = Para%Wavenumber*d(i, j)*v(i, j)*w_fluc(i, j)*inv_y &
                                        +(d_fluc(i, j)*v(i, j)**2+ &
                                          2.0d0*d(i, j)*v(i, j)*v_fluc(i, j))*inv_y

                Source_Invis(i, j, 4) = -Para%Wavenumber*p_fluc(i, j)*inv_y &
                                        +(2.0d0*d(i, j)*v(i, j)*w_fluc(i, j))*inv_y

                Source_Invis(i, j, 5) = Para%Wavenumber*(rhoe(i, j)+p(i, j))*w_fluc(i, j)*inv_y &
                                        +(v_fluc(i, j)*(rhoe(i, j)+p(i, j))+ &
                                          v(i, j)*(rhoe_fluc(i, j)+p_fluc(i, j)))*inv_y
            else
                Source_Invis(i, j, 1) = 0.0d0
                Source_Invis(i, j, 2) = 0.0d0
                Source_Invis(i, j, 3) = 0.0d0
                Source_Invis(i, j, 4) = 0.0d0
                Source_Invis(i, j, 5) = 0.0d0
            end if
        end do
        end do

        ! axis treatment
        if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
            if (mod(m_int, 2) == 0) then
                do j = 1, ny
                    Source_Invis(1, j, 1) = (56.0d0*Source_Invis(2, j, 1)- &
                                             28.0d0*Source_Invis(3, j, 1)+ &
                                             8.0d0*Source_Invis(4, j, 1)- &
                                             Source_Invis(5, j, 1))/35.0d0
                    Source_Invis(1, j, 2) = (56.0d0*Source_Invis(2, j, 2)- &
                                             28.0d0*Source_Invis(3, j, 2)+ &
                                             8.0d0*Source_Invis(4, j, 2)- &
                                             Source_Invis(5, j, 2))/35.0d0
                    Source_Invis(1, j, 3) = 0.0d0
                    Source_Invis(1, j, 4) = 0.0d0
                    Source_Invis(1, j, 5) = (56.0d0*Source_Invis(2, j, 5)- &
                                             28.0d0*Source_Invis(3, j, 5)+ &
                                             8.0d0*Source_Invis(4, j, 5)- &
                                             Source_Invis(5, j, 5))/35.0d0
                end do
            else
                do j = 1, ny
                    Source_Invis(1, j, 1) = 0.0d0
                    Source_Invis(1, j, 2) = 0.0d0
                    Source_Invis(1, j, 3) = (56.0d0*Source_Invis(2, j, 3)- &
                                             28.0d0*Source_Invis(3, j, 3)+ &
                                             8.0d0*Source_Invis(4, j, 3)- &
                                             Source_Invis(5, j, 3))/35.0d0
                    Source_Invis(1, j, 4) = (56.0d0*Source_Invis(2, j, 4)- &
                                             28.0d0*Source_Invis(3, j, 4)+ &
                                             8.0d0*Source_Invis(4, j, 4)- &
                                             Source_Invis(5, j, 4))/35.0d0
                    Source_Invis(1, j, 5) = 0.0d0
                end do
            end if
        end if

        if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
            if (mod(m_int, 2) == 0) then
                do j = 1, ny
                    Source_Invis(nx, j, 1) = (56.0d0*Source_Invis(nx-1, j, 1)- &
                                              28.0d0*Source_Invis(nx-2, j, 1)+ &
                                              8.0d0*Source_Invis(nx-3, j, 1)- &
                                              Source_Invis(nx-4, j, 1))/35.0d0
                    Source_Invis(nx, j, 2) = (56.0d0*Source_Invis(nx-1, j, 2)- &
                                              28.0d0*Source_Invis(nx-2, j, 2)+ &
                                              8.0d0*Source_Invis(nx-3, j, 2)- &
                                              Source_Invis(nx-4, j, 2))/35.0d0
                    Source_Invis(nx, j, 3) = 0.0d0
                    Source_Invis(nx, j, 4) = 0.0d0
                    Source_Invis(nx, j, 5) = (56.0d0*Source_Invis(nx-1, j, 5)- &
                                              28.0d0*Source_Invis(nx-2, j, 5)+ &
                                              8.0d0*Source_Invis(nx-3, j, 5)- &
                                              Source_Invis(nx-4, j, 5))/35.0d0
                end do
            else
                do j = 1, ny
                    Source_Invis(nx, j, 1) = 0.0d0
                    Source_Invis(nx, j, 2) = 0.0d0
                    Source_Invis(nx, j, 3) = (56.0d0*Source_Invis(nx-1, j, 3)- &
                                              28.0d0*Source_Invis(nx-2, j, 3)+ &
                                              8.0d0*Source_Invis(nx-3, j, 3)- &
                                              Source_Invis(nx-4, j, 3))/35.0d0
                    Source_Invis(nx, j, 4) = (56.0d0*Source_Invis(nx-1, j, 4)- &
                                              28.0d0*Source_Invis(nx-2, j, 4)+ &
                                              8.0d0*Source_Invis(nx-3, j, 4)- &
                                              Source_Invis(nx-4, j, 4))/35.0d0
                    Source_Invis(nx, j, 5) = 0.0d0
                end do
            end if
        end if

        if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
            if (mod(m_int, 2) == 0) then
                do i = 1, nx
                    Source_Invis(i, 1, 1) = (56.0d0*Source_Invis(i, 2, 1)- &
                                             28.0d0*Source_Invis(i, 3, 1)+ &
                                             8.0d0*Source_Invis(i, 4, 1)- &
                                             Source_Invis(i, 5, 1))/35.0d0
                    Source_Invis(i, 1, 2) = (56.0d0*Source_Invis(i, 2, 2)- &
                                             28.0d0*Source_Invis(i, 3, 2)+ &
                                             8.0d0*Source_Invis(i, 4, 2)- &
                                             Source_Invis(i, 5, 2))/35.0d0
                    Source_Invis(i, 1, 3) = 0.0d0
                    Source_Invis(i, 1, 4) = 0.0d0
                    Source_Invis(i, 1, 5) = (56.0d0*Source_Invis(i, 2, 5)- &
                                             28.0d0*Source_Invis(i, 3, 5)+ &
                                             8.0d0*Source_Invis(i, 4, 5)- &
                                             Source_Invis(i, 5, 5))/35.0d0
                end do
            else
                do i = 1, nx
                    Source_Invis(i, 1, 1) = 0.0d0
                    Source_Invis(i, 1, 2) = 0.0d0
                    Source_Invis(i, 1, 3) = (56.0d0*Source_Invis(i, 2, 3)- &
                                             28.0d0*Source_Invis(i, 3, 3)+ &
                                             8.0d0*Source_Invis(i, 4, 3)- &
                                             Source_Invis(i, 5, 3))/35.0d0
                    Source_Invis(i, 1, 4) = (56.0d0*Source_Invis(i, 2, 4)- &
                                             28.0d0*Source_Invis(i, 3, 4)+ &
                                             8.0d0*Source_Invis(i, 4, 4)- &
                                             Source_Invis(i, 5, 4))/35.0d0
                    Source_Invis(i, 1, 5) = 0.0d0
                end do
            end if
        end if

        if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
            if (mod(m_int, 2) == 0) then
                do i = 1, nx
                    Source_Invis(i, ny, 1) = (56.0d0*Source_Invis(i, ny-1, 1)- &
                                              28.0d0*Source_Invis(i, ny-2, 1)+ &
                                              8.0d0*Source_Invis(i, ny-3, 1)- &
                                              Source_Invis(i, ny-4, 1))/35.0d0
                    Source_Invis(i, ny, 2) = (56.0d0*Source_Invis(i, ny-1, 2)- &
                                              28.0d0*Source_Invis(i, ny-2, 2)+ &
                                              8.0d0*Source_Invis(i, ny-3, 2)- &
                                              Source_Invis(i, ny-4, 2))/35.0d0
                    Source_Invis(i, ny, 3) = 0.0d0
                    Source_Invis(i, ny, 4) = 0.0d0
                    Source_Invis(i, ny, 5) = (56.0d0*Source_Invis(i, ny-1, 5)- &
                                              28.0d0*Source_Invis(i, ny-2, 5)+ &
                                              8.0d0*Source_Invis(i, ny-3, 5)- &
                                              Source_Invis(i, ny-4, 5))/35.0d0
                end do
            else
                do i = 1, nx
                    Source_Invis(i, ny, 1) = 0.0d0
                    Source_Invis(i, ny, 2) = 0.0d0
                    Source_Invis(i, ny, 3) = (56.0d0*Source_Invis(i, ny-1, 3)- &
                                              28.0d0*Source_Invis(i, ny-2, 3)+ &
                                              8.0d0*Source_Invis(i, ny-3, 3)- &
                                              Source_Invis(i, ny-4, 3))/35.0d0
                    Source_Invis(i, ny, 4) = (56.0d0*Source_Invis(i, ny-1, 4)- &
                                              28.0d0*Source_Invis(i, ny-2, 4)+ &
                                              8.0d0*Source_Invis(i, ny-3, 4)- &
                                              Source_Invis(i, ny-4, 4))/35.0d0
                    Source_Invis(i, ny, 5) = 0.0d0
                end do
            end if
        end if

    else
        do j = 1, ny
        do i = 1, nx
            ! Cartesian coordinate z-direction wavenumber source term
            Source_Invis(i, j, 1) = Para%Wavenumber*d(i, j)*w_fluc(i, j)
            Source_Invis(i, j, 2) = Para%Wavenumber*d(i, j)*u(i, j)*w_fluc(i, j)
            Source_Invis(i, j, 3) = Para%Wavenumber*d(i, j)*v(i, j)*w_fluc(i, j)
            Source_Invis(i, j, 4) = -Para%Wavenumber*p_fluc(i, j)
            Source_Invis(i, j, 5) = Para%Wavenumber*(rhoe(i, j)+p(i, j))*w_fluc(i, j)
        end do
        end do
    end if

    do j = 1, ny
    do i = 1, nx
        du(i, j, 1) = du(i, j, 1)-Source_Invis(i, j, 1)
        du(i, j, 2) = du(i, j, 2)-Source_Invis(i, j, 2)
        du(i, j, 3) = du(i, j, 3)-Source_Invis(i, j, 3)
        du(i, j, 4) = du(i, j, 4)-Source_Invis(i, j, 4)
        du(i, j, 5) = du(i, j, 5)-Source_Invis(i, j, 5)
    end do
    end do

end subroutine flux_inviscid2d_d
