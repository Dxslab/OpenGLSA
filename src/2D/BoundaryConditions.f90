! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine ocfd2d_bc_Flat_Plate
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j

    real(dp) :: x0, Tw

    x0 = Para%BC_para(1)
    Tw = Para%BC_para(2)

    if (npx == 0) then
        do j = 1, ny
            d(1, j) = 1.0d0
            u(1, j) = 1.0d0
            v(1, j) = 0.0d0
            T(1, j) = 1.0d0
        end do
    end if

    if (npy == npy0-1) then
        do i = 1, nx
            d(i, ny) = 1.0d0
            u(i, ny) = 1.0d0
            v(i, ny) = 0.0d0
            T(i, ny) = 1.0d0
        end do
    end if

    if (npy == 0) then
        do i = 1, nx
            if (Axx(i, 1) < x0) then
                d(i, 1) = (4.0d0*d(i, 2)-d(i, 3))/3.0d0
                u(i, 1) = (4.0d0*u(i, 2)-u(i, 3))/3.0d0
                v(i, 1) = 0.0d0
                T(i, 1) = (4.0d0*T(i, 2)-T(i, 3))/3.0d0
            else
                u(i, 1) = 0.0d0
                v(i, 1) = 0.0d0
                if (Tw <= 0) then
                    d(i, 1) = (4.0d0*d(i, 2)-d(i, 3))/3.0d0
                    T(i, 1) = (4.0d0*T(i, 2)-T(i, 3))/3.0d0
                else
                    T(i, 1) = Tw
                    d(i, 1) = (4.0d0*d(i, 2)*T(i, 2)-d(i, 3)*T(i, 3))/(3.0d0*Tw)
                end if
            end if
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            d(nx, j) = 2.0d0*d(nx-1, j)-d(nx-2, j)
            u(nx, j) = 2.0d0*u(nx-1, j)-u(nx-2, j)
            v(nx, j) = 2.0d0*v(nx-1, j)-v(nx-2, j)
            T(nx, j) = 2.0d0*T(nx-1, j)-T(nx-2, j)
        end do
    end if

end subroutine ocfd2d_bc_Flat_Plate

subroutine ocfd2d_bc_Flat_Plate_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    integer :: i, j

    real(dp) :: x0, Tw, p1_fluc, p2_fluc, pb_fluc

    x0 = Para%BC_para(1)
    Tw = Para%BC_para(2)

    if (npx == 0) then
        do j = 1, ny
            d_fluc(1, j) = 0.0d0
            u_fluc(1, j) = 0.0d0
            v_fluc(1, j) = 0.0d0
            w_fluc(1, j) = 0.0d0
            T_fluc(1, j) = 0.0d0
        end do
    end if

    if (npy == npy0-1) then
        do i = 1, nx
            d_fluc(i, ny) = 0.0d0
            u_fluc(i, ny) = 0.0d0
            v_fluc(i, ny) = 0.0d0
            w_fluc(i, ny) = 0.0d0
            T_fluc(i, ny) = 0.0d0
        end do
    end if

    if (npy == 0) then
        do i = 1, nx
            if (Axx(i, 1) < x0) then
                d_fluc(i, 1) = (4.0d0*d_fluc(i, 2)-d_fluc(i, 3))/3.0d0
                u_fluc(i, 1) = (4.0d0*u_fluc(i, 2)-u_fluc(i, 3))/3.0d0
                v_fluc(i, 1) = 0.0d0
                w_fluc(i, 1) = 0.0d0
                T_fluc(i, 1) = (4.0d0*T_fluc(i, 2)-T_fluc(i, 3))/3.0d0
            else
                u_fluc(i, 1) = 0.0d0
                v_fluc(i, 1) = 0.0d0
                w_fluc(i, 1) = 0.0d0
                if (Tw <= 0) then
                    d_fluc(i, 1) = (4.0d0*d_fluc(i, 2)-d_fluc(i, 3))/3.0d0
                    T_fluc(i, 1) = (4.0d0*T_fluc(i, 2)-T_fluc(i, 3))/3.0d0
                else
                    T_fluc(i, 1) = 0.0d0
                    p1_fluc = d_fluc(i, 2)*T(i, 2)+d(i, 2)*T_fluc(i, 2)
                    p2_fluc = d_fluc(i, 3)*T(i, 3)+d(i, 3)*T_fluc(i, 3)
                    pb_fluc = (4.0d0*p1_fluc-p2_fluc)/3.0d0
                    d_fluc(i, 1) = pb_fluc/Tw
                end if
            end if
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            d_fluc(nx, j) = 2.0d0*d_fluc(nx-1, j)-d_fluc(nx-2, j)
            u_fluc(nx, j) = 2.0d0*u_fluc(nx-1, j)-u_fluc(nx-2, j)
            v_fluc(nx, j) = 2.0d0*v_fluc(nx-1, j)-v_fluc(nx-2, j)
            w_fluc(nx, j) = 2.0d0*w_fluc(nx-1, j)-w_fluc(nx-2, j)
            T_fluc(nx, j) = 2.0d0*T_fluc(nx-1, j)-T_fluc(nx-2, j)
        end do
    end if

end subroutine ocfd2d_bc_Flat_Plate_d

subroutine ocfd2d_bc_Lid_Driven_Cavity
    use mesh_data
    use flow_data
    use Centrals
    use Flux_Viscous, only: diff_mu
    implicit none

    integer :: i, j
    real(dp) :: Rgas, R_inv, muw
    real(dp) :: pb_bot(nx), pb_top(nx), pb_lft(ny), pb_rgt(ny)

    real(dp), parameter :: twothird = 2.0d0/3.0d0, onethird = 1.0d0/3.0d0, fourthird = 4.0d0/3.0d0
    real(dp), parameter :: c_os4 = 12.0d0/25.0d0

    muw = 1.0d0/Para%Re
    Rgas = 1.0d0/(Para%Ma**2*Para%gamma)
    R_inv = Para%Ma**2*Para%gamma

    if (npy == npy0-1) then
        do i = 1, nx
            u(i, ny) = (1.0d0-(2.0d0*Axx(i, ny)-1.0d0)**18)**2
            v(i, ny) = 0.0d0
            T(i, ny) = 1.0d0
        end do
    end if

    if (npx == 0) then
        do j = 1, ny
            u(1, j) = 0.0d0
            v(1, j) = 0.0d0
            T(1, j) = 1.0d0
        end do
    end if

    if (npy == 0) then
        do i = 1, nx
            u(i, 1) = 0.0d0
            v(i, 1) = 0.0d0
            T(i, 1) = 1.0d0
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            u(nx, j) = 0.0d0
            v(nx, j) = 0.0d0
            T(nx, j) = 1.0d0
        end do
    end if

    call mpi_exchange_xy2d(u)
    call mpi_exchange_xy2d(v)
    call mpi_exchange_xy2d(T)

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
    call OCFD2d_dy2(u, d2udi2, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy2(v, d2vdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy2(T, d2Tdi2, nx, ny, Scheme%Scheme_Vis, BcJ)

    if (npy == 0) then
        do i = 1, nx
            pb_bot(i) = pwall_eta(i, 1, 1)
        end do
    end if

    if (npy == npy0-1) then
        do i = 1, nx
            pb_top(i) = pwall_eta(i, ny, -1)
        end do
    end if

    if (npx == 0) then
        do j = 1, ny
            pb_lft(j) = pwall_xi(1, j, 1)
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            pb_rgt(j) = pwall_xi(nx, j, -1)
        end do
    end if

    if (npx == 0 .and. npy == 0) call corner_mean(pb_bot(1), pb_lft(1))
    if (npx == 0 .and. npy == npy0-1) call corner_mean(pb_top(1), pb_lft(ny))
    if (npx == npx0-1 .and. npy == 0) call corner_mean(pb_bot(nx), pb_rgt(1))
    if (npx == npx0-1 .and. npy == npy0-1) call corner_mean(pb_top(nx), pb_rgt(ny))

    if (npy == 0) d(1:nx, 1) = pb_bot*R_inv
    if (npy == npy0-1) d(1:nx, ny) = pb_top*R_inv
    if (npx == 0) d(1, 1:ny) = pb_lft*R_inv
    if (npx == npx0-1) d(nx, 1:ny) = pb_rgt*R_inv

contains

    real(dp) function pwall_eta(ii, jj, inward)
        integer, intent(IN) :: ii, jj, inward
        real(dp) :: p1, p2, p3, p4, dmuwdk, dmuwdi, mu_term, dmudk_term, dmudi_term

        p1 = d(ii, jj+inward)*T(ii, jj+inward)*Rgas
        p2 = d(ii, jj+2*inward)*T(ii, jj+2*inward)*Rgas
        p3 = d(ii, jj+3*inward)*T(ii, jj+3*inward)*Rgas
        p4 = d(ii, jj+4*inward)*T(ii, jj+4*inward)*Rgas

        dmuwdk = diff_mu(dTdk(ii, jj), T(ii, jj))
        dmuwdi = diff_mu(dTdi(ii, jj), T(ii, jj))

        dmudk_term = Akx(ii, jj)*dudi(ii, jj)+Akx(ii, jj)**2*dvdk(ii, jj)/Aiy(ii, jj)
        dmudi_term = Aiy(ii, jj)*dvdi(ii, jj)*fourthird-Akx(ii, jj)*dudk(ii, jj)*twothird
        mu_term = Aiy(ii, jj)*d2vdi2(ii, jj)*fourthird+Aiyi(ii, jj)*dvdi(ii, jj)*fourthird+ &
                  Akx(ii, jj)*d2udkdi(ii, jj)*onethird-Akxi(ii, jj)*dudk(ii, jj)*twothird+ &
                  Aiyk(ii, jj)*Akx(ii, jj)*dudi(ii, jj)/Aiy(ii, jj)+Akx(ii, jj)**2*d2vdk2(ii, jj)/Aiy(ii, jj)+ &
                  Akx(ii, jj)*Akxk(ii, jj)*dvdk(ii, jj)/Aiy(ii, jj)

        pwall_eta = (48.0d0*p1-36.0d0*p2+16.0d0*p3-3.0d0*p4)/25.0d0 &
                    -inward*c_os4*(muw*mu_term+dmuwdk*dmudk_term+dmuwdi*dmudi_term)
    end function pwall_eta

    real(dp) function pwall_xi(ii, jj, inward)
        integer, intent(IN) :: ii, jj, inward
        real(dp) :: p1, p2, p3, p4, dmuwdk, dmuwdi, mu_term, dmudk_term, dmudi_term

        p1 = d(ii+inward, jj)*T(ii+inward, jj)*Rgas
        p2 = d(ii+2*inward, jj)*T(ii+2*inward, jj)*Rgas
        p3 = d(ii+3*inward, jj)*T(ii+3*inward, jj)*Rgas
        p4 = d(ii+4*inward, jj)*T(ii+4*inward, jj)*Rgas

        dmuwdk = diff_mu(dTdk(ii, jj), T(ii, jj))
        dmuwdi = diff_mu(dTdi(ii, jj), T(ii, jj))

        dmudk_term = -Aiy(ii, jj)*dvdi(ii, jj)*twothird+Akx(ii, jj)*dudk(ii, jj)*fourthird
        dmudi_term = Aiy(ii, jj)**2*dudi(ii, jj)/Akx(ii, jj)+Aiy(ii, jj)*dvdk(ii, jj)
        mu_term = Aiy(ii, jj)**2*d2udi2(ii, jj)/Akx(ii, jj)+Aiy(ii, jj)*Aiyi(ii, jj)*dudi(ii, jj)/Akx(ii, jj)+ &
                  Aiy(ii, jj)*d2vdkdi(ii, jj)*onethird+Aiy(ii, jj)*Akxi(ii, jj)*dvdk(ii, jj)/Akx(ii, jj)- &
                  Aiyk(ii, jj)*dvdi(ii, jj)*twothird+Akx(ii, jj)*d2udk2(ii, jj)*fourthird+Akxk(ii, jj)*dudk(ii, jj)*fourthird

        pwall_xi = (48.0d0*p1-36.0d0*p2+16.0d0*p3-3.0d0*p4)/25.0d0 &
                   -inward*c_os4*(muw*mu_term+dmuwdk*dmudk_term+dmuwdi*dmudi_term)
    end function pwall_xi

    subroutine corner_mean(pa, pb)
        real(dp), intent(INOUT) :: pa, pb

        pa = 0.5d0*(pa+pb)
        pb = pa
    end subroutine corner_mean

end subroutine ocfd2d_bc_Lid_Driven_Cavity

subroutine ocfd2d_bc_Lid_Driven_Cavity_d
    use mesh_data
    use flow_data
    use fluc_data
    use Centrals
    use Flux_Viscous_d, only: diff_mu_d
    implicit none

    integer :: i, j
    real(dp) :: Rgas, R_inv, muw
    real(dp) :: pb_bot(nx), pb_top(nx), pb_lft(ny), pb_rgt(ny)

    real(dp), parameter :: twothird = 2.0d0/3.0d0, onethird = 1.0d0/3.0d0, fourthird = 4.0d0/3.0d0
    real(dp), parameter :: c_os4 = 12.0d0/25.0d0

    muw = 1.0d0/Para%Re
    Rgas = 1.0d0/(Para%Ma**2*Para%gamma)
    R_inv = Para%Ma**2*Para%gamma

    if (npy == npy0-1) then
        do i = 1, nx
            u_fluc(i, ny) = 0.0d0
            v_fluc(i, ny) = 0.0d0
            w_fluc(i, ny) = 0.0d0
            T_fluc(i, ny) = 0.0d0
        end do
    end if

    if (npx == 0) then
        do j = 1, ny
            u_fluc(1, j) = 0.0d0
            v_fluc(1, j) = 0.0d0
            w_fluc(1, j) = 0.0d0
            T_fluc(1, j) = 0.0d0
        end do
    end if

    if (npy == 0) then
        do i = 1, nx
            u_fluc(i, 1) = 0.0d0
            v_fluc(i, 1) = 0.0d0
            w_fluc(i, 1) = 0.0d0
            T_fluc(i, 1) = 0.0d0
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            u_fluc(nx, j) = 0.0d0
            v_fluc(nx, j) = 0.0d0
            w_fluc(nx, j) = 0.0d0
            T_fluc(nx, j) = 0.0d0
        end do
    end if

    call mpi_exchange_xy2d(u_fluc)
    call mpi_exchange_xy2d(v_fluc)
    call mpi_exchange_xy2d(w_fluc)
    call mpi_exchange_xy2d(T_fluc)

    call OCFD2d_dx0(u_fluc, dufdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx0(v_fluc, dvfdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx0(w_fluc, dwfdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx0(T_fluc, dTfdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)

    call OCFD2d_dy0(u_fluc, dufdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(v_fluc, dvfdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(w_fluc, dwfdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(T_fluc, dTfdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)

    call mpi_exchange_y2d(dufdk)
    call mpi_exchange_y2d(dvfdk)
    call mpi_exchange_y2d(dwfdk)
    call mpi_exchange_y2d(dTfdk)

    call OCFD2d_dy0(dufdk, d2ufdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(dvfdk, d2vfdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(dwfdk, d2wfdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(dTfdk, d2Tfdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)

    call OCFD2d_dx2(u_fluc, d2ufdk2, nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx2(v_fluc, d2vfdk2, nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx2(w_fluc, d2wfdk2, nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx2(T_fluc, d2Tfdk2, nx, ny, Scheme%Scheme_Vis, BcI)

    call OCFD2d_dy2(u_fluc, d2ufdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy2(v_fluc, d2vfdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy2(w_fluc, d2wfdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy2(T_fluc, d2Tfdi2, nx, ny, Scheme%Scheme_Vis, BcJ)

    if (npy == 0) then
        do i = 1, nx
            pb_bot(i) = pwall_eta_d(i, 1, 1)
        end do
    end if

    if (npy == npy0-1) then
        do i = 1, nx
            pb_top(i) = pwall_eta_d(i, ny, -1)
        end do
    end if

    if (npx == 0) then
        do j = 1, ny
            pb_lft(j) = pwall_xi_d(1, j, 1)
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            pb_rgt(j) = pwall_xi_d(nx, j, -1)
        end do
    end if

    if (npx == 0 .and. npy == 0) call corner_mean(pb_bot(1), pb_lft(1))
    if (npx == 0 .and. npy == npy0-1) call corner_mean(pb_top(1), pb_lft(ny))
    if (npx == npx0-1 .and. npy == 0) call corner_mean(pb_bot(nx), pb_rgt(1))
    if (npx == npx0-1 .and. npy == npy0-1) call corner_mean(pb_top(nx), pb_rgt(ny))

    if (npy == 0) d_fluc(1:nx, 1) = pb_bot*R_inv-d(1:nx, 1)*T_fluc(1:nx, 1)
    if (npy == npy0-1) d_fluc(1:nx, ny) = pb_top*R_inv-d(1:nx, ny)*T_fluc(1:nx, ny)
    if (npx == 0) d_fluc(1, 1:ny) = pb_lft*R_inv-d(1, 1:ny)*T_fluc(1, 1:ny)
    if (npx == npx0-1) d_fluc(nx, 1:ny) = pb_rgt*R_inv-d(nx, 1:ny)*T_fluc(nx, 1:ny)

contains

    real(dp) function pwall_eta_d(ii, jj, inward)
        integer, intent(IN) :: ii, jj, inward
        real(dp) :: p1_fluc, p2_fluc, p3_fluc, p4_fluc, dmufwdk, dmufwdi
        real(dp) :: mu_term, dmudk_term, dmudi_term, dmufdk_term, dmufdi_term

        dmufwdk = diff_mu_d(dTfdk(ii, jj), T_fluc(ii, jj), dTdk(ii, jj), T(ii, jj))
        dmufwdi = diff_mu_d(dTfdi(ii, jj), T_fluc(ii, jj), dTdi(ii, jj), T(ii, jj))

        dmufdk_term = Akx(ii, jj)*dudi(ii, jj)+Akx(ii, jj)**2*dvdk(ii, jj)/Aiy(ii, jj)
        dmufdi_term = Aiy(ii, jj)*dvdi(ii, jj)*fourthird-Akx(ii, jj)*dudk(ii, jj)*twothird
        dmudk_term = Akx(ii, jj)*dufdi(ii, jj)+Akx(ii, jj)**2*dvfdk(ii, jj)/Aiy(ii, jj)
        dmudi_term = Aiy(ii, jj)*dvfdi(ii, jj)*fourthird-Akx(ii, jj)*dufdk(ii, jj)*twothird- &
                     twothird*w_fluc(ii, jj)*Para%wavenumber
        mu_term = Aiy(ii, jj)*d2vfdi2(ii, jj)*fourthird+Aiyi(ii, jj)*dvfdi(ii, jj)*fourthird+ &
                  Akx(ii, jj)*d2ufdkdi(ii, jj)*onethird-Akxi(ii, jj)*dufdk(ii, jj)*twothird+ &
                  dwfdi(ii, jj)*onethird*Para%wavenumber+Aiyk(ii, jj)*Akx(ii, jj)*dufdi(ii, jj)/Aiy(ii, jj)+ &
                  Akx(ii, jj)**2*d2vfdk2(ii, jj)/Aiy(ii, jj)+Akx(ii, jj)*Akxk(ii, jj)*dvfdk(ii, jj)/Aiy(ii, jj)- &
                  v_fluc(ii, jj)*Para%wavenumber**2/Aiy(ii, jj)

        p1_fluc = (d_fluc(ii, jj+inward)*T(ii, jj+inward)+d(ii, jj+inward)*T_fluc(ii, jj+inward))*Rgas
        p2_fluc = (d_fluc(ii, jj+2*inward)*T(ii, jj+2*inward)+d(ii, jj+2*inward)*T_fluc(ii, jj+2*inward))*Rgas
        p3_fluc = (d_fluc(ii, jj+3*inward)*T(ii, jj+3*inward)+d(ii, jj+3*inward)*T_fluc(ii, jj+3*inward))*Rgas
        p4_fluc = (d_fluc(ii, jj+4*inward)*T(ii, jj+4*inward)+d(ii, jj+4*inward)*T_fluc(ii, jj+4*inward))*Rgas
        pwall_eta_d = (48.0d0*p1_fluc-36.0d0*p2_fluc+16.0d0*p3_fluc-3.0d0*p4_fluc)/25.0d0- &
                      inward*c_os4*(muw*mu_term+dmudk(ii, jj)*dmudk_term+dmudi(ii, jj)*dmudi_term+ &
                                    dmufwdk*dmufdk_term+dmufwdi*dmufdi_term)
    end function pwall_eta_d

    real(dp) function pwall_xi_d(ii, jj, inward)
        integer, intent(IN) :: ii, jj, inward
        real(dp) :: p1_fluc, p2_fluc, p3_fluc, p4_fluc, dmufwdk, dmufwdi
        real(dp) :: mu_term, dmudk_term, dmudi_term, dmufdk_term, dmufdi_term

        dmufwdk = diff_mu_d(dTfdk(ii, jj), T_fluc(ii, jj), dTdk(ii, jj), T(ii, jj))
        dmufwdi = diff_mu_d(dTfdi(ii, jj), T_fluc(ii, jj), dTdi(ii, jj), T(ii, jj))

        dmufdk_term = -Aiy(ii, jj)*dvdi(ii, jj)*twothird+Akx(ii, jj)*dudk(ii, jj)*fourthird
        dmufdi_term = Aiy(ii, jj)**2*dudi(ii, jj)/Akx(ii, jj)+Aiy(ii, jj)*dvdk(ii, jj)
        dmudk_term = -Aiy(ii, jj)*dvfdi(ii, jj)*twothird+Akx(ii, jj)*dufdk(ii, jj)*fourthird- &
                     twothird*w_fluc(ii, jj)*Para%wavenumber
        dmudi_term = Aiy(ii, jj)**2*dufdi(ii, jj)/Akx(ii, jj)+Aiy(ii, jj)*dvfdk(ii, jj)
        mu_term = Aiy(ii, jj)**2*d2ufdi2(ii, jj)/Akx(ii, jj)+Aiy(ii, jj)*Aiyi(ii, jj)*dufdi(ii, jj)/Akx(ii, jj)+ &
                  Aiy(ii, jj)*d2vfdkdi(ii, jj)*onethird+Aiy(ii, jj)*Akxi(ii, jj)*dvfdk(ii, jj)/Akx(ii, jj)- &
                  Aiyk(ii, jj)*dvfdi(ii, jj)*twothird+Akx(ii, jj)*d2ufdk2(ii, jj)*fourthird+ &
                  Akxk(ii, jj)*dufdk(ii, jj)*fourthird+dwfdk(ii, jj)*onethird*Para%wavenumber- &
                  u_fluc(ii, jj)*Para%wavenumber**2/Akx(ii, jj)

        p1_fluc = (d_fluc(ii+inward, jj)*T(ii+inward, jj)+d(ii+inward, jj)*T_fluc(ii+inward, jj))*Rgas
        p2_fluc = (d_fluc(ii+2*inward, jj)*T(ii+2*inward, jj)+d(ii+2*inward, jj)*T_fluc(ii+2*inward, jj))*Rgas
        p3_fluc = (d_fluc(ii+3*inward, jj)*T(ii+3*inward, jj)+d(ii+3*inward, jj)*T_fluc(ii+3*inward, jj))*Rgas
        p4_fluc = (d_fluc(ii+4*inward, jj)*T(ii+4*inward, jj)+d(ii+4*inward, jj)*T_fluc(ii+4*inward, jj))*Rgas
        pwall_xi_d = (48.0d0*p1_fluc-36.0d0*p2_fluc+16.0d0*p3_fluc-3.0d0*p4_fluc)/25.0d0- &
                     inward*c_os4*(muw*mu_term+dmudk(ii, jj)*dmudk_term+dmudi(ii, jj)*dmudi_term+ &
                                   dmufwdk*dmufdk_term+dmufwdi*dmufdi_term)
    end function pwall_xi_d

    subroutine corner_mean(pa, pb)
        real(dp), intent(INOUT) :: pa, pb

        pa = 0.5d0*(pa+pb)
        pb = pa
    end subroutine corner_mean

end subroutine ocfd2d_bc_Lid_Driven_Cavity_d

subroutine farfield_top
    use mesh_data
    use flow_data
    implicit none

    integer :: i

    real(dp) :: n1, n2, t1, t2, normn, dot_pos
    real(dp) :: gamma, rgp, eps
    real(dp) :: d_inf, u_inf, v_inf, T_inf, p_inf, un_inf, c_inf
    real(dp) :: d_e, u_e, v_e, T_e, p_e, c_e, un_e, Ma_b, Zimp
    real(dp) :: d_b, u_b, v_b, T_b, p_b

    gamma = Para%gamma
    rgp = gamma*Para%Ma*Para%Ma
    eps = 1.0d-12

    d_inf = 1.0d0
    u_inf = 1.0d0
    v_inf = 0.0d0
    T_inf = 1.0d0
    p_inf = 1.0d0/rgp
    c_inf = sqrt(gamma*p_inf/d_inf)

    do i = 1, nx
        t1 = Axx(i+1, ny)-Axx(i-1, ny)
        t2 = Ayy(i+1, ny)-Ayy(i-1, ny)
        n1 = t2; n2 = -t1
        normn = sqrt(n1*n1+n2*n2)+eps
        n1 = n1/normn; n2 = n2/normn

        dot_pos = n1*(Axx(i, ny-1)-Axx(i, ny))+n2*(Ayy(i, ny-1)-Ayy(i, ny))
        if (dot_pos > 0.0d0) then
            n1 = -n1; n2 = -n2
        end if

        d_e = 2.0d0*d(i, ny-1)-d(i, ny-2)
        u_e = 2.0d0*u(i, ny-1)-u(i, ny-2)
        v_e = 2.0d0*v(i, ny-1)-v(i, ny-2)
        T_e = 2.0d0*T(i, ny-1)-T(i, ny-2)
        p_e = d_e*T_e/rgp
        c_e = sqrt(gamma*p_e/d_e)

        Zimp = d_e*c_e

        un_e = n1*u_e+n2*v_e
        un_inf = n1*u_inf+n2*v_inf
        Ma_b = un_e/c_e

        if (Ma_b <= -1.0d0) then
            ! supersonic inflow
            d_b = d_inf; u_b = u_inf; v_b = v_inf; T_b = T_inf; p_b = p_inf
        else if (Ma_b >= 1.0d0) then
            ! supersonic outflow
            d_b = d_e; u_b = u_e; v_b = v_e; T_b = T_e; p_b = p_e
        else if (Ma_b <= 0.0d0) then
            ! subsonic inflow: entropy and tangential velocity from freestream
            p_b = 0.5d0*(p_e+p_inf)-0.5d0*Zimp*(un_inf-un_e)
            d_b = d_inf+(p_b-p_inf)/(c_e**2)
            u_b = u_inf+n1*(p_b-p_inf)/Zimp
            v_b = v_inf+n2*(p_b-p_inf)/Zimp
            T_b = p_b*rgp/d_b
        else
            ! subsonic outflow: entropy and tangential velocity convected out
            p_b = 0.5d0*(p_e+p_inf)-0.5d0*Zimp*(un_inf-un_e)
            d_b = d_e+(p_b-p_e)/(c_e**2)
            u_b = u_e-n1*(p_b-p_e)/Zimp
            v_b = v_e-n2*(p_b-p_e)/Zimp
            T_b = p_b*rgp/d_b
        end if

        d(i, ny) = d_b
        u(i, ny) = u_b
        v(i, ny) = v_b
        T(i, ny) = T_b
    end do

end subroutine farfield_top

subroutine farfield_top_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    integer :: i

    real(dp) :: n1, n2, t1, t2, normn, dot_pos
    real(dp) :: gamma, rgp, eps
    real(dp) :: d_inf, u_inf, v_inf, p_inf, un_inf, c_inf
    real(dp) :: d_e, u_e, v_e, T_e, p_e, c_e, un_e, Ma_b, Zimp, Zimpd
    real(dp) :: d_ed, u_ed, v_ed, w_ed, T_ed, p_ed, un_ed, c_ed
    real(dp) :: d_b, p_b, dp_b
    real(dp) :: d_bd, u_bd, v_bd, w_bd, T_bd, p_bd, dp_bd

    gamma = Para%gamma
    rgp = gamma*Para%Ma*Para%Ma
    eps = 1.0d-12

    d_inf = 1.0d0
    u_inf = 1.0d0
    v_inf = 0.0d0
    p_inf = 1.0d0/rgp
    c_inf = sqrt(gamma*p_inf/d_inf)

    d_bd = 0.0d0; u_bd = 0.0d0; v_bd = 0.0d0; w_bd = 0.0d0; T_bd = 0.0d0

    do i = 1, nx
        t1 = Axx(i+1, ny)-Axx(i-1, ny)
        t2 = Ayy(i+1, ny)-Ayy(i-1, ny)
        n1 = t2; n2 = -t1
        normn = sqrt(n1*n1+n2*n2)+eps
        n1 = n1/normn; n2 = n2/normn

        dot_pos = n1*(Axx(i, ny-1)-Axx(i, ny))+n2*(Ayy(i, ny-1)-Ayy(i, ny))
        if (dot_pos > 0.0d0) then
            n1 = -n1; n2 = -n2
        end if

        d_e = 2.0d0*d(i, ny-1)-d(i, ny-2)
        u_e = 2.0d0*u(i, ny-1)-u(i, ny-2)
        v_e = 2.0d0*v(i, ny-1)-v(i, ny-2)
        T_e = 2.0d0*T(i, ny-1)-T(i, ny-2)
        p_e = d_e*T_e/rgp
        c_e = sqrt(gamma*p_e/d_e)

        d_ed = 2.0d0*d_fluc(i, ny-1)-d_fluc(i, ny-2)
        u_ed = 2.0d0*u_fluc(i, ny-1)-u_fluc(i, ny-2)
        v_ed = 2.0d0*v_fluc(i, ny-1)-v_fluc(i, ny-2)
        w_ed = 2.0d0*w_fluc(i, ny-1)-w_fluc(i, ny-2)
        T_ed = 2.0d0*T_fluc(i, ny-1)-T_fluc(i, ny-2)
        p_ed = (T_e*d_ed+d_e*T_ed)/rgp
        c_ed = 0.5d0*gamma*T_ed/(rgp*c_e)

        Zimp = d_e*c_e
        Zimpd = d_ed*c_e+d_e*c_ed

        un_e = n1*u_e+n2*v_e
        un_inf = n1*u_inf+n2*v_inf
        Ma_b = un_e/c_e

        un_ed = n1*u_ed+n2*v_ed

        if (Ma_b <= -1.0d0) then
            ! supersonic inflow
            d_bd = 0.0d0; u_bd = 0.0d0; v_bd = 0.0d0; w_bd = 0.0d0; T_bd = 0.0d0
        else if (Ma_b >= 1.0d0) then
            ! supersonic outflow
            d_bd = d_ed; u_bd = u_ed; v_bd = v_ed; w_bd = w_ed; T_bd = T_ed
        else if (Ma_b <= 0.0d0) then
            ! subsonic inflow.  Baseflow:
            p_b = 0.5d0*(p_e+p_inf)-0.5d0*Zimp*(un_inf-un_e)
            p_bd = 0.5d0*p_ed-0.5d0*Zimpd*(un_inf-un_e)+0.5d0*Zimp*un_ed
            dp_b = p_b-p_inf
            dp_bd = p_bd
            d_b = d_inf+dp_b/c_e**2
            d_bd = dp_bd/c_e**2-2.0d0*dp_b*c_ed/c_e**3
            u_bd = n1*(dp_bd/Zimp-dp_b*Zimpd/Zimp**2)
            v_bd = n2*(dp_bd/Zimp-dp_b*Zimpd/Zimp**2)
            w_bd = 0.0d0
            T_bd = rgp*(p_bd-p_b*d_bd/d_b)/d_b
        else
            ! subsonic outflow.  Baseflow:
            p_b = 0.5d0*(p_e+p_inf)-0.5d0*Zimp*(un_inf-un_e)
            p_bd = 0.5d0*p_ed-0.5d0*Zimpd*(un_inf-un_e)+0.5d0*Zimp*un_ed
            dp_b = p_b-p_e
            dp_bd = p_bd-p_ed
            d_b = d_e+dp_b/c_e**2
            d_bd = d_ed+dp_bd/c_e**2-2.0d0*dp_b*c_ed/c_e**3
            u_bd = u_ed-n1*(dp_bd/Zimp-dp_b*Zimpd/Zimp**2)
            v_bd = v_ed-n2*(dp_bd/Zimp-dp_b*Zimpd/Zimp**2)
            w_bd = w_ed
            T_bd = rgp*(p_bd-p_b*d_bd/d_b)/d_b
        end if

        d_fluc(i, ny) = d_bd
        u_fluc(i, ny) = u_bd
        v_fluc(i, ny) = v_bd
        w_fluc(i, ny) = w_bd
        T_fluc(i, ny) = T_bd
    end do

end subroutine farfield_top_d

subroutine ocfd2d_bc_Sphere
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, i0

    if (npy == 0) then
        do i = 1, nx
            u(i, 1) = 0.0d0
            v(i, 1) = 0.0d0
            T(i, 1) = (4.0d0*T(i, 2)-T(i, 3))/3.0d0
            d(i, 1) = (4.0d0*d(i, 2)*T(i, 2)-d(i, 3)*T(i, 3))/(3.0d0*T(i, 1))
        end do
    end if

    if (npy == npy0-1) call farfield_top

    if (npx == 0) then
        do j = 1, ny
            v(1, j) = 0.0d0
            do i0 = 1-LAP, 0
                d(i0, j) = d(2-i0, j)
                u(i0, j) = u(2-i0, j)
                v(i0, j) = -v(2-i0, j)
                T(i0, j) = T(2-i0, j)
            end do
        end do
    end if

    if (npx == npx0-1) then
        do j = 1, ny
            v(nx, j) = 0.0d0
            do i0 = nx+1, nx+LAP
                d(i0, j) = d(2*nx-i0, j)
                u(i0, j) = u(2*nx-i0, j)
                v(i0, j) = -v(2*nx-i0, j)
                T(i0, j) = T(2*nx-i0, j)
            end do
        end do
    end if

end subroutine ocfd2d_bc_Sphere

subroutine ocfd2d_bc_Sphere_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    integer :: i, j, m_int, i0
    real(dp) :: p1_fluc, p2_fluc, pb_fluc

    m_int = nint(Para%wavenumber)

    if (npy == 0) then
        do i = 1, nx
            u_fluc(i, 1) = 0.0d0
            v_fluc(i, 1) = 0.0d0
            w_fluc(i, 1) = 0.0d0
            T_fluc(i, 1) = (4.0d0*T_fluc(i, 2)-T_fluc(i, 3))/3.0d0
            p1_fluc = d_fluc(i, 2)*T(i, 2)+d(i, 2)*T_fluc(i, 2)
            p2_fluc = d_fluc(i, 3)*T(i, 3)+d(i, 3)*T_fluc(i, 3)
            pb_fluc = (4.0d0*p1_fluc-p2_fluc)/3.0d0
            d_fluc(i, 1) = (pb_fluc-d(i, 1)*T_fluc(i, 1))/T(i, 1)
        end do
    end if

    if (npy == npy0-1) call farfield_top_d

    if (mod(m_int, 2) == 0) then

        if (npx == 0) then
            do j = 1, ny
                v_fluc(1, j) = 0.0d0
                w_fluc(1, j) = 0.0d0
                do i0 = 1-LAP, 0
                    d_fluc(i0, j) = d_fluc(2-i0, j)
                    u_fluc(i0, j) = u_fluc(2-i0, j)
                    v_fluc(i0, j) = -v_fluc(2-i0, j)
                    w_fluc(i0, j) = -w_fluc(2-i0, j)
                    T_fluc(i0, j) = T_fluc(2-i0, j)
                end do
            end do
        end if

        if (npx == npx0-1) then
            do j = 1, ny
                v_fluc(nx, j) = 0.0d0
                w_fluc(nx, j) = 0.0d0
                do i0 = nx+1, nx+LAP
                    d_fluc(i0, j) = d_fluc(2*nx-i0, j)
                    u_fluc(i0, j) = u_fluc(2*nx-i0, j)
                    v_fluc(i0, j) = -v_fluc(2*nx-i0, j)
                    w_fluc(i0, j) = -w_fluc(2*nx-i0, j)
                    T_fluc(i0, j) = T_fluc(2*nx-i0, j)
                end do
            end do
        end if

    else

        if (npx == 0) then
            do j = 1, ny
                d_fluc(1, j) = 0.0d0
                u_fluc(1, j) = 0.0d0
                T_fluc(1, j) = 0.0d0
                do i0 = 1-LAP, 0
                    d_fluc(i0, j) = -d_fluc(2-i0, j)
                    u_fluc(i0, j) = -u_fluc(2-i0, j)
                    v_fluc(i0, j) = v_fluc(2-i0, j)
                    w_fluc(i0, j) = w_fluc(2-i0, j)
                    T_fluc(i0, j) = -T_fluc(2-i0, j)
                end do
            end do
        end if

        if (npx == npx0-1) then
            do j = 1, ny
                d_fluc(nx, j) = 0.0d0
                u_fluc(nx, j) = 0.0d0
                T_fluc(nx, j) = 0.0d0
                do i0 = nx+1, nx+LAP
                    d_fluc(i0, j) = -d_fluc(2*nx-i0, j)
                    u_fluc(i0, j) = -u_fluc(2*nx-i0, j)
                    v_fluc(i0, j) = v_fluc(2*nx-i0, j)
                    w_fluc(i0, j) = w_fluc(2*nx-i0, j)
                    T_fluc(i0, j) = -T_fluc(2*nx-i0, j)
                end do
            end do
        end if

    end if

end subroutine ocfd2d_bc_Sphere_d

subroutine ocfd2d_bc_Cylinder
    use mesh_data
    use flow_data
    implicit none

    integer :: i

    if (npy == 0) then
        do i = 1, nx
            u(i, 1) = 0.0d0
            v(i, 1) = 0.0d0
            T(i, 1) = (4.0d0*T(i, 2)-T(i, 3))/3.0d0
            d(i, 1) = (4.0d0*d(i, 2)*T(i, 2)-d(i, 3)*T(i, 3))/(3.0d0*T(i, 1))
        end do
    end if

    ! far field
    if (npy == npy0-1) call farfield_top

end subroutine ocfd2d_bc_Cylinder

subroutine ocfd2d_bc_Cylinder_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    integer :: i
    real(dp) :: p1_fluc, p2_fluc, pb_fluc

    if (npy == 0) then
        do i = 1, nx
            u_fluc(i, 1) = 0.0d0
            v_fluc(i, 1) = 0.0d0
            w_fluc(i, 1) = 0.0d0
            T_fluc(i, 1) = (4.0d0*T_fluc(i, 2)-T_fluc(i, 3))/3.0d0
            p1_fluc = d_fluc(i, 2)*T(i, 2)+d(i, 2)*T_fluc(i, 2)
            p2_fluc = d_fluc(i, 3)*T(i, 3)+d(i, 3)*T_fluc(i, 3)
            pb_fluc = (4.0d0*p1_fluc-p2_fluc)/3.0d0
            d_fluc(i, 1) = (pb_fluc-d(i, 1)*T_fluc(i, 1))/T(i, 1)
        end do
    end if

    ! far field
    if (npy == npy0-1) call farfield_top_d

end subroutine ocfd2d_bc_Cylinder_d

subroutine ocfd2d_bc_User_Def
    use mesh_data
    use flow_data
    implicit none

end subroutine ocfd2d_bc_User_Def

subroutine ocfd2d_bc_User_Def_d
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

end subroutine ocfd2d_bc_User_Def_d
