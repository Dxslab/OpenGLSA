! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! CopyRight by Li Xinliang  Email: lixl@imech.ac.cn
! Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS, lixl@imech.ac.cn
! Copyright by Li Xinliang
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine Allocate_arrays
    implicit none

    call Allocate_mesh_data
    call Allocate_flow_data
    call Allocate_invis_data
    call Allocate_vis_data

end subroutine Allocate_arrays

subroutine Deallocate_arrays
    implicit none

    call Deallocate_mesh_data
    call Deallocate_flow_data
    call Deallocate_invis_data
    call Deallocate_vis_data

end subroutine Deallocate_arrays

subroutine Allocate_mesh_data
    use mesh_data
    implicit none

    allocate (Axx(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Ayy(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Azz(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))

    allocate (Akx(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aix(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Asx(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aky(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aiy(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Asy(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Akz(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aiz(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Asz(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Ajac(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))

    allocate (Akx1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aix1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Asx1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aky1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aiy1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Asy1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Akz1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Aiz1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Asz1(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))

    allocate (BcI(ny, nz, 2)); BcI = 0
    allocate (BcJ(nx, nz, 2)); BcJ = 0
    allocate (BcK(nx, ny, 2)); BcK = 0

    if (Para%Iperiodic_X == 0) then
        if (npx == 0 .and. Para%Isymmetric_X(1) == 0) BcI(:, :, 1) = 1
        if (npx == npx0-1 .and. Para%Isymmetric_X(2) == 0) BcI(:, :, 2) = 1
    end if
    if (Para%Iperiodic_Y == 0) then
        if (npy == 0 .and. Para%Isymmetric_Y(1) == 0) BcJ(:, :, 1) = 1
        if (npy == npy0-1 .and. Para%Isymmetric_Y(2) == 0) BcJ(:, :, 2) = 1
    end if
    if (Para%Iperiodic_Z == 0) then
        if (npz == 0 .and. Para%Isymmetric_Z(1) == 0) BcK(:, :, 1) = 1
        if (npz == npz0-1 .and. Para%Isymmetric_Z(2) == 0) BcK(:, :, 2) = 1
    end if

    if (Para%IF_Viscous == 1) then
        allocate (Akxk(1:nx, 1:ny, 1:nz))
        allocate (Akxi(1:nx, 1:ny, 1:nz))
        allocate (Akxs(1:nx, 1:ny, 1:nz))
        allocate (Aixk(1:nx, 1:ny, 1:nz))
        allocate (Aixi(1:nx, 1:ny, 1:nz))
        allocate (Aixs(1:nx, 1:ny, 1:nz))
        allocate (Asxk(1:nx, 1:ny, 1:nz))
        allocate (Asxi(1:nx, 1:ny, 1:nz))
        allocate (Asxs(1:nx, 1:ny, 1:nz))

        allocate (Akyk(1:nx, 1:ny, 1:nz))
        allocate (Akyi(1:nx, 1:ny, 1:nz))
        allocate (Akys(1:nx, 1:ny, 1:nz))
        allocate (Aiyk(1:nx, 1:ny, 1:nz))
        allocate (Aiyi(1:nx, 1:ny, 1:nz))
        allocate (Aiys(1:nx, 1:ny, 1:nz))
        allocate (Asyk(1:nx, 1:ny, 1:nz))
        allocate (Asyi(1:nx, 1:ny, 1:nz))
        allocate (Asys(1:nx, 1:ny, 1:nz))

        allocate (Akzk(1:nx, 1:ny, 1:nz))
        allocate (Akzi(1:nx, 1:ny, 1:nz))
        allocate (Akzs(1:nx, 1:ny, 1:nz))
        allocate (Aizk(1:nx, 1:ny, 1:nz))
        allocate (Aizi(1:nx, 1:ny, 1:nz))
        allocate (Aizs(1:nx, 1:ny, 1:nz))
        allocate (Aszk(1:nx, 1:ny, 1:nz))
        allocate (Aszi(1:nx, 1:ny, 1:nz))
        allocate (Aszs(1:nx, 1:ny, 1:nz))
    end if
end subroutine Allocate_mesh_data

subroutine Deallocate_mesh_data
    use mesh_data
    implicit none

    deallocate (Axx)
    deallocate (Ayy)
    deallocate (Azz)
    deallocate (Akx)
    deallocate (Aix)
    deallocate (Asx)
    deallocate (Aky)
    deallocate (Aiy)
    deallocate (Asy)
    deallocate (Akz)
    deallocate (Aiz)
    deallocate (Asz)
    deallocate (Ajac)
    deallocate (Akx1)
    deallocate (Aix1)
    deallocate (Asx1)
    deallocate (Aky1)
    deallocate (Aiy1)
    deallocate (Asy1)
    deallocate (Akz1)
    deallocate (Aiz1)
    deallocate (Asz1)

    deallocate (BcI)
    deallocate (BcJ)
    deallocate (BcK)

    if (Para%IF_Viscous == 1) then
        deallocate (Akxk)
        deallocate (Akxi)
        deallocate (Akxs)
        deallocate (Aixk)
        deallocate (Aixi)
        deallocate (Aixs)
        deallocate (Asxk)
        deallocate (Asxi)
        deallocate (Asxs)

        deallocate (Akyk)
        deallocate (Akyi)
        deallocate (Akys)
        deallocate (Aiyk)
        deallocate (Aiyi)
        deallocate (Aiys)
        deallocate (Asyk)
        deallocate (Asyi)
        deallocate (Asys)

        deallocate (Akzk)
        deallocate (Akzi)
        deallocate (Akzs)
        deallocate (Aizk)
        deallocate (Aizi)
        deallocate (Aizs)
        deallocate (Aszk)
        deallocate (Aszi)
        deallocate (Aszs)
    end if
end subroutine Deallocate_mesh_data

subroutine Allocate_flow_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    allocate (d(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (u(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (v(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (w(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (T(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (p(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (c(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (Hc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (rhoe(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))

    allocate (f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP, 5))
    allocate (fn(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP, 5))
    allocate (du(1:nx, 1:ny, 1:nz, 5))

    if (Para%runMode == 0) then
        if (Para%IF_SFD == 1) then
            allocate (fsfd(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP, 5))
            fsfd = 0.0d0
        end if
    else
        allocate (d_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (u_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (v_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (w_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (T_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (p_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (rhoe_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    end if

end subroutine Allocate_flow_data

subroutine Deallocate_flow_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    deallocate (d)
    deallocate (u)
    deallocate (v)
    deallocate (w)
    deallocate (T)
    deallocate (p)
    deallocate (c)
    deallocate (Hc)
    deallocate (rhoe)

    deallocate (f)
    deallocate (fn)
    deallocate (du)

    if (Para%runMode == 0) then
        if (Para%IF_SFD == 1) then
            deallocate (fsfd)
        end if
    else
        deallocate (d_fluc)
        deallocate (u_fluc)
        deallocate (v_fluc)
        deallocate (w_fluc)
        deallocate (T_fluc)
        deallocate (p_fluc)
        deallocate (rhoe_fluc)
    end if

end subroutine Deallocate_flow_data

subroutine Allocate_invis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    allocate (hhx(0:nx, 5))
    allocate (hhy(0:ny, 5))
    allocate (hhz(0:nz, 5))
    allocate (Hc1x(1-LAP:nx+LAP), c1x(1-LAP:nx+LAP), d1x(1-LAP:nx+LAP), u1x(1-LAP:nx+LAP), v1x(1-LAP:nx+LAP), w1x(1-LAP:nx+LAP))
    allocate (Hc1y(1-LAP:ny+LAP), c1y(1-LAP:ny+LAP), d1y(1-LAP:ny+LAP), u1y(1-LAP:ny+LAP), v1y(1-LAP:ny+LAP), w1y(1-LAP:ny+LAP))
    allocate (Hc1z(1-LAP:nz+LAP), c1z(1-LAP:nz+LAP), d1z(1-LAP:nz+LAP), u1z(1-LAP:nz+LAP), v1z(1-LAP:nz+LAP), w1z(1-LAP:nz+LAP))
    allocate (A1x(1-LAP:nx+LAP), A2x(1-LAP:nx+LAP), A3x(1-LAP:nx+LAP), hhx1(0:nx), hhx2(0:nx))
    allocate (A1y(1-LAP:ny+LAP), A2y(1-LAP:ny+LAP), A3y(1-LAP:ny+LAP), hhy1(0:ny), hhy2(0:ny))
    allocate (A1z(1-LAP:nz+LAP), A2z(1-LAP:nz+LAP), A3z(1-LAP:nz+LAP), hhz1(0:nz), hhz2(0:nz))
    allocate (fpx(1-LAP:nx+LAP, 5), fmx(1-LAP:nx+LAP, 5))
    allocate (fpy(1-LAP:ny+LAP, 5), fmy(1-LAP:ny+LAP, 5))
    allocate (fpz(1-LAP:nz+LAP, 5), fmz(1-LAP:nz+LAP, 5))

    if (Para%runMode /= 0) then
        allocate (p2x(1-LAP:nx+LAP), d2x(1-LAP:nx+LAP), u2x(1-LAP:nx+LAP), v2x(1-LAP:nx+LAP), w2x(1-LAP:nx+LAP))
        allocate (p2y(1-LAP:ny+LAP), d2y(1-LAP:ny+LAP), u2y(1-LAP:ny+LAP), v2y(1-LAP:ny+LAP), w2y(1-LAP:ny+LAP))
        allocate (p2z(1-LAP:nz+LAP), d2z(1-LAP:nz+LAP), u2z(1-LAP:nz+LAP), v2z(1-LAP:nz+LAP), w2z(1-LAP:nz+LAP))
    end if

    !> Activation map of the positivity-preserving limiter (only used by IF_Character = 2)
    if (Para%runMode == 0 .and. Para%IF_Character == 2) then
        allocate (pp_mark_rho_x(1:nx, 1:ny, 1:nz), pp_mark_p_x(1:nx, 1:ny, 1:nz))
        allocate (pp_mark_rho_y(1:nx, 1:ny, 1:nz), pp_mark_p_y(1:nx, 1:ny, 1:nz))
        allocate (pp_mark_rho_z(1:nx, 1:ny, 1:nz), pp_mark_p_z(1:nx, 1:ny, 1:nz))
        pp_mark_rho_x = 0
        pp_mark_p_x = 0
        pp_mark_rho_y = 0
        pp_mark_p_y = 0
        pp_mark_rho_z = 0
        pp_mark_p_z = 0
    end if

end subroutine Allocate_invis_data

subroutine Deallocate_invis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    deallocate (hhx, hhy, hhz)
    deallocate (Hc1x, c1x, d1x, u1x, v1x, w1x)
    deallocate (Hc1y, c1y, d1y, u1y, v1y, w1y)
    deallocate (Hc1z, c1z, d1z, u1z, v1z, w1z)
    deallocate (A1x, A2x, A3x, hhx1, hhx2)
    deallocate (A1y, A2y, A3y, hhy1, hhy2)
    deallocate (A1z, A2z, A3z, hhz1, hhz2)
    deallocate (fpx, fmx)
    deallocate (fpy, fmy)
    deallocate (fpz, fmz)

    if (Para%runMode /= 0) then
        deallocate (p2x, d2x, u2x, v2x, w2x)
        deallocate (p2y, d2y, u2y, v2y, w2y)
        deallocate (p2z, d2z, u2z, v2z, w2z)
    end if

    if (allocated(pp_mark_rho_x)) then
        deallocate (pp_mark_rho_x, pp_mark_p_x)
        deallocate (pp_mark_rho_y, pp_mark_p_y)
        deallocate (pp_mark_rho_z, pp_mark_p_z)
    end if

end subroutine Deallocate_invis_data

subroutine Allocate_vis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    allocate (mu(1:nx, 1:ny, 1:nz))
    allocate (dmudk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dmudi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dmuds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dudk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dvdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dwdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dTdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dudi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dvdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dwdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dTdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (duds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dvds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dwds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
    allocate (dTds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))

    allocate (d2udkdi(1:nx, 1:ny, 1:nz), d2vdkdi(1:nx, 1:ny, 1:nz), d2wdkdi(1:nx, 1:ny, 1:nz), d2Tdkdi(1:nx, 1:ny, 1:nz))
    allocate (d2udkds(1:nx, 1:ny, 1:nz), d2vdkds(1:nx, 1:ny, 1:nz), d2wdkds(1:nx, 1:ny, 1:nz), d2Tdkds(1:nx, 1:ny, 1:nz))
    allocate (d2udids(1:nx, 1:ny, 1:nz), d2vdids(1:nx, 1:ny, 1:nz), d2wdids(1:nx, 1:ny, 1:nz), d2Tdids(1:nx, 1:ny, 1:nz))
    allocate (d2udk2(1:nx, 1:ny, 1:nz), d2vdk2(1:nx, 1:ny, 1:nz), d2wdk2(1:nx, 1:ny, 1:nz), d2Tdk2(1:nx, 1:ny, 1:nz))
    allocate (d2udi2(1:nx, 1:ny, 1:nz), d2vdi2(1:nx, 1:ny, 1:nz), d2wdi2(1:nx, 1:ny, 1:nz), d2Tdi2(1:nx, 1:ny, 1:nz))
    allocate (d2uds2(1:nx, 1:ny, 1:nz), d2vds2(1:nx, 1:ny, 1:nz), d2wds2(1:nx, 1:ny, 1:nz), d2Tds2(1:nx, 1:ny, 1:nz))

    if (Para%runMode /= 0) then
        allocate (mu_fluc(1:nx, 1:ny, 1:nz))
        allocate (dmufdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dmufdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dmufds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dufdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dvfdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dwfdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dTfdk(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dufdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dvfdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dwfdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dTfdi(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dufds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dvfds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dwfds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))
        allocate (dTfds(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP))

        allocate (d2ufdkdi(1:nx, 1:ny, 1:nz), d2vfdkdi(1:nx, 1:ny, 1:nz), d2wfdkdi(1:nx, 1:ny, 1:nz), d2Tfdkdi(1:nx, 1:ny, 1:nz))
        allocate (d2ufdkds(1:nx, 1:ny, 1:nz), d2vfdkds(1:nx, 1:ny, 1:nz), d2wfdkds(1:nx, 1:ny, 1:nz), d2Tfdkds(1:nx, 1:ny, 1:nz))
        allocate (d2ufdids(1:nx, 1:ny, 1:nz), d2vfdids(1:nx, 1:ny, 1:nz), d2wfdids(1:nx, 1:ny, 1:nz), d2Tfdids(1:nx, 1:ny, 1:nz))
        allocate (d2ufdk2(1:nx, 1:ny, 1:nz), d2vfdk2(1:nx, 1:ny, 1:nz), d2wfdk2(1:nx, 1:ny, 1:nz), d2Tfdk2(1:nx, 1:ny, 1:nz))
        allocate (d2ufdi2(1:nx, 1:ny, 1:nz), d2vfdi2(1:nx, 1:ny, 1:nz), d2wfdi2(1:nx, 1:ny, 1:nz), d2Tfdi2(1:nx, 1:ny, 1:nz))
        allocate (d2ufds2(1:nx, 1:ny, 1:nz), d2vfds2(1:nx, 1:ny, 1:nz), d2wfds2(1:nx, 1:ny, 1:nz), d2Tfds2(1:nx, 1:ny, 1:nz))
    end if

end subroutine Allocate_vis_data

subroutine Deallocate_vis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    deallocate (mu)
    deallocate (dmudk, dmudi, dmuds)
    deallocate (dudk, dvdk, dwdk, dTdk)
    deallocate (dudi, dvdi, dwdi, dTdi)
    deallocate (duds, dvds, dwds, dTds)
    deallocate (d2udkdi, d2vdkdi, d2wdkdi, d2Tdkdi)
    deallocate (d2udkds, d2vdkds, d2wdkds, d2Tdkds)
    deallocate (d2udids, d2vdids, d2wdids, d2Tdids)
    deallocate (d2udk2, d2vdk2, d2wdk2, d2Tdk2)
    deallocate (d2udi2, d2vdi2, d2wdi2, d2Tdi2)
    deallocate (d2uds2, d2vds2, d2wds2, d2Tds2)

    if (Para%runMode /= 0) then
        deallocate (mu_fluc)
        deallocate (dmufdk, dmufdi, dmufds)
        deallocate (dufdk, dvfdk, dwfdk, dTfdk)
        deallocate (dufdi, dvfdi, dwfdi, dTfdi)
        deallocate (dufds, dvfds, dwfds, dTfds)
        deallocate (d2ufdkdi, d2vfdkdi, d2wfdkdi, d2Tfdkdi)
        deallocate (d2ufdkds, d2vfdkds, d2wfdkds, d2Tfdkds)
        deallocate (d2ufdids, d2vfdids, d2wfdids, d2Tfdids)
        deallocate (d2ufdk2, d2vfdk2, d2wfdk2, d2Tfdk2)
        deallocate (d2ufdi2, d2vfdi2, d2wfdi2, d2Tfdi2)
        deallocate (d2ufds2, d2vfds2, d2wfds2, d2Tfds2)
    end if

end subroutine Deallocate_vis_data
