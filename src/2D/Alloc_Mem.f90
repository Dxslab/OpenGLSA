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

    allocate (Axx(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Ayy(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Akx(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Aix(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Aky(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Aiy(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Ajac(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Akx1(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Aix1(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Aky1(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Aiy1(1-LAP:nx+LAP, 1-LAP:ny+LAP))

    allocate (BcI(ny, 2)); BcI = 0
    allocate (BcJ(nx, 2)); BcJ = 0
    if (Para%Iperiodic_X == 0) then
        if (npx == 0 .and. Para%Isymmetric_X(1) == 0) BcI(:, 1) = 1
        if (npx == npx0-1 .and. Para%Isymmetric_X(2) == 0) BcI(:, 2) = 1
    end if
    if (Para%Iperiodic_Y == 0) then
        if (npy == 0 .and. Para%Isymmetric_Y(1) == 0) BcJ(:, 1) = 1
        if (npy == npy0-1 .and. Para%Isymmetric_Y(2) == 0) BcJ(:, 2) = 1
    end if

    if (Para%IF_Viscous == 1) then
        allocate (Akxk(1:nx, 1:ny))
        allocate (Akxi(1:nx, 1:ny))
        allocate (Aixk(1:nx, 1:ny))
        allocate (Aixi(1:nx, 1:ny))
        allocate (Akyk(1:nx, 1:ny))
        allocate (Akyi(1:nx, 1:ny))
        allocate (Aiyk(1:nx, 1:ny))
        allocate (Aiyi(1:nx, 1:ny))
        if (Para%IF_Cylindrical_Coordinate == 1) then
            allocate (dydk(1:nx, 1:ny))
            allocate (dydi(1:nx, 1:ny))
        end if
    end if
end subroutine Allocate_mesh_data

subroutine Deallocate_mesh_data
    use mesh_data
    implicit none

    deallocate (Axx)
    deallocate (Ayy)
    deallocate (Akx)
    deallocate (Aix)
    deallocate (Aky)
    deallocate (Aiy)
    deallocate (Ajac)
    deallocate (Akx1)
    deallocate (Aix1)
    deallocate (Aky1)
    deallocate (Aiy1)

    deallocate (BcI)
    deallocate (BcJ)

    if (Para%IF_Viscous == 1) then
        deallocate (Akxk)
        deallocate (Akxi)
        deallocate (Aixk)
        deallocate (Aixi)
        deallocate (Akyk)
        deallocate (Akyi)
        deallocate (Aiyk)
        deallocate (Aiyi)
        if (Para%IF_Cylindrical_Coordinate == 1) then
            deallocate (dydk)
            deallocate (dydi)
        end if
    end if
end subroutine Deallocate_mesh_data

subroutine Allocate_flow_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    allocate (d(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (u(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (v(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (T(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (p(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (c(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (Hc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (rhoe(1-LAP:nx+LAP, 1-LAP:ny+LAP))

    if (Para%runMode == 0) then
        allocate (f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 4))
        allocate (fn(1-LAP:nx+LAP, 1-LAP:ny+LAP, 4))
        allocate (du(1:nx, 1:ny, 4))
        if (Para%IF_SFD == 1) then
            allocate (fsfd(1-LAP:nx+LAP, 1-LAP:ny+LAP, 4))
            fsfd = 0.0d0
        end if
    else
        allocate (d_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (u_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (v_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (w_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (T_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (p_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (rhoe_fluc(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 5))
        allocate (fn(1-LAP:nx+LAP, 1-LAP:ny+LAP, 5))
        allocate (du(1:nx, 1:ny, 5))
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

    allocate (hhx(0:nx, 4))
    allocate (hhy(0:ny, 4))
    allocate (Hc1x(1-LAP:nx+LAP), c1x(1-LAP:nx+LAP), d1x(1-LAP:nx+LAP), u1x(1-LAP:nx+LAP), v1x(1-LAP:nx+LAP))
    allocate (Hc1y(1-LAP:ny+LAP), c1y(1-LAP:ny+LAP), d1y(1-LAP:ny+LAP), u1y(1-LAP:ny+LAP), v1y(1-LAP:ny+LAP))
    allocate (A1x(1-LAP:nx+LAP), A2x(1-LAP:nx+LAP), hhx1(0:nx), hhx2(0:nx))
    allocate (A1y(1-LAP:ny+LAP), A2y(1-LAP:ny+LAP), hhy1(0:ny), hhy2(0:ny))

    if (Para%runMode == 0) then
        allocate (fpx(1-LAP:nx+LAP, 4), fmx(1-LAP:nx+LAP, 4))
        allocate (fpy(1-LAP:ny+LAP, 4), fmy(1-LAP:ny+LAP, 4))
        allocate (Source_Invis(1:nx, 1:ny, 4))
    else
        allocate (fpx(1-LAP:nx+LAP, 5), fmx(1-LAP:nx+LAP, 5))
        allocate (fpy(1-LAP:ny+LAP, 5), fmy(1-LAP:ny+LAP, 5))
        allocate (p2x(1-LAP:nx+LAP), d2x(1-LAP:nx+LAP), u2x(1-LAP:nx+LAP), v2x(1-LAP:nx+LAP), w2x(1-LAP:nx+LAP))
        allocate (p2y(1-LAP:ny+LAP), d2y(1-LAP:ny+LAP), u2y(1-LAP:ny+LAP), v2y(1-LAP:ny+LAP), w2y(1-LAP:ny+LAP))
        allocate (Source_Invis(1:nx, 1:ny, 5))
    end if

    Source_Invis = 0.0d0

    if (Para%runMode == 0 .and. Para%IF_Character == 2) then
        allocate (pp_mark_rho_x(1:nx, 1:ny), pp_mark_p_x(1:nx, 1:ny))
        allocate (pp_mark_rho_y(1:nx, 1:ny), pp_mark_p_y(1:nx, 1:ny))
        pp_mark_rho_x = 0
        pp_mark_p_x = 0
        pp_mark_rho_y = 0
        pp_mark_p_y = 0
    end if

end subroutine Allocate_invis_data

subroutine Deallocate_invis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    deallocate (hhx, hhy)
    deallocate (Hc1x, c1x, d1x, u1x, v1x)
    deallocate (Hc1y, c1y, d1y, u1y, v1y)
    deallocate (A1x, A2x, hhx1, hhx2)
    deallocate (A1y, A2y, hhy1, hhy2)
    deallocate (fpx, fmx)
    deallocate (fpy, fmy)
    deallocate (Source_Invis)

    if (Para%runMode /= 0) then
        deallocate (p2x, d2x, u2x, v2x, w2x)
        deallocate (p2y, d2y, u2y, v2y, w2y)
    end if

    if (allocated(pp_mark_rho_x)) then
        deallocate (pp_mark_rho_x, pp_mark_p_x)
        deallocate (pp_mark_rho_y, pp_mark_p_y)
    end if

end subroutine Deallocate_invis_data

subroutine Allocate_vis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    allocate (mu(1:nx, 1:ny))
    allocate (dmudk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dmudi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dudk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dvdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dTdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dudi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dvdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (dTdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
    allocate (d2udkdi(1:nx, 1:ny), d2vdkdi(1:nx, 1:ny), d2Tdkdi(1:nx, 1:ny))
    allocate (d2udk2(1:nx, 1:ny), d2vdk2(1:nx, 1:ny), d2Tdk2(1:nx, 1:ny))
    allocate (d2udi2(1:nx, 1:ny), d2vdi2(1:nx, 1:ny), d2Tdi2(1:nx, 1:ny))


    if (Para%runMode == 0) then
        allocate (Source_Vis(1:nx, 1:ny, 4))
    else
        allocate (mu_fluc(1:nx, 1:ny))
        allocate (dmufdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dmufdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dufdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dvfdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dwfdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dTfdk(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dufdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dvfdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dwfdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (dTfdi(1-LAP:nx+LAP, 1-LAP:ny+LAP))
        allocate (d2ufdkdi(1:nx, 1:ny), d2vfdkdi(1:nx, 1:ny), d2wfdkdi(1:nx, 1:ny), d2Tfdkdi(1:nx, 1:ny))
        allocate (d2ufdk2(1:nx, 1:ny), d2vfdk2(1:nx, 1:ny), d2wfdk2(1:nx, 1:ny), d2Tfdk2(1:nx, 1:ny))
        allocate (d2ufdi2(1:nx, 1:ny), d2vfdi2(1:nx, 1:ny), d2wfdi2(1:nx, 1:ny), d2Tfdi2(1:nx, 1:ny))
        allocate (Source_Vis(1:nx, 1:ny, 5))
    end if

    Source_Vis = 0.0d0

end subroutine Allocate_vis_data

subroutine Deallocate_vis_data
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    deallocate (mu)
    deallocate (dmudk, dmudi)
    deallocate (dudk, dvdk, dTdk)
    deallocate (dudi, dvdi, dTdi)
    deallocate (d2udkdi, d2vdkdi, d2Tdkdi)
    deallocate (d2udk2, d2vdk2, d2Tdk2)
    deallocate (d2udi2, d2vdi2, d2Tdi2)
    deallocate (Source_Vis)

    if (Para%runMode /= 0) then
        deallocate (mu_fluc)
        deallocate (dmufdk, dmufdi)
        deallocate (dufdk, dvfdk, dwfdk, dTfdk)
        deallocate (dufdi, dvfdi, dwfdi, dTfdi)
        deallocate (d2ufdkdi, d2vfdkdi, d2wfdkdi, d2Tfdkdi)
        deallocate (d2ufdk2, d2vfdk2, d2wfdk2, d2Tfdk2)
        deallocate (d2ufdi2, d2vfdi2, d2wfdi2, d2Tfdi2)
    end if

end subroutine Deallocate_vis_data
