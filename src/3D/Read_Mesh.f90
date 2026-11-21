! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine read_mesh
    use mesh_data
    use IO
    use Centrals
    implicit none

    integer :: i, j, k

    if (my_id == 0) open (66, FILE="grid3d.dat", FORM="UNFORMATTED")

    call read_data(66, Axx)
    call read_data(66, Ayy)
    call read_data(66, Azz)

    if (my_id == 0) close (66)

    call mpi_exchange_xyz3d(Axx)
    call mpi_exchange_xyz3d(Ayy)
    call mpi_exchange_xyz3d(Azz)
    call update_xyz_periodic
    call update_xyz_symmetric

    call compute_Jocabian
    call mpi_exchange_xyz3d(Akx)
    call mpi_exchange_xyz3d(Aky)
    call mpi_exchange_xyz3d(Akz)
    call mpi_exchange_xyz3d(Aix)
    call mpi_exchange_xyz3d(Aiy)
    call mpi_exchange_xyz3d(Aiz)
    call mpi_exchange_xyz3d(Asx)
    call mpi_exchange_xyz3d(Asy)
    call mpi_exchange_xyz3d(Asz)
    call mpi_exchange_xyz3d(Ajac)
    call update_metrics_symmetric

    if (Para%IF_Viscous == 1) then
        !
        call OCFD_dx0(Akx, Akxk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Aky, Akyk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Akz, Akzk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Aix, Aixk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Aiy, Aiyk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Aiz, Aizk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Asx, Asxk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Asy, Asyk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(Asz, Aszk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        !
        call OCFD_dy0(Akx, Akxi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Aky, Akyi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Akz, Akzi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Aix, Aixi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Aiy, Aiyi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Aiz, Aizi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Asx, Asxi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Asy, Asyi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(Asz, Aszi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        !
        call OCFD_dz0(Akx, Akxs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Aky, Akys, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Akz, Akzs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Aix, Aixs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Aiy, Aiys, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Aiz, Aizs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Asx, Asxs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Asy, Asys, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(Asz, Aszs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        !
    end if

    do k = 1-LAP, nz+LAP
    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        Akx1(i, j, k) = Akx(i, j, k)/Ajac(i, j, k)
        Aky1(i, j, k) = Aky(i, j, k)/Ajac(i, j, k)
        Akz1(i, j, k) = Akz(i, j, k)/Ajac(i, j, k)
        Aix1(i, j, k) = Aix(i, j, k)/Ajac(i, j, k)
        Aiy1(i, j, k) = Aiy(i, j, k)/Ajac(i, j, k)
        Aiz1(i, j, k) = Aiz(i, j, k)/Ajac(i, j, k)
        Asx1(i, j, k) = Asx(i, j, k)/Ajac(i, j, k)
        Asy1(i, j, k) = Asy(i, j, k)/Ajac(i, j, k)
        Asz1(i, j, k) = Asz(i, j, k)/Ajac(i, j, k)
    end do
    end do
    end do

    if (my_id == 0) print *, "Read grid OK"

end subroutine read_mesh

subroutine update_xyz_periodic
    use mesh_data
    use flow_para
    implicit none

    integer :: i, j, k

    if (Para%Iperiodic_X == 1 .and. npx == 0) then
        do k = 1, nz
        do j = 1, ny
        do i = 1-LAP, 0
            Axx(i, j, k) = Axx(i, j, k)-Para%Periodic_ISpan(1)
            Ayy(i, j, k) = Ayy(i, j, k)-Para%Periodic_ISpan(2)
            Azz(i, j, k) = Azz(i, j, k)-Para%Periodic_ISpan(3)
        end do
        end do
        end do
    end if

    if (Para%Iperiodic_X == 1 .and. npx == npx0-1) then
        do k = 1, nz
        do j = 1, ny
        do i = nx+1, nx+LAP
            Axx(i, j, k) = Axx(i, j, k)+Para%Periodic_ISpan(1)
            Ayy(i, j, k) = Ayy(i, j, k)+Para%Periodic_ISpan(2)
            Azz(i, j, k) = Azz(i, j, k)+Para%Periodic_ISpan(3)
        end do
        end do
        end do
    end if

    if (Para%Iperiodic_Y == 1 .and. npy == 0) then
        do k = 1, nz
        do j = 1-LAP, 0
        do i = 1, nx
            Axx(i, j, k) = Axx(i, j, k)-Para%Periodic_JSpan(1)
            Ayy(i, j, k) = Ayy(i, j, k)-Para%Periodic_JSpan(2)
            Azz(i, j, k) = Azz(i, j, k)-Para%Periodic_JSpan(3)
        end do
        end do
        end do
    end if

    if (Para%Iperiodic_Y == 1 .and. npy == npy0-1) then
        do k = 1, nz
        do j = ny+1, ny+LAP
        do i = 1, nx
            Axx(i, j, k) = Axx(i, j, k)+Para%Periodic_JSpan(1)
            Ayy(i, j, k) = Ayy(i, j, k)+Para%Periodic_JSpan(2)
            Azz(i, j, k) = Azz(i, j, k)+Para%Periodic_JSpan(3)
        end do
        end do
        end do
    end if

    if (Para%Iperiodic_Z == 1 .and. npz == 0) then
        do k = 1-LAP, 0
        do j = 1, ny
        do i = 1, nx
            Axx(i, j, k) = Axx(i, j, k)-Para%Periodic_KSpan(1)
            Ayy(i, j, k) = Ayy(i, j, k)-Para%Periodic_KSpan(2)
            Azz(i, j, k) = Azz(i, j, k)-Para%Periodic_KSpan(3)
        end do
        end do
        end do
    end if

    if (Para%Iperiodic_Z == 1 .and. npz == npz0-1) then
        do k = nz+1, nz+LAP
        do j = 1, ny
        do i = 1, nx
            Axx(i, j, k) = Axx(i, j, k)+Para%Periodic_KSpan(1)
            Ayy(i, j, k) = Ayy(i, j, k)+Para%Periodic_KSpan(2)
            Azz(i, j, k) = Azz(i, j, k)+Para%Periodic_KSpan(3)
        end do
        end do
        end do
    end if

end subroutine update_xyz_periodic

subroutine update_xyz_symmetric
    use mesh_data
    use flow_para
    use precision
    implicit none

    integer :: i, j, k, ii, jj, kk
    real(dp) :: normal_x, normal_y, normal_z, vx, vy, vz, dot2

    if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
        do k = 1, nz
        do j = 1, ny
            call get_symmetric_normal(1, 1, j, k, normal_x, normal_y, normal_z)
            do i = 1-LAP, 0
                ii = 2-i
                vx = Axx(ii, j, k)-Axx(1, j, k)
                vy = Ayy(ii, j, k)-Ayy(1, j, k)
                vz = Azz(ii, j, k)-Azz(1, j, k)
                dot2 = 2.0d0*(vx*normal_x+vy*normal_y+vz*normal_z)
                Axx(i, j, k) = Axx(1, j, k)+vx-dot2*normal_x
                Ayy(i, j, k) = Ayy(1, j, k)+vy-dot2*normal_y
                Azz(i, j, k) = Azz(1, j, k)+vz-dot2*normal_z
            end do
        end do
        end do
    end if

    if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
        do k = 1, nz
        do j = 1, ny
            call get_symmetric_normal(1, nx, j, k, normal_x, normal_y, normal_z)
            do i = nx+1, nx+LAP
                ii = 2*nx-i
                vx = Axx(ii, j, k)-Axx(nx, j, k)
                vy = Ayy(ii, j, k)-Ayy(nx, j, k)
                vz = Azz(ii, j, k)-Azz(nx, j, k)
                dot2 = 2.0d0*(vx*normal_x+vy*normal_y+vz*normal_z)
                Axx(i, j, k) = Axx(nx, j, k)+vx-dot2*normal_x
                Ayy(i, j, k) = Ayy(nx, j, k)+vy-dot2*normal_y
                Azz(i, j, k) = Azz(nx, j, k)+vz-dot2*normal_z
            end do
        end do
        end do
    end if

    if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
        do k = 1, nz
        do j = 1-LAP, 0
        do i = 1, nx
            jj = 2-j
            call get_symmetric_normal(2, i, 1, k, normal_x, normal_y, normal_z)
            vx = Axx(i, jj, k)-Axx(i, 1, k)
            vy = Ayy(i, jj, k)-Ayy(i, 1, k)
            vz = Azz(i, jj, k)-Azz(i, 1, k)
            dot2 = 2.0d0*(vx*normal_x+vy*normal_y+vz*normal_z)
            Axx(i, j, k) = Axx(i, 1, k)+vx-dot2*normal_x
            Ayy(i, j, k) = Ayy(i, 1, k)+vy-dot2*normal_y
            Azz(i, j, k) = Azz(i, 1, k)+vz-dot2*normal_z
        end do
        end do
        end do
    end if

    if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
        do k = 1, nz
        do j = ny+1, ny+LAP
        do i = 1, nx
            jj = 2*ny-j
            call get_symmetric_normal(2, i, ny, k, normal_x, normal_y, normal_z)
            vx = Axx(i, jj, k)-Axx(i, ny, k)
            vy = Ayy(i, jj, k)-Ayy(i, ny, k)
            vz = Azz(i, jj, k)-Azz(i, ny, k)
            dot2 = 2.0d0*(vx*normal_x+vy*normal_y+vz*normal_z)
            Axx(i, j, k) = Axx(i, ny, k)+vx-dot2*normal_x
            Ayy(i, j, k) = Ayy(i, ny, k)+vy-dot2*normal_y
            Azz(i, j, k) = Azz(i, ny, k)+vz-dot2*normal_z
        end do
        end do
        end do
    end if

    if (Para%Isymmetric_Z(1) == 1 .and. npz == 0) then
        do k = 1-LAP, 0
        do j = 1, ny
        do i = 1, nx
            kk = 2-k
            call get_symmetric_normal(3, i, j, 1, normal_x, normal_y, normal_z)
            vx = Axx(i, j, kk)-Axx(i, j, 1)
            vy = Ayy(i, j, kk)-Ayy(i, j, 1)
            vz = Azz(i, j, kk)-Azz(i, j, 1)
            dot2 = 2.0d0*(vx*normal_x+vy*normal_y+vz*normal_z)
            Axx(i, j, k) = Axx(i, j, 1)+vx-dot2*normal_x
            Ayy(i, j, k) = Ayy(i, j, 1)+vy-dot2*normal_y
            Azz(i, j, k) = Azz(i, j, 1)+vz-dot2*normal_z
        end do
        end do
        end do
    end if

    if (Para%Isymmetric_Z(2) == 1 .and. npz == npz0-1) then
        do k = nz+1, nz+LAP
        do j = 1, ny
        do i = 1, nx
            kk = 2*nz-k
            call get_symmetric_normal(3, i, j, nz, normal_x, normal_y, normal_z)
            vx = Axx(i, j, kk)-Axx(i, j, nz)
            vy = Ayy(i, j, kk)-Ayy(i, j, nz)
            vz = Azz(i, j, kk)-Azz(i, j, nz)
            dot2 = 2.0d0*(vx*normal_x+vy*normal_y+vz*normal_z)
            Axx(i, j, k) = Axx(i, j, nz)+vx-dot2*normal_x
            Ayy(i, j, k) = Ayy(i, j, nz)+vy-dot2*normal_y
            Azz(i, j, k) = Azz(i, j, nz)+vz-dot2*normal_z
        end do
        end do
        end do
    end if

end subroutine update_xyz_symmetric

subroutine update_metrics_symmetric
    use mesh_data
    use flow_para
    use precision
    implicit none

    integer :: i, j, k, ii, jj, kk
    real(dp) :: normal_x, normal_y, normal_z

    if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
        do k = 1, nz
        do j = 1, ny
            call get_symmetric_normal(1, 1, j, k, normal_x, normal_y, normal_z)
            do i = 1-LAP, 0
                ii = 2-i
                call set_reflected_metrics(i, j, k, ii, j, k, normal_x, normal_y, normal_z, -1.0d0, 1.0d0, 1.0d0)
            end do
        end do
        end do
    end if

    if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
        do k = 1, nz
        do j = 1, ny
            call get_symmetric_normal(1, nx, j, k, normal_x, normal_y, normal_z)
            do i = nx+1, nx+LAP
                ii = 2*nx-i
                call set_reflected_metrics(i, j, k, ii, j, k, normal_x, normal_y, normal_z, -1.0d0, 1.0d0, 1.0d0)
            end do
        end do
        end do
    end if

    if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
        do k = 1, nz
        do j = 1-LAP, 0
        do i = 1, nx
            jj = 2-j
            call get_symmetric_normal(2, i, 1, k, normal_x, normal_y, normal_z)
            call set_reflected_metrics(i, j, k, i, jj, k, normal_x, normal_y, normal_z, 1.0d0, -1.0d0, 1.0d0)
        end do
        end do
        end do
    end if

    if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
        do k = 1, nz
        do j = ny+1, ny+LAP
        do i = 1, nx
            jj = 2*ny-j
            call get_symmetric_normal(2, i, ny, k, normal_x, normal_y, normal_z)
            call set_reflected_metrics(i, j, k, i, jj, k, normal_x, normal_y, normal_z, 1.0d0, -1.0d0, 1.0d0)
        end do
        end do
        end do
    end if

    if (Para%Isymmetric_Z(1) == 1 .and. npz == 0) then
        do k = 1-LAP, 0
        do j = 1, ny
        do i = 1, nx
            kk = 2-k
            call get_symmetric_normal(3, i, j, 1, normal_x, normal_y, normal_z)
            call set_reflected_metrics(i, j, k, i, j, kk, normal_x, normal_y, normal_z, 1.0d0, 1.0d0, -1.0d0)
        end do
        end do
        end do
    end if

    if (Para%Isymmetric_Z(2) == 1 .and. npz == npz0-1) then
        do k = nz+1, nz+LAP
        do j = 1, ny
        do i = 1, nx
            kk = 2*nz-k
            call get_symmetric_normal(3, i, j, nz, normal_x, normal_y, normal_z)
            call set_reflected_metrics(i, j, k, i, j, kk, normal_x, normal_y, normal_z, 1.0d0, 1.0d0, -1.0d0)
        end do
        end do
        end do
    end if

contains

    subroutine set_reflected_metrics(i0, j0, k0, is, js, ks, nnormal_x, nnormal_y, nnormal_z, scale_k, scale_i, scale_s)
        implicit none

        integer, intent(IN) :: i0, j0, k0, is, js, ks
        real(dp), intent(IN) :: nnormal_x, nnormal_y, nnormal_z, scale_k, scale_i, scale_s

        call reflect_metric_vector(Akx(is, js, ks), Aky(is, js, ks), Akz(is, js, ks), &
                                   nnormal_x, nnormal_y, nnormal_z, scale_k, Akx(i0, j0, k0), Aky(i0, j0, k0), Akz(i0, j0, k0))
        call reflect_metric_vector(Aix(is, js, ks), Aiy(is, js, ks), Aiz(is, js, ks), &
                                   nnormal_x, nnormal_y, nnormal_z, scale_i, Aix(i0, j0, k0), Aiy(i0, j0, k0), Aiz(i0, j0, k0))
        call reflect_metric_vector(Asx(is, js, ks), Asy(is, js, ks), Asz(is, js, ks), &
                                   nnormal_x, nnormal_y, nnormal_z, scale_s, Asx(i0, j0, k0), Asy(i0, j0, k0), Asz(i0, j0, k0))

        Ajac(i0, j0, k0) = Ajac(is, js, ks)
    end subroutine set_reflected_metrics

    subroutine reflect_metric_vector(ax, ay, az, nnormal_x, nnormal_y, nnormal_z, scale, bx, by, bz)
        implicit none

        real(dp), intent(IN) :: ax, ay, az, nnormal_x, nnormal_y, nnormal_z, scale
        real(dp), intent(OUT) :: bx, by, bz
        real(dp) :: dot2

        dot2 = 2.0d0*(ax*nnormal_x+ay*nnormal_y+az*nnormal_z)
        bx = scale*(ax-dot2*nnormal_x)
        by = scale*(ay-dot2*nnormal_y)
        bz = scale*(az-dot2*nnormal_z)
    end subroutine reflect_metric_vector

end subroutine update_metrics_symmetric

subroutine get_symmetric_normal(face, ib, jb, kb, normal_x, normal_y, normal_z)
    use mesh_data
    use precision
    implicit none

    integer, intent(IN) :: face, ib, jb, kb
    real(dp), intent(OUT) :: normal_x, normal_y, normal_z
    integer :: IERR
    real(dp) :: t1x, t1y, t1z, t2x, t2y, t2z, norm_n
    real(dp), parameter :: eps_n = 1.0d-14

    select case (face)
    case (1)
        call tangent_j(ib, jb, kb, t1x, t1y, t1z)
        call tangent_k(ib, jb, kb, t2x, t2y, t2z)
    case (2)
        call tangent_k(ib, jb, kb, t1x, t1y, t1z)
        call tangent_i(ib, jb, kb, t2x, t2y, t2z)
    case (3)
        call tangent_i(ib, jb, kb, t1x, t1y, t1z)
        call tangent_j(ib, jb, kb, t2x, t2y, t2z)
    case DEFAULT
        print *, "*** Error: invalid symmetric face in get_symmetric_normal ***"
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select

    normal_x = t1y*t2z-t1z*t2y
    normal_y = t1z*t2x-t1x*t2z
    normal_z = t1x*t2y-t1y*t2x
    norm_n = sqrt(normal_x*normal_x+normal_y*normal_y+normal_z*normal_z)

    if (norm_n > eps_n) then
        normal_x = normal_x/norm_n
        normal_y = normal_y/norm_n
        normal_z = normal_z/norm_n
    else
        print *, "*** Error: zero surface normal in symmetric mesh update ***"
        print *, "Index:", ib, jb, kb
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end if

contains

    subroutine tangent_i(i0, j0, k0, tx, ty, tz)
        implicit none

        integer, intent(IN) :: i0, j0, k0
        real(dp), intent(OUT) :: tx, ty, tz

        if (nx < 2) then
            print *, "*** Error: nx must be greater than 1 for symmetric mesh update ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        else if (i0 == 1) then
            tx = Axx(2, j0, k0)-Axx(1, j0, k0)
            ty = Ayy(2, j0, k0)-Ayy(1, j0, k0)
            tz = Azz(2, j0, k0)-Azz(1, j0, k0)
        else if (i0 == nx) then
            tx = Axx(nx, j0, k0)-Axx(nx-1, j0, k0)
            ty = Ayy(nx, j0, k0)-Ayy(nx-1, j0, k0)
            tz = Azz(nx, j0, k0)-Azz(nx-1, j0, k0)
        else
            tx = 0.5d0*(Axx(i0+1, j0, k0)-Axx(i0-1, j0, k0))
            ty = 0.5d0*(Ayy(i0+1, j0, k0)-Ayy(i0-1, j0, k0))
            tz = 0.5d0*(Azz(i0+1, j0, k0)-Azz(i0-1, j0, k0))
        end if
    end subroutine tangent_i

    subroutine tangent_j(i0, j0, k0, tx, ty, tz)
        implicit none

        integer, intent(IN) :: i0, j0, k0
        real(dp), intent(OUT) :: tx, ty, tz

        if (ny < 2) then
            print *, "*** Error: ny must be greater than 1 for symmetric mesh update ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        else if (j0 == 1) then
            tx = Axx(i0, 2, k0)-Axx(i0, 1, k0)
            ty = Ayy(i0, 2, k0)-Ayy(i0, 1, k0)
            tz = Azz(i0, 2, k0)-Azz(i0, 1, k0)
        else if (j0 == ny) then
            tx = Axx(i0, ny, k0)-Axx(i0, ny-1, k0)
            ty = Ayy(i0, ny, k0)-Ayy(i0, ny-1, k0)
            tz = Azz(i0, ny, k0)-Azz(i0, ny-1, k0)
        else
            tx = 0.5d0*(Axx(i0, j0+1, k0)-Axx(i0, j0-1, k0))
            ty = 0.5d0*(Ayy(i0, j0+1, k0)-Ayy(i0, j0-1, k0))
            tz = 0.5d0*(Azz(i0, j0+1, k0)-Azz(i0, j0-1, k0))
        end if
    end subroutine tangent_j

    subroutine tangent_k(i0, j0, k0, tx, ty, tz)
        implicit none

        integer, intent(IN) :: i0, j0, k0
        real(dp), intent(OUT) :: tx, ty, tz

        if (nz < 2) then
            print *, "*** Error: nz must be greater than 1 for symmetric mesh update ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        else if (k0 == 1) then
            tx = Axx(i0, j0, 2)-Axx(i0, j0, 1)
            ty = Ayy(i0, j0, 2)-Ayy(i0, j0, 1)
            tz = Azz(i0, j0, 2)-Azz(i0, j0, 1)
        else if (k0 == nz) then
            tx = Axx(i0, j0, nz)-Axx(i0, j0, nz-1)
            ty = Ayy(i0, j0, nz)-Ayy(i0, j0, nz-1)
            tz = Azz(i0, j0, nz)-Azz(i0, j0, nz-1)
        else
            tx = 0.5d0*(Axx(i0, j0, k0+1)-Axx(i0, j0, k0-1))
            ty = 0.5d0*(Ayy(i0, j0, k0+1)-Ayy(i0, j0, k0-1))
            tz = 0.5d0*(Azz(i0, j0, k0+1)-Azz(i0, j0, k0-1))
        end if
    end subroutine tangent_k

end subroutine get_symmetric_normal

subroutine compute_Jocabian
    use mesh_data
    use Centrals
    implicit none

    integer :: i, j, k, IERR
    real(dp), allocatable, dimension(:, :, :) :: xk, xi, xs, yk, yi, ys, zk, zi, zs
    real(dp) :: xi1, xj1, xk1, yi1, yj1, yk1, zi1, zj1, zk1, Jac1

    allocate (xk(nx, ny, nz), xi(nx, ny, nz), xs(nx, ny, nz))
    allocate (yk(nx, ny, nz), yi(nx, ny, nz), ys(nx, ny, nz))
    allocate (zk(nx, ny, nz), zi(nx, ny, nz), zs(nx, ny, nz))

    call OCFD_dx0(Axx, xk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
    call OCFD_dx0(Ayy, yk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
    call OCFD_dx0(Azz, zk, nx, ny, nz, Scheme%Scheme_Vis, BcI)
    call OCFD_dy0(Axx, xi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
    call OCFD_dy0(Ayy, yi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
    call OCFD_dy0(Azz, zi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
    call OCFD_dz0(Axx, xs, nx, ny, nz, Scheme%Scheme_Vis, BcK)
    call OCFD_dz0(Ayy, ys, nx, ny, nz, Scheme%Scheme_Vis, BcK)
    call OCFD_dz0(Azz, zs, nx, ny, nz, Scheme%Scheme_Vis, BcK)

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
        xi1 = xk(i, j, k); xj1 = xi(i, j, k); xk1 = xs(i, j, k)
        yi1 = yk(i, j, k); yj1 = yi(i, j, k); yk1 = ys(i, j, k)
        zi1 = zk(i, j, k); zj1 = zi(i, j, k); zk1 = zs(i, j, k)
        Jac1 = 1.d0/(xi1*yj1*zk1+yi1*zj1*xk1+zi1*xj1*yk1-zi1*yj1*xk1-yi1*xj1*zk1-xi1*zj1*yk1)   ! 1./Jocabian = d(x,y,z)/d(i,j,k)
        Ajac(i, j, k) = Jac1
        Akx(i, j, k) = Jac1*(yj1*zk1-zj1*yk1)
        Aky(i, j, k) = Jac1*(zj1*xk1-xj1*zk1)
        Akz(i, j, k) = Jac1*(xj1*yk1-yj1*xk1)
        Aix(i, j, k) = Jac1*(yk1*zi1-zk1*yi1)
        Aiy(i, j, k) = Jac1*(zk1*xi1-xk1*zi1)
        Aiz(i, j, k) = Jac1*(xk1*yi1-yk1*xi1)
        Asx(i, j, k) = Jac1*(yi1*zj1-zi1*yj1)
        Asy(i, j, k) = Jac1*(zi1*xj1-xi1*zj1)
        Asz(i, j, k) = Jac1*(xi1*yj1-yi1*xj1)

        if (Jac1 <= 0.0d0) then
            print *, "*** Error: negative Jacobian ***"
            print *, "Please Check Grid Orientation or Dimensions! "
            print *, i_offset(npx)+i-1, j_offset(npy)+j-1, k_offset(npz)+k-1
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        end if
    end do
    end do
    end do

    deallocate (xk, xi, xs, yk, yi, ys, zk, zi, zs)
    if (my_id == 0) print *, "Compute Jacobian OK"

end subroutine compute_Jocabian
