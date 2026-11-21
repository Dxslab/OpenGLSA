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

    integer :: i, j

    if (my_id == 0) open (66, FILE="grid.dat", FORM="UNFORMATTED")

    call read_data(66, Axx)
    call read_data(66, Ayy)

    if (my_id == 0) close (66)

    call mpi_exchange_xy2d(Axx)
    call mpi_exchange_xy2d(Ayy)
    call update_xy_periodic
    call update_xy_symmetric

    call compute_Jocabian
    call mpi_exchange_xy2d(Akx)
    call mpi_exchange_xy2d(Aky)
    call mpi_exchange_xy2d(Aix)
    call mpi_exchange_xy2d(Aiy)
    call mpi_exchange_xy2d(Ajac)
    call update_metrics_symmetric

    if (Para%IF_Viscous == 1) then
        !
        call OCFD2d_dx0(Akx, Akxk, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(Aky, Akyk, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(Aix, Aixk, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(Aiy, Aiyk, nx, ny, Scheme%Scheme_Vis, BcI)
        !
        call OCFD2d_dy0(Akx, Akxi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(Aky, Akyi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(Aix, Aixi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(Aiy, Aiyi, nx, ny, Scheme%Scheme_Vis, BcJ)
        !
        if (Para%IF_Cylindrical_Coordinate == 1) then
            call OCFD2d_dx0(Ayy, dydk, nx, ny, Scheme%Scheme_Vis, BcI)
            call OCFD2d_dy0(Ayy, dydi, nx, ny, Scheme%Scheme_Vis, BcJ)
        end if
    end if

    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        Akx1(i, j) = Akx(i, j)/Ajac(i, j)
        Aky1(i, j) = Aky(i, j)/Ajac(i, j)
        Aix1(i, j) = Aix(i, j)/Ajac(i, j)
        Aiy1(i, j) = Aiy(i, j)/Ajac(i, j)
    end do
    end do

    if (my_id == 0) print *, "Read grid OK"

end subroutine read_mesh

subroutine update_xy_periodic
    use mesh_data
    use flow_para
    implicit none

    integer :: i, j

    if (Para%Iperiodic_X == 1 .and. npx == 0) then
        do j = 1, ny
        do i = 1-LAP, 0
            Axx(i, j) = Axx(i, j)-Para%Periodic_ISpan(1)
            Ayy(i, j) = Ayy(i, j)-Para%Periodic_ISpan(2)
        end do
        end do
    end if

    if (Para%Iperiodic_X == 1 .and. npx == npx0-1) then
        do j = 1, ny
        do i = nx+1, nx+LAP
            Axx(i, j) = Axx(i, j)+Para%Periodic_ISpan(1)
            Ayy(i, j) = Ayy(i, j)+Para%Periodic_ISpan(2)
        end do
        end do
    end if

    if (Para%Iperiodic_Y == 1 .and. npy == 0) then
        do j = 1-LAP, 0
        do i = 1, nx
            Axx(i, j) = Axx(i, j)-Para%Periodic_JSpan(1)
            Ayy(i, j) = Ayy(i, j)-Para%Periodic_JSpan(2)
        end do
        end do
    end if

    if (Para%Iperiodic_Y == 1 .and. npy == npy0-1) then
        do j = ny+1, ny+LAP
        do i = 1, nx
            Axx(i, j) = Axx(i, j)+Para%Periodic_JSpan(1)
            Ayy(i, j) = Ayy(i, j)+Para%Periodic_JSpan(2)
        end do
        end do
    end if

end subroutine update_xy_periodic

subroutine update_xy_symmetric
    use mesh_data
    use flow_para
    use precision
    implicit none

    integer :: i, j, ii, jj, IERR
    real(dp) :: tx, ty, vx, vy, dot2
    real(dp) :: norm_t, eps_t

    eps_t = 1.0d-14

    if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
        do j = 1, ny
            if (j == 1) then
                tx = Axx(1, 2)-Axx(1, 1)
                ty = Ayy(1, 2)-Ayy(1, 1)
            else if (j == ny) then
                tx = Axx(1, ny)-Axx(1, ny-1)
                ty = Ayy(1, ny)-Ayy(1, ny-1)
            else
                tx = 0.5d0*(Axx(1, j+1)-Axx(1, j-1))
                ty = 0.5d0*(Ayy(1, j+1)-Ayy(1, j-1))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_xy_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            do i = 1-LAP, 0
                ii = 2-i
                vx = Axx(ii, j)-Axx(1, j)
                vy = Ayy(ii, j)-Ayy(1, j)
                dot2 = 2.0d0*(vx*tx+vy*ty)
                Axx(i, j) = Axx(1, j)+dot2*tx-vx
                Ayy(i, j) = Ayy(1, j)+dot2*ty-vy
            end do
        end do
    end if

    if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
        do j = 1, ny
            if (j == 1) then
                tx = Axx(nx, 2)-Axx(nx, 1)
                ty = Ayy(nx, 2)-Ayy(nx, 1)
            else if (j == ny) then
                tx = Axx(nx, ny)-Axx(nx, ny-1)
                ty = Ayy(nx, ny)-Ayy(nx, ny-1)
            else
                tx = 0.5d0*(Axx(nx, j+1)-Axx(nx, j-1))
                ty = 0.5d0*(Ayy(nx, j+1)-Ayy(nx, j-1))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_xy_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            do i = nx+1, nx+LAP
                ii = 2*nx-i
                vx = Axx(ii, j)-Axx(nx, j)
                vy = Ayy(ii, j)-Ayy(nx, j)
                dot2 = 2.0d0*(vx*tx+vy*ty)
                Axx(i, j) = Axx(nx, j)+dot2*tx-vx
                Ayy(i, j) = Ayy(nx, j)+dot2*ty-vy
            end do
        end do
    end if

    if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
        do j = 1-LAP, 0
        do i = 1, nx
            jj = 2-j
            if (i == 1) then
                tx = Axx(2, 1)-Axx(1, 1)
                ty = Ayy(2, 1)-Ayy(1, 1)
            else if (i == nx) then
                tx = Axx(nx, 1)-Axx(nx-1, 1)
                ty = Ayy(nx, 1)-Ayy(nx-1, 1)
            else
                tx = 0.5d0*(Axx(i+1, 1)-Axx(i-1, 1))
                ty = 0.5d0*(Ayy(i+1, 1)-Ayy(i-1, 1))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_xy_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            vx = Axx(i, jj)-Axx(i, 1)
            vy = Ayy(i, jj)-Ayy(i, 1)
            dot2 = 2.0d0*(vx*tx+vy*ty)
            Axx(i, j) = Axx(i, 1)+dot2*tx-vx
            Ayy(i, j) = Ayy(i, 1)+dot2*ty-vy
        end do
        end do
    end if

    if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
        do j = ny+1, ny+LAP
        do i = 1, nx
            jj = 2*ny-j
            if (i == 1) then
                tx = Axx(2, ny)-Axx(1, ny)
                ty = Ayy(2, ny)-Ayy(1, ny)
            else if (i == nx) then
                tx = Axx(nx, ny)-Axx(nx-1, ny)
                ty = Ayy(nx, ny)-Ayy(nx-1, ny)
            else
                tx = 0.5d0*(Axx(i+1, ny)-Axx(i-1, ny))
                ty = 0.5d0*(Ayy(i+1, ny)-Ayy(i-1, ny))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_xy_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            vx = Axx(i, jj)-Axx(i, ny)
            vy = Ayy(i, jj)-Ayy(i, ny)
            dot2 = 2.0d0*(vx*tx+vy*ty)
            Axx(i, j) = Axx(i, ny)+dot2*tx-vx
            Ayy(i, j) = Ayy(i, ny)+dot2*ty-vy
        end do
        end do
    end if

end subroutine update_xy_symmetric

subroutine update_metrics_symmetric
    use mesh_data
    use flow_para
    use precision
    implicit none

    integer :: i, j, ii, jj, IERR
    real(dp) :: tx, ty, dot2
    real(dp) :: ax, ay, bx, by
    real(dp) :: norm_t, eps_t

    eps_t = 1.0d-14

    if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
        do j = 1, ny
            if (j == 1) then
                tx = Axx(1, 2)-Axx(1, 1)
                ty = Ayy(1, 2)-Ayy(1, 1)
            else if (j == ny) then
                tx = Axx(1, ny)-Axx(1, ny-1)
                ty = Ayy(1, ny)-Ayy(1, ny-1)
            else
                tx = 0.5d0*(Axx(1, j+1)-Axx(1, j-1))
                ty = 0.5d0*(Ayy(1, j+1)-Ayy(1, j-1))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_metrics_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            do i = 1-LAP, 0
                ii = 2-i

                ax = Akx(ii, j)
                ay = Aky(ii, j)
                dot2 = 2.0d0*(ax*tx+ay*ty)
                Akx(i, j) = -(dot2*tx-ax)
                Aky(i, j) = -(dot2*ty-ay)

                bx = Aix(ii, j)
                by = Aiy(ii, j)
                dot2 = 2.0d0*(bx*tx+by*ty)
                Aix(i, j) = dot2*tx-bx
                Aiy(i, j) = dot2*ty-by

                Ajac(i, j) = Ajac(ii, j)
            end do
        end do
    end if

    if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
        do j = 1, ny
            if (j == 1) then
                tx = Axx(nx, 2)-Axx(nx, 1)
                ty = Ayy(nx, 2)-Ayy(nx, 1)
            else if (j == ny) then
                tx = Axx(nx, ny)-Axx(nx, ny-1)
                ty = Ayy(nx, ny)-Ayy(nx, ny-1)
            else
                tx = 0.5d0*(Axx(nx, j+1)-Axx(nx, j-1))
                ty = 0.5d0*(Ayy(nx, j+1)-Ayy(nx, j-1))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_metrics_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            do i = nx+1, nx+LAP
                ii = 2*nx-i

                ax = Akx(ii, j)
                ay = Aky(ii, j)
                dot2 = 2.0d0*(ax*tx+ay*ty)
                Akx(i, j) = -(dot2*tx-ax)
                Aky(i, j) = -(dot2*ty-ay)

                bx = Aix(ii, j)
                by = Aiy(ii, j)
                dot2 = 2.0d0*(bx*tx+by*ty)
                Aix(i, j) = dot2*tx-bx
                Aiy(i, j) = dot2*ty-by

                Ajac(i, j) = Ajac(ii, j)
            end do
        end do
    end if

    if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
        do j = 1-LAP, 0
        do i = 1, nx
            jj = 2-j
            if (i == 1) then
                tx = Axx(2, 1)-Axx(1, 1)
                ty = Ayy(2, 1)-Ayy(1, 1)
            else if (i == nx) then
                tx = Axx(nx, 1)-Axx(nx-1, 1)
                ty = Ayy(nx, 1)-Ayy(nx-1, 1)
            else
                tx = 0.5d0*(Axx(i+1, 1)-Axx(i-1, 1))
                ty = 0.5d0*(Ayy(i+1, 1)-Ayy(i-1, 1))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_metrics_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            ax = Akx(i, jj)
            ay = Aky(i, jj)
            dot2 = 2.0d0*(ax*tx+ay*ty)
            Akx(i, j) = dot2*tx-ax
            Aky(i, j) = dot2*ty-ay

            bx = Aix(i, jj)
            by = Aiy(i, jj)
            dot2 = 2.0d0*(bx*tx+by*ty)
            Aix(i, j) = -(dot2*tx-bx)
            Aiy(i, j) = -(dot2*ty-by)

            Ajac(i, j) = Ajac(i, jj)
        end do
        end do
    end if

    if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
        do j = ny+1, ny+LAP
        do i = 1, nx
            jj = 2*ny-j
            if (i == 1) then
                tx = Axx(2, ny)-Axx(1, ny)
                ty = Ayy(2, ny)-Ayy(1, ny)
            else if (i == nx) then
                tx = Axx(nx, ny)-Axx(nx-1, ny)
                ty = Ayy(nx, ny)-Ayy(nx-1, ny)
            else
                tx = 0.5d0*(Axx(i+1, ny)-Axx(i-1, ny))
                ty = 0.5d0*(Ayy(i+1, ny)-Ayy(i-1, ny))
            end if

            norm_t = sqrt(tx*tx+ty*ty)
            if (norm_t > eps_t) then
                tx = tx/norm_t
                ty = ty/norm_t
            else
                print *, "*** Error: zero tangent vector in update_metrics_symmetric ***"
                call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
            end if

            ax = Akx(i, jj)
            ay = Aky(i, jj)
            dot2 = 2.0d0*(ax*tx+ay*ty)
            Akx(i, j) = dot2*tx-ax
            Aky(i, j) = dot2*ty-ay

            bx = Aix(i, jj)
            by = Aiy(i, jj)
            dot2 = 2.0d0*(bx*tx+by*ty)
            Aix(i, j) = -(dot2*tx-bx)
            Aiy(i, j) = -(dot2*ty-by)

            Ajac(i, j) = Ajac(i, jj)
        end do
        end do
    end if

end subroutine update_metrics_symmetric

subroutine compute_Jocabian
    use mesh_data
    use Centrals
    implicit none

    integer :: i, j, IERR
    real(dp), allocatable, dimension(:, :) :: xk, xi, yk, yi
    allocate (xk(nx, ny), xi(nx, ny), yk(nx, ny), yi(nx, ny))

    call OCFD2d_dx0(Axx, xk, nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dx0(Ayy, yk, nx, ny, Scheme%Scheme_Vis, BcI)
    call OCFD2d_dy0(Axx, xi, nx, ny, Scheme%Scheme_Vis, BcJ)
    call OCFD2d_dy0(Ayy, yi, nx, ny, Scheme%Scheme_Vis, BcJ)

    do j = 1, ny
    do i = 1, nx
        Ajac(i, j) = 1.0d0/(xk(i, j)*yi(i, j)-xi(i, j)*yk(i, j))
        Akx(i, j) = Ajac(i, j)*yi(i, j)
        Aky(i, j) = -Ajac(i, j)*xi(i, j)
        Aix(i, j) = -Ajac(i, j)*yk(i, j)
        Aiy(i, j) = Ajac(i, j)*xk(i, j)
        if (Ajac(i, j) <= 0.0d0) then
            print *, "*** Error: negative Jacobian ***"
            print *, "Please Check Grid Orientation or Dimensions! "
            print *, i_offset(npx)+i-1, j_offset(npy)+j-1
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        end if
    end do
    end do

    deallocate (xk, xi, yk, yi)
    if (my_id == 0) print *, "Compute Jacobian OK"

end subroutine compute_Jocabian
