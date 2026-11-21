! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid
    use flow_data
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
        c(i, j, k) = sqrt(T(i, j, k))/Para%Ma
    end do
    end do
    end do

    !> x direction
    do k = 1, nz
    do j = 1, ny
        do i = 1-LAP, nx+LAP
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
        !
        Bc_type = BcI(j, k, :)
        do m = 1, 5
            call OCFD_flux1(fpx(1-LAP, m), hhx1, nx, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD_flux2(fmx(1-LAP, m), hhx2, nx, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)

            do i = 1, nx
                du(i, j, k, m) = (hhx1(i)-hhx1(i-1)+hhx2(i)-hhx2(i-1))
            end do
        end do
    end do
    end do
    !
    !> y direction
    do k = 1, nz
    do i = 1, nx
        do j = 1-LAP, ny+LAP
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
        !
        Bc_type = BcJ(i, k, :)
        do m = 1, 5
            call OCFD_flux1(fpy(1-LAP, m), hhy1, ny, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD_flux2(fmy(1-LAP, m), hhy2, ny, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)

            do j = 1, ny
                du(i, j, k, m) = du(i, j, k, m)+hhy1(j)-hhy1(j-1)+hhy2(j)-hhy2(j-1)
            end do
        end do
    end do
    end do
    !
    !> z direction
    do j = 1, ny
    do i = 1, nx
        do k = 1-LAP, nz+LAP
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
        !
        Bc_type = BcK(i, j, :)
        do m = 1, 5
            call OCFD_flux1(fpz(1-LAP, m), hhz1, nz, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD_flux2(fmz(1-LAP, m), hhz2, nz, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)

            do k = 1, nz
                du(i, j, k, m) = -(du(i, j, k, m)+hhz1(k)-hhz1(k-1)+hhz2(k)-hhz2(k-1))*Ajac(i, j, k)
            end do
        end do
    end do
    end do

end subroutine flux_inviscid

subroutine OCFD_flux1(f1d, hh, nx1, LAP1, Ka, Kb, Scheme_Invis, Bc_type)
    use Constants
    use WENO_Reconstruction
    implicit none
    integer ::  nx1, LAP1, Ka, Kb, Scheme_Invis, Bc_type(2), i, IERR
    real(dp) :: f1d(1-LAP1:nx1+LAP1), hh(0:nx1)

    select case (Scheme_Invis)
    case (OCFD_Scheme_WENOZ5)
        call OCFD_WENOZ5P(f1d, hh, nx1, LAP1)
    case (OCFD_Scheme_LINEAR5)
        call OCFD_LINEAR5P(f1d, hh, nx1, LAP1)
    case (OCFD_Scheme_WENOZQ5)
        call OCFD_WENOZQ5P(f1d, hh, nx1, LAP1)
    case (OCFD_Scheme_WENOMR5)
        call OCFD_WENOMR5P(f1d, hh, nx1, LAP1)
    case DEFAULT
        print *, '*** Error: Unknown scheme number in flux_x1 ***'
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select
    !> Boundary scheme
    if (Bc_type(1) == 1) then
        do i = 0, -Ka
            if (i == 0 .or. i == 1) then
                hh(i) = (2.0d0*f1d(1)+5.0d0*f1d(2)-f1d(3))/6.0d0
            else if (i == 2) then
                if (Scheme_Invis == OCFD_Scheme_LINEAR5) then
                    call hh_LINEAR5P_boundary(Ka, Kb, f1d(i+Ka:i+Kb), hh(i), del_left)
                else
                    call hh_WENOMR3P(Ka, Kb, f1d(i+Ka:i+Kb), hh(i))
                end if
            else
                call hh_WENOZQ5P(Ka, Kb, f1d(i+Ka:i+Kb), hh(i))
            end if
        end do
    end if

    if (Bc_type(2) == 1) then
        do i = nx1-Kb+1, nx1
            if (i == nx1) then
                hh(i) = (2.0d0*f1d(nx1-2)-7.0d0*f1d(nx1-1)+11.0d0*f1d(nx1))/6.0d0
            else if (i == nx1-1) then
                if (Scheme_Invis == OCFD_Scheme_LINEAR5) then
                    call hh_LINEAR5P_boundary(Ka, Kb, f1d(i+Ka:i+Kb), hh(i), del_right)
                else
                    call hh_WENOMR3P(Ka, Kb, f1d(i+Ka:i+Kb), hh(i))
                end if
            else
                call hh_WENOZQ5P(Ka, Kb, f1d(i+Ka:i+Kb), hh(i))
            end if
        end do
    end if

end subroutine OCFD_flux1

subroutine OCFD_flux2(f1d, hh, ny1, LAP1, Ka, Kb, Scheme_Invis, Bc_type)
    use Constants
    use WENO_Reconstruction
    implicit none
    integer ::  ny1, LAP1, Ka, Kb, Scheme_Invis, Bc_type(2), j, IERR
    real(dp) :: f1d(1-LAP1:ny1+LAP1), hh(0:ny1)

    select case (Scheme_Invis)
    case (OCFD_Scheme_WENOZ5)
        call OCFD_WENOZ5M(f1d, hh, ny1, LAP1)
    case (OCFD_Scheme_LINEAR5)
        call OCFD_LINEAR5M(f1d, hh, ny1, LAP1)
    case (OCFD_Scheme_WENOZQ5)
        call OCFD_WENOZQ5M(f1d, hh, ny1, LAP1)
    case (OCFD_Scheme_WENOMR5)
        call OCFD_WENOMR5M(f1d, hh, ny1, LAP1)
    case DEFAULT
        print *, '*** Error: Unknown scheme number in flux_x2 ***'
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select
    !> Boundary scheme
    if (Bc_type(1) == 1) then
        do j = 0, -Ka
            if (j == 0) then
                hh(j) = (2.0d0*f1d(3)-7.0d0*f1d(2)+11.0d0*f1d(1))/6.0d0
            else if (j == 1) then
                if (Scheme_Invis == OCFD_Scheme_LINEAR5) then
                    call hh_LINEAR5M_boundary(Ka, Kb, f1d(j+Ka:j+Kb), hh(j), del_left)
                else
                    call hh_WENOMR3M(Ka, Kb, f1d(j+Ka:j+Kb), hh(j))
                end if
            else
                call hh_WENOMR5M(Ka, Kb, f1d(j+Ka:j+Kb), hh(j))
            end if
        end do
    end if

    if (Bc_type(2) == 1) then
        do j = ny1-Kb+1, ny1
            if (j == ny1 .or. j == ny1-1) then
                hh(j) = (2.0d0*f1d(ny1)+5.0d0*f1d(ny1-1)-f1d(ny1-2))/6.0d0
            else if (j == ny1-2) then
                if (Scheme_Invis == OCFD_Scheme_LINEAR5) then
                    call hh_LINEAR5M_boundary(Ka, Kb, f1d(j+Ka:j+Kb), hh(j), del_right)
                else
                    call hh_WENOMR3M(Ka, Kb, f1d(j+Ka:j+Kb), hh(j))
                end if
            else
                call hh_WENOMR5M(Ka, Kb, f1d(j+Ka:j+Kb), hh(j))
            end if
        end do
    end if

end subroutine OCFD_flux2
