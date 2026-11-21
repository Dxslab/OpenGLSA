! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine flux_inviscid2d
    use flow_data
    use mesh_data
    implicit none

    integer :: i, j, m
    integer :: Ka1, Kb1, Ka2, Kb2, Bc_type(2)

    Ka1 = Scheme%Ka1
    Kb1 = Scheme%Kb1
    Ka2 = Scheme%Ka2
    Kb2 = Scheme%Kb2

    do j = 1-LAP, ny+LAP
    do i = 1-LAP, nx+LAP
        c(i, j) = sqrt(T(i, j))/Para%Ma
    end do
    end do

    !> x direction
    do j = 1, ny
        do i = 1-LAP, nx+LAP
            c1x(i) = c(i, j)
            d1x(i) = d(i, j)
            u1x(i) = u(i, j)
            v1x(i) = v(i, j)
            A1x(i) = Akx1(i, j)
            A2x(i) = Aky1(i, j)
        end do
        !
        call split2d_jac(nx, LAP, c1x, d1x, u1x, v1x, A1x, A2x, fpx, fmx)
        !
        Bc_type = BcI(j, :)
        do m = 1, 4
            call OCFD2d_flux1(fpx(1-LAP, m), hhx1, nx, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD2d_flux2(fmx(1-LAP, m), hhx2, nx, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)

            do i = 1, nx
                du(i, j, m) = (hhx1(i)-hhx1(i-1)+hhx2(i)-hhx2(i-1))
            end do
        end do
    end do
    !
    !> y direction
    do i = 1, nx
        do j = 1-LAP, ny+LAP
            c1y(j) = c(i, j)
            d1y(j) = d(i, j)
            u1y(j) = u(i, j)
            v1y(j) = v(i, j)
            A1y(j) = Aix1(i, j)
            A2y(j) = Aiy1(i, j)
        end do
        !
        call split2d_jac(ny, LAP, c1y, d1y, u1y, v1y, A1y, A2y, fpy, fmy)
        !
        Bc_type = BcJ(i, :)
        do m = 1, 4
            call OCFD2d_flux1(fpy(1-LAP, m), hhy1, ny, LAP, Ka1, Kb1, Scheme%Scheme_Invis, Bc_type)
            call OCFD2d_flux2(fmy(1-LAP, m), hhy2, ny, LAP, Ka2, Kb2, Scheme%Scheme_Invis, Bc_type)

            do j = 1, ny
                du(i, j, m) = -(du(i, j, m)+hhy1(j)-hhy1(j-1)+hhy2(j)-hhy2(j-1))*Ajac(i, j)
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
                Source_Invis(i, j, 4) = d(i, j)*v(i, j)*(c(i, j)**2/(Para%gamma-1.d0)+0.5d0*(u(i, j)**2+v(i, j)**2))/Ayy(i, j)
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

end subroutine flux_inviscid2d

subroutine OCFD2d_flux1(f1d, hh, nx1, LAP1, Ka, Kb, Scheme_Invis, Bc_type)
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

end subroutine OCFD2d_flux1

subroutine OCFD2d_flux2(f1d, hh, ny1, LAP1, Ka, Kb, Scheme_Invis, Bc_type)
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

end subroutine OCFD2d_flux2
