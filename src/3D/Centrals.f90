! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

!> Finite difference for viscous terms  (Central Schemes)
module Centrals
    use Constants
    use Para_mpi, only: LAP
    implicit none

contains

    subroutine OCFD_dx0(f, fx, nx, ny, nz, Num_Scheme, Bc)
        integer :: nx, ny, nz, Num_Scheme, i, j, k, ierr
        integer :: Bc(ny, nz, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp) :: fx(:, :, :)
        real(dp) :: a1, a2, b1, b2, b3
        a1 = 8.d0/12.d0            ! 4th central scheme
        a2 = -1.d0/12.d0

        b1 = 1.d0/60.d0             ! 6th central scheme
        b2 = -3.d0/20.d0
        b3 = 3.d0/4.d0

        !> Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fx(i, j, k) = a1*(f(i+1, j, k)-f(i-1, j, k))+ &
                              a2*(f(i+2, j, k)-f(i-2, j, k))
            end do
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fx(i, j, k) = b1*(f(i+3, j, k)-f(i-3, j, k))+ &
                              b2*(f(i+2, j, k)-f(i-2, j, k))+ &
                              b3*(f(i+1, j, k)-f(i-1, j, k))
            end do
            end do
            end do
        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- i- boundary  ---------------------------
        do k = 1, nz
        do j = 1, ny
            if (Bc(j, k, 1) == 1) then
                fx(1, j, k) = (-25.0d0*f(1, j, k)+48.0d0*f(2, j, k)-36.0d0*f(3, j, k)+ &
                               16.0d0*f(4, j, k)-3.0d0*f(5, j, k))/12.0d0
                fx(2, j, k) = (-3.0d0*f(1, j, k)-10.0d0*f(2, j, k)+18.0d0*f(3, j, k)- &
                               6.0d0*f(4, j, k)+f(5, j, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(3, j, k) = a1*(f(4, j, k)-f(2, j, k))+ &
                                  a2*(f(5, j, k)-f(1, j, k))  ! 4th central scheme
                end if
            end if
        end do
        end do
        !--------- i+ boundary ------------------------
        do k = 1, nz
        do j = 1, ny
            if (Bc(j, k, 2) == 1) then
                fx(nx, j, k) = (3.0d0*f(nx-4, j, k)-16.0d0*f(nx-3, j, k)+36.0d0*f(nx-2, j, k)- &
                                48.0d0*f(nx-1, j, k)+25.0d0*f(nx, j, k))/12.0d0
                fx(nx-1, j, k) = (-f(nx-4, j, k)+6.0d0*f(nx-3, j, k)-18.0d0*f(nx-2, j, k)+ &
                                  10.0d0*f(nx-1, j, k)+3.0d0*f(nx, j, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(nx-2, j, k) = a1*(f(nx-1, j, k)-f(nx-3, j, k))+ &
                                     a2*(f(nx, j, k)-f(nx-4, j, k))   ! 4th central scheme
                end if
            end if
        end do
        end do

    end subroutine OCFD_dx0

    subroutine OCFD_dy0(f, fy, nx, ny, nz, Num_Scheme, Bc)
        integer :: nx, ny, nz, Num_Scheme, i, j, k, ierr
        integer :: Bc(nx, nz, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp) :: fy(:, :, :)
        real(dp) :: a1, a2, b1, b2, b3

        a1 = 8.d0/12.d0             ! 4th central scheme
        a2 = -1.d0/12.d0

        b1 = 1.d0/60.d0
        b2 = -3.d0/20.d0
        b3 = 3.d0/4.d0

        !------Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fy(i, j, k) = a1*(f(i, j+1, k)-f(i, j-1, k))+ &
                              a2*(f(i, j+2, k)-f(i, j-2, k))
            end do
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fy(i, j, k) = b1*(f(i, j+3, k)-f(i, j-3, k))+ &
                              b2*(f(i, j+2, k)-f(i, j-2, k))+ &
                              b3*(f(i, j+1, k)-f(i, j-1, k))
            end do
            end do
            end do
        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- j- boundary  ---------------------------
        do k = 1, nz
        do i = 1, nx
            if (Bc(i, k, 1) == 1) then
                fy(i, 1, k) = (-25.0d0*f(i, 1, k)+48.0d0*f(i, 2, k)-36.0d0*f(i, 3, k)+ &
                               16.0d0*f(i, 4, k)-3.0d0*f(i, 5, k))/12.0d0
                fy(i, 2, k) = (-3.0d0*f(i, 1, k)-10.0d0*f(i, 2, k)+18.0d0*f(i, 3, k)- &
                               6.0d0*f(i, 4, k)+f(i, 5, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, 3, k) = a1*(f(i, 4, k)-f(i, 2, k))+ &
                                  a2*(f(i, 5, k)-f(i, 1, k))  ! 4th central scheme
                end if
            end if
        end do
        end do

        !--------- j+ boundary ------------------------
        do k = 1, nz
        do i = 1, nx
            if (Bc(i, k, 2) == 1) then
                fy(i, ny, k) = (3.0d0*f(i, ny-4, k)-16.0d0*f(i, ny-3, k)+36.0d0*f(i, ny-2, k)- &
                                48.0d0*f(i, ny-1, k)+25.0d0*f(i, ny, k))/12.0d0
                fy(i, ny-1, k) = (-f(i, ny-4, k)+6.0d0*f(i, ny-3, k)-18.0d0*f(i, ny-2, k)+ &
                                  10.0d0*f(i, ny-1, k)+3.0d0*f(i, ny, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, ny-2, k) = a1*(f(i, ny-1, k)-f(i, ny-3, k))+ &
                                     a2*(f(i, ny, k)-f(i, ny-4, k))   ! 4th central scheme
                end if
            end if
        end do
        end do

    end subroutine OCFD_dy0

    subroutine OCFD_dz0(f, fz, nx, ny, nz, Num_Scheme, Bc)
        integer :: nx, ny, nz, Num_Scheme, i, j, k, ierr
        integer :: Bc(nx, ny, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp) :: fz(:, :, :)
        real(dp) :: a1, a2, b1, b2, b3

        a1 = 8.d0/12.d0             ! 4th central scheme
        a2 = -1.d0/12.d0

        b1 = 1.d0/60.d0
        b2 = -3.d0/20.d0
        b3 = 3.d0/4.d0

        !------Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fz(i, j, k) = a1*(f(i, j, k+1)-f(i, j, k-1))+ &
                              a2*(f(i, j, k+2)-f(i, j, k-2))
            end do
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fz(i, j, k) = b1*(f(i, j, k+3)-f(i, j, k-3))+ &
                              b2*(f(i, j, k+2)-f(i, j, k-2))+ &
                              b3*(f(i, j, k+1)-f(i, j, k-1))
            end do
            end do
            end do
        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- j- boundary  ---------------------------
        do j = 1, ny
        do i = 1, nx
            if (Bc(i, j, 1) == 1) then
                fz(i, j, 1) = (-25.0d0*f(i, j, 1)+48.0d0*f(i, j, 2)-36.0d0*f(i, j, 3)+ &
                               16.0d0*f(i, j, 4)-3.0d0*f(i, j, 5))/12.0d0
                fz(i, j, 2) = (-3.0d0*f(i, j, 1)-10.0d0*f(i, j, 2)+18.0d0*f(i, j, 3)- &
                               6.0d0*f(i, j, 4)+f(i, j, 5))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fz(i, j, 3) = a1*(f(i, j, 4)-f(i, j, 2))+ &
                                  a2*(f(i, j, 5)-f(i, j, 1))  ! 4th central scheme
                end if
            end if
        end do
        end do

        !--------- j+ boundary ------------------------
        do j = 1, ny
        do i = 1, nx
            if (Bc(i, j, 2) == 1) then
                fz(i, j, nz) = (3.0d0*f(i, j, nz-4)-16.0d0*f(i, j, nz-3)+36.0d0*f(i, j, nz-2)- &
                                48.0d0*f(i, j, nz-1)+25.0d0*f(i, j, nz))/12.0d0
                fz(i, j, nz-1) = (-f(i, j, nz-4)+6.0d0*f(i, j, nz-3)-18.0d0*f(i, j, nz-2)+ &
                                  10.0d0*f(i, j, nz-1)+3.0d0*f(i, j, nz))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fz(i, j, nz-2) = a1*(f(i, j, nz-1)-f(i, j, nz-3))+ &
                                     a2*(f(i, j, nz)-f(i, j, nz-4))   ! 4th central scheme
                end if
            end if
        end do
        end do

    end subroutine OCFD_dz0

! Finite difference for viscous terms  (second derivative)
    subroutine OCFD_dx2(f, fx, nx, ny, nz, Num_Scheme, Bc)
        integer :: nx, ny, nz, Num_Scheme, i, j, k, ierr
        integer :: Bc(ny, nz, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp) :: fx(:, :, :)
        real(dp) :: a1, a2, a3, b1, b2, b3, b4

        a1 = -30.d0/12.d0             ! 4th central scheme
        a2 = 16.d0/12.d0
        a3 = -1.d0/12.d0

        b1 = -490.d0/180.d0             ! 6th central scheme
        b2 = 270.d0/180.d0
        b3 = -27.d0/180.d0
        b4 = 2.d0/180.d0

        !------Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fx(i, j, k) = a1*f(i, j, k)+ &
                              a2*(f(i+1, j, k)+f(i-1, j, k))+ &
                              a3*(f(i+2, j, k)+f(i-2, j, k))
            end do
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fx(i, j, k) = b1*f(i, j, k)+ &
                              b2*(f(i+1, j, k)+f(i-1, j, k))+ &
                              b3*(f(i+2, j, k)+f(i-2, j, k))+ &
                              b4*(f(i+3, j, k)+f(i-3, j, k))
            end do
            end do
            end do
        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- i- boundary  ---------------------------
        do k = 1, nz
        do j = 1, ny
            if (Bc(j, k, 1) == 1) then
                fx(1, j, k) = (35.0d0*f(1, j, k)-104.0d0*f(2, j, k)+114.0d0*f(3, j, k)- &
                               56.0d0*f(4, j, k)+11.0d0*f(5, j, k))/12.0d0
                fx(2, j, k) = (11.0d0*f(1, j, k)-20.0d0*f(2, j, k)+6.0d0*f(3, j, k)+ &
                               4.0d0*f(4, j, k)-f(5, j, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(3, j, k) = a1*f(3, j, k)+ &
                                  a2*(f(4, j, k)+f(2, j, k))+ &
                                  a3*(f(5, j, k)+f(1, j, k))
                end if
            end if
        end do
        end do
        !--------- i+ boundary ------------------------
        do k = 1, nz
        do j = 1, ny
            if (Bc(j, k, 2) == 1) then
                fx(nx, j, k) = (11.0d0*f(nx-4, j, k)-56.0d0*f(nx-3, j, k)+114.0d0*f(nx-2, j, k)- &
                                104.0d0*f(nx-1, j, k)+35.0d0*f(nx, j, k))/12.0d0
                fx(nx-1, j, k) = (-f(nx-4, j, k)+4.0d0*f(nx-3, j, k)+6.0d0*f(nx-2, j, k)- &
                                  20.0d0*f(nx-1, j, k)+11.0d0*f(nx, j, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(nx-2, j, k) = a1*f(nx-2, j, k)+ &
                                     a2*(f(nx-3, j, k)+f(nx-1, j, k))+ &
                                     a3*(f(nx-4, j, k)+f(nx, j, k))
                end if
            end if
        end do
        end do

    end subroutine OCFD_dx2

    subroutine OCFD_dy2(f, fy, nx, ny, nz, Num_Scheme, Bc)
        integer :: nx, ny, nz, Num_Scheme, i, j, k, ierr
        integer :: Bc(nx, nz, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp) :: fy(:, :, :)
        real(dp) :: a1, a2, a3, b1, b2, b3, b4

        a1 = -30.d0/12.d0             ! 4th central scheme
        a2 = 16.d0/12.d0
        a3 = -1.d0/12.d0

        b1 = -490.d0/180.d0             ! 6th central scheme
        b2 = 270.d0/180.d0
        b3 = -27.d0/180.d0
        b4 = 2.d0/180.d0

        !------Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fy(i, j, k) = a1*f(i, j, k)+ &
                              a2*(f(i, j+1, k)+f(i, j-1, k))+ &
                              a3*(f(i, j+2, k)+f(i, j-2, k))
            end do
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fy(i, j, k) = b1*f(i, j, k)+ &
                              b2*(f(i, j+1, k)+f(i, j-1, k))+ &
                              b3*(f(i, j+2, k)+f(i, j-2, k))+ &
                              b4*(f(i, j+3, k)+f(i, j-3, k))
            end do
            end do
            end do

        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- j- boundary  ---------------------------
        do k = 1, nz
        do i = 1, nx
            if (Bc(i, k, 1) == 1) then
                fy(i, 1, k) = (35.0d0*f(i, 1, k)-104.0d0*f(i, 2, k)+114.0d0*f(i, 3, k)- &
                               56.0d0*f(i, 4, k)+11.0d0*f(i, 5, k))/12.0d0
                fy(i, 2, k) = (11.0d0*f(i, 1, k)-20.0d0*f(i, 2, k)+6.0d0*f(i, 3, k)+ &
                               4.0d0*f(i, 4, k)-f(i, 5, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, 3, k) = a1*f(i, 3, k)+ &
                                  a2*(f(i, 4, k)+f(i, 2, k))+ &
                                  a3*(f(i, 5, k)+f(i, 1, k))
                end if
            end if
        end do
        end do
        !--------- j+ boundary ------------------------
        do k = 1, nz
        do i = 1, nx
            if (Bc(i, k, 2) == 1) then
                fy(i, ny, k) = (11.0d0*f(i, ny-4, k)-56.0d0*f(i, ny-3, k)+114.0d0*f(i, ny-2, k)- &
                                104.0d0*f(i, ny-1, k)+35.0d0*f(i, ny, k))/12.0d0
                fy(i, ny-1, k) = (-f(i, ny-4, k)+4.0d0*f(i, ny-3, k)+6.0d0*f(i, ny-2, k)- &
                                  20.0d0*f(i, ny-1, k)+11.0d0*f(i, ny, k))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, ny-2, k) = a1*f(i, ny-2, k)+ &
                                     a2*(f(i, ny-3, k)+f(i, ny-1, k))+ &
                                     a3*(f(i, ny-4, k)+f(i, ny, k))
                end if
            end if
        end do
        end do

    end subroutine OCFD_dy2

    subroutine OCFD_dz2(f, fz, nx, ny, nz, Num_Scheme, Bc)
        integer :: nx, ny, nz, Num_Scheme, i, j, k, ierr
        integer :: Bc(nx, ny, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1-LAP:nz+LAP)
        real(dp) :: fz(:, :, :)
        real(dp) :: a1, a2, a3, b1, b2, b3, b4

        a1 = -30.d0/12.d0             ! 4th central scheme
        a2 = 16.d0/12.d0
        a3 = -1.d0/12.d0

        b1 = -490.d0/180.d0             ! 6th central scheme
        b2 = 270.d0/180.d0
        b3 = -27.d0/180.d0
        b4 = 2.d0/180.d0

        !------Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fz(i, j, k) = a1*f(i, j, k)+ &
                              a2*(f(i, j, k+1)+f(i, j, k-1))+ &
                              a3*(f(i, j, k+2)+f(i, j, k-2))
            end do
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                fz(i, j, k) = b1*f(i, j, k)+ &
                              b2*(f(i, j, k+1)+f(i, j, k-1))+ &
                              b3*(f(i, j, k+2)+f(i, j, k-2))+ &
                              b4*(f(i, j, k+3)+f(i, j, k-3))
            end do
            end do
            end do
        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- j- boundary  ---------------------------
        do j = 1, ny
        do i = 1, nx
            if (Bc(i, j, 1) == 1) then
                fz(i, j, 1) = (35.0d0*f(i, j, 1)-104.0d0*f(i, j, 2)+114.0d0*f(i, j, 3)- &
                               56.0d0*f(i, j, 4)+11.0d0*f(i, j, 5))/12.0d0
                fz(i, j, 2) = (11.0d0*f(i, j, 1)-20.0d0*f(i, j, 2)+6.0d0*f(i, j, 3)+ &
                               4.0d0*f(i, j, 4)-f(i, j, 5))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fz(i, j, 3) = a1*f(i, j, 3)+ &
                                  a2*(f(i, j, 4)+f(i, j, 2))+ &
                                  a3*(f(i, j, 5)+f(i, j, 1))
                end if
            end if
        end do
        end do

        !--------- j+ boundary ------------------------
        do j = 1, ny
        do i = 1, nx
            if (Bc(i, j, 2) == 1) then
                fz(i, j, nz) = (11.0d0*f(i, j, nz-4)-56.0d0*f(i, j, nz-3)+114.0d0*f(i, j, nz-2)- &
                                104.0d0*f(i, j, nz-1)+35.0d0*f(i, j, nz))/12.0d0
                fz(i, j, nz-1) = (-f(i, j, nz-4)+4.0d0*f(i, j, nz-3)+6.0d0*f(i, j, nz-2)- &
                                  20.0d0*f(i, j, nz-1)+11.0d0*f(i, j, nz))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fz(i, j, nz-2) = a1*f(i, j, nz-2)+ &
                                     a2*(f(i, j, nz-3)+f(i, j, nz-1))+ &
                                     a3*(f(i, j, nz-4)+f(i, j, nz))
                end if
            end if
        end do
        end do

    end subroutine OCFD_dz2

end module Centrals
