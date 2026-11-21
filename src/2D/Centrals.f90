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

    subroutine OCFD2d_dx0(f, fx, nx, ny, Num_Scheme, Bc)
        integer :: nx, ny, Num_Scheme, i, j, ierr
        integer :: Bc(ny, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)
        real(dp) :: fx(:, :)
        real(dp) :: a1, a2, b1, b2, b3
        a1 = 8.d0/12.d0            ! 4th central scheme
        a2 = -1.d0/12.d0

        b1 = 1.d0/60.d0             ! 6th central scheme
        b2 = -3.d0/20.d0
        b3 = 3.d0/4.d0

        !> Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do j = 1, ny
            do i = 1, nx
                fx(i, j) = a1*(f(i+1, j)-f(i-1, j))+ &
                           a2*(f(i+2, j)-f(i-2, j))
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do j = 1, ny
            do i = 1, nx
                fx(i, j) = b1*(f(i+3, j)-f(i-3, j))+ &
                           b2*(f(i+2, j)-f(i-2, j))+ &
                           b3*(f(i+1, j)-f(i-1, j))
            end do
            end do

        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- i- boundary  ---------------------------
        do j = 1, ny
            if (Bc(j, 1) == 1) then
                fx(1, j) = (-25.0d0*f(1, j)+48.0d0*f(2, j)-36.0d0*f(3, j)+ &
                            16.0d0*f(4, j)-3.0d0*f(5, j))/12.0d0
                fx(2, j) = (-3.0d0*f(1, j)-10.0d0*f(2, j)+18.0d0*f(3, j)- &
                            6.0d0*f(4, j)+f(5, j))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(3, j) = a1*(f(4, j)-f(2, j))+ &
                               a2*(f(5, j)-f(1, j))  ! 4th central scheme
                end if
            end if
        end do

        !--------- i+ boundary ------------------------
        do j = 1, ny
            if (Bc(j, 2) == 1) then
                fx(nx, j) = (3.0d0*f(nx-4, j)-16.0d0*f(nx-3, j)+36.0d0*f(nx-2, j)- &
                             48.0d0*f(nx-1, j)+25.0d0*f(nx, j))/12.0d0
                fx(nx-1, j) = (-f(nx-4, j)+6.0d0*f(nx-3, j)-18.0d0*f(nx-2, j)+ &
                               10.0d0*f(nx-1, j)+3.0d0*f(nx, j))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(nx-2, j) = a1*(f(nx-1, j)-f(nx-3, j))+ &
                                  a2*(f(nx, j)-f(nx-4, j))   ! 4th central scheme
                end if
            end if
        end do

    end subroutine OCFD2d_dx0

    subroutine OCFD2d_dy0(f, fy, nx, ny, Num_Scheme, Bc)
        integer :: nx, ny, Num_Scheme, i, j, ierr
        integer :: Bc(nx, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)
        real(dp) :: fy(:, :)
        real(dp) :: a1, a2, b1, b2, b3

        a1 = 8.d0/12.d0             ! 4th central scheme
        a2 = -1.d0/12.d0

        b1 = 1.d0/60.d0
        b2 = -3.d0/20.d0
        b3 = 3.d0/4.d0

        !------Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            do j = 1, ny
            do i = 1, nx
                fy(i, j) = a1*(f(i, j+1)-f(i, j-1))+ &
                           a2*(f(i, j+2)-f(i, j-2))
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do j = 1, ny
            do i = 1, nx
                fy(i, j) = b1*(f(i, j+3)-f(i, j-3))+ &
                           b2*(f(i, j+2)-f(i, j-2))+ &
                           b3*(f(i, j+1)-f(i, j-1))
            end do
            end do

        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- j- boundary  ---------------------------
        do i = 1, nx
            if (Bc(i, 1) == 1) then
                fy(i, 1) = (-25.0d0*f(i, 1)+48.0d0*f(i, 2)-36.0d0*f(i, 3)+ &
                            16.0d0*f(i, 4)-3.0d0*f(i, 5))/12.0d0
                fy(i, 2) = (-3.0d0*f(i, 1)-10.0d0*f(i, 2)+18.0d0*f(i, 3)- &
                            6.0d0*f(i, 4)+f(i, 5))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, 3) = a1*(f(i, 4)-f(i, 2))+ &
                               a2*(f(i, 5)-f(i, 1))  ! 4th central scheme
                end if
            end if
        end do

        !--------- j+ boundary ------------------------
        do i = 1, nx
            if (Bc(i, 2) == 1) then
                fy(i, ny) = (3.0d0*f(i, ny-4)-16.0d0*f(i, ny-3)+36.0d0*f(i, ny-2)- &
                             48.0d0*f(i, ny-1)+25.0d0*f(i, ny))/12.0d0
                fy(i, ny-1) = (-f(i, ny-4)+6.0d0*f(i, ny-3)-18.0d0*f(i, ny-2)+ &
                               10.0d0*f(i, ny-1)+3.0d0*f(i, ny))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, ny-2) = a1*(f(i, ny-1)-f(i, ny-3))+ &
                                  a2*(f(i, ny)-f(i, ny-4))   ! 4th central scheme
                end if
            end if
        end do

    end subroutine OCFD2d_dy0

! Finite difference for viscous terms  (second derivative)
    subroutine OCFD2d_dx2(f, fx, nx, ny, Num_Scheme, Bc)
        integer :: nx, ny, Num_Scheme, i, j, ierr
        integer :: Bc(ny, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)
        real(dp) :: fx(:, :)
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
            do j = 1, ny
            do i = 1, nx
                fx(i, j) = a1*f(i, j)+ &
                           a2*(f(i+1, j)+f(i-1, j))+ &
                           a3*(f(i+2, j)+f(i-2, j))
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do j = 1, ny
            do i = 1, nx
                fx(i, j) = b1*f(i, j)+ &
                           b2*(f(i+1, j)+f(i-1, j))+ &
                           b3*(f(i+2, j)+f(i-2, j))+ &
                           b4*(f(i+3, j)+f(i-3, j))
            end do
            end do

        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- i- boundary  ---------------------------
        do j = 1, ny
            if (Bc(j, 1) == 1) then
                fx(1, j) = (35.0d0*f(1, j)-104.0d0*f(2, j)+114.0d0*f(3, j)- &
                            56.0d0*f(4, j)+11.0d0*f(5, j))/12.0d0
                fx(2, j) = (11.0d0*f(1, j)-20.0d0*f(2, j)+6.0d0*f(3, j)+ &
                            4.0d0*f(4, j)-f(5, j))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(3, j) = a1*f(3, j)+ &
                               a2*(f(4, j)+f(2, j))+ &
                               a3*(f(5, j)+f(1, j))
                end if
            end if
        end do

        !--------- i+ boundary ------------------------
        do j = 1, ny
            if (Bc(j, 2) == 1) then
                fx(nx, j) = (11.0d0*f(nx-4, j)-56.0d0*f(nx-3, j)+114.0d0*f(nx-2, j)- &
                             104.0d0*f(nx-1, j)+35.0d0*f(nx, j))/12.0d0
                fx(nx-1, j) = (-f(nx-4, j)+4.0d0*f(nx-3, j)+6.0d0*f(nx-2, j)- &
                               20.0d0*f(nx-1, j)+11.0d0*f(nx, j))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fx(nx-2, j) = a1*f(nx-2, j)+ &
                                  a2*(f(nx-3, j)+f(nx-1, j))+ &
                                  a3*(f(nx-4, j)+f(nx, j))
                end if
            end if
        end do

    end subroutine OCFD2d_dx2

    subroutine OCFD2d_dy2(f, fy, nx, ny, Num_Scheme, Bc)
        integer :: nx, ny, Num_Scheme, i, j, ierr
        integer :: Bc(nx, 2)
        real(dp) :: f(1-LAP:nx+LAP, 1-LAP:ny+LAP)
        real(dp) :: fy(:, :)
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
            do j = 1, ny
            do i = 1, nx
                fy(i, j) = a1*f(i, j)+ &
                           a2*(f(i, j+1)+f(i, j-1))+ &
                           a3*(f(i, j+2)+f(i, j-2))
            end do
            end do
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            do j = 1, ny
            do i = 1, nx
                fy(i, j) = b1*f(i, j)+ &
                           b2*(f(i, j+1)+f(i, j-1))+ &
                           b3*(f(i, j+2)+f(i, j-2))+ &
                           b4*(f(i, j+3)+f(i, j-3))
            end do
            end do

        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        !---------Boundary Scheme ------- (low-order scheme)----
        !---------- j- boundary  ---------------------------
        do i = 1, nx
            if (Bc(i, 1) == 1) then
                fy(i, 1) = (35.0d0*f(i, 1)-104.0d0*f(i, 2)+114.0d0*f(i, 3)- &
                            56.0d0*f(i, 4)+11.0d0*f(i, 5))/12.0d0
                fy(i, 2) = (11.0d0*f(i, 1)-20.0d0*f(i, 2)+6.0d0*f(i, 3)+ &
                            4.0d0*f(i, 4)-f(i, 5))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, 3) = a1*f(i, 3)+ &
                               a2*(f(i, 4)+f(i, 2))+ &
                               a3*(f(i, 5)+f(i, 1))
                end if
            end if
        end do

        !--------- j+ boundary ------------------------
        do i = 1, nx
            if (Bc(i, 2) == 1) then
                fy(i, ny) = (11.0d0*f(i, ny-4)-56.0d0*f(i, ny-3)+114.0d0*f(i, ny-2)- &
                             104.0d0*f(i, ny-1)+35.0d0*f(i, ny))/12.0d0
                fy(i, ny-1) = (-f(i, ny-4)+4.0d0*f(i, ny-3)+6.0d0*f(i, ny-2)- &
                               20.0d0*f(i, ny-1)+11.0d0*f(i, ny))/12.0d0
                if (Num_Scheme == OCFD_Scheme_CD6) then
                    fy(i, ny-2) = a1*f(i, ny-2)+ &
                                  a2*(f(i, ny-3)+f(i, ny-1))+ &
                                  a3*(f(i, ny-4)+f(i, ny))
                end if
            end if
        end do

    end subroutine OCFD2d_dy2

    subroutine OCFD2d_dx0_single_node(f, df, Num_Scheme)
        integer :: Num_Scheme, ierr
        real(dp) :: f(-LAP:LAP), df
        real(dp) :: a1, a2, b1, b2, b3

        a1 = 8.d0/12.d0            ! 4th central scheme
        a2 = -1.d0/12.d0

        b1 = 1.d0/60.d0             ! 6th central scheme
        b2 = -3.d0/20.d0
        b3 = 3.d0/4.d0

        !> Scheme for inner points  ( CD4, CD6, CD8)
        if (Num_Scheme == OCFD_Scheme_CD4) then
            df = a1*(f(1)-f(-1))+a2*(f(2)-f(-2))
        else if (Num_Scheme == OCFD_Scheme_CD6) then
            df = b1*(f(3)-f(-3))+b2*(f(2)-f(-2))+b3*(f(1)-f(-1))
        else
            print *, '*** This Numerical Scheme is not supported in viscous terms ! ***'
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine OCFD2d_dx0_single_node

end module Centrals
