! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine apply_bc
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, k, IERR

    select case (Para%IBC)
    case (BC_None)
    case (BC_Flat_Plate)
        call ocfd_bc_Flat_Plate ! flat plate boundary
    case (BC_Lid_Driven_Cavity)
        call ocfd_bc_Lid_Driven_Cavity ! lid-driven cavity flow
    case (BC_User_Def)
        call ocfd_bc_User_Def ! user-define boundaries
    case DEFAULT
        if (my_id == 0) print *, "This boundary condition is not supported!!!"
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select

    if (Para%runMode == 0) then
        !
        if (npx == 0 .and. Para%Iperiodic_X == 0) then
            do k = 1, nz
            do j = 1, ny
                f(1, j, k, 1) = d(1, j, k)
                f(1, j, k, 2) = d(1, j, k)*u(1, j, k)
                f(1, j, k, 3) = d(1, j, k)*v(1, j, k)
                f(1, j, k, 4) = d(1, j, k)*w(1, j, k)
                f(1, j, k, 5) = d(1, j, k)*(Para%Cv*T(1, j, k)+(u(1, j, k)**2+v(1, j, k)**2+w(1, j, k)**2)*0.5d0)
            end do
            end do
        end if

        if (npx == npx0-1 .and. Para%Iperiodic_X == 0) then
            do k = 1, nz
            do j = 1, ny
                f(nx, j, k, 1) = d(nx, j, k)
                f(nx, j, k, 2) = d(nx, j, k)*u(nx, j, k)
                f(nx, j, k, 3) = d(nx, j, k)*v(nx, j, k)
                f(nx, j, k, 4) = d(nx, j, k)*w(nx, j, k)
                f(nx, j, k, 5) = d(nx, j, k)*(Para%Cv*T(nx, j, k)+(u(nx, j, k)**2+v(nx, j, k)**2+w(nx, j, k)**2)*0.5d0)
            end do
            end do
        end if

        if (npy == 0 .and. Para%Iperiodic_Y == 0) then
            do k = 1, nz
            do i = 1, nx
                f(i, 1, k, 1) = d(i, 1, k)
                f(i, 1, k, 2) = d(i, 1, k)*u(i, 1, k)
                f(i, 1, k, 3) = d(i, 1, k)*v(i, 1, k)
                f(i, 1, k, 4) = d(i, 1, k)*w(i, 1, k)
                f(i, 1, k, 5) = d(i, 1, k)*(Para%Cv*T(i, 1, k)+(u(i, 1, k)**2+v(i, 1, k)**2+w(i, 1, k)**2)*0.5d0)
            end do
            end do
        end if

        if (npy == npy0-1 .and. Para%Iperiodic_Y == 0) then
            do k = 1, nz
            do i = 1, nx
                f(i, ny, k, 1) = d(i, ny, k)
                f(i, ny, k, 2) = d(i, ny, k)*u(i, ny, k)
                f(i, ny, k, 3) = d(i, ny, k)*v(i, ny, k)
                f(i, ny, k, 4) = d(i, ny, k)*w(i, ny, k)
                f(i, ny, k, 5) = d(i, ny, k)*(Para%Cv*T(i, ny, k)+(u(i, ny, k)**2+v(i, ny, k)**2+w(i, ny, k)**2)*0.5d0)
            end do
            end do
        end if

        if (npz == 0 .and. Para%Iperiodic_Z == 0) then
            do j = 1, ny
            do i = 1, nx
                f(i, j, 1, 1) = d(i, j, 1)
                f(i, j, 1, 2) = d(i, j, 1)*u(i, j, 1)
                f(i, j, 1, 3) = d(i, j, 1)*v(i, j, 1)
                f(i, j, 1, 4) = d(i, j, 1)*w(i, j, 1)
                f(i, j, 1, 5) = d(i, j, 1)*(Para%Cv*T(i, j, 1)+(u(i, j, 1)**2+v(i, j, 1)**2+w(i, j, 1)**2)*0.5d0)
            end do
            end do
        end if

        if (npz == npz0-1 .and. Para%Iperiodic_Z == 0) then
            do j = 1, ny
            do i = 1, nx
                f(i, j, nz, 1) = d(i, j, nz)
                f(i, j, nz, 2) = d(i, j, nz)*u(i, j, nz)
                f(i, j, nz, 3) = d(i, j, nz)*v(i, j, nz)
                f(i, j, nz, 4) = d(i, j, nz)*w(i, j, nz)
                f(i, j, nz, 5) = d(i, j, nz)*(Para%Cv*T(i, j, nz)+(u(i, j, nz)**2+v(i, j, nz)**2+w(i, j, nz)**2)*0.5d0)
            end do
            end do
        end if
        !
    end if

end subroutine apply_bc

subroutine apply_bc_d
    use mesh_data
    use fluc_data
    use flow_data
    implicit none

    integer :: i, j, k, IERR

    select case (Para%IBC)
    case (BC_None)
    case (BC_Flat_Plate)
        call ocfd_bc_Flat_Plate_d ! flat plate boundary
    case (BC_Lid_Driven_Cavity)
        call ocfd_bc_Lid_Driven_Cavity_d ! lid-driven cavity flow
    case (BC_User_Def)
        call ocfd_bc_User_Def_d ! user-define boundaries
    case DEFAULT
        if (my_id == 0) print *, "This boundary condition is not supported!!!"
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select

    if (npx == 0 .and. Para%Iperiodic_X == 0) then
        do k = 1, nz
        do j = 1, ny
            f(1, j, k, 1) = d_fluc(1, j, k)
            f(1, j, k, 2) = d_fluc(1, j, k)*u(1, j, k)+d(1, j, k)*u_fluc(1, j, k)
            f(1, j, k, 3) = d_fluc(1, j, k)*v(1, j, k)+d(1, j, k)*v_fluc(1, j, k)
            f(1, j, k, 4) = d_fluc(1, j, k)*w(1, j, k)+d(1, j, k)*w_fluc(1, j, k)
            f(1, j, k, 5) = d_fluc(1, j, k)*(Para%Cv*T(1, j, k)+(u(1, j, k)**2+v(1, j, k)**2+w(1, j, k)**2)*0.5d0)+ &
                            d(1, j, k)*(Para%Cv*T_fluc(1, j, k)+ &
                                        (u(1, j, k)*u_fluc(1, j, k)+ &
                                         v(1, j, k)*v_fluc(1, j, k)+ &
                                         w(1, j, k)*w_fluc(1, j, k)))
        end do
        end do
    end if

    if (npx == npx0-1 .and. Para%Iperiodic_X == 0) then
        do k = 1, nz
        do j = 1, ny
            f(nx, j, k, 1) = d_fluc(nx, j, k)
            f(nx, j, k, 2) = d_fluc(nx, j, k)*u(nx, j, k)+d(nx, j, k)*u_fluc(nx, j, k)
            f(nx, j, k, 3) = d_fluc(nx, j, k)*v(nx, j, k)+d(nx, j, k)*v_fluc(nx, j, k)
            f(nx, j, k, 4) = d_fluc(nx, j, k)*w(nx, j, k)+d(nx, j, k)*w_fluc(nx, j, k)
            f(nx, j, k, 5) = d_fluc(nx, j, k)*(Para%Cv*T(nx, j, k)+(u(nx, j, k)**2+v(nx, j, k)**2+w(nx, j, k)**2)*0.5d0)+ &
                             d(nx, j, k)*(Para%Cv*T_fluc(nx, j, k)+ &
                                          (u(nx, j, k)*u_fluc(nx, j, k)+ &
                                           v(nx, j, k)*v_fluc(nx, j, k)+ &
                                           w(nx, j, k)*w_fluc(nx, j, k)))
        end do
        end do
    end if

    if (npy == 0 .and. Para%Iperiodic_Y == 0) then
        do k = 1, nz
        do i = 1, nx
            f(i, 1, k, 1) = d_fluc(i, 1, k)
            f(i, 1, k, 2) = d_fluc(i, 1, k)*u(i, 1, k)+d(i, 1, k)*u_fluc(i, 1, k)
            f(i, 1, k, 3) = d_fluc(i, 1, k)*v(i, 1, k)+d(i, 1, k)*v_fluc(i, 1, k)
            f(i, 1, k, 4) = d_fluc(i, 1, k)*w(i, 1, k)+d(i, 1, k)*w_fluc(i, 1, k)
            f(i, 1, k, 5) = d_fluc(i, 1, k)*(Para%Cv*T(i, 1, k)+(u(i, 1, k)**2+v(i, 1, k)**2+w(i, 1, k)**2)*0.5d0)+ &
                            d(i, 1, k)*(Para%Cv*T_fluc(i, 1, k)+ &
                                        (u(i, 1, k)*u_fluc(i, 1, k)+ &
                                         v(i, 1, k)*v_fluc(i, 1, k)+ &
                                         w(i, 1, k)*w_fluc(i, 1, k)))
        end do
        end do
    end if

    if (npy == npy0-1 .and. Para%Iperiodic_Y == 0) then
        do k = 1, nz
        do i = 1, nx
            f(i, ny, k, 1) = d_fluc(i, ny, k)
            f(i, ny, k, 2) = d_fluc(i, ny, k)*u(i, ny, k)+d(i, ny, k)*u_fluc(i, ny, k)
            f(i, ny, k, 3) = d_fluc(i, ny, k)*v(i, ny, k)+d(i, ny, k)*v_fluc(i, ny, k)
            f(i, ny, k, 4) = d_fluc(i, ny, k)*w(i, ny, k)+d(i, ny, k)*w_fluc(i, ny, k)
            f(i, ny, k, 5) = d_fluc(i, ny, k)*(Para%Cv*T(i, ny, k)+(u(i, ny, k)**2+v(i, ny, k)**2+w(i, ny, k)**2)*0.5d0)+ &
                             d(i, ny, k)*(Para%Cv*T_fluc(i, ny, k)+ &
                                          (u(i, ny, k)*u_fluc(i, ny, k)+ &
                                           v(i, ny, k)*v_fluc(i, ny, k)+ &
                                           w(i, ny, k)*w_fluc(i, ny, k)))
        end do
        end do
    end if

    if (npz == 0 .and. Para%Iperiodic_Z == 0) then
        do j = 1, ny
        do i = 1, nx
            f(i, j, 1, 1) = d_fluc(i, j, 1)
            f(i, j, 1, 2) = d_fluc(i, j, 1)*u(i, j, 1)+d(i, j, 1)*u_fluc(i, j, 1)
            f(i, j, 1, 3) = d_fluc(i, j, 1)*v(i, j, 1)+d(i, j, 1)*v_fluc(i, j, 1)
            f(i, j, 1, 4) = d_fluc(i, j, 1)*w(i, j, 1)+d(i, j, 1)*w_fluc(i, j, 1)
            f(i, j, 1, 5) = d_fluc(i, j, 1)*(Para%Cv*T(i, j, 1)+(u(i, j, 1)**2+v(i, j, 1)**2+w(i, j, 1)**2)*0.5d0)+ &
                            d(i, j, 1)*(Para%Cv*T_fluc(i, j, 1)+ &
                                        (u(i, j, 1)*u_fluc(i, j, 1)+ &
                                         v(i, j, 1)*v_fluc(i, j, 1)+ &
                                         w(i, j, 1)*w_fluc(i, j, 1)))
        end do
        end do
    end if

    if (npz == npz0-1 .and. Para%Iperiodic_Z == 0) then
        do j = 1, ny
        do i = 1, nx
            f(i, j, nz, 1) = d_fluc(i, j, nz)
            f(i, j, nz, 2) = d_fluc(i, j, nz)*u(i, j, nz)+d(i, j, nz)*u_fluc(i, j, nz)
            f(i, j, nz, 3) = d_fluc(i, j, nz)*v(i, j, nz)+d(i, j, nz)*v_fluc(i, j, nz)
            f(i, j, nz, 4) = d_fluc(i, j, nz)*w(i, j, nz)+d(i, j, nz)*w_fluc(i, j, nz)
            f(i, j, nz, 5) = d_fluc(i, j, nz)*(Para%Cv*T(i, j, nz)+(u(i, j, nz)**2+v(i, j, nz)**2+w(i, j, nz)**2)*0.5d0)+ &
                             d(i, j, nz)*(Para%Cv*T_fluc(i, j, nz)+ &
                                          (u(i, j, nz)*u_fluc(i, j, nz)+ &
                                           v(i, j, nz)*v_fluc(i, j, nz)+ &
                                           w(i, j, nz)*w_fluc(i, j, nz)))
        end do
        end do
    end if

end subroutine apply_bc_d
