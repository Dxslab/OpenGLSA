! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine apply_bc
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, IERR

    select case (Para%IBC)
    case (BC_None)
        !> No boundary Condition (such as periodic boundary in both x- and y- )
    case (BC_Flat_Plate)
        call ocfd2d_bc_Flat_Plate ! flat plate boundary
    case (BC_Lid_Driven_Cavity)
        call ocfd2d_bc_Lid_Driven_Cavity ! lid-driven cavity flow
    case (BC_Sphere)
        call ocfd2d_bc_Sphere ! flow over sphere
    case (BC_Cylinder)
        call ocfd2d_bc_Cylinder ! flow over cylinder
    case (BC_User_Def)
        call ocfd2d_bc_User_Def ! user-define boundaries
    case DEFAULT
        if (my_id == 0) print *, "This boundary condition is not supported!!!"
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select

    if (Para%runMode == 0) then
        if (npx == 0 .and. Para%Iperiodic_X == 0) then
            do j = 1, ny
                f(1, j, 1) = d(1, j)
                f(1, j, 2) = d(1, j)*u(1, j)
                f(1, j, 3) = d(1, j)*v(1, j)
                f(1, j, 4) = d(1, j)*(Para%Cv*T(1, j)+(u(1, j)*u(1, j)+v(1, j)*v(1, j))*0.5d0)
            end do
        end if

        if (npx == npx0-1 .and. Para%Iperiodic_X == 0) then
            do j = 1, ny
                f(nx, j, 1) = d(nx, j)
                f(nx, j, 2) = d(nx, j)*u(nx, j)
                f(nx, j, 3) = d(nx, j)*v(nx, j)
                f(nx, j, 4) = d(nx, j)*(Para%Cv*T(nx, j)+(u(nx, j)*u(nx, j)+v(nx, j)*v(nx, j))*0.5d0)
            end do
        end if

        if (npy == 0 .and. Para%Iperiodic_Y == 0) then
            do i = 1, nx
                f(i, 1, 1) = d(i, 1)
                f(i, 1, 2) = d(i, 1)*u(i, 1)
                f(i, 1, 3) = d(i, 1)*v(i, 1)
                f(i, 1, 4) = d(i, 1)*(Para%Cv*T(i, 1)+(u(i, 1)*u(i, 1)+v(i, 1)*v(i, 1))*0.5d0)
            end do
        end if

        if (npy == npy0-1 .and. Para%Iperiodic_Y == 0) then
            do i = 1, nx
                f(i, ny, 1) = d(i, ny)
                f(i, ny, 2) = d(i, ny)*u(i, ny)
                f(i, ny, 3) = d(i, ny)*v(i, ny)
                f(i, ny, 4) = d(i, ny)*(Para%Cv*T(i, ny)+(u(i, ny)*u(i, ny)+v(i, ny)*v(i, ny))*0.5d0)
            end do
        end if
    end if

end subroutine apply_bc

subroutine apply_bc_d
    use mesh_data
    use fluc_data
    use flow_data
    implicit none

    integer :: i, j, IERR

    select case (Para%IBC)
    case (BC_None)
        !> No boundary Condition (such as periodic boundary in both x- and y- )
    case (BC_Flat_Plate)
        call ocfd2d_bc_Flat_Plate_d ! flat plate boundary
    case (BC_Lid_Driven_Cavity)
        call ocfd2d_bc_Lid_Driven_Cavity_d ! lid-driven cavity flow
    case (BC_Sphere)
        call ocfd2d_bc_Sphere_d ! flow over sphere
    case (BC_Cylinder)
        call ocfd2d_bc_Cylinder_d ! flow over cylinder
    case (BC_User_Def)
        call ocfd2d_bc_User_Def_d ! user-define boundaries
    case DEFAULT
        if (my_id == 0) print *, "This boundary condition is not supported!!!"
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select

    if (npx == 0 .and. Para%Iperiodic_X == 0) then
        do j = 1, ny
            f(1, j, 1) = d_fluc(1, j)
            f(1, j, 2) = d_fluc(1, j)*u(1, j)+d(1, j)*u_fluc(1, j)
            f(1, j, 3) = d_fluc(1, j)*v(1, j)+d(1, j)*v_fluc(1, j)
            f(1, j, 4) = d(1, j)*w_fluc(1, j)
            f(1, j, 5) = d(1, j)*(Para%Cv*T_fluc(1, j)+u(1, j)*u_fluc(1, j)+v(1, j)*v_fluc(1, j))+ &
                         d_fluc(1, j)*(Para%Cv*T(1, j)+(u(1, j)**2+v(1, j)**2)*0.5d0)
        end do
    end if

    if (npx == npx0-1 .and. Para%Iperiodic_X == 0) then
        do j = 1, ny
            f(nx, j, 1) = d_fluc(nx, j)
            f(nx, j, 2) = d_fluc(nx, j)*u(nx, j)+d(nx, j)*u_fluc(nx, j)
            f(nx, j, 3) = d_fluc(nx, j)*v(nx, j)+d(nx, j)*v_fluc(nx, j)
            f(nx, j, 4) = d(nx, j)*w_fluc(nx, j)
            f(nx, j, 5) = d(nx, j)*(Para%Cv*T_fluc(nx, j)+u(nx, j)*u_fluc(nx, j)+v(nx, j)*v_fluc(nx, j))+ &
                          d_fluc(nx, j)*(Para%Cv*T(nx, j)+(u(nx, j)**2+v(nx, j)**2)*0.5d0)
        end do
    end if

    if (npy == 0 .and. Para%Iperiodic_Y == 0) then
        do i = 1, nx
            f(i, 1, 1) = d_fluc(i, 1)
            f(i, 1, 2) = d_fluc(i, 1)*u(i, 1)+d(i, 1)*u_fluc(i, 1)
            f(i, 1, 3) = d_fluc(i, 1)*v(i, 1)+d(i, 1)*v_fluc(i, 1)
            f(i, 1, 4) = d(i, 1)*w_fluc(i, 1)
            f(i, 1, 5) = d(i, 1)*(Para%Cv*T_fluc(i, 1)+u(i, 1)*u_fluc(i, 1)+v(i, 1)*v_fluc(i, 1))+ &
                         d_fluc(i, 1)*(Para%Cv*T(i, 1)+(u(i, 1)**2+v(i, 1)**2)*0.5d0)
        end do
    end if

    if (npy == npy0-1 .and. Para%Iperiodic_Y == 0) then
        do i = 1, nx
            f(i, ny, 1) = d_fluc(i, ny)
            f(i, ny, 2) = d_fluc(i, ny)*u(i, ny)+d(i, ny)*u_fluc(i, ny)
            f(i, ny, 3) = d_fluc(i, ny)*v(i, ny)+d(i, ny)*v_fluc(i, ny)
            f(i, ny, 4) = d(i, ny)*w_fluc(i, ny)
            f(i, ny, 5) = d(i, ny)*(Para%Cv*T_fluc(i, ny)+u(i, ny)*u_fluc(i, ny)+v(i, ny)*v_fluc(i, ny))+ &
                          d_fluc(i, ny)*(Para%Cv*T(i, ny)+(u(i, ny)**2+v(i, ny)**2)*0.5d0)
        end do
    end if

end subroutine apply_bc_d
