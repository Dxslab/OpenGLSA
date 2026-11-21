! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

subroutine LNS2d_Solver
    use mesh_data
    use flow_data
    use fluc_data
    use Flux_Viscous_d
    use Flux_Viscous_cyl_d
    use IO
    implicit none

    integer :: KRK
    real(dp) :: wall_time
    
    if (Para%runMode == 2) then
        call init_flow
    end if

    call apply_bc_d
    if (Para%IBC == BC_Lid_Driven_Cavity) call remove_mean_density_fluc

    wall_time = MPI_Wtime()
    do while (tt < Para%End_time-Para%dt*1.d-4)

        fn = f
        du = 0.0d0

        do KRK = 1, 3
            !
            call mpi_exchange_xy2d(d_fluc)
            if (Para%IBC /= BC_Lid_Driven_Cavity) then
                call mpi_exchange_xy2d(u_fluc)
                call mpi_exchange_xy2d(v_fluc)
                call mpi_exchange_xy2d(w_fluc)
                call mpi_exchange_xy2d(T_fluc)
            end if
            !
            call flux_inviscid2d_d
            !
            if (Para%IF_Viscous == 1 .and. Para%IF_Cylindrical_Coordinate == 0) then
                call flux_viscous2d_d
            else if (Para%IF_Viscous == 1 .and. Para%IF_Cylindrical_Coordinate == 1) then
                call flux_viscous2d_cyl_d
            end if
            !
            call RK3(KRK, f, fn, du, 5)
            !
            call compute_duvwT
            call apply_bc_d
            if (Para%IBC == BC_Lid_Driven_Cavity) call remove_mean_density_fluc
            !
        end do  ! KRK

        Istep = Istep+1
        tt = tt+Para%dt

        if (mod(Istep, Para%Istep_Show) == 0) call compute_Residual_d
        if (mod(Istep, Para%Istep_Show*10) == 0) call show_CPU_time(wall_time)
        if (mod(Istep, Para%Istep_Save) == 0 .and. Para%runMode == 2) call write_flow

    end do

    if (Para%runMode == 2) then
        call write_flow
    end if

end subroutine LNS2d_Solver

subroutine compute_duvwT
    use mesh_data
    use flow_data
    use fluc_data
    implicit none

    integer :: i, j, ierr

    do j = 1, ny
    do i = 1, nx
        d_fluc(i, j) = f(i, j, 1)
        u_fluc(i, j) = (f(i, j, 2)-d_fluc(i, j)*u(i, j))/d(i, j)
        v_fluc(i, j) = (f(i, j, 3)-d_fluc(i, j)*v(i, j))/d(i, j)
        w_fluc(i, j) = f(i, j, 4)/d(i, j)
        T_fluc(i, j) = (f(i, j, 5)-d_fluc(i, j)*(Para%Cv*T(i, j)+(u(i, j)**2+v(i, j)**2)*0.5d0)- &
                        d(i, j)*u(i, j)*u_fluc(i, j)-d(i, j)*v(i, j)*v_fluc(i, j))/(d(i, j)*Para%Cv)

   if (ISNAN(d_fluc(i, j)) .or. ISNAN(u_fluc(i, j)) .or. ISNAN(v_fluc(i, j)) .or. ISNAN(w_fluc(i, j)) .or. ISNAN(T_fluc(i, j))) then
            print *, "*** Error: NAN encountered ***"
            print *, "x= ", Axx(i, j), " y= ", Ayy(i, j)
            print *, "d= ", d_fluc(i, j)
            print *, "u= ", u_fluc(i, j)
            print *, "v= ", v_fluc(i, j)
            print *, "w= ", w_fluc(i, j)
            print *, "T= ", T_fluc(i, j)
            print *, "*** Error: NAN encountered ***"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if (abs(d_fluc(i, j)/d(i, j)) > 10.d0 .and. Para%runMode == 1) then
            print *, "*** Error: DIVERGED! Try smaller time step ***"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end do
    end do

end subroutine compute_duvwT

subroutine compute_Residual_d
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, m, ierr
    real(dp) :: resn2(5), resn2s(5), resn2total
    real(dp) :: resid

    resn2 = 0.0d0

    do m = 1, 5
        do j = 1, ny
        do i = 1, nx
            if ((i == 1 .and. BcI(j, 1) == 1) .or. &
                (i == nx .and. BcI(j, 2) == 1) .or. &
                (j == 1 .and. BcJ(i, 1) == 1) .or. &
                (j == ny .and. BcJ(i, 2) == 1)) cycle
            resid = abs(f(i, j, m)-fn(i, j, m))/Para%dt
            resn2(m) = resn2(m)+resid**2
        end do
        end do
    end do

    call MPI_REDUCE(resn2, resn2s, 5, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (my_id == 0) then
        do m = 1, 5
            resn2s(m) = sqrt(resn2s(m)/real(nx_global*ny_global, dp))
        end do
        resn2total = sqrt(sum(resn2s**2)/5.0d0)

        write (*, '(A12,6ES14.6)') "Residual: ", resn2s(1), resn2s(2), resn2s(3), resn2s(4), resn2s(5), resn2total
        if (Para%runMode == 2) then
            open (78, FILE='resid.history', POSITION='APPEND')
            write (78, '(7ES14.6)') tt, resn2s(1), resn2s(2), resn2s(3), resn2s(4), resn2s(5), resn2total
            close (78)
        end if
    end if

end subroutine compute_Residual_d

subroutine remove_mean_density_fluc
    ! Closed cavity: keep sum(rho'/Ajac) = 0 by removing the uniform-density (mass) mode
    ! rho' = c, u' = v' = w' = T' = 0, which otherwise appears as a spurious near-zero eigenvalue
    use mesh_data
    use flow_data
    use fluc_data
    implicit none
    real(dp) :: sums_loc(2), sums(2), rho_mean
    integer :: i, j, ierr

    sums_loc = 0.0d0
    do j = 1, ny
    do i = 1, nx
        sums_loc(1) = sums_loc(1)+d_fluc(i, j)/Ajac(i, j)
        sums_loc(2) = sums_loc(2)+1.0d0/Ajac(i, j)
    end do
    end do
    call MPI_ALLREDUCE(sums_loc, sums, 2, OCFD_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    rho_mean = sums(1)/sums(2)

    do j = 1, ny
    do i = 1, nx
        d_fluc(i, j) = d_fluc(i, j)-rho_mean
        f(i, j, 1) = f(i, j, 1)-rho_mean
        f(i, j, 2) = f(i, j, 2)-rho_mean*u(i, j)
        f(i, j, 3) = f(i, j, 3)-rho_mean*v(i, j)
        f(i, j, 5) = f(i, j, 5)-rho_mean*(Para%Cv*T(i, j)+0.5d0*(u(i, j)**2+v(i, j)**2))
    end do
    end do

end subroutine remove_mean_density_fluc
