! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! CopyRight by Li Xinliang  Email: lixl@imech.ac.cn
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine NS_Solver
    use mesh_data
    use flow_data
    use Flux_Viscous
    use SFD
    use IO
    implicit none

    real(dp) :: wall_time, mass0
    integer :: KRK

    call init_flow
    call apply_bc

    if (my_id == 0) print *, 'Init OK'
    wall_time = MPI_wtime()

    if (Para%IF_SFD == 1) fsfd = f

    call get_CFL

    ! Closed cavity: the wall closure is not mass conservative, fix the total mass to its initial value
    if (Para%IBC == BC_Lid_Driven_Cavity) call total_mass(mass0)

    do while (tt < Para%End_time-Para%dt*1.d-4)
        !
        fn = f
        du = 0.0d0
        !
        do KRK = 1, 3
            !
            call mpi_exchange_xyz3d(d)
            call mpi_exchange_xyz3d(u)
            call mpi_exchange_xyz3d(v)
            call mpi_exchange_xyz3d(w)
            call mpi_exchange_xyz3d(T)
            !
            if (Para%IF_Character == 1) then
                call flux_inviscid_character
            else if (Para%IF_Character == 2) then
                call flux_inviscid_character_pp
            else
                call flux_inviscid
            end if
            !
            if (Para%IF_Viscous == 1) then
                call flux_viscous3d
            end if
            !
            ! CALL apply_sponge_layer
            !
            call RK3(KRK, f, fn, du, 5)
            !
            call compute_duvwTp
            call apply_bc
            !
        end do
        !
        Istep = Istep+1
        tt = tt+Para%dt
        !
        if (Para%IBC == BC_Lid_Driven_Cavity) call fix_total_mass(mass0)
        !
        ! NS residual (f_RK - fn)/dt, evaluated before SFD filtering
        if (mod(Istep, Para%Istep_Show) == 0) call compute_Residual
        !
        if (Para%IF_SFD == 1) then
            call apply_SFD
            call compute_duvwTp
            call apply_bc
        end if
        if (mod(Istep, Para%Istep_Show*10) == 0) then
            call show_CPU_time(wall_time)
            call get_CFL
        end if
        !
        if (mod(Istep, Para%Istep_Save) == 0) call write_flow
        !
    end do
    !
    call write_flow
    !
    if (my_id == 0) print *, "Simulation completed. XD"
    !
end subroutine NS_Solver

subroutine compute_duvwTp
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, k, ierr
    real(dp) :: tmp0

    tmp0 = 1.0d0/(Para%gamma*Para%Ma*Para%Ma)

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
        d(i, j, k) = f(i, j, k, 1)
        u(i, j, k) = f(i, j, k, 2)/d(i, j, k)
        v(i, j, k) = f(i, j, k, 3)/d(i, j, k)
        w(i, j, k) = f(i, j, k, 4)/d(i, j, k)
        T(i, j, k) = (f(i, j, k, 5)-0.5d0*d(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2))/(d(i, j, k)*Para%Cv)
        p(i, j, k) = d(i, j, k)*T(i, j, k)*tmp0

        if (ISNAN(d(i, j, k)) .or. ISNAN(u(i, j, k)) .or. ISNAN(v(i, j, k)) .or. ISNAN(w(i, j, k)) .or. ISNAN(T(i, j, k))) then
            print *, "*** Error: NaN detected in compute_duvwTp! ***"
            print *, " x=", Axx(i, j, k), " y=", Ayy(i, j, k), " z=", Azz(i, j, k)
            print *, "*** Error: NaN detected in compute_duvwTp! ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        if (T(i, j, k) < 0.0d0) then
            print *, "*** Error: Negative temperature detected ***"
            print *, " x=", Axx(i, j, k), " y=", Ayy(i, j, k), " z=", Azz(i, j, k)
            print *, "  f(1)=rho     =", f(i, j, k, 1)
            print *, "  f(2)=rho*u   =", f(i, j, k, 2)
            print *, "  f(3)=rho*v   =", f(i, j, k, 3)
            print *, "  f(4)=rho*w   =", f(i, j, k, 4)
            print *, "  f(5)=rho*E   =", f(i, j, k, 5)
            print *, "  d=", d(i, j, k)
            print *, "  u=", u(i, j, k)
            print *, "  v=", v(i, j, k)
            print *, "  w=", w(i, j, k)
            print *, "  KE = 0.5*rho*(u^2+v^2+w^2) =", 0.5d0*d(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)
            print *, "  Internal energy = f(5)-KE =", f(i, j, k, 5)-0.5d0*d(i, j, k)*(u(i, j, k)**2+v(i, j, k)**2+w(i, j, k)**2)
            print *, "  T=", T(i, j, k)
            print *, "  P=", p(i, j, k)
            print *, "*** Error: Negative temperature detected ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
    end do
    end do
    end do

end subroutine compute_duvwTp

subroutine compute_Residual
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, k, m, ierr
    real(dp) :: resn2(5), resn2s(5), resn2total
    real(dp) :: resid

    resn2 = 0.0d0

    do m = 1, 5
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            if ((i == 1 .and. BcI(j, k, 1) == 1) .or. &
                (i == nx .and. BcI(j, k, 2) == 1) .or. &
                (j == 1 .and. BcJ(i, k, 1) == 1) .or. &
                (j == ny .and. BcJ(i, k, 2) == 1) .or. &
                (k == 1 .and. BcK(i, j, 1) == 1) .or. &
                (k == nz .and. BcK(i, j, 2) == 1)) cycle
            resid = abs(f(i, j, k, m)-fn(i, j, k, m))/Para%dt
            resn2(m) = resn2(m)+resid
        end do
        end do
        end do
    end do

    call MPI_REDUCE(resn2, resn2s, 5, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (my_id == 0) then
        do m = 1, 5
            resn2s(m) = resn2s(m)/real(nx_global*ny_global*nz_global, dp)
        end do
        resn2total = sum(resn2s)/5.0d0

        write (*, '(A12,6ES14.6)') "Residual: ", resn2s(1), resn2s(2), resn2s(3), resn2s(4), resn2s(5), resn2total
        open (78, FILE='resid.history', POSITION='APPEND')
        write (78, '(7ES14.6)') tt, resn2s(1), resn2s(2), resn2s(3), resn2s(4), resn2s(5), resn2total
        close (78)
    end if

end subroutine compute_Residual

subroutine show_CPU_time(wall_time)
    use precision
    use flow_para
    use flow_data, only: limiter_count_rho, limiter_count_pressure, limiter_count_total
    implicit none

    real(dp):: wall_time, wtmp
    integer :: limiter_counts(3), limiter_counts_global(3), ierr

    if (Para%IF_Character == 2 .and. Para%RunMode == 0) then
        limiter_counts(1) = limiter_count_rho
        limiter_counts(2) = limiter_count_pressure
        limiter_counts(3) = limiter_count_total
        call MPI_REDUCE(limiter_counts, limiter_counts_global, 3, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if

    wtmp = wall_time
    wall_time = MPI_wtime()

    if (my_id == 0) then
        if (Para%IF_Character == 2 .and. Para%RunMode == 0 .and. Para%pp_limiter /= 0) then
            write (*, '(A,I0,A,I0)') 'Limiter: rho=', limiter_counts_global(1), ' | p=', limiter_counts_global(2)
        end if
        write (*, '(A)') '--------------------------------------------------------------------------------'
        write (*, '(A,I0,A,F10.4,A,F10.6)') 'Iteration: ', Istep, ' | Simulation time: ', tt, &
            ' | Time/iteration: ', (wall_time-wtmp)/(real(Para%Istep_show*10, dp))
        write (*, '(A)') '--------------------------------------------------------------------------------'
    end if

end subroutine show_CPU_time

subroutine get_CFL
    use mesh_data
    use flow_data
    implicit none

    integer :: i, j, k, ierr

    real(dp) :: freq_local, freq_global
    real(dp) :: wave_xi, wave_eta, wave_zeta
    real(dp) :: wave_xi_vis, wave_eta_vis, wave_zeta_vis
    real(dp) :: mu_local, rho_local, nu_max
    real(dp) :: uvel, vvel, wvel, ssi, ssj, ssk
    real(dp) :: dt_cfl, cfl

    freq_local = 0.0d0

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
        uvel = u(i, j, k)*Akx(i, j, k)+v(i, j, k)*Aky(i, j, k)+w(i, j, k)*Akz(i, j, k)
        ssi = sqrt(Akx(i, j, k)**2+Aky(i, j, k)**2+Akz(i, j, k)**2)
        wave_xi = abs(uvel)+c(i, j, k)*ssi

        vvel = u(i, j, k)*Aix(i, j, k)+v(i, j, k)*Aiy(i, j, k)+w(i, j, k)*Aiz(i, j, k)
        ssj = sqrt(Aix(i, j, k)**2+Aiy(i, j, k)**2+Aiz(i, j, k)**2)
        wave_eta = abs(vvel)+c(i, j, k)*ssj

        wvel = u(i, j, k)*Asx(i, j, k)+v(i, j, k)*Asy(i, j, k)+w(i, j, k)*Asz(i, j, k)
        ssk = sqrt(Asx(i, j, k)**2+Asy(i, j, k)**2+Asz(i, j, k)**2)
        wave_zeta = abs(wvel)+c(i, j, k)*ssk

        if (Para%IF_Viscous == 1) then
            mu_local = mu(i, j, k)*Para%Re
            rho_local = d(i, j, k)
            nu_max = max((4.0d0/3.0d0), Para%gamma/Para%Pr)*(mu_local/rho_local)

            wave_xi_vis = nu_max*ssi**2
            wave_eta_vis = nu_max*ssj**2
            wave_zeta_vis = nu_max*ssk**2
        else
            wave_xi_vis = 0.0d0
            wave_eta_vis = 0.0d0
            wave_zeta_vis = 0.0d0
        end if

        freq_local = max(freq_local, wave_xi+wave_eta+wave_zeta+2.0d0*(wave_xi_vis+wave_eta_vis+wave_zeta_vis)/Para%Re)
    end do
    end do
    end do

    call MPI_ALLREDUCE(freq_local, freq_global, 1, OCFD_DATA_TYPE, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (freq_global > 1.0d-14) then
        dt_cfl = 1.0d0/freq_global
        cfl = Para%dt*freq_global
    else
        dt_cfl = 9999.0d0
        cfl = 0.0d0
    end if

    if (my_id == 0) then
        write (*, '(A)') ' '
        write (*, '(A,F10.6,A,ES14.6)') 'Current CFL number: ', cfl, ' | Time step size for cfl = 1.0: ', dt_cfl
        write (*, '(A)') ' '
    end if

end subroutine get_CFL

subroutine total_mass(mass)
    ! Total mass sum(rho/Ajac) over all nodes, including the boundary nodes
    use mesh_data
    use flow_data
    implicit none
    real(dp), intent(out) :: mass
    real(dp) :: mass_loc
    integer :: i, j, k, ierr

    mass_loc = 0.0d0
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
        mass_loc = mass_loc+d(i, j, k)/Ajac(i, j, k)
    end do
    end do
    end do
    call MPI_ALLREDUCE(mass_loc, mass, 1, OCFD_DATA_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine total_mass

subroutine fix_total_mass(mass0)
    ! Rescale the density (u, v, w, T unchanged) so that the total mass stays at mass0
    use mesh_data
    use flow_data
    implicit none
    real(dp), intent(in) :: mass0
    real(dp) :: mass, ratio

    call total_mass(mass)
    ratio = mass0/mass
    f(1:nx, 1:ny, 1:nz, :) = f(1:nx, 1:ny, 1:nz, :)*ratio
    d(1:nx, 1:ny, 1:nz) = d(1:nx, 1:ny, 1:nz)*ratio
    p(1:nx, 1:ny, 1:nz) = p(1:nx, 1:ny, 1:nz)*ratio

end subroutine fix_total_mass
