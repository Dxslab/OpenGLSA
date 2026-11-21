! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0
!
! Compile: gfortran -O2 -o init_vector init_vector.f90
! Run in the case directory: ./init_vector [nx ny nz [amplitude]]
! Defaults: 101 101 101 1.0d0. Dimensions must match control.nml.
! Input: grid3d.dat. Output: fluc3d.dat (replaced if it exists).
! Arnoldi: RunMode=1; set IF_init=1 in control_gsa.
! Use the same unformatted-record conventions as the solver.
! Requires Fortran 2003; no IEEE module needed.
program generate_ldc_initial_vector
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    integer :: dims(3), nx, ny, nz, i, j, k, m, unit, ios
    real(dp) :: amplitude, envelope, x, y, z, lo, length, peak
    real(dp) :: phase, mode1, mode2
    real(dp), allocatable :: xyz(:, :, :, :)
    real(dp), allocatable :: velocity(:, :, :, :), zero(:, :)
    character(len=128) :: arg

    unit = 20
    dims = 101
    amplitude = 1.0_dp
    select case (command_argument_count())
    case (0, 3, 4)
    case default
        call fail('Usage: ./init_vector [nx ny nz [amplitude]]')
    end select
    if (command_argument_count() >= 3) then
        do m = 1, 3
            call get_command_argument(m, arg)
            read (arg, *, iostat=ios) dims(m)
            if (ios /= 0) call fail('Invalid grid dimension')
        end do
    end if
    if (any(dims < 3)) call fail('Grid dimensions must be >= 3')
    if (command_argument_count() == 4) then
        call get_command_argument(4, arg)
        read (arg, *, iostat=ios) amplitude
        if (ios /= 0) call fail('Invalid amplitude')
    end if
    if (.not. (abs(amplitude) <= huge(amplitude))) then
        call fail('Amplitude must be finite')
    end if
    if (amplitude <= 0.0_dp) call fail('Amplitude must be positive')

    nx = dims(1); ny = dims(2); nz = dims(3)
    allocate (xyz(nx, ny, nz, 3), velocity(nx, ny, nz, 3), zero(nx, ny))
    open (unit=unit, file='grid3d.dat', form='unformatted', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('Cannot open grid3d.dat')
    do m = 1, 3
        do k = 1, nz
            read (unit, iostat=ios) xyz(:, :, k, m)
            if (ios /= 0) then
                call fail('Grid read failed: check dimensions/format')
            end if
        end do
    end do
    close (unit)

    ! Normalize physical coordinates, including stretched grids.
    do m = 1, 3
        if (.not. all(abs(xyz(:, :, :, m)) <= huge(1.0_dp))) then
            call fail('Nonfinite grid coordinates')
        end if
        lo = minval(xyz(:, :, :, m))
        length = maxval(xyz(:, :, :, m))-lo
        if (length <= 0.0_dp) call fail('Grid has zero extent')
        xyz(:, :, :, m) = (xyz(:, :, :, m)-lo)/length
    end do

    ! Deterministic mixed modes perturb all velocity components.
    ! This is an Arnoldi seed; no divergence-free projection is applied.
    velocity = 0.0_dp
    do k = 2, nz-1
    do j = 2, ny-1
    do i = 2, nx-1
        x = xyz(i, j, k, 1); y = xyz(i, j, k, 2); z = xyz(i, j, k, 3)
        envelope = (sin(pi*x)*sin(pi*y)*sin(pi*z))**2
        phase = 2.0_dp*pi*(x+0.37_dp*y+0.61_dp*z)
        mode1 = sin(phase)
        mode2 = 0.5_dp*cos(4.0_dp*pi*y+0.3_dp)
        velocity(i, j, k, 1) = envelope*(mode1+mode2)
        phase = 2.0_dp*pi*(0.53_dp*x+y+0.29_dp*z)+0.7_dp
        mode1 = cos(phase)
        mode2 = 0.5_dp*sin(4.0_dp*pi*z+0.2_dp)
        velocity(i, j, k, 2) = envelope*(mode1+mode2)
        phase = 2.0_dp*pi*(0.41_dp*x+0.73_dp*y+z)+1.1_dp
        mode1 = sin(phase)
        mode2 = 0.5_dp*cos(4.0_dp*pi*x+0.9_dp)
        velocity(i, j, k, 3) = envelope*(mode1+mode2)
    end do
    end do
    end do
    peak = maxval(abs(velocity))
    if (peak <= tiny(1.0_dp)) call fail('Velocity vector is zero')
    velocity = (velocity/peak)*amplitude
    zero = 0.0_dp

    open (unit=unit, file='fluc3d.dat', form='unformatted', &
          status='replace', action='write', iostat=ios)
    if (ios /= 0) call fail('Cannot create fluc3d.dat')
    write (unit) 0, 0.0_dp
    do k = 1, nz
        write (unit) zero                  ! Density perturbation
    end do
    do m = 1, 3
        do k = 1, nz
            write (unit) velocity(:, :, k, m)
        end do
    end do
    do k = 1, nz
        write (unit) zero                  ! Temperature perturbation
    end do
    close (unit)
    write (*, '(a,3(i0,1x))') 'Grid dimensions: ', dims
    peak = maxval(abs(velocity))
    write (*, '(a,es14.6)') 'Maximum velocity component: ', peak
    write (*, '(a)') 'Created fluc3d.dat.'
    write (*, '(a)') 'Density and temperature perturbations: zero.'
contains

    subroutine fail(message)
        character(len=*), intent(in) :: message
        write (*, '(a)') 'ERROR: '//message
        stop 1
    end subroutine fail

end program generate_ldc_initial_vector
