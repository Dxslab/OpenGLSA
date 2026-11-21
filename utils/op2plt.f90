! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
program op2plt
    use iso_fortran_env, only: error_unit, int64
    use iso_c_binding, only: c_int32_t, c_double, c_null_char
    implicit none

    include "tecio.f90"

    integer :: nx, ny, nx_out, Istep, nvar, k, ios, unit
    integer :: slash, dot
    real(c_double) :: tt
    real(c_double), allocatable :: xx(:, :), yy(:, :), fields(:, :, :)
    character(len=1024) :: data_filename, output_filename, tmp, base
    character(len=128) :: var_list
    ! TecIO 142 uses 32-bit counts and nodal, unshared variables.
    integer(c_int32_t) :: count, imax, jmax, kmax, status
    integer(c_int32_t) :: passive(8), location(8), sharing(8)
    integer(int64) :: points
    real(c_double), allocatable :: block(:, :)
    logical :: is_5var, omesh

    if (command_argument_count() < 3 .or. command_argument_count() > 4) then
        write (error_unit, '(A)') 'Usage: op2plt <data_file> <nx> <ny> [omesh]'
        stop 1
    end if
    call get_command_argument(1, data_filename)
    call get_command_argument(2, tmp)
    read (tmp, *, iostat=ios) nx
    call check_io(ios, 'nx must be a positive integer')
    call get_command_argument(3, tmp)
    read (tmp, *, iostat=ios) ny
    call check_io(ios, 'ny must be a positive integer')
    if (nx <= 0 .or. ny <= 0) then
        write (error_unit, '(A)') 'nx and ny must be positive'
        stop 1
    end if

    omesh = .false.
    if (command_argument_count() == 4) then
        call get_command_argument(4, tmp)
        if (trim(tmp) /= 'omesh' .and. trim(tmp) /= 'OMESH') then
            write (error_unit, '(A)') 'Unknown option: '//trim(tmp)
            stop 1
        end if
        omesh = .true.
    end if

    slash = scan(trim(data_filename), '/', back=.true.)
    base = data_filename(slash+1:)
    is_5var = (base(1:4) == 'fluc' .or. base(1:4) == 'FLUC' &
               .or. base(1:4) == 'Mode' .or. base(1:4) == 'MODE')
    if (is_5var) then
        nvar = 5
        var_list = 'X Y Density U V W Temperature'
    else
        nvar = 4
        var_list = 'X Y Density U V Temperature'
    end if

    ! Replace the final extension; never overwrite the binary input itself.
    dot = scan(trim(data_filename), '.', back=.true.)
    if (dot > slash+1) then
        output_filename = data_filename(:dot-1)//'.plt'
    else
        output_filename = trim(data_filename)//'.plt'
    end if
    if (output_filename == data_filename) output_filename = trim(data_filename)//'.plt'

    nx_out = nx
    if (omesh) then
        if (nx == huge(nx)) then
            write (error_unit, '(A)') 'nx is too large to close the O-mesh'
            stop 1
        end if
        nx_out = nx+1
    end if
    points = int(nx_out, int64)*int(ny, int64)
    if (points > int(huge(count), int64)) then
        write (error_unit, '(A)') 'Grid exceeds the TecIO 142 point-count limit'
        stop 1
    end if

    allocate (xx(nx, ny), yy(nx, ny), fields(nx, ny, nvar))
    open (newunit=unit, file='grid.dat', form='unformatted', status='old', action='read', iostat=ios)
    call check_io(ios, 'Cannot open grid.dat')
    read (unit, iostat=ios) xx
    call check_io(ios, 'Cannot read X coordinates from grid.dat')
    read (unit, iostat=ios) yy
    call check_io(ios, 'Cannot read Y coordinates from grid.dat')
    close (unit)

    open (newunit=unit, file=trim(data_filename), form='unformatted', status='old', action='read', iostat=ios)
    call check_io(ios, 'Cannot open '//trim(data_filename))
    read (unit, iostat=ios) Istep, tt
    call check_io(ios, 'Cannot read snapshot header')
    do k = 1, nvar
        read (unit, iostat=ios) fields(:, :, k)
        call check_io(ios, 'Cannot read snapshot field; check dimensions and filename prefix')
    end do
    close (unit)

    count = int(points, c_int32_t)
    imax = nx_out; jmax = ny; kmax = 1
    passive = 0; location = 1; sharing = 0
    call check_io(TECINI142('OpenGLSA solution'//c_null_char, trim(var_list)//c_null_char, &
        trim(output_filename)//c_null_char, '.'//c_null_char, &
        0_c_int32_t, 0_c_int32_t, 0_c_int32_t, 1_c_int32_t), 'TECINI142: cannot create output')
    status = TECZNE142('OpenGLSA solution'//c_null_char, &
        0_c_int32_t, imax, jmax, kmax, 0_c_int32_t, 0_c_int32_t, 0_c_int32_t, &
        tt, 1_c_int32_t, 0_c_int32_t, 1_c_int32_t, 0_c_int32_t, 0_c_int32_t, &
        0_c_int32_t, 0_c_int32_t, 0_c_int32_t, passive, location, sharing, 0_c_int32_t)
    call check_io(status, 'TECZNE142: cannot create zone')
    allocate (block(nx_out, ny))
    block(1:nx, :) = xx
    if (omesh) block(nx_out, :) = xx(1, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write xx')
    block(1:nx, :) = yy
    if (omesh) block(nx_out, :) = yy(1, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write yy')
    do k = 1, nvar
        block(1:nx, :) = fields(:, :, k)
        if (omesh) block(nx_out, :) = fields(1, :, k)
        call check_io(TECDATD142(count, block), 'TECDATD142: cannot write solution field')
    end do
    call check_io(TECEND142(), 'TECEND142: cannot finish output')
    print *, 'Wrote Tecplot binary: ', trim(output_filename)

contains

    subroutine check_io(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message
        if (status /= 0) then
            write (error_unit, '(A)') trim(message)
            stop 1
        end if
    end subroutine check_io
end program op2plt
