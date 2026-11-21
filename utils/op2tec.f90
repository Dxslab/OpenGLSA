! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
program op2tec
    use iso_fortran_env, only: real64, error_unit
    implicit none

    integer :: nx, ny, nx_out, Istep, nvar, i, j, k, ii, ios, unit
    integer :: slash, dot
    real(real64) :: tt
    real(real64), allocatable :: xx(:, :), yy(:, :), fields(:, :, :)
    character(len=1024) :: data_filename, output_filename, tmp, base
    character(len=128) :: var_list
    logical :: is_5var, omesh

    if (command_argument_count() < 3 .or. command_argument_count() > 4) then
        write (error_unit, '(A)') 'Usage: op2tec <data_file> <nx> <ny> [omesh]'
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
        var_list = '"X" "Y" "Density" "U" "V" "W" "Temperature"'
    else
        nvar = 4
        var_list = '"X" "Y" "Density" "U" "V" "Temperature"'
    end if

    ! Replace the final extension; never overwrite the binary input itself.
    dot = scan(trim(data_filename), '.', back=.true.)
    if (dot > slash+1) then
        output_filename = data_filename(:dot-1)//'.tec'
    else
        output_filename = trim(data_filename)//'.tec'
    end if
    if (output_filename == data_filename) output_filename = trim(data_filename)//'.tec'

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

    nx_out = nx
    if (omesh) nx_out = nx+1
    open (newunit=unit, file=trim(output_filename), form='formatted', status='replace', action='write', iostat=ios)
    call check_io(ios, 'Cannot create '//trim(output_filename))
    write (unit, '(A)') 'TITLE = "OpenGLSA solution"'
    write (unit, '(A)') 'VARIABLES = '//trim(var_list)
    write (unit, '(A,I0,A,I0,A,I0,A,ES25.16E3)') &
        'ZONE T="Istep=', Istep, '", I=', nx_out, ', J=', ny, &
        ', K=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=', tt
    do j = 1, ny
        do i = 1, nx_out
            ii = i
            if (i > nx) ii = 1
            write (unit, '(*(ES25.16E3,1X))', iostat=ios) xx(ii, j), yy(ii, j), fields(ii, j, :)
            call check_io(ios, 'Cannot write '//trim(output_filename))
        end do
    end do
    close (unit, iostat=ios)
    call check_io(ios, 'Cannot close '//trim(output_filename))
    print *, 'Wrote Tecplot ASCII: ', trim(output_filename)
    print *, 'Istep, tt=', Istep, tt

contains

    subroutine check_io(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message
        if (status /= 0) then
            write (error_unit, '(A)') trim(message)
            stop 1
        end if
    end subroutine check_io
end program op2tec
