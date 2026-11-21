! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

program op2tec3d
    use iso_fortran_env, only: error_unit
    implicit none

    integer :: nx, ny, nz, Istep, i, j, k, ii, unit, ios, slash, dot
    real(kind=8) :: tt
    real(kind=8), dimension(:, :, :), allocatable :: xx, yy, zz, d, u, v, w, T

    character(len=1024) :: data_filename, output_filename, tmp
    character(len=64)  :: var_list
    logical :: omesh
    integer :: nx_out

    if (command_argument_count() < 4 .or. command_argument_count() > 5) then
        print *, "Usage: op2tec3d <data_file> <nx> <ny> <nz> [omesh]"
        print *, "  All 3D files (flow3d/fluc3d/Mode) are read as 5-field data (d,u,v,w,T)."
        stop 1
    end if
    call get_command_argument(1, data_filename)
    call get_command_argument(2, tmp)
    read (tmp, *, iostat=ios) nx
    call check_io(ios, 'nx must be a positive integer')
    call get_command_argument(3, tmp)
    read (tmp, *, iostat=ios) ny
    call check_io(ios, 'ny must be a positive integer')
    call get_command_argument(4, tmp)
    read (tmp, *, iostat=ios) nz
    call check_io(ios, 'nz must be a positive integer')
    if (nx <= 0 .or. ny <= 0 .or. nz <= 0) then
        write (error_unit, '(A)') 'nx, ny and nz must be positive'
        stop 1
    end if

    omesh = .false.
    if (command_argument_count() >= 5) then
        call get_command_argument(5, tmp)
        if (trim(tmp) /= 'omesh' .and. trim(tmp) /= 'OMESH') then
            write (error_unit, '(A)') 'Unknown option: '//trim(tmp)
            stop 1
        end if
        omesh = .true.
    end if

    var_list = '"X" "Y" "Z" "Density" "U" "V" "W" "Temperature"'

    ! Replace the final extension; never overwrite the binary input itself.
    slash = scan(trim(data_filename), '/', back=.true.)
    dot = scan(trim(data_filename), '.', back=.true.)
    if (dot > slash+1) then
        output_filename = data_filename(:dot-1)//'.tec'
    else
        output_filename = trim(data_filename)//'.tec'
    end if
    if (output_filename == data_filename) output_filename = trim(data_filename)//'.tec'

    allocate (xx(nx, ny, nz), yy(nx, ny, nz), zz(nx, ny, nz))
    allocate (d(nx, ny, nz), u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz), T(nx, ny, nz))

    open (55, file='grid3d.dat', form='unformatted', status='old', action='read', iostat=ios)
    call check_io(ios, 'Cannot open grid3d.dat')
    do k = 1, nz
        read (55, iostat=ios) xx(:, :, k)
        call check_io(ios, 'Cannot read xx coordinates from grid3d.dat')
    end do
    do k = 1, nz
        read (55, iostat=ios) yy(:, :, k)
        call check_io(ios, 'Cannot read yy coordinates from grid3d.dat')
    end do
    do k = 1, nz
        read (55, iostat=ios) zz(:, :, k)
        call check_io(ios, 'Cannot read zz coordinates from grid3d.dat')
    end do
    close (55)

    open (44, file=trim(data_filename), form='unformatted', status='old', action='read', iostat=ios)
    call check_io(ios, 'Cannot open '//trim(data_filename))
    read (44, iostat=ios) Istep, tt
    call check_io(ios, 'Cannot read snapshot header')
    do k = 1, nz
        read (44, iostat=ios) d(:, :, k)
        call check_io(ios, 'Cannot read snapshot field d; check dimensions')
    end do
    do k = 1, nz
        read (44, iostat=ios) u(:, :, k)
        call check_io(ios, 'Cannot read snapshot field u; check dimensions')
    end do
    do k = 1, nz
        read (44, iostat=ios) v(:, :, k)
        call check_io(ios, 'Cannot read snapshot field v; check dimensions')
    end do
    do k = 1, nz
        read (44, iostat=ios) w(:, :, k)
        call check_io(ios, 'Cannot read snapshot field w; check dimensions')
    end do
    do k = 1, nz
        read (44, iostat=ios) T(:, :, k)
        call check_io(ios, 'Cannot read snapshot field T; check dimensions')
    end do
    close (44)
    print *, 'Istep, tt=', Istep, tt

    nx_out = nx
    if (omesh) nx_out = nx+1
    open (newunit=unit, file=trim(output_filename), form='formatted', status='replace', action='write', iostat=ios)
    call check_io(ios, 'Cannot create '//trim(output_filename))
    write (unit, '(A)') 'TITLE = "OpenGLSA solution"'
    write (unit, '(A)') 'VARIABLES = '//trim(var_list)
    write (unit, '(A,I0,A,I0,A,I0,A,I0,A,ES25.16E3)') &
        'ZONE T="Istep=', Istep, '", I=', nx_out, ', J=', ny, ', K=', nz, &
        ', DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=', tt
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx_out
                ii = i
                if (i > nx) ii = 1
                write (unit, '(*(ES25.16E3,1X))', iostat=ios) &
                    xx(ii, j, k), yy(ii, j, k), zz(ii, j, k), &
                    d(ii, j, k), u(ii, j, k), v(ii, j, k), w(ii, j, k), T(ii, j, k)
                call check_io(ios, 'Cannot write '//trim(output_filename))
            end do
        end do
    end do
    close (unit, iostat=ios)
    call check_io(ios, 'Cannot close '//trim(output_filename))
    print *, 'Wrote Tecplot ASCII: ', trim(output_filename)

    deallocate (xx, yy, zz, d, u, v, w, T)

contains

    subroutine check_io(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message
        if (status /= 0) then
            write (error_unit, '(A)') trim(message)
            stop 1
        end if
    end subroutine check_io
end program op2tec3d
