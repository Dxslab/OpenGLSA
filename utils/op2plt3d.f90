! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

program op2plt3d
    use iso_fortran_env, only: error_unit, int64
    use iso_c_binding, only: c_int32_t, c_double, c_null_char
    implicit none

    include "tecio.f90"

    integer :: nx, ny, nz, Istep, k, ios, slash, dot
    real(c_double) :: tt
    real(c_double), dimension(:, :, :), allocatable :: xx, yy, zz, d, u, v, w, T

    character(len=1024) :: data_filename, output_filename, tmp
    character(len=64)  :: var_list
    ! TecIO 142 uses 32-bit counts and nodal, unshared variables.
    integer(c_int32_t) :: count, imax, jmax, kmax, status
    integer(c_int32_t) :: passive(8), location(8), sharing(8)
    integer(int64) :: points
    real(c_double), allocatable :: block(:, :, :)
    logical :: omesh
    integer :: nx_out

    if (command_argument_count() < 4 .or. command_argument_count() > 5) then
        print *, "Usage: op2plt3d <data_file> <nx> <ny> <nz> [omesh]"
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

    var_list = 'X Y Z Density U V W Temperature'

    ! Replace the final extension; never overwrite the binary input itself.
    slash = scan(trim(data_filename), '/', back=.true.)
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
    if (points > int(huge(count), int64)/int(nz, int64)) then
        write (error_unit, '(A)') 'Grid exceeds the TecIO 142 point-count limit'
        stop 1
    end if
    points = points*int(nz, int64)

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

    count = int(points, c_int32_t)
    imax = nx_out; jmax = ny; kmax = nz
    passive = 0; location = 1; sharing = 0
    call check_io(TECINI142('OpenGLSA solution'//c_null_char, trim(var_list)//c_null_char, &
        trim(output_filename)//c_null_char, '.'//c_null_char, &
        0_c_int32_t, 0_c_int32_t, 0_c_int32_t, 1_c_int32_t), 'TECINI142: cannot create output')
    status = TECZNE142('OpenGLSA solution'//c_null_char, &
        0_c_int32_t, imax, jmax, kmax, 0_c_int32_t, 0_c_int32_t, 0_c_int32_t, &
        tt, 1_c_int32_t, 0_c_int32_t, 1_c_int32_t, 0_c_int32_t, 0_c_int32_t, &
        0_c_int32_t, 0_c_int32_t, 0_c_int32_t, passive, location, sharing, 0_c_int32_t)
    call check_io(status, 'TECZNE142: cannot create zone')
    allocate (block(nx_out, ny, nz))
    block(1:nx, :, :) = xx
    if (omesh) block(nx_out, :, :) = xx(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write xx')
    block(1:nx, :, :) = yy
    if (omesh) block(nx_out, :, :) = yy(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write yy')
    block(1:nx, :, :) = zz
    if (omesh) block(nx_out, :, :) = zz(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write zz')
    block(1:nx, :, :) = d
    if (omesh) block(nx_out, :, :) = d(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write d')
    block(1:nx, :, :) = u
    if (omesh) block(nx_out, :, :) = u(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write u')
    block(1:nx, :, :) = v
    if (omesh) block(nx_out, :, :) = v(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write v')
    block(1:nx, :, :) = w
    if (omesh) block(nx_out, :, :) = w(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write w')
    block(1:nx, :, :) = T
    if (omesh) block(nx_out, :, :) = T(1, :, :)
    call check_io(TECDATD142(count, block), 'TECDATD142: cannot write T')
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
end program op2plt3d
