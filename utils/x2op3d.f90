! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

program X2OP3D
    implicit none

    integer :: i, j, k
    integer :: nx, ny, nz, nx_out, Istep
    logical :: omesh

    real(kind=8) :: tt
    real(kind=8), allocatable, dimension(:, :, :) :: xx, yy, zz, d, u, v, w, T
    character(len=200) :: filename, tmp

    if (command_argument_count() == 0) then
        write (*, *) 'ERROR, PLEASE INPUT FILENAME!'
        stop
    elseif (command_argument_count() .gt. 2) then
        write (*, *) 'ERROR, TOO MANY ARGUMENTS!'
        stop
    end if

    call get_command_argument(1, filename)

    omesh = .false.
    if (command_argument_count() >= 2) then
        call get_command_argument(2, tmp)
        if (trim(tmp) == 'omesh' .or. trim(tmp) == 'OMESH') omesh = .true.
    end if

    open (77, file=trim(filename))
    read (77, *) nx, ny, nz

    allocate (xx(nx, ny, nz), yy(nx, ny, nz), zz(nx, ny, nz))

    read (77, *) (((xx(i, j, k), i=1, nx), j=1, ny), k=1, nz), &
        (((yy(i, j, k), i=1, nx), j=1, ny), k=1, nz), &
        (((zz(i, j, k), i=1, nx), j=1, ny), k=1, nz)
    close (77)

    nx_out = nx
    if (omesh) nx_out = nx-1
    allocate (d(nx_out, ny, nz), u(nx_out, ny, nz), v(nx_out, ny, nz), &
              w(nx_out, ny, nz), T(nx_out, ny, nz))

    write (*, *) "Generating grid3d.dat ..."
    open (78, file='grid3d.dat', form='unformatted')
    do k = 1, nz
        write (78) xx(1:nx_out, :, k)
    end do
    do k = 1, nz
        write (78) yy(1:nx_out, :, k)
    end do
    do k = 1, nz
        write (78) zz(1:nx_out, :, k)
    end do
    close (78)

    write (*, *) "Generating flow3d.dat ..."
    Istep = 0; tt = 0.d0
    d = 1.d0; u = 1.d0; v = 0.d0; w = 0.d0; T = 1.d0

    open (44, file='flow3d.dat', form='unformatted')
    write (44) Istep, tt
    do k = 1, nz
        write (44) d(:, :, k)
    end do
    do k = 1, nz
        write (44) u(:, :, k)
    end do
    do k = 1, nz
        write (44) v(:, :, k)
    end do
    do k = 1, nz
        write (44) w(:, :, k)
    end do
    do k = 1, nz
        write (44) T(:, :, k)
    end do
    close (44)

    open (79, file='grid3d.tec')
    write (79, *) 'variables=x,y,z,d,u,v,w,T'
    write (79, *) 'zone i=', nx_out, 'j=', ny, 'k=', nz
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx_out
        write (79, '(8f20.12)') xx(i, j, k), yy(i, j, k), zz(i, j, k), &
            d(i, j, k), u(i, j, k), v(i, j, k), w(i, j, k), T(i, j, k)
    end do
    end do
    end do
    close (79)

    write (*, *) "Done!:>"
    deallocate (xx, yy, zz)
    deallocate (d, u, v, w, T)
end program X2OP3D
