! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

program X2OP
    implicit none

    integer :: i, j
    integer :: nx, ny, nx_out, Istep
    logical :: omesh

    real(kind=8) :: tt
    real(kind=8), allocatable, dimension(:, :) :: xx, yy, d, u, v, T
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
    read (77, *) nx, ny

    allocate (xx(nx, ny), yy(nx, ny))

    read (77, *) ((xx(i, j), i=1, nx), j=1, ny), ((yy(i, j), i=1, nx), j=1, ny)
    close (77)

    nx_out = nx
    if (omesh) nx_out = nx-1
    allocate (d(nx_out, ny), u(nx_out, ny), v(nx_out, ny), T(nx_out, ny))

    write (*, *) "Generating grid.dat ..."
    open (78, file='grid.dat', form='unformatted')
    write (78) xx(1:nx_out, :)
    write (78) yy(1:nx_out, :)
    close (78)

    ! write (*, *) "Generating flow-00000000.dat ..."
    ! Istep = 0; tt = 0.d0
    ! d = 1.d0; u = 1.d0; v = 0.d0; T = 1.d0

    ! open (44, file='flow-00000000.dat', form='unformatted')
    ! write (44) Istep, tt
    ! write (44) d
    ! write (44) u
    ! write (44) v
    ! write (44) T
    ! close (44)

    ! Freestream values are used only in the ASCII grid preview.
    d = 1.0d0; u = 1.0d0; v = 0.0d0; T = 1.0d0

    open (79, file='grid2d.tec')
    write (79, *) 'variables=x,y,d,u,v,T'
    write (79, *) 'zone i=', nx_out, 'j=', ny
    do j = 1, ny
        do i = 1, nx_out
            write (79, '(6f20.12)') xx(i, j), yy(i, j), d(i, j), u(i, j), v(i, j), T(i, j)
        end do
    end do
    close (79)

    write (*, *) "Done!:>"
    deallocate (xx, yy)
    deallocate (d, u, v, T)
end program X2OP
