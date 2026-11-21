! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

program RESET_STEP
    implicit none

    integer, parameter :: HEADER_BYTES = 4+8
    integer :: iarg, nargs, ios, marker_head, marker_tail, Istep
    real(kind=8) :: tt
    character(len=512) :: filename

    nargs = command_argument_count()
    if (nargs == 0) then
        write (*, *) 'Usage: reset_step <flow_file> [more files ...]'
        stop 1
    end if

    do iarg = 1, nargs
        call get_command_argument(iarg, filename)

        open (66, FILE=trim(filename), ACCESS='STREAM', FORM='UNFORMATTED', &
              STATUS='OLD', ACTION='READWRITE', IOSTAT=ios)
        if (ios /= 0) then
            write (*, *) 'ERROR: cannot open ', trim(filename)
            cycle
        end if

        read (66, POS=1, IOSTAT=ios) marker_head, Istep, tt, marker_tail
        if (ios /= 0 .or. marker_head /= HEADER_BYTES .or. marker_tail /= HEADER_BYTES) then
            write (*, *) 'ERROR: unexpected header layout, skipped: ', trim(filename)
            close (66)
            cycle
        end if

        if (Istep == 0 .and. tt == 0.0d0) then
            write (*, '(A,A)') 'Already zero: ', trim(filename)
        else
            write (*, '(A,A,A,I10,A,ES20.12,A)') 'Reset ', trim(filename), &
                ' (Istep =', Istep, ', tt =', tt, ') -> 0'
            write (66, POS=5) 0, 0.0d0
        end if
        close (66)
    end do
end program RESET_STEP
