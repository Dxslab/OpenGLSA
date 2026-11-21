! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

subroutine IRAM_Solver
    use flow_data
    use fluc_data
    use mesh_data
    use Arnoldi
    use IO
    implicit none

    integer :: j, cyc, INFO_NAUPD

    if (my_id == 0) then
        print *, "Start IRAM solver"
    end if

    NEV = Para%nev
    NCV = max(30, 2*NEV+1, Para%ncv)
    NLOC = 5*nx*ny

    BMAT = 'I'
    WHICH = 'LM'

    LWORKL = 3*NCV*NCV+6*NCV
    TOL = Para%tol
    IDO = 0
    INFO = 0
    LDV = NLOC

    IPARAM(1) = 1
    IPARAM(3) = Para%max_iter
    IPARAM(7) = 1

    NDIGIT = 12
    LOGFIL = 6
    MNAUPD = 0
    MNAUP2 = 0

    if (my_id == 0) then
        MNAUPD = 3
        MNAUP2 = 2
    end if

    allocate (select(NCV), EIVAL(NCV, 2), RESID(NLOC), EIVEC(LDV, NCV), WORKD(3*NLOC), WORKEV(3*NCV), WORKL(3*NCV*NCV+6*NCV))
    if (my_id == 0) then
        print *, " "
        write (*, "(A, I5)") " Required number of eigenvalues:       ", NEV
        write (*, "(A, I5)") " Number of Lanczos vectors:            ", NCV
        write (*, "(A, I5)") " Maximum iteration cycles:             ", IPARAM(3)
        write (*, "(A, ES12.5)") " Arnoldi convergence tolerance:        ", TOL
        print *, " "
    end if

    call init_flow

    if (Para%IF_init == 1) then
        call ReadFlucToArnoldiArray(RESID)
        INFO = 1
    end if

    cyc = 0
    do while (ido /= 99)
        call PDNAUPD(MPI_COMM_WORLD, IDO, BMAT, NLOC, WHICH, NEV, TOL, RESID, NCV, &
                     EIVEC, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO)
        cyc = cyc+1
        if (IDO == -1 .or. IDO == 1) then
            call CopyArnoldiArrayToField(WORKD(IPNTR(1)))
            call LNS2d_Solver
            call CopyFieldToArnoldiArray(WORKD(IPNTR(2)))
        end if
    end do

    INFO_NAUPD = INFO

    if (INFO_NAUPD < 0) then
        if (my_id == 0) print *, "*** Error with PDNAUPD, INFO = ", INFO, " ***"
    else
        RVEC = .true.
        call PDNEUPD(MPI_COMM_WORLD, RVEC, 'A', select, EIVAL, EIVAL(1, 2), EIVEC, LDV, &
                     SIGMAR, SIGMAI, WORKEV, BMAT, NLOC, WHICH, NEV, TOL, &
                     RESID, NCV, EIVEC, LDV, IPARAM, IPNTR, WORKD, WORKL, &
                     LWORKL, INFO)

        if (INFO /= 0) then
            if (my_id == 0) print *, "*** Error with PDNEUPD, INFO = ", INFO, " ***"
        else
            NCONV = IPARAM(5)
            call PDMOUT(MPI_COMM_WORLD, 6, NCONV, 2, EIVAL, NCV, -6, "Ritz values (Real,Imag)")

            if (my_id == 0) call writeEVS

            do j = 1, NCONV
                call CopyArnoldiArrayToField(EIVEC(1, j))
                call writeMode(j)
            end do

        end if

        if (my_id == 0) then
            if (INFO_NAUPD == 1) then
                print *, '** Maximum number of Arnoldi cycles reached **'
            else if (INFO_NAUPD == 3) then
                print *, '** No shifts could be applied during implicit Arnoldi update, try increasing NCV **'
            end if

            print *, " "
            write (*, "(A, A)") " What portion of the spectrum:                     ", WHICH
            write (*, "(A, I5)") " The number of converged Ritz values is:           ", NCONV
            write (*, "(A, I5)") " The number of Implicit Arnoldi update iterations: ", IPARAM(3)
            write (*, "(A, I5)") " The number of OP*x is:                            ", IPARAM(9)
            write (*, "(A, ES12.5)") " The convergence criterion is:                     ", TOL
            print *, " "
        end if

    end if

    deallocate (select, EIVAL, RESID, EIVEC, WORKD, WORKEV, WORKL)

end subroutine IRAM_Solver

subroutine CopyArnoldiArrayToField(vector)
    use flow_data
    use mesh_data
    use fluc_data
    use Arnoldi
    implicit none

    integer :: i, j, idx_offset
    real(dp) :: vector(NLOC)

    idx_offset = 0
    do j = 1, ny
    do i = 1, nx
        d_fluc(i, j) = vector(idx_offset+1+(i-1)*5+(j-1)*5*nx)
        u_fluc(i, j) = vector(idx_offset+2+(i-1)*5+(j-1)*5*nx)
        v_fluc(i, j) = vector(idx_offset+3+(i-1)*5+(j-1)*5*nx)
        w_fluc(i, j) = vector(idx_offset+4+(i-1)*5+(j-1)*5*nx)
        T_fluc(i, j) = vector(idx_offset+5+(i-1)*5+(j-1)*5*nx)

        f(i, j, 1) = d_fluc(i, j)
        f(i, j, 2) = d(i, j)*u_fluc(i, j)+u(i, j)*d_fluc(i, j)
        f(i, j, 3) = d(i, j)*v_fluc(i, j)+v(i, j)*d_fluc(i, j)
        f(i, j, 4) = d(i, j)*w_fluc(i, j)
        f(i, j, 5) = d_fluc(i, j)*((u(i, j)*u(i, j)+v(i, j)*v(i, j))*0.5d0+ &
                                   Para%Cv*T(i, j))+d(i, j)*(u(i, j)*u_fluc(i, j)+v(i, j)*v_fluc(i, j)+Para%Cv*T_fluc(i, j))
    end do
    end do
    idx_offset = idx_offset+5*nx*ny

    call apply_bc_d
    tt = 0d0; Istep = 0

end subroutine CopyArnoldiArrayToField

subroutine CopyFieldToArnoldiArray(vector)
    use mesh_data
    use flow_data
    use fluc_data
    use Arnoldi
    implicit none

    integer :: i, j, idx_offset
    real(dp) :: vector(NLOC)

    idx_offset = 0
    do j = 1, ny
    do i = 1, nx
        vector(idx_offset+1+(i-1)*5+(j-1)*5*nx) = d_fluc(i, j)
        vector(idx_offset+2+(i-1)*5+(j-1)*5*nx) = u_fluc(i, j)
        vector(idx_offset+3+(i-1)*5+(j-1)*5*nx) = v_fluc(i, j)
        vector(idx_offset+4+(i-1)*5+(j-1)*5*nx) = w_fluc(i, j)
        vector(idx_offset+5+(i-1)*5+(j-1)*5*nx) = T_fluc(i, j)
    end do
    end do
    idx_offset = idx_offset+5*nx*ny

end subroutine CopyFieldToArnoldiArray

subroutine ReadFlucToArnoldiArray(vector)
    use mesh_data
    use flow_data
    use fluc_data
    use Arnoldi
    use IO
    implicit none

    integer :: IERR
    logical :: ext
    real(dp) :: vector(NLOC)

    if (my_id == 0) then
        inquire (FILE="fluc.dat", EXIST=ext)
        if (.not. ext) then
            print *, "*** ERROR: fluc.dat does not exist. ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        end if
        print *, "Reading initial Arnoldi vector from fluc.dat"
        open (55, FILE='fluc.dat', FORM='UNFORMATTED')
        read (55) Istep, tt
    end if

    call read_data(55, d_fluc)
    call read_data(55, u_fluc)
    call read_data(55, v_fluc)
    call read_data(55, w_fluc)
    call read_data(55, T_fluc)

    if (my_id == 0) then
        print *, "Read fluc field ok!"
        close (55)
    end if
    Istep = 0; tt = 0.0d0

    call CopyFieldToArnoldiArray(vector)

end subroutine ReadFlucToArnoldiArray

subroutine writeEVS
    use mesh_data
    use flow_para
    use Arnoldi
    implicit none

    integer :: i
    real(dp) :: abs_ev(NCONV), ang_ev(NCONV), delta_t
    real(dp), parameter :: PI = 4.0_dp*atan(1.0_dp)

    delta_t = real(Istep, dp)*Para%dt
    if (delta_t <= 0.0d0) delta_t = Para%End_time

    if (my_id == 0) then
        write (*, '(A, ES20.10)') " Target End_time: ", Para%End_time
        write (*, '(A, ES20.10)') " Actual dt*Nstep: ", delta_t
    end if

    open (102, FILE="Eigenvalues.dat")
    write (102, '(A4,5A25)') "No.", "Magnitude", "Angle", "Growth", "Frequency", "Strouhal Number"
    write (*, '(A4,5A25)') "No.", "Magnitude", "Angle", "Growth", "Frequency", "Strouhal Number"

    do i = 1, NCONV

        abs_ev(i) = hypot(EIVAL(i, 1), EIVAL(i, 2))
        ang_ev(i) = atan2(EIVAL(i, 2), EIVAL(i, 1))

      write (102, '(I4.2,5F25.16)') i, abs_ev(i), ang_ev(i), log(abs_ev(i))/delta_t, ang_ev(i)/delta_t, ang_ev(i)/delta_t/(2.0d0*PI)
        write (*, '(I4.2,5F25.16)') i, abs_ev(i), ang_ev(i), log(abs_ev(i))/delta_t, ang_ev(i)/delta_t, ang_ev(i)/delta_t/(2.0d0*PI)

    end do

    close (102)

end subroutine writeEVS

subroutine writeMode(kk)
    use fluc_data
    use IO
    implicit none

    integer :: kk
    character(LEN=128) :: filename

    if (my_id == 0) then
        write (filename, "('Mode-',I4.4,'.dat')") kk
        open (34, FILE=trim(filename), FORM="unformatted")
        write (34) kk, 0.0d0
    end if

    call write_data(34, d_fluc)
    call write_data(34, u_fluc)
    call write_data(34, v_fluc)
    call write_data(34, w_fluc)
    call write_data(34, T_fluc)

    if (my_id == 0) then
        close (34)
    end if

end subroutine writeMode
