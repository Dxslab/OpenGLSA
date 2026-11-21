! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine read_parameters
    use mesh_data
    use flow_para
    implicit none

    integer :: nparameters(100), IERR
    logical :: ext
    integer :: RunMode, IF_Viscous, IF_Cylindrical_Coordinate, IF_init, IF_SFD, IF_Character, pp_limiter
    integer :: Iperiodic_X, Iperiodic_Y, Isymmetric_X(2), Isymmetric_Y(2)
    integer :: Istep_Show, Istep_Save

    real(dp) :: rparameters(100)
    real(dp) :: Ma, Re, gamma, Pr, T_inf
    real(dp) :: dt, End_time
    real(dp) :: Periodic_ISpan(2), Periodic_JSpan(2), BC_Para(100)
    real(dp) :: Chi, Delta

    integer :: nev, ncv, max_iter
    real(dp) :: tol, wavenumber

    character(LEN=50) :: Scheme_Invis, Scheme_Vis, BoundaryCondition

    namelist /control_common/ &
        RunMode, &
        IF_Character, &
        pp_limiter, &
        IF_Viscous, IF_Cylindrical_Coordinate, &
        Scheme_Invis, Scheme_Vis, &
        nx_global, ny_global, LAP, &
        Ma, Re, gamma, Pr, T_inf, &
        Istep_Show, Istep_Save, &
        Isymmetric_X, Isymmetric_Y, &
        Iperiodic_X, Iperiodic_Y, &
        Periodic_ISpan, Periodic_JSpan, &
        BoundaryCondition, BC_Para

    namelist /control_cfd2d/ &
        IF_init, &
        dt, End_time, &
        IF_SFD, Chi, Delta

    namelist /control_gsa2d/ &
        IF_init, &
        dt, End_time, &
        tol, nev, ncv, max_iter, wavenumber

    !> set default values for parameters
    RunMode = 0       ! 0: NS solver, 1: global stability analysis solver, 2: LNS solver
    IF_Viscous = 1    ! 0: inviscid flow, 1: viscous flow
    IF_Cylindrical_Coordinate = 0  ! 0: Cartesian coordinate, 1: cylindrical coordinate
    IF_Character = 0  ! 0: conserved reconstruction, 1: physical-space splitting, 2: positivity-preserving characteristic LLF
    pp_limiter = 0      ! 0: no positivity-preserving limiter, 1: use positivity-preserving limiter
    Scheme_Invis = "LINEAR5"
    Scheme_Vis = "CD6"
    nx_global = 256
    ny_global = 256
    LAP = 4          ! overlap size
    Ma = 0.1d0
    Re = 100.0d0
    gamma = 1.4d0
    Pr = 0.72d0
    T_inf = 288.0d0
    Istep_Show = 100 ! show flow field every Istep_Show steps
    Istep_Save = 100 ! save flow field every Istep_Save steps
    Isymmetric_X = (/0, 0/)  ! 0: non-symmetric, 1: symmetric
    Isymmetric_Y = (/0, 0/)
    Iperiodic_X = 0  ! 0: non-periodic in x direction, 1: periodic in x direction
    Iperiodic_Y = 0  ! 0: non-periodic in y direction
    Periodic_ISpan = (/0.0d0, 1.0d0/)  ! periodic span in x direction (only used when Iperiodic_X=1)
    Periodic_JSpan = (/0.0d0, 1.0d0/)  ! periodic span in y direction (only used when Iperiodic_Y=1)
    BoundaryCondition = "BC_None"
    BC_Para = 0.0d0  ! parameters for boundary conditions
    !> parameters for different solvers
    IF_init = 0       ! 0: uniform flow (u = 1), 1: restart file, 2: quiescent flow (u = 0, e.g. closed cavities)
    dt = 0.001d0
    End_time = 10.0d0
    !> parameters for cfd2d solver
    IF_SFD = 0       ! 0: no selective frequency damping, 1: use selective frequency damping
    Chi = 0.0d0      ! control coefficient for selective frequency damping
    Delta = 1.0d0    ! filter width for selective frequency damping
    !> parameters for global stability analysis solver
    tol = 1.0d-6     ! convergence tolerance for eigenvalue solver
    nev = 10         ! number of eigenvalues to compute
    ncv = 30         ! number of Lanczos vectors to use in the eigenvalue
    max_iter = 10    ! maximum number of iterations for the eigenvalue solver
    wavenumber = 0.0d0  ! wavenumber for global stability analysis (only used when RunMode=2)

    if (my_id == 0) then
        inquire (FILE="control.nml", EXIST=ext)
        if (.not. ext) then
            print *, "*** ERROR: control.nml does not exist. ***"
            call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
        end if
        open (66, FILE="control.nml")
        read (66, NML=control_common)
        if (RunMode == 0) then
            read (66, NML=control_cfd2d)
        else if (RunMode == 1 .or. RunMode == 2) then
            read (66, NML=control_gsa2d)
        else
            print *, "*** ERROR: Invalid RunMode in control.nml. ***"
            call MPI_ABORT(MPI_COMM_WORLD, 2, IERR)
        end if
        close (66)

        if (RunMode /= 0) then
            IF_Character = 0             ! LNS/GSA: Steger-Warming splitting; comment out to use other splitting (IF_Character)
            Scheme_Invis = "LINEAR5"     ! LNS/GSA: linear 5th-order upwind; comment out to use other schemes (Scheme_Invis)
        end if

        print *, "====================== Parameters read from control.nml ======================="
        if (RunMode == 0) then
            print *, "*** Solver: Navier-Stokes time-stepping solver ***"
            if (IF_SFD == 1) then
                print *, "Selective Frequency Damping activated"
                print *, "Control coefficient Chi: ", Chi
                print *, "Filter width Delta: ", Delta
            end if
            print *, " "
        else if (RunMode == 1) then
            print *, "*** Solver: Global stability analysis solver ***"
            print *, "Convergence tolerance for eigenvalue solver: ", tol
            print *, "Number of eigenvalues to compute: ", nev
            print *, "Number of Lanczos vectors to use in the eigenvalue solver: ", ncv
            print *, "Maximum number of iterations for the eigenvalue solver: ", max_iter
            print *, "Wavenumber for global stability analysis: ", wavenumber
            print *, " "
        else if (RunMode == 2) then
            print *, "*** Solver: Linearized Navier-Stokes time-stepping solver ***"
            print *, " "
        end if
        if (IF_Viscous == 0) then
            print *, "Flow Type: Inviscid flow"
        else if (IF_Viscous == 1) then
            print *, "Flow Type: Viscous flow"
        end if
        if (IF_Cylindrical_Coordinate == 0) then
            print *, "Coordinate System: Cartesian coordinate"
        else if (IF_Cylindrical_Coordinate == 1) then
            print *, "Coordinate System: Cylindrical coordinate"
        end if

        print *, "Mach number: ", Ma
        print *, "Reynolds number: ", Re
        print *, "Specific heat ratio: ", gamma
        print *, "Prandtl number: ", Pr
        print *, "Reference temperature: ", T_inf
        print *, "Time step size: ", dt
        print *, "End time: ", End_time
        print *, "Inviscid scheme: ", Scheme_Invis
        print *, "Viscous scheme: ", Scheme_Vis
        print *, "Boundary condition: ", BoundaryCondition
        print *, "==============================================================================="

        call capitalize(Scheme_Invis)
        call capitalize(Scheme_Vis)
        call capitalize(BoundaryCondition)
        call find_scheme(Scheme_Invis, Scheme%Scheme_Invis)
        call find_scheme(Scheme_Vis, Scheme%Scheme_Vis)
        call find_bc(BoundaryCondition, Para%IBC)

        nparameters = 0
        rparameters = 0.0d0

        nparameters(1) = RunMode
        nparameters(2) = IF_Character
        nparameters(3) = pp_limiter
        nparameters(4) = IF_Viscous
        nparameters(6) = IF_Cylindrical_Coordinate
        nparameters(7) = Scheme%Scheme_Invis
        nparameters(8) = Scheme%Scheme_Vis
        nparameters(9) = nx_global
        nparameters(10) = ny_global
        nparameters(11) = LAP
        nparameters(12) = Istep_Show
        nparameters(13) = Istep_Save
        nparameters(14:15) = Isymmetric_X(1:2)
        nparameters(16:17) = Isymmetric_Y(1:2)
        nparameters(18) = Iperiodic_X
        nparameters(19) = Iperiodic_Y
        nparameters(20) = Para%IBC

        nparameters(21) = IF_init
        nparameters(22) = IF_SFD
        nparameters(23) = nev
        nparameters(24) = ncv
        nparameters(25) = max_iter

        rparameters(1) = Ma
        rparameters(2) = Re
        rparameters(3) = gamma
        rparameters(4) = Pr
        rparameters(5) = T_inf
        rparameters(6) = dt
        rparameters(7) = End_time
        rparameters(8:9) = Periodic_ISpan
        rparameters(10:11) = Periodic_JSpan
        rparameters(12) = Chi
        rparameters(13) = Delta
        rparameters(14) = tol
        rparameters(15) = wavenumber

    end if

    call MPI_BCAST(nparameters, 100, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    call MPI_BCAST(rparameters, 100, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, IERR)
    call MPI_BCAST(BC_Para, 100, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, IERR)

    Para%RunMode = nparameters(1)
    Para%IF_Character = nparameters(2)
    Para%pp_limiter = nparameters(3)
    Para%IF_Viscous = nparameters(4)
    Para%IF_Cylindrical_Coordinate = nparameters(6)
    Scheme%Scheme_Invis = nparameters(7)
    Scheme%Scheme_Vis = nparameters(8)
    nx_global = nparameters(9)
    ny_global = nparameters(10)
    LAP = nparameters(11)
    Para%Istep_Show = nparameters(12)
    Para%Istep_Save = nparameters(13)
    Para%Isymmetric_X(1:2) = nparameters(14:15)
    Para%Isymmetric_Y(1:2) = nparameters(16:17)
    Para%Iperiodic_X = nparameters(18)
    Para%Iperiodic_Y = nparameters(19)
    Para%IBC = nparameters(20)

    Para%IF_init = nparameters(21)
    Para%IF_SFD = nparameters(22)
    Para%nev = nparameters(23)
    Para%ncv = nparameters(24)
    Para%max_iter = nparameters(25)

    Para%Ma = rparameters(1)
    Para%Re = rparameters(2)
    Para%gamma = rparameters(3)
    Para%Pr = rparameters(4)
    Para%T_inf = rparameters(5)
    Para%dt = rparameters(6)
    Para%End_time = rparameters(7)
    Para%Periodic_ISpan(1:2) = rparameters(8:9)
    Para%Periodic_JSpan(1:2) = rparameters(10:11)
    Para%Chi = rparameters(12)
    Para%Delta = rparameters(13)
    Para%tol = rparameters(14)
    Para%wavenumber = rparameters(15)

    Para%BC_Para(1:100) = BC_Para(1:100)

    Para%Cv = 1.0d0/(Para%gamma*(Para%gamma-1.0d0)*Para%Ma**2)
    Para%Cp = Para%gamma*Para%Cv

    select case (Scheme%Scheme_Invis)
    case (OCFD_Scheme_WENOZ5, OCFD_Scheme_WENOZQ5, OCFD_Scheme_WENOMR5, OCFD_Scheme_LINEAR5)
        Scheme%Ka1 = -2; Scheme%Kb1 = 2
        Scheme%Ka2 = -1; Scheme%Kb2 = 3
    case DEFAULT
        print *, "*** Error: The Inviscous Scheme is not supported! ***"
        call MPI_ABORT(MPI_COMM_WORLD, 4, ierr)
    end select

    if (Para%runMode == 1) then
        Para%Istep_Show = nint(Para%End_time/Para%dt/10.0d0)
        if (my_id == 0) then
            print *, "runMode = 1 (Arnoldi Method), Istep_Show is reset to ", Para%Istep_Show
            print *, " "
        end if
    end if

end subroutine read_parameters

subroutine find_scheme(Schm1, scheme_type)
    use Constants
    implicit none
    character(LEN=50) :: Schm1
    integer :: scheme_type, ierr

    select case (trim(Schm1))

    case ("WENOZ5")
        scheme_type = OCFD_Scheme_WENOZ5
    case ("WENOZQ5")
        scheme_type = OCFD_Scheme_WENOZQ5
    case ("WENOMR5")
        scheme_type = OCFD_Scheme_WENOMR5
    case ("LINEAR5")
        scheme_type = OCFD_Scheme_LINEAR5
    case ("CD4")
        scheme_type = OCFD_Scheme_CD4
    case ("CD6")
        scheme_type = OCFD_Scheme_CD6
    case DEFAULT
        print *, "*** Error: This Scheme is not supported! ", trim(Schm1), " ***"
        print *, "*** Scheme_Invis can be WENOZ5/WENOZQ5/WENOMR5/WENOALWMR5/WENOMZQ3 ***"
        print *, "*** Scheme_Vis can be CD4/CD6 (4th/6th Central Scheme) ***"
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end select
end subroutine find_scheme

subroutine find_bc(bc, ibc)
    use Constants
    implicit none

    character(LEN=50):: bc
    integer :: ibc, ierr

    select case (trim(bc))
    case ("BC_NONE")
        ibc = BC_None
    case ("BC_FLAT_PLATE")
        ibc = BC_Flat_Plate
    case ("BC_LID_DRIVEN_CAVITY")
        ibc = BC_Lid_Driven_Cavity
    case ("BC_SPHERE")
        ibc = BC_Sphere
    case ("BC_CYLINDER")
        ibc = BC_Cylinder
    case ("BC_USER_DEF")
        ibc = BC_User_Def
    case DEFAULT
        print *, "*** ERROR: Invalid Boundary Condition ", trim(bc), " ***"
        call MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
    end select

end subroutine find_bc

subroutine capitalize(str)
    implicit none
    character(LEN=*):: str
    integer :: i, length

    length = len(str)

    do i = 1, length
        if (lge(str(i:i), 'a') .and. lle(str(i:i), 'z')) then
            str(i:i) = achar(iachar(str(i:i))-32)
        end if
    end do

end subroutine capitalize
