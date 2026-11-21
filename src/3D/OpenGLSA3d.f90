! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Copyright by Li Xinliang, LHD, Institute of Mechanics.   Email: lixl@imech.ac.cn
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

program OpenGLSA3d
    use flow_para
    implicit none
    !
    integer :: IERR
    !
    call MPI_INIT(IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np_size, IERR)
    !
    call read_parameters
    call partition_mpi
    call Allocate_arrays
    call read_mesh
    !
    select case (Para%RunMode)
    case (0)
        call NS_Solver
    case (1)
        call IRAM_Solver
    case (2)
        call LNS_Solver
    case DEFAULT
        if (my_id == 0) then
            print *, "Invalid RunMode! Please check the input file."
        end if
        call MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    end select
    !
    call Deallocate_arrays
    call MPI_FINALIZE(IERR)

end program OpenGLSA3d
