! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

module SFD
    use mesh_data
    use flow_data
    implicit none
    private
    public :: apply_SFD

contains

    subroutine apply_SFD
        integer :: i, j, k, m, IERR
        real(dp) :: csfd1, csfd2, Lsfd(4)
        real(dp) :: f_tmp(5)
        real(dp) :: resid(5), resid_sfd(5), resid_sfd_total

        csfd1 = 1.0d0/(1.0d0+Para%Chi*Para%Delta)
        csfd2 = exp(-1.0d0*(Para%chi+1.0d0/Para%Delta)*Para%dt)
        Lsfd(1) = (1.0d0+Para%chi*Para%Delta*csfd2)*csfd1
        Lsfd(2) = (Para%chi*Para%Delta*(1.0d0-csfd2))*csfd1
        Lsfd(3) = (1.0d0-csfd2)*csfd1
        Lsfd(4) = (Para%chi*Para%Delta+csfd2)*csfd1

        if (my_id == 0 .and. Istep == 1) then
            print *, "SFD COEFFICIENTS: "
            print *, Lsfd(1), Lsfd(2)
            print *, Lsfd(3), Lsfd(4)
        end if

        if (mod(Istep, Para%Istep_Show) == 0) then
            resid = 0.0d0
            do m = 1, 5
                do k = 1, nz
                do j = 1, ny
                do i = 1, nx
                    resid(m) = resid(m)+abs(f(i, j, k, m)-fsfd(i, j, k, m))
                end do
                end do
                end do
            end do

            call MPI_REDUCE(resid, resid_sfd, 5, OCFD_DATA_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, IERR)

            if (my_id == 0) then
                do m = 1, 5
                    resid_sfd(m) = resid_sfd(m)/real(nx_global*ny_global*nz_global, dp)
                end do
                resid_sfd_total = sum(resid_sfd)/5.0d0

                write (*, '(A15,6ES14.6)') "Residual(SFD): ", &
                    resid_sfd(1), resid_sfd(2), resid_sfd(3), resid_sfd(4), resid_sfd(5), resid_sfd_total
                open (78, FILE='resid-sfd.history', POSITION='APPEND')
                write (78, '(7ES14.6)') tt, resid_sfd(1), resid_sfd(2), resid_sfd(3), resid_sfd(4), resid_sfd(5), resid_sfd_total
                close (78)
            end if
        end if

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            f_tmp(1) = Lsfd(1)*f(i, j, k, 1)+Lsfd(2)*fsfd(i, j, k, 1)
            f_tmp(2) = Lsfd(1)*f(i, j, k, 2)+Lsfd(2)*fsfd(i, j, k, 2)
            f_tmp(3) = Lsfd(1)*f(i, j, k, 3)+Lsfd(2)*fsfd(i, j, k, 3)
            f_tmp(4) = Lsfd(1)*f(i, j, k, 4)+Lsfd(2)*fsfd(i, j, k, 4)
            f_tmp(5) = Lsfd(1)*f(i, j, k, 5)+Lsfd(2)*fsfd(i, j, k, 5)
            !
            fsfd(i, j, k, 1) = Lsfd(3)*f(i, j, k, 1)+Lsfd(4)*fsfd(i, j, k, 1)
            fsfd(i, j, k, 2) = Lsfd(3)*f(i, j, k, 2)+Lsfd(4)*fsfd(i, j, k, 2)
            fsfd(i, j, k, 3) = Lsfd(3)*f(i, j, k, 3)+Lsfd(4)*fsfd(i, j, k, 3)
            fsfd(i, j, k, 4) = Lsfd(3)*f(i, j, k, 4)+Lsfd(4)*fsfd(i, j, k, 4)
            fsfd(i, j, k, 5) = Lsfd(3)*f(i, j, k, 5)+Lsfd(4)*fsfd(i, j, k, 5)
            !
            f(i, j, k, 1) = f_tmp(1)
            f(i, j, k, 2) = f_tmp(2)
            f(i, j, k, 3) = f_tmp(3)
            f(i, j, k, 4) = f_tmp(4)
            f(i, j, k, 5) = f_tmp(5)
        end do
        end do
        end do

    end subroutine apply_SFD

end module SFD
