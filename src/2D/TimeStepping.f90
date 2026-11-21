! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine rk3(KRK, f, fn, du, neq)
    use mesh_data
    use flow_para
    implicit none

    integer :: KRK, m, i, j, neq
    real(dp), parameter:: Ralfa(3) = [1.0d0, 3.d0/4.d0, 1.0d0/3.0d0]
    real(dp), parameter:: Rbeta(3) = [1.0d0, 1.d0/4.d0, 2.0d0/3.0d0]

    real(dp), dimension(1:nx, 1:ny, 1:neq), intent(IN) :: du
    real(dp), dimension(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1:neq), intent(IN) :: fn
    real(dp), dimension(1-LAP:nx+LAP, 1-LAP:ny+LAP, 1:neq), intent(INOUT) :: f

    if (KRK == 1) then
        do m = 1, neq
            do j = 1, ny
            do i = 1, nx
                f(i, j, m) = fn(i, j, m)+Para%dt*du(i, j, m)
            end do
            end do
        end do
    else
        do m = 1, neq
            do j = 1, ny
            do i = 1, nx
                f(i, j, m) = Ralfa(KRK)*fn(i, j, m)+Rbeta(KRK)*(f(i, j, m)+Para%dt*du(i, j, m))
            end do
            end do
        end do
    end if

end subroutine rk3
