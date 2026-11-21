! Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! SPDX-License-Identifier: Apache-2.0

program generate_vortex
  implicit none

  integer, parameter:: nx = 128
  integer, parameter:: ny = 128
  real(8), parameter:: x_min = 0.0d0, x_max = 10.0d0
  real(8), parameter:: y_min = 0.0d0, y_max = 10.0d0
  real(8), parameter:: gamma = 1.4d0, Minf = sqrt(1.d0/gamma)
  real(8), parameter:: pi = 4.0d0*atan(1.0d0)
  real(8), parameter :: Sx = 6.0d0, Sy = 6.0d0, Ax = 0.04d0, Ay = 0.04d0, omega = 0.25d0, Lx = 2.0d0, Ly = 2.0d0
  real(8), parameter :: dx = (x_max-x_min)/dble(nx), dy = (y_max-y_min)/dble(ny)
  ! ...existing code...

  real(8) :: x(nx, ny), y(nx, ny)
  real(8) :: d(nx, ny), u(nx, ny), v(nx, ny), T(nx, ny), p(nx, ny)
  real(8) :: rr
  integer :: i, j
  integer, parameter :: iounit = 10

  do j = 1, ny
  do i = 1, nx
    x(i, j) = dble(i-1)*dx+Ax*sin(2.0d0*pi*omega)*sin(Sx*pi*dble(j-1)*dy/Ly)
    y(i, j) = dble(j-1)*dy+Ay*sin(2.0d0*pi*omega)*sin(Sy*pi*dble(i-1)*dx/Lx)
  end do
  end do

  do j = 1, ny
  do i = 1, nx
    rr = (x(i, j)-5.d0)**2.d0+(y(i, j)-5.d0)**2.d0
    u(i, j) = 1.d0-(y(i, j)-5.d0)*5.d0/(2.d0*pi)*exp(0.5d0*(1.d0-rr))
    v(i, j) = 1.d0+(x(i, j)-5.d0)*5.d0/(2.d0*pi)*exp(0.5d0*(1.d0-rr))
    T(i, j) = 1.d0-(5.d0**2.d0*(gamma-1.d0)/(8.d0*gamma*pi**2.d0))*exp(1.d0-rr)
    d(i, j) = T(i, j)**(1.d0/(gamma-1.d0))
    p(i, j) = d(i, j)**gamma
  end do
  end do

  open (unit=iounit, file='vortex.tec')
  write (iounit, '(A)') 'TITLE = "Isentropic Vortex"'
  write (iounit, '(A)') 'VARIABLES = "X", "Y", "Rho", "RhoU", "RhoV", "RhoE"'
  write (iounit, '(A, I5, A, I5, A)') 'ZONE T="Vortex", I=', nx, ', J=', ny, ', F=POINT'
  do j = 1, ny
  do i = 1, nx
    write(iounit, '(6F12.6)') x(i,j), y(i,j), d(i, j), d(i, j)*u(i, j), d(i, j)*v(i, j), p(i, j)/(gamma-1.d0)+0.5d0*d(i, j)*(u(i, j)**2.d0+v(i, j)**2.d0)
  end do
  end do
  close (iounit)

  open (unit=iounit+1, file='grid.dat', form="unformatted")
  write (iounit+1) x
  write (iounit+1) y
  close (iounit+1)

  open (unit=iounit+2, file='flow.dat', form="unformatted")
  write (iounit+2) 0, 0.d0
  write (iounit+2) d
  write (iounit+2) u
  write (iounit+2) v
  write (iounit+2) T
  close (iounit+2)

end program generate_vortex
