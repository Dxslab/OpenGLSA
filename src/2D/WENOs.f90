! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module WENO_Reconstruction
    use precision
    implicit none

    real(dp), parameter:: ep_wenoz = 1.d-20
    real(dp), parameter:: ep_wenozq = 1.d-6
    real(dp), parameter:: ep_wenomr = 1.d-6

contains

    subroutine hh_WENO5ZP(Ka, Kb, v, hh)          ! Ka=-2, Kb=2

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: S0, S1, S2, a0, a1, a2, am, q03, q13, q23, tau5
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        S0 = 13.d0/12.d0*(v(0)-2.d0*v(1)+v(2))**2+1.d0/4.d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
        S1 = 13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+1.d0/4.d0*(v(-1)-v(1))**2
        S2 = 13.d0/12.d0*(v(-2)-2.d0*v(-1)+v(0))**2+1.d0/4.d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
        tau5 = abs(S0-S2)
        a0 = C03*(1.d0+(tau5/(ep_wenoz+S0))**2)
        a1 = C13*(1.d0+(tau5/(ep_wenoz+S1))**2)
        a2 = C23*(1.d0+(tau5/(ep_wenoz+S2))**2)
        am = a0+a1+a2
        q03 = 1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
        q13 = -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)
        q23 = 1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
        hh = (a0*q03+a1*q13+a2*q23)/am

    end subroutine hh_WENO5ZP

    subroutine hh_WENO5ZM(Ka, Kb, v, hh)           ! Ka=-1,  Kb=3

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: S0, S1, S2, a0, a1, a2, am, q03, q13, q23, tau5
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        S0 = 13.d0/12.d0*(v(1)-2.d0*v(0)+v(-1))**2+1.d0/4.d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
        S1 = 13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+1.d0/4.d0*(v(2)-v(0))**2
        S2 = 13.d0/12.d0*(v(3)-2.d0*v(2)+v(1))**2+1.d0/4.d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
        tau5 = abs(S0-S2)
        a0 = C03*(1.d0+(tau5/(ep_wenoz+S0))**2)
        a1 = C13*(1.d0+(tau5/(ep_wenoz+S1))**2)
        a2 = C23*(1.d0+(tau5/(ep_wenoz+S2))**2)
        am = a0+a1+a2
        q03 = 1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)
        q13 = -1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)
        q23 = 1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)
        hh = (a0*q03+a1*q13+a2*q23)/am

    end subroutine hh_WENO5ZM

    subroutine hh_WENOZQ5P(Ka, Kb, v, hh)

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3
        real(dp), parameter :: d1 = 0.98d0, d2 = 0.01d0, d3 = 0.01d0

        beta1 = (v(-2)-8.d0*v(-1)+8.d0*v(1)-v(2))**2/144.d0+ &
                (-11.d0*v(-2)+174.d0*v(-1)-326.d0*v(0)+174.d0*v(1)-11.d0*v(2))**2/15600.d0+ &
                781.d0*(-v(-2)+2.d0*v(-1)-2.d0*v(1)+v(2))**2/2880.d0+ &
                1421461.d0*(v(-2)-4.d0*v(-1)+6.d0*v(0)-4.d0*v(1)+v(2))**2/1310400.d0
        beta2 = (v(-1)-v(0))**2
        beta3 = (v(0)-v(1))**2

        tau = ((abs(beta1-beta2)+abs(beta1-beta3))/2.d0)**2
        alpha1 = d1*(1.d0+tau/(ep_wenozq+beta1))
        alpha2 = d2*(1.d0+tau/(ep_wenozq+beta2))
        alpha3 = d3*(1.d0+tau/(ep_wenozq+beta3))
        w1 = alpha1/(alpha1+alpha2+alpha3)
        w2 = alpha2/(alpha1+alpha2+alpha3)
        w3 = alpha3/(alpha1+alpha2+alpha3)

        p1 = (2.d0*v(-2)-13.d0*v(-1)+47.d0*v(0)+27.d0*v(1)-3.d0*v(2))/60.d0
        p2 = (-v(-1)+3.d0*v(0))/2.d0
        p3 = (v(0)+v(1))/2.d0

        hh = (w1/d1)*(p1-d2*p2-d3*p3)+w2*p2+w3*p3

    end subroutine hh_WENOZQ5P

    subroutine hh_WENOZQ5M(Ka, Kb, v, hh)

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3
        real(dp), parameter :: d1 = 0.98d0, d2 = 0.01d0, d3 = 0.01d0

        beta1 = (v(3)-8.d0*v(2)+8.d0*v(0)-v(-1))**2/144.d0+ &
                (-11.d0*v(3)+174.d0*v(2)-326.d0*v(1)+174.d0*v(0)-11.d0*v(-1))**2/15600.d0+ &
                781.d0*(-v(3)+2.d0*v(2)-2.d0*v(0)+v(-1))**2/2880.d0+ &
                1421461.d0*(v(3)-4.d0*v(2)+6.d0*v(1)-4.d0*v(0)+v(-1))**2/1310400.d0
        beta2 = (v(2)-v(1))**2
        beta3 = (v(1)-v(0))**2

        tau = ((abs(beta1-beta2)+abs(beta1-beta3))/2.d0)**2
        alpha1 = d1*(1.d0+tau/(ep_wenozq+beta1))
        alpha2 = d2*(1.d0+tau/(ep_wenozq+beta2))
        alpha3 = d3*(1.d0+tau/(ep_wenozq+beta3))
        w1 = alpha1/(alpha1+alpha2+alpha3)
        w2 = alpha2/(alpha1+alpha2+alpha3)
        w3 = alpha3/(alpha1+alpha2+alpha3)

        p1 = (2.d0*v(3)-13.d0*v(2)+47.d0*v(1)+27.d0*v(0)-3.d0*v(-1))/60.d0
        p2 = (-v(2)+3.d0*v(1))/2.d0
        p3 = (v(1)+v(0))/2.d0

        hh = (w1/d1)*(p1-d2*p2-d3*p3)+w2*p2+w3*p3

    end subroutine hh_WENOZQ5M

    subroutine hh_WENOMR5P(Ka, Kb, v, hh)

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3, q1, q2, q3
        real(dp) :: s0, s1, g01s, g11s, g01, g11, si0, si1, si
        real(dp), parameter :: d12 = 1.d0/11.d0, d22 = 10.d0/11.d0
        real(dp), parameter :: d13 = 1.d0/111.d0, d23 = 10.d0/111.d0, d33 = 100.d0/111.d0

        s0 = (v(0)-v(-1))**2; s1 = (v(1)-v(0))**2
        if (s0 >= s1) then
            g01s = 1.d0
        else
            g01s = 10.d0
        end if
        g11s = 11.d0-g01s
        g01 = g01s/11.d0
        g11 = 1.d0-g01
        si0 = g01*(1.d0+(s0-s1)**2/(s0+ep_wenomr))
        si1 = g11*(1.d0+(s0-s1)**2/(s1+ep_wenomr))
        si = si0+si1

        beta1 = (si0*(v(0)-v(-1))+si1*(v(1)-v(0)))**2/si**2
        beta2 = 13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+1.d0/4.d0*(v(-1)-v(1))**2
        beta3 = (v(-2)-8.d0*v(-1)+8.d0*v(1)-v(2))**2/144.d0+ &
                (-11.d0*v(-2)+174.d0*v(-1)-326.d0*v(0)+174.d0*v(1)-11.d0*v(2))**2/15600.d0+ &
                781.d0*(-v(-2)+2.d0*v(-1)-2.d0*v(1)+v(2))**2/2880.d0+ &
                1421461.d0*(v(-2)-4.d0*v(-1)+6.d0*v(0)-4.d0*v(1)+v(2))**2/1310400.d0

        tau = ((abs(beta3-beta1)+abs(beta3-beta2))/2.d0)**2
        alpha1 = d13*(1.d0+tau/(ep_wenomr+beta1))
        alpha2 = d23*(1.d0+tau/(ep_wenomr+beta2))
        alpha3 = d33*(1.d0+tau/(ep_wenomr+beta3))
        w1 = alpha1/(alpha1+alpha2+alpha3)
        w2 = alpha2/(alpha1+alpha2+alpha3)
        w3 = alpha3/(alpha1+alpha2+alpha3)

        q1 = v(0)
        q2 = -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)
        q3 = (2.d0*v(-2)-13.d0*v(-1)+47.d0*v(0)+27.d0*v(1)-3.d0*v(2))/60.d0

        p1 = q1
        p2 = (1.d0/d22)*q2-(d12/d22)*p1
        p3 = (1.d0/d33)*q3-(d13/d33)*p1-(d23/d33)*p2

        hh = w1*p1+w2*p2+w3*p3

    end subroutine hh_WENOMR5P

    subroutine hh_WENOMR5M(Ka, Kb, v, hh)

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3, q1, q2, q3
        real(dp) :: s0, s1, g01s, g11s, g01, g11, si0, si1, si
        real(dp), parameter :: d12 = 1.d0/11.d0, d22 = 10.d0/11.d0
        real(dp), parameter :: d13 = 1.d0/111.d0, d23 = 10.d0/111.d0, d33 = 100.d0/111.d0

        s0 = (v(1)-v(2))**2; s1 = (v(0)-v(1))**2
        if (s0 >= s1) then
            g01s = 1.d0
        else
            g01s = 10.d0
        end if
        g11s = 11.d0-g01s
        g01 = g01s/11.d0
        g11 = 1.d0-g01
        si0 = g01*(1.d0+(s0-s1)**2/(s0+ep_wenomr))
        si1 = g11*(1.d0+(s0-s1)**2/(s1+ep_wenomr))
        si = si0+si1

        beta1 = (si0*(v(1)-v(2))+si1*(v(0)-v(1)))**2/si**2
        beta2 = 13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+1.d0/4.d0*(v(2)-v(0))**2
        beta3 = (v(3)-8.d0*v(2)+8.d0*v(0)-v(-1))**2/144.d0+ &
                (-11.d0*v(3)+174.d0*v(2)-326.d0*v(1)+174.d0*v(0)-11.d0*v(-1))**2/15600.d0+ &
                781.d0*(-v(3)+2.d0*v(2)-2.d0*v(0)+v(-1))**2/2880.d0+ &
                1421461.d0*(v(3)-4.d0*v(2)+6.d0*v(1)-4.d0*v(0)+v(-1))**2/1310400.d0

        tau = ((abs(beta3-beta1)+abs(beta3-beta2))/2.d0)**2
        alpha1 = d13*(1.d0+tau/(ep_wenomr+beta1))
        alpha2 = d23*(1.d0+tau/(ep_wenomr+beta2))
        alpha3 = d33*(1.d0+tau/(ep_wenomr+beta3))
        w1 = alpha1/(alpha1+alpha2+alpha3)
        w2 = alpha2/(alpha1+alpha2+alpha3)
        w3 = alpha3/(alpha1+alpha2+alpha3)

        q1 = v(1)
        q2 = -1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)
        q3 = (2.d0*v(3)-13.d0*v(2)+47.d0*v(1)+27.d0*v(0)-3.d0*v(-1))/60.d0

        p1 = q1
        p2 = (1.d0/d22)*q2-(d12/d22)*p1
        p3 = (1.d0/d33)*q3-(d13/d33)*p1-(d23/d33)*p2

        hh = w1*p1+w2*p2+w3*p3

    end subroutine hh_WENOMR5M

    subroutine hh_WENOMR3P(Ka, Kb, v, hh)

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: beta1, beta2, alpha1, alpha2, w1, w2, tau
        real(dp) :: p1, p2, q1, q2
        real(dp) :: s0, s1, g01s, g11s, g01, g11, si0, si1, si
        real(dp), parameter :: d12 = 1.d0/11.d0, d22 = 10.d0/11.d0

        s0 = (v(0)-v(-1))**2; s1 = (v(1)-v(0))**2
        if (s0 >= s1) then
            g01s = 1.d0
        else
            g01s = 10.d0
        end if
        g11s = 11.d0-g01s
        g01 = g01s/11.d0
        g11 = 1.d0-g01
        si0 = g01*(1.d0+abs(s0-s1)/(s0+ep_wenomr))
        si1 = g11*(1.d0+abs(s0-s1)/(s1+ep_wenomr))
        si = si0+si1

        beta1 = (si0*(v(0)-v(-1))+si1*(v(1)-v(0)))**2/si**2
        beta2 = 13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+1.d0/4.d0*(v(-1)-v(1))**2

        tau = abs(beta2-beta1)
        alpha1 = d12*(1.d0+tau/(ep_wenomr+beta1))
        alpha2 = d22*(1.d0+tau/(ep_wenomr+beta2))
        w1 = alpha1/(alpha1+alpha2)
        w2 = alpha2/(alpha1+alpha2)

        q1 = v(0)
        q2 = -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)

        p1 = q1
        p2 = (1.d0/d22)*q2-(d12/d22)*p1

        hh = w1*p1+w2*p2

    end subroutine hh_WENOMR3P

    subroutine hh_WENOMR3M(Ka, Kb, v, hh)

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: beta1, beta2, alpha1, alpha2, w1, w2, tau
        real(dp) :: p1, p2, q1, q2
        real(dp) :: s0, s1, g01s, g11s, g01, g11, si0, si1, si
        real(dp), parameter :: d12 = 1.d0/11.d0, d22 = 10.d0/11.d0

        s0 = (v(1)-v(2))**2; s1 = (v(0)-v(1))**2
        if (s0 >= s1) then
            g01s = 1.d0
        else
            g01s = 10.d0
        end if
        g11s = 11.d0-g01s
        g01 = g01s/11.d0
        g11 = 1.d0-g01
        si0 = g01*(1.d0+abs(s0-s1)/(s0+ep_wenomr))
        si1 = g11*(1.d0+abs(s0-s1)/(s1+ep_wenomr))
        si = si0+si1

        beta1 = (si0*(v(1)-v(2))+si1*(v(0)-v(1)))**2/si**2
        beta2 = 13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+1.d0/4.d0*(v(2)-v(0))**2

        tau = abs(beta2-beta1)
        alpha1 = d12*(1.d0+tau/(ep_wenomr+beta1))
        alpha2 = d22*(1.d0+tau/(ep_wenomr+beta2))
        w1 = alpha1/(alpha1+alpha2)
        w2 = alpha2/(alpha1+alpha2)

        q1 = v(1)
        q2 = -1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)

        p1 = q1
        p2 = (1.d0/d22)*q2-(d12/d22)*p1

        hh = w1*p1+w2*p2

    end subroutine hh_WENOMR3M

    subroutine hh_LINEAR5P(Ka, Kb, v, hh)          ! Ka=-2, Kb=2

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        q03 = 1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
        q13 = -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)
        q23 = 1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)

        hh = C03*q03+C13*q13+C23*q23

    end subroutine hh_LINEAR5P

    subroutine hh_LINEAR5M(Ka, Kb, v, hh)           ! Ka=-1, Kb=3

        integer :: Ka, Kb
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        q03 = 1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)
        q13 = -1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)
        q23 = 1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)

        hh = C03*q03+C13*q13+C23*q23

    end subroutine hh_LINEAR5M

    subroutine hh_WENO5P_boundary(Ka, Kb, v, hh, del_stencil)
        use Constants, only: del_left, del_right
        implicit none

        integer :: Ka, Kb, del_stencil
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: S0, S1, S2, a0, a1, a2, am, q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        a0 = 0.0d0; q03 = 0.0d0
        a2 = 0.0d0; q23 = 0.0d0

        if (del_stencil /= del_right) then
            S0 = 13.d0/12.d0*(v(0)-2.d0*v(1)+v(2))**2+1.d0/4.d0*(3.d0*v(0)-4.d0*v(1)+v(2))**2
            a0 = C03/((ep_wenozq+S0)**2)
            q03 = 1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
        end if

        S1 = 13.d0/12.d0*(v(-1)-2.d0*v(0)+v(1))**2+1.d0/4.d0*(v(-1)-v(1))**2
        a1 = C13/((ep_wenozq+S1)**2)
        q13 = -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)

        if (del_stencil /= del_left) then
            S2 = 13.d0/12.d0*(v(-2)-2.d0*v(-1)+v(0))**2+1.d0/4.d0*(v(-2)-4.d0*v(-1)+3.d0*v(0))**2
            a2 = C23/((ep_wenozq+S2)**2)
            q23 = 1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
        end if

        am = a0+a1+a2
        hh = (a0*q03+a1*q13+a2*q23)/am

    end subroutine hh_WENO5P_boundary

    subroutine hh_WENO5M_boundary(Ka, Kb, v, hh, del_stencil)
        use Constants, only: del_left, del_right
        implicit none

        integer :: Ka, Kb, del_stencil
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: S0, S1, S2, a0, a1, a2, am, q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        a0 = 0.0d0; q03 = 0.0d0
        a2 = 0.0d0; q23 = 0.0d0

        if (del_stencil /= del_left) then
            S0 = 13.d0/12.d0*(v(1)-2.d0*v(0)+v(-1))**2+1.d0/4.d0*(3.d0*v(1)-4.d0*v(0)+v(-1))**2
            a0 = C03/((ep_wenozq+S0)**2)
            q03 = 1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)
        end if

        S1 = 13.d0/12.d0*(v(2)-2.d0*v(1)+v(0))**2+1.d0/4.d0*(v(2)-v(0))**2
        a1 = C13/((ep_wenozq+S1)**2)
        q13 = -1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)

        if (del_stencil /= del_right) then
            S2 = 13.d0/12.d0*(v(3)-2.d0*v(2)+v(1))**2+1.d0/4.d0*(v(3)-4.d0*v(2)+3.d0*v(1))**2
            a2 = C23/((ep_wenozq+S2)**2)
            q23 = 1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)
        end if

        am = a0+a1+a2
        hh = (a0*q03+a1*q13+a2*q23)/am

    end subroutine hh_WENO5M_boundary

    subroutine hh_LINEAR5P_boundary(Ka, Kb, v, hh, del_stencil)
        use Constants, only: del_left, del_right
        implicit none

        integer :: Ka, Kb, del_stencil
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: a0, a1, a2, am, q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        a0 = 0.0d0; q03 = 0.0d0
        a2 = 0.0d0; q23 = 0.0d0

        if (del_stencil /= del_right) then
            a0 = C03
            q03 = 1.d0/3.d0*v(0)+5.d0/6.d0*v(1)-1.d0/6.d0*v(2)
        end if

        a1 = C13
        q13 = -1.d0/6.d0*v(-1)+5.d0/6.d0*v(0)+1.d0/3.d0*v(1)

        if (del_stencil /= del_left) then
            a2 = C23
            q23 = 1.d0/3.d0*v(-2)-7.d0/6.d0*v(-1)+11.d0/6.d0*v(0)
        end if

        am = a0+a1+a2
        hh = (a0*q03+a1*q13+a2*q23)/am

    end subroutine hh_LINEAR5P_boundary

    subroutine hh_LINEAR5M_boundary(Ka, Kb, v, hh, del_stencil)
        use Constants, only: del_left, del_right
        implicit none

        integer :: Ka, Kb, del_stencil
        real(dp) :: v(Ka:Kb), hh
        real(dp) :: a0, a1, a2, am, q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        a0 = 0.0d0; q03 = 0.0d0
        a2 = 0.0d0; q23 = 0.0d0

        if (del_stencil /= del_left) then
            a0 = C03
            q03 = 1.d0/3.d0*v(1)+5.d0/6.d0*v(0)-1.d0/6.d0*v(-1)
        end if

        a1 = C13
        q13 = -1.d0/6.d0*v(2)+5.d0/6.d0*v(1)+1.d0/3.d0*v(0)

        if (del_stencil /= del_right) then
            a2 = C23
            q23 = 1.d0/3.d0*v(3)-7.d0/6.d0*v(2)+11.d0/6.d0*v(1)
        end if

        am = a0+a1+a2
        hh = (a0*q03+a1*q13+a2*q23)/am

    end subroutine hh_LINEAR5M_boundary

    subroutine OCFD_WENOZ5P(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) ::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23, tau5
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        do i = 0, nx
            S0 = 13.d0/12.d0*(v(i)-2.d0*v(i+1)+v(i+2))**2+1.d0/4.d0*(3.d0*v(i)-4.d0*v(i+1)+v(i+2))**2
            S1 = 13.d0/12.d0*(v(i-1)-2.d0*v(i)+v(i+1))**2+1.d0/4.d0*(v(i-1)-v(i+1))**2
            S2 = 13.d0/12.d0*(v(i-2)-2.d0*v(i-1)+v(i))**2+1.d0/4.d0*(v(i-2)-4.d0*v(i-1)+3.d0*v(i))**2
            tau5 = abs(S0-S2)
            a0 = C03*(1.d0+(tau5/(ep_wenoz+S0))**2)
            a1 = C13*(1.d0+(tau5/(ep_wenoz+S1))**2)
            a2 = C23*(1.d0+(tau5/(ep_wenoz+S2))**2)
            am = a0+a1+a2
            q03 = 1.d0/3.d0*v(i)+5.d0/6.d0*v(i+1)-1.d0/6.d0*v(i+2)
            q13 = -1.d0/6.d0*v(i-1)+5.d0/6.d0*v(i)+1.d0/3.d0*v(i+1)
            q23 = 1.d0/3.d0*v(i-2)-7.d0/6.d0*v(i-1)+11.d0/6.d0*v(i)
            hh(i) = (a0*q03+a1*q13+a2*q23)/am
        end do

    end subroutine OCFD_WENOZ5P

    subroutine OCFD_WENOZ5M(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) ::  S0, S1, S2, a0, a1, a2, am, q03, q13, q23, tau5
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        do i = 1, nx+1
            S0 = 13.d0/12.d0*(v(i)-2.d0*v(i-1)+v(i-2))**2+1.d0/4.d0*(3.d0*v(i)-4.d0*v(i-1)+v(i-2))**2
            S1 = 13.d0/12.d0*(v(i+1)-2.d0*v(i)+v(i-1))**2+1.d0/4.d0*(v(i+1)-v(i-1))**2
            S2 = 13.d0/12.d0*(v(i+2)-2.d0*v(i+1)+v(i))**2+1.d0/4.d0*(v(i+2)-4.d0*v(i+1)+3.d0*v(i))**2
            tau5 = abs(S0-S2)
            a0 = C03*(1.d0+(tau5/(ep_wenoz+S0))**2)
            a1 = C13*(1.d0+(tau5/(ep_wenoz+S1))**2)
            a2 = C23*(1.d0+(tau5/(ep_wenoz+S2))**2)
            am = a0+a1+a2
            q03 = 1.d0/3.d0*v(i)+5.d0/6.d0*v(i-1)-1.d0/6.d0*v(i-2)
            q13 = -1.d0/6.d0*v(i+1)+5.d0/6.d0*v(i)+1.d0/3.d0*v(i-1)
            q23 = 1.d0/3.d0*v(i+2)-7.d0/6.d0*v(i+1)+11.d0/6.d0*v(i)
            hh(i-1) = (a0*q03+a1*q13+a2*q23)/am
        end do

    end subroutine OCFD_WENOZ5M

    subroutine OCFD_WENOZQ5P(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3
        real(dp), parameter :: d1 = 0.98d0, d2 = 0.01d0, d3 = 0.01d0

        do i = 0, nx
            beta1 = (v(i-2)-8.d0*v(i-1)+8.d0*v(i+1)-v(i+2))**2/144.d0+ &
                    (-11.d0*v(i-2)+174.d0*v(i-1)-326.d0*v(i)+174.d0*v(i+1)-11.d0*v(i+2))**2/15600.d0+ &
                    781.d0*(-v(i-2)+2.d0*v(i-1)-2.d0*v(i+1)+v(i+2))**2/2880.d0+ &
                    1421461.d0*(v(i-2)-4.d0*v(i-1)+6.d0*v(i)-4.d0*v(i+1)+v(i+2))**2/1310400.d0
            beta2 = (v(i-1)-v(i))**2
            beta3 = (v(i)-v(i+1))**2

            tau = ((abs(beta1-beta2)+abs(beta1-beta3))/2.d0)**2
            alpha1 = d1*(1.d0+tau/(ep_wenozq+beta1))
            alpha2 = d2*(1.d0+tau/(ep_wenozq+beta2))
            alpha3 = d3*(1.d0+tau/(ep_wenozq+beta3))
            w1 = alpha1/(alpha1+alpha2+alpha3)
            w2 = alpha2/(alpha1+alpha2+alpha3)
            w3 = alpha3/(alpha1+alpha2+alpha3)

            p1 = (2.d0*v(i-2)-13.d0*v(i-1)+47.d0*v(i)+27.d0*v(i+1)-3.d0*v(i+2))/60.d0
            p2 = (-v(i-1)+3.d0*v(i))/2.d0
            p3 = (v(i)+v(i+1))/2.d0

            hh(i) = (w1/d1)*(p1-d2*p2-d3*p3)+w2*p2+w3*p3
        end do
    end subroutine OCFD_WENOZQ5P

    subroutine OCFD_WENOZQ5M(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3
        real(dp), parameter :: d1 = 0.98d0, d2 = 0.01d0, d3 = 0.01d0

        do i = 1, nx+1
            beta1 = (v(i+2)-8.d0*v(i+1)+8.d0*v(i-1)-v(i-2))**2/144.d0+ &
                    (-11.d0*v(i+2)+174.d0*v(i+1)-326.d0*v(i)+174.d0*v(i-1)-11.d0*v(i-2))**2/15600.d0+ &
                    781.d0*(-v(i+2)+2.d0*v(i+1)-2.d0*v(i-1)+v(i-2))**2/2880.d0+ &
                    1421461.d0*(v(i+2)-4.d0*v(i+1)+6.d0*v(i)-4.d0*v(i-1)+v(i-2))**2/1310400.d0
            beta2 = (v(i+1)-v(i))**2
            beta3 = (v(i)-v(i-1))**2

            tau = ((abs(beta1-beta2)+abs(beta1-beta3))/2.d0)**2
            alpha1 = d1*(1.d0+tau/(ep_wenozq+beta1))
            alpha2 = d2*(1.d0+tau/(ep_wenozq+beta2))
            alpha3 = d3*(1.d0+tau/(ep_wenozq+beta3))
            w1 = alpha1/(alpha1+alpha2+alpha3)
            w2 = alpha2/(alpha1+alpha2+alpha3)
            w3 = alpha3/(alpha1+alpha2+alpha3)

            p1 = (2.d0*v(i+2)-13.d0*v(i+1)+47.d0*v(i)+27.d0*v(i-1)-3.d0*v(i-2))/60.d0
            p2 = (-v(i+1)+3.d0*v(i))/2.d0
            p3 = (v(i)+v(i-1))/2.d0

            hh(i-1) = (w1/d1)*(p1-d2*p2-d3*p3)+w2*p2+w3*p3
        end do
    end subroutine OCFD_WENOZQ5M

    subroutine OCFD_WENOMR5P(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3, q1, q2, q3
        real(dp) :: s0, s1, g01s, g11s, g01, g11, si0, si1, si
        real(dp), parameter :: d12 = 1.d0/11.d0, d22 = 10.d0/11.d0
        real(dp), parameter :: d13 = 1.d0/111.d0, d23 = 10.d0/111.d0, d33 = 100.d0/111.d0

        do i = 0, nx
            s0 = (v(i)-v(i-1))**2; s1 = (v(i+1)-v(i))**2
            if (s0 >= s1) then
                g01s = 1.d0
            else
                g01s = 10.d0
            end if
            g11s = 11.d0-g01s
            g01 = g01s/11.d0
            g11 = 1.d0-g01
            si0 = g01*(1.d0+(s0-s1)**2/(s0+ep_wenomr))
            si1 = g11*(1.d0+(s0-s1)**2/(s1+ep_wenomr))
            si = si0+si1

            beta1 = (si0*(v(i)-v(i-1))+si1*(v(i+1)-v(i)))**2/si**2
            beta2 = 13.d0/12.d0*(v(i-1)-2.d0*v(i)+v(i+1))**2+1.d0/4.d0*(v(i-1)-v(i+1))**2
            beta3 = (v(i-2)-8.d0*v(i-1)+8.d0*v(i+1)-v(i+2))**2/144.d0+ &
                    (-11.d0*v(i-2)+174.d0*v(i-1)-326.d0*v(i)+174.d0*v(i+1)-11.d0*v(i+2))**2/15600.d0+ &
                    781.d0*(-v(i-2)+2.d0*v(i-1)-2.d0*v(i+1)+v(i+2))**2/2880.d0+ &
                    1421461.d0*(v(i-2)-4.d0*v(i-1)+6.d0*v(i)-4.d0*v(i+1)+v(i+2))**2/1310400.d0

            tau = ((abs(beta3-beta1)+abs(beta3-beta2))/2.d0)**2
            alpha1 = d13*(1.d0+tau/(ep_wenomr+beta1))
            alpha2 = d23*(1.d0+tau/(ep_wenomr+beta2))
            alpha3 = d33*(1.d0+tau/(ep_wenomr+beta3))
            w1 = alpha1/(alpha1+alpha2+alpha3)
            w2 = alpha2/(alpha1+alpha2+alpha3)
            w3 = alpha3/(alpha1+alpha2+alpha3)

            q1 = v(i)
            q2 = -1.d0/6.d0*v(i-1)+5.d0/6.d0*v(i)+1.d0/3.d0*v(i+1)
            q3 = (2.d0*v(i-2)-13.d0*v(i-1)+47.d0*v(i)+27.d0*v(i+1)-3.d0*v(i+2))/60.d0

            p1 = q1
            p2 = (1.d0/d22)*q2-(d12/d22)*p1
            p3 = (1.d0/d33)*q3-(d13/d33)*p1-(d23/d33)*p2

            hh(i) = w1*p1+w2*p2+w3*p3
        end do
    end subroutine OCFD_WENOMR5P

    subroutine OCFD_WENOMR5M(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) :: beta1, beta2, beta3, alpha1, alpha2, alpha3, w1, w2, w3, tau
        real(dp) :: p1, p2, p3, q1, q2, q3
        real(dp) :: s0, s1, g01s, g11s, g01, g11, si0, si1, si
        real(dp), parameter :: d12 = 1.d0/11.d0, d22 = 10.d0/11.d0
        real(dp), parameter :: d13 = 1.d0/111.d0, d23 = 10.d0/111.d0, d33 = 100.d0/111.d0

        do i = 1, nx+1
            s0 = (v(i)-v(i+1))**2; s1 = (v(i-1)-v(i))**2
            if (s0 >= s1) then
                g01s = 1.d0
            else
                g01s = 10.d0
            end if
            g11s = 11.d0-g01s
            g01 = g01s/11.d0
            g11 = 1.d0-g01
            si0 = g01*(1.d0+(s0-s1)**2/(s0+ep_wenomr))
            si1 = g11*(1.d0+(s0-s1)**2/(s1+ep_wenomr))
            si = si0+si1

            beta1 = (si0*(v(i)-v(i+1))+si1*(v(i-1)-v(i)))**2/si**2
            beta2 = 13.d0/12.d0*(v(i+1)-2.d0*v(i)+v(i-1))**2+1.d0/4.d0*(v(i+1)-v(i-1))**2
            beta3 = (v(i+2)-8.d0*v(i+1)+8.d0*v(i-1)-v(i-2))**2/144.d0+ &
                    (-11.d0*v(i+2)+174.d0*v(i+1)-326.d0*v(i)+174.d0*v(i-1)-11.d0*v(i-2))**2/15600.d0+ &
                    781.d0*(-v(i+2)+2.d0*v(i+1)-2.d0*v(i-1)+v(i-2))**2/2880.d0+ &
                    1421461.d0*(v(i+2)-4.d0*v(i+1)+6.d0*v(i)-4.d0*v(i-1)+v(i-2))**2/1310400.d0

            tau = ((abs(beta3-beta1)+abs(beta3-beta2))/2.d0)**2
            alpha1 = d13*(1.d0+tau/(ep_wenomr+beta1))
            alpha2 = d23*(1.d0+tau/(ep_wenomr+beta2))
            alpha3 = d33*(1.d0+tau/(ep_wenomr+beta3))
            w1 = alpha1/(alpha1+alpha2+alpha3)
            w2 = alpha2/(alpha1+alpha2+alpha3)
            w3 = alpha3/(alpha1+alpha2+alpha3)

            q1 = v(i)
            q2 = -1.d0/6.d0*v(i+1)+5.d0/6.d0*v(i)+1.d0/3.d0*v(i-1)
            q3 = (2.d0*v(i+2)-13.d0*v(i+1)+47.d0*v(i)+27.d0*v(i-1)-3.d0*v(i-2))/60.d0

            p1 = q1
            p2 = (1.d0/d22)*q2-(d12/d22)*p1
            p3 = (1.d0/d33)*q3-(d13/d33)*p1-(d23/d33)*p2

            hh(i-1) = w1*p1+w2*p2+w3*p3
        end do
    end subroutine OCFD_WENOMR5M

    subroutine OCFD_LINEAR5P(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) :: q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        do i = 0, nx
            q03 = 1.d0/3.d0*v(i)+5.d0/6.d0*v(i+1)-1.d0/6.d0*v(i+2)
            q13 = -1.d0/6.d0*v(i-1)+5.d0/6.d0*v(i)+1.d0/3.d0*v(i+1)
            q23 = 1.d0/3.d0*v(i-2)-7.d0/6.d0*v(i-1)+11.d0/6.d0*v(i)

            hh(i) = C03*q03+C13*q13+C23*q23
        end do
    end subroutine OCFD_LINEAR5P

    subroutine OCFD_LINEAR5M(v, hh, nx, LAP)
        integer :: nx, LAP, i
        real(dp) :: v(1-LAP:nx+LAP), hh(0:nx)
        real(dp) :: q03, q13, q23
        real(dp), parameter :: C03 = 3.d0/10.d0, C13 = 3.d0/5.d0, C23 = 1.d0/10.d0

        do i = 1, nx+1
            q03 = 1.d0/3.d0*v(i)+5.d0/6.d0*v(i-1)-1.d0/6.d0*v(i-2)
            q13 = -1.d0/6.d0*v(i+1)+5.d0/6.d0*v(i)+1.d0/3.d0*v(i-1)
            q23 = 1.d0/3.d0*v(i+2)-7.d0/6.d0*v(i+1)+11.d0/6.d0*v(i)

            hh(i-1) = C03*q03+C13*q13+C23*q23
        end do
    end subroutine OCFD_LINEAR5M

end module WENO_Reconstruction
