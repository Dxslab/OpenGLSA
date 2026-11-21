! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

subroutine split2d_jac(nn, LAP1, c1, d1, u1, v1, A1, A2, fp, fm)
    use flow_para
    implicit none

    integer :: nn, LAP1, i
    real(dp), dimension(1-LAP1:nn+LAP1) :: c1, d1, u1, v1, A1, A2
    real(dp), dimension(1-LAP1:nn+LAP1, 4) :: fp, fm
    real(dp) :: Ak1, Ak2, tmp0, tmp1, tmp2, ss
    real(dp) :: E1, E2, E3, E1P, E2P, E3P, E1M, E2M, E3M
    real(dp) :: vs, uc1, uc2, vc1, vc2, vvc1, vvc2, vv, W2
    real(dp) :: eps, eps2
    real(dp), parameter :: delta_ef = 0.01d0
    real(dp), parameter :: eps_floor = 1.0d-12

    tmp1 = 2.0d0*(Para%gamma-1.0d0)
    tmp2 = (3.0d0-Para%gamma)/(2.0d0*(Para%gamma-1.0d0))

    do i = 1-LAP1, nn+LAP1
        vs = A1(i)*u1(i)+A2(i)*v1(i)
        ss = sqrt(A1(i)*A1(i)+A2(i)*A2(i))
        ak1 = A1(i)/ss
        ak2 = A2(i)/ss

        E1 = vs
        E2 = vs-c1(i)*ss
        E3 = vs+c1(i)*ss

        eps = max(eps_floor, delta_ef*(abs(vs)+c1(i)*ss))
        eps2 = eps*eps

        E1P = 0.5d0*(E1+sqrt(E1*E1+eps2))
        E2P = 0.5d0*(E2+sqrt(E2*E2+eps2))
        E3P = 0.5d0*(E3+sqrt(E3*E3+eps2))
        E1M = E1-E1P
        E2M = E2-E2P
        E3M = E3-E3P

        tmp0 = d1(i)/(2.d0*Para%gamma)
        uc1 = u1(i)-c1(i)*ak1
        uc2 = u1(i)+c1(i)*ak1
        vc1 = v1(i)-c1(i)*ak2
        vc2 = v1(i)+c1(i)*ak2
        vvc1 = (uc1*uc1+vc1*vc1)*0.5d0
        vvc2 = (uc2*uc2+vc2*vc2)*0.5d0
        vv = (Para%gamma-1.d0)*(u1(i)*u1(i)+v1(i)*v1(i))
        W2 = tmp2*c1(i)*c1(i)

        fp(i, 1) = tmp0*(tmp1*E1P+E2P+E3P)
        fp(i, 2) = tmp0*(tmp1*E1P*u1(i)+E2P*uc1+E3P*uc2)
        fp(i, 3) = tmp0*(tmp1*E1P*v1(i)+E2P*vc1+E3P*vc2)
        fp(i, 4) = tmp0*(E1P*vv+E2p*vvc1+E3P*vvc2+W2*(E2P+E3P))

        fm(i, 1) = tmp0*(tmp1*E1M+E2M+E3M)
        fm(i, 2) = tmp0*(tmp1*E1M*u1(i)+E2M*uc1+E3M*uc2)
        fm(i, 3) = tmp0*(tmp1*E1M*v1(i)+E2M*vc1+E3M*vc2)
        fm(i, 4) = tmp0*(E1M*vv+E2M*vvc1+E3M*vvc2+W2*(E2M+E3M))
    end do

end subroutine split2d_jac

subroutine split2d_jac_d(nn, LAP1, c1, d1, u1, v1, p2, d2, u2, v2, w2, A1, A2, fp, fm)
    use flow_para
    implicit none

    integer :: nn, LAP1, i
    real(dp), dimension(1-LAP1:nn+LAP1) :: c1, d1, u1, v1, A1, A2
    real(dp), dimension(1-LAP1:nn+LAP1) :: p2, d2, u2, v2, w2
    real(dp), dimension(1-LAP1:nn+LAP1, 5) :: fp, fm
    real(dp):: Ak1, AK2, tmp0, tmp1, tmp2, ss, vs, &
               E1, E2, E3, E1P, E2P, E3P, E1M, E2M, E3M, &
               q1, uu1
    real(dp):: q2, uu2, f1tmp1, ppdcuu, pmdcuu, &
               f2tmp1, f2tmp2, f2tmp3, &
               f3tmp1, f3tmp2, f3tmp3, &
               f5tmp1, f5tmp2, f5tmp3
    real(dp) :: eps
    real(dp), parameter :: delta_ef = 0.01d0
    real(dp), parameter :: eps_floor = 1.0d-12

    tmp2 = 1.0d0/(Para%gamma-1.0d0)

    do i = 1-LAP1, nn+LAP1
        vs = A1(i)*u1(i)+A2(i)*v1(i)
        ss = sqrt(A1(i)*A1(i)+A2(i)*A2(i))
        ak1 = A1(i)/ss
        ak2 = A2(i)/ss

        E1 = vs
        E2 = vs-c1(i)*ss
        E3 = vs+c1(i)*ss

        eps = max(eps_floor, delta_ef*(abs(vs)+c1(i)*ss))

        E1P = 0.5d0*(E1+sqrt(E1*E1+eps*eps))
        E2P = 0.5d0*(E2+sqrt(E2*E2+eps*eps))
        E3P = 0.5d0*(E3+sqrt(E3*E3+eps*eps))
        E1M = E1-E1P
        E2M = E2-E2P
        E3M = E3-E3P

        tmp0 = 1.0d0/(c1(i)*c1(i))
        tmp1 = c1(i)*c1(i)*tmp2
        q1 = (u1(i)*u1(i)+v1(i)*v1(i))*0.5d0
        q2 = u1(i)*u2(i)+v1(i)*v2(i)
        uu1 = ak1*u1(i)+ak2*v1(i)
        uu2 = ak1*u2(i)+ak2*v2(i)

        f1tmp1 = d2(i)-p2(i)*tmp0
        ppdcuu = (p2(i)+d1(i)*c1(i)*uu2)*tmp0*0.5d0
        pmdcuu = (p2(i)-d1(i)*c1(i)*uu2)*tmp0*0.5d0

        f2tmp1 = u1(i)*f1tmp1+d1(i)*(u2(i)-ak1*uu2)
        f2tmp2 = (u1(i)-c1(i)*ak1)*pmdcuu
        f2tmp3 = (u1(i)+c1(i)*ak1)*ppdcuu

        f3tmp1 = v1(i)*f1tmp1+d1(i)*(v2(i)-ak2*uu2)
        f3tmp2 = (v1(i)-c1(i)*ak2)*pmdcuu
        f3tmp3 = (v1(i)+c1(i)*ak2)*ppdcuu

        f5tmp1 = q1*f1tmp1+d1(i)*(q2-uu1*uu2)
        f5tmp2 = (tmp1+q1-c1(i)*uu1)*pmdcuu
        f5tmp3 = (tmp1+q1+c1(i)*uu1)*ppdcuu

        fp(i, 1) = f1tmp1*E1P+pmdcuu*E2P+ppdcuu*E3P
        fp(i, 2) = f2tmp1*E1P+f2tmp2*E2P+f2tmp3*E3P
        fp(i, 3) = f3tmp1*E1P+f3tmp2*E2P+f3tmp3*E3P
        fp(i, 4) = d1(i)*w2(i)*E1P
        fp(i, 5) = f5tmp1*E1P+f5tmp2*E2P+f5tmp3*E3P

        fm(i, 1) = f1tmp1*E1M+pmdcuu*E2M+ppdcuu*E3M
        fm(i, 2) = f2tmp1*E1M+f2tmp2*E2M+f2tmp3*E3M
        fm(i, 3) = f3tmp1*E1M+f3tmp2*E2M+f3tmp3*E3M
        fm(i, 4) = d1(i)*w2(i)*E1M
        fm(i, 5) = f5tmp1*E1M+f5tmp2*E2M+f5tmp3*E3M
    end do

end subroutine split2d_jac_d

subroutine split2d_jac_llf_d(nn, LAP1, c1, d1, u1, v1, p2, d2, u2, v2, w2, A1, A2, fp, fm)
    use flow_para
    implicit none

    integer :: nn, LAP1, i, m
    real(dp), dimension(1-LAP1:nn+LAP1) :: c1, d1, u1, v1, A1, A2
    real(dp), dimension(1-LAP1:nn+LAP1) :: p2, d2, u2, v2, w2
    real(dp), dimension(1-LAP1:nn+LAP1, 5) :: fp, fm
    real(dp) :: ss, vsb, vsp, Hb, lambda, tmp2
    real(dp) :: q(5), F(5)
    real(dp), parameter :: eps_floor = 1.0d-12

    tmp2 = 1.0d0/(Para%gamma-1.0d0)

    do i = 1-LAP1, nn+LAP1
        ss = sqrt(A1(i)*A1(i)+A2(i)*A2(i))
        vsb = A1(i)*u1(i)+A2(i)*v1(i)          ! baseflow contravariant velocity
        vsp = A1(i)*u2(i)+A2(i)*v2(i)          ! contravariant velocity perturbation
        Hb = c1(i)*c1(i)*tmp2+0.5d0*(u1(i)*u1(i)+v1(i)*v1(i))  ! baseflow total enthalpy

        q(1) = d2(i)
        q(2) = d1(i)*u2(i)+d2(i)*u1(i)
        q(3) = d1(i)*v2(i)+d2(i)*v1(i)
        q(4) = d1(i)*w2(i)
        q(5) = p2(i)*tmp2+0.5d0*d2(i)*(u1(i)*u1(i)+v1(i)*v1(i))+d1(i)*(u1(i)*u2(i)+v1(i)*v2(i))

        F(1) = q(1)*vsb+d1(i)*vsp
        F(2) = q(2)*vsb+d1(i)*u1(i)*vsp+A1(i)*p2(i)
        F(3) = q(3)*vsb+d1(i)*v1(i)*vsp+A2(i)*p2(i)
        F(4) = q(4)*vsb
        F(5) = (q(5)+p2(i))*vsb+d1(i)*Hb*vsp

        lambda = max(abs(vsb)+c1(i)*ss, eps_floor)

        do m = 1, 5
            fp(i, m) = 0.5d0*(F(m)+lambda*q(m))
            fm(i, m) = 0.5d0*(F(m)-lambda*q(m))
        end do
    end do

end subroutine split2d_jac_llf_d
