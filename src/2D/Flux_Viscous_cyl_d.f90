! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Copyright by Li Xinliang
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module Flux_Viscous_cyl_d
    implicit none

contains

    subroutine flux_viscous2d_cyl_d
        use mesh_data
        use flow_data
        use fluc_data
        use Centrals
        use Flux_Viscous_d, only: comput_mu_d, diff_mu_d
        implicit none

        integer :: i, j, m_int
        real(dp) :: inv_y
        real(dp) :: dmufdk_term(4), dmufdi_term(4), muf_term(4)
        real(dp) :: dmudk_term(4), dmudi_term(4), mu_term(4)

        real(dp), parameter :: onethird = 1.0d0/3.0d0
        real(dp), parameter :: twothird = 2.0d0/3.0d0
        real(dp), parameter :: fourthird = 4.0d0/3.0d0
        real(dp), parameter :: eighthird = 8.0d0/3.0d0
        real(dp), parameter :: tenthird = 10.0d0/3.0d0

        m_int = abs(nint(Para%Wavenumber))

        call comput_mu_d

        call OCFD2d_dx0(u_fluc, dufdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(v_fluc, dvfdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(w_fluc, dwfdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx0(T_fluc, dTfdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)

        call OCFD2d_dy0(u_fluc, dufdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(v_fluc, dvfdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(w_fluc, dwfdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(T_fluc, dTfdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)

        ! Mixed k/i derivatives require only the y face halos.
        call mpi_exchange_y2d(dufdk)
        call mpi_exchange_y2d(dvfdk)
        call mpi_exchange_y2d(dwfdk)
        call mpi_exchange_y2d(dTfdk)

        call OCFD2d_dy0(dufdk, d2ufdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(dvfdk, d2vfdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(dwfdk, d2wfdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy0(dTfdk, d2Tfdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)

        call OCFD2d_dx2(u_fluc, d2ufdk2, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx2(v_fluc, d2vfdk2, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx2(w_fluc, d2wfdk2, nx, ny, Scheme%Scheme_Vis, BcI)
        call OCFD2d_dx2(T_fluc, d2Tfdk2, nx, ny, Scheme%Scheme_Vis, BcI)

        call OCFD2d_dy2(u_fluc, d2ufdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy2(v_fluc, d2vfdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy2(w_fluc, d2wfdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        call OCFD2d_dy2(T_fluc, d2Tfdi2, nx, ny, Scheme%Scheme_Vis, BcJ)

        dmufdk = diff_mu_d(dTfdk, T_fluc, dTdk, T)
        dmufdi = diff_mu_d(dTfdi, T_fluc, dTdi, T)
        do j = 1, ny
        do i = 1, nx
            dmufdi_term(1) = Aix(i, j)**2*dudi(i, j)*fourthird+Aix(i, j)*Aiy(i, j)*dvdi(i, j)*onethird+ &
                             Aix(i, j)*Akx(i, j)*dudk(i, j)*fourthird-Aix(i, j)*Aky(i, j)*dvdk(i, j)*twothird+ &
                             Aiy(i, j)**2*dudi(i, j)+Aiy(i, j)*Akx(i, j)*dvdk(i, j)+Aiy(i, j)*Aky(i, j)*dudk(i, j)

            dmufdk_term(1) = Aix(i, j)*Akx(i, j)*dudi(i, j)*fourthird+Aix(i, j)*Aky(i, j)*dvdi(i, j)- &
                             Aiy(i, j)*Akx(i, j)*dvdi(i, j)*twothird+Aiy(i, j)*Aky(i, j)*dudi(i, j)+ &
                             Akx(i, j)**2*dudk(i, j)*fourthird+Akx(i, j)*Aky(i, j)*dvdk(i, j)*onethird+Aky(i, j)**2*dudk(i, j)

            muf_term(1) = Aix(i, j)**2*d2udi2(i, j)*fourthird+ &
                          Aix(i, j)*Aixi(i, j)*dudi(i, j)*fourthird+Aix(i, j)*Aiy(i, j)*d2vdi2(i, j)*onethird- &
                          Aix(i, j)*Aiyi(i, j)*dvdi(i, j)*twothird+Aix(i, j)*Akx(i, j)*d2udkdi(i, j)*eighthird+ &
                          Aix(i, j)*Akxi(i, j)*dudk(i, j)*fourthird+Aix(i, j)*Aky(i, j)*d2vdkdi(i, j)*onethird- &
                          Aix(i, j)*Akyi(i, j)*dvdk(i, j)*twothird+Aixi(i, j)*Aiy(i, j)*dvdi(i, j)+ &
                          Aixk(i, j)*Akx(i, j)*dudi(i, j)*fourthird+Aixk(i, j)*Aky(i, j)*dvdi(i, j)+ &
                          Aiy(i, j)**2*d2udi2(i, j)+Aiy(i, j)*Aiyi(i, j)*dudi(i, j)+ &
                          Aiy(i, j)*Akx(i, j)*d2vdkdi(i, j)*onethird+Aiy(i, j)*Akxi(i, j)*dvdk(i, j)+ &
                          2.0d0*Aiy(i, j)*Aky(i, j)*d2udkdi(i, j)+Aiy(i, j)*Akyi(i, j)*dudk(i, j)- &
                          Aiyk(i, j)*Akx(i, j)*dvdi(i, j)*twothird+Aiyk(i, j)*Aky(i, j)*dudi(i, j)+ &
                          Akx(i, j)**2*d2udk2(i, j)*fourthird+Akx(i, j)*Akxk(i, j)*dudk(i, j)*fourthird+ &
                          Akx(i, j)*Aky(i, j)*d2vdk2(i, j)*onethird-Akx(i, j)*Akyk(i, j)*dvdk(i, j)*twothird+ &
                          Akxk(i, j)*Aky(i, j)*dvdk(i, j)+Aky(i, j)**2*d2udk2(i, j)+Aky(i, j)*Akyk(i, j)*dudk(i, j)

            dmudi_term(1) = Aix(i, j)**2*dufdi(i, j)*fourthird+ &
                            Aix(i, j)*Aiy(i, j)*dvfdi(i, j)*onethird+Aix(i, j)*Akx(i, j)*dufdk(i, j)*fourthird- &
                            Aix(i, j)*Aky(i, j)*dvfdk(i, j)*twothird+Aiy(i, j)**2*dufdi(i, j)+Aiy(i, j)*Akx(i, j)*dvfdk(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dufdk(i, j)

            dmudk_term(1) = Aix(i, j)*Akx(i, j)*dufdi(i, j)*fourthird+Aix(i, j)*Aky(i, j)*dvfdi(i, j)- &
                            Aiy(i, j)*Akx(i, j)*dvfdi(i, j)*twothird+Aiy(i, j)*Aky(i, j)*dufdi(i, j)+ &
                            Akx(i, j)**2*dufdk(i, j)*fourthird+Akx(i, j)*Aky(i, j)*dvfdk(i, j)*onethird+ &
                            Aky(i, j)**2*dufdk(i, j)

            mu_term(1) = Aix(i, j)**2*d2ufdi2(i, j)*fourthird+ &
                         Aix(i, j)*Aixi(i, j)*dufdi(i, j)*fourthird+Aix(i, j)*Aiy(i, j)*d2vfdi2(i, j)*onethird- &
                         Aix(i, j)*Aiyi(i, j)*dvfdi(i, j)*twothird+Aix(i, j)*Akx(i, j)*d2ufdkdi(i, j)*eighthird+ &
                         Aix(i, j)*Akxi(i, j)*dufdk(i, j)*fourthird+Aix(i, j)*Aky(i, j)*d2vfdkdi(i, j)*onethird- &
                         Aix(i, j)*Akyi(i, j)*dvfdk(i, j)*twothird+Aixi(i, j)*Aiy(i, j)*dvfdi(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*dufdi(i, j)*fourthird+Aixk(i, j)*Aky(i, j)*dvfdi(i, j)+ &
                         Aiy(i, j)**2*d2ufdi2(i, j)+Aiy(i, j)*Aiyi(i, j)*dufdi(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2vfdkdi(i, j)*onethird+Aiy(i, j)*Akxi(i, j)*dvfdk(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*d2ufdkdi(i, j)+Aiy(i, j)*Akyi(i, j)*dufdk(i, j)- &
                         Aiyk(i, j)*Akx(i, j)*dvfdi(i, j)*twothird+Aiyk(i, j)*Aky(i, j)*dufdi(i, j)+ &
                         Akx(i, j)**2*d2ufdk2(i, j)*fourthird+Akx(i, j)*Akxk(i, j)*dufdk(i, j)*fourthird+ &
                         Akx(i, j)*Aky(i, j)*d2vfdk2(i, j)*onethird-Akx(i, j)*Akyk(i, j)*dvfdk(i, j)*twothird+ &
                         Akxk(i, j)*Aky(i, j)*dvfdk(i, j)+Aky(i, j)**2*d2ufdk2(i, j)+Aky(i, j)*Akyk(i, j)*dufdk(i, j)

            dmufdi_term(2) = Aix(i, j)**2*dvdi(i, j)+Aix(i, j)*Aiy(i, j)*dudi(i, j)*onethird+ &
                             Aix(i, j)*Akx(i, j)*dvdk(i, j)+Aix(i, j)*Aky(i, j)*dudk(i, j)+Aiy(i, j)**2*dvdi(i, j)*fourthird- &
                             Aiy(i, j)*Akx(i, j)*dudk(i, j)*twothird+Aiy(i, j)*Aky(i, j)*dvdk(i, j)*fourthird

            dmufdk_term(2) = Aix(i, j)*Akx(i, j)*dvdi(i, j)-Aix(i, j)*Aky(i, j)*dudi(i, j)*twothird+ &
                             Aiy(i, j)*Akx(i, j)*dudi(i, j)+Aiy(i, j)*Aky(i, j)*dvdi(i, j)*fourthird+Akx(i, j)**2*dvdk(i, j)+ &
                             Akx(i, j)*Aky(i, j)*dudk(i, j)*onethird+Aky(i, j)**2*dvdk(i, j)*fourthird

            muf_term(2) = Aix(i, j)**2*d2vdi2(i, j)+Aix(i, j)*Aixi(i, j)*dvdi(i, j)+ &
                          Aix(i, j)*Aiy(i, j)*d2udi2(i, j)*onethird+Aix(i, j)*Aiyi(i, j)*dudi(i, j)+ &
                          2.0d0*Aix(i, j)*Akx(i, j)*d2vdkdi(i, j)+Aix(i, j)*Akxi(i, j)*dvdk(i, j)+ &
                          Aix(i, j)*Aky(i, j)*d2udkdi(i, j)*onethird+Aix(i, j)*Akyi(i, j)*dudk(i, j)- &
                          Aixi(i, j)*Aiy(i, j)*dudi(i, j)*twothird+Aixk(i, j)*Akx(i, j)*dvdi(i, j)- &
                          Aixk(i, j)*Aky(i, j)*dudi(i, j)*twothird+Aiy(i, j)**2*d2vdi2(i, j)*fourthird+ &
                          Aiy(i, j)*Aiyi(i, j)*dvdi(i, j)*fourthird+Aiy(i, j)*Akx(i, j)*d2udkdi(i, j)*onethird- &
                          Aiy(i, j)*Akxi(i, j)*dudk(i, j)*twothird+Aiy(i, j)*Aky(i, j)*d2vdkdi(i, j)*eighthird+ &
                          Aiy(i, j)*Akyi(i, j)*dvdk(i, j)*fourthird+Aiyk(i, j)*Akx(i, j)*dudi(i, j)+ &
                          Aiyk(i, j)*Aky(i, j)*dvdi(i, j)*fourthird+Akx(i, j)**2*d2vdk2(i, j)+ &
                          Akx(i, j)*Akxk(i, j)*dvdk(i, j)+Akx(i, j)*Aky(i, j)*d2udk2(i, j)*onethird+ &
                          Akx(i, j)*Akyk(i, j)*dudk(i, j)-Akxk(i, j)*Aky(i, j)*dudk(i, j)*twothird+ &
                          Aky(i, j)**2*d2vdk2(i, j)*fourthird+Aky(i, j)*Akyk(i, j)*dvdk(i, j)*fourthird

            dmudi_term(2) = Aix(i, j)**2*dvfdi(i, j)+Aix(i, j)*Aiy(i, j)*dufdi(i, j)*onethird+ &
                           Aix(i, j)*Akx(i, j)*dvfdk(i, j)+Aix(i, j)*Aky(i, j)*dufdk(i, j)+Aiy(i, j)**2*dvfdi(i, j)*fourthird- &
                            Aiy(i, j)*Akx(i, j)*dufdk(i, j)*twothird+Aiy(i, j)*Aky(i, j)*dvfdk(i, j)*fourthird

            dmudk_term(2) = Aix(i, j)*Akx(i, j)*dvfdi(i, j)-Aix(i, j)*Aky(i, j)*dufdi(i, j)*twothird+ &
                           Aiy(i, j)*Akx(i, j)*dufdi(i, j)+Aiy(i, j)*Aky(i, j)*dvfdi(i, j)*fourthird+Akx(i, j)**2*dvfdk(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dufdk(i, j)*onethird+Aky(i, j)**2*dvfdk(i, j)*fourthird

            mu_term(2) = Aix(i, j)**2*d2vfdi2(i, j)+Aix(i, j)*Aixi(i, j)*dvfdi(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2ufdi2(i, j)*onethird+Aix(i, j)*Aiyi(i, j)*dufdi(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*d2vfdkdi(i, j)+Aix(i, j)*Akxi(i, j)*dvfdk(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2ufdkdi(i, j)*onethird+Aix(i, j)*Akyi(i, j)*dufdk(i, j)- &
                         Aixi(i, j)*Aiy(i, j)*dufdi(i, j)*twothird+Aixk(i, j)*Akx(i, j)*dvfdi(i, j)- &
                         Aixk(i, j)*Aky(i, j)*dufdi(i, j)*twothird+Aiy(i, j)**2*d2vfdi2(i, j)*fourthird+ &
                         Aiy(i, j)*Aiyi(i, j)*dvfdi(i, j)*fourthird+Aiy(i, j)*Akx(i, j)*d2ufdkdi(i, j)*onethird- &
                         Aiy(i, j)*Akxi(i, j)*dufdk(i, j)*twothird+Aiy(i, j)*Aky(i, j)*d2vfdkdi(i, j)*eighthird+ &
                         Aiy(i, j)*Akyi(i, j)*dvfdk(i, j)*fourthird+Aiyk(i, j)*Akx(i, j)*dufdi(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*dvfdi(i, j)*fourthird+Akx(i, j)**2*d2vfdk2(i, j)+ &
                         Akx(i, j)*Akxk(i, j)*dvfdk(i, j)+Akx(i, j)*Aky(i, j)*d2ufdk2(i, j)*onethird+ &
                         Akx(i, j)*Akyk(i, j)*dufdk(i, j)-Akxk(i, j)*Aky(i, j)*dufdk(i, j)*twothird+ &
                         Aky(i, j)**2*d2vfdk2(i, j)*fourthird+Aky(i, j)*Akyk(i, j)*dvfdk(i, j)*fourthird

            dmufdi_term(3) = 0.0d0

            dmufdk_term(3) = 0.0d0

            muf_term(3) = 0.0d0

            dmudi_term(3) = Aix(i, j)**2*dwfdi(i, j)+Aix(i, j)*Akx(i, j)*dwfdk(i, j)+ &
                            Aiy(i, j)**2*dwfdi(i, j)+Aiy(i, j)*Aky(i, j)*dwfdk(i, j)

            dmudk_term(3) = Aix(i, j)*Akx(i, j)*dwfdi(i, j)+Aiy(i, j)*Aky(i, j)*dwfdi(i, j)+ &
                            Akx(i, j)**2*dwfdk(i, j)+Aky(i, j)**2*dwfdk(i, j)

            mu_term(3) = Aix(i, j)**2*d2wfdi2(i, j)+Aix(i, j)*Aixi(i, j)*dwfdi(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*d2wfdkdi(i, j)+Aix(i, j)*Akxi(i, j)*dwfdk(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*dwfdi(i, j)+Aiy(i, j)**2*d2wfdi2(i, j)+Aiy(i, j)*Aiyi(i, j)*dwfdi(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*d2wfdkdi(i, j)+Aiy(i, j)*Akyi(i, j)*dwfdk(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*dwfdi(i, j)+Akx(i, j)**2*d2wfdk2(i, j)+Akx(i, j)*Akxk(i, j)*dwfdk(i, j)+ &
                         Aky(i, j)**2*d2wfdk2(i, j)+Aky(i, j)*Akyk(i, j)*dwfdk(i, j)

            dmufdi_term(4) = Aix(i, j)**2*Para%Cp*dTdi(i, j)/Para%Pr+ &
                             Aix(i, j)**2*dudi(i, j)*fourthird*u(i, j)+Aix(i, j)**2*dvdi(i, j)*v(i, j)+ &
                             Aix(i, j)*Aiy(i, j)*dudi(i, j)*onethird*v(i, j)+Aix(i, j)*Aiy(i, j)*dvdi(i, j)*onethird*u(i, j)+ &
                             Aix(i, j)*Akx(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+Aix(i, j)*Akx(i, j)*dudk(i, j)*fourthird*u(i, j)+ &
                             Aix(i, j)*Akx(i, j)*dvdk(i, j)*v(i, j)+Aix(i, j)*Aky(i, j)*dudk(i, j)*v(i, j)- &
                             Aix(i, j)*Aky(i, j)*dvdk(i, j)*twothird*u(i, j)+Aiy(i, j)**2*Para%Cp*dTdi(i, j)/Para%Pr+ &
                             Aiy(i, j)**2*dudi(i, j)*u(i, j)+Aiy(i, j)**2*dvdi(i, j)*fourthird*v(i, j)- &
                             Aiy(i, j)*Akx(i, j)*dudk(i, j)*twothird*v(i, j)+Aiy(i, j)*Akx(i, j)*dvdk(i, j)*u(i, j)+ &
                             Aiy(i, j)*Aky(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+Aiy(i, j)*Aky(i, j)*dudk(i, j)*u(i, j)+ &
                             Aiy(i, j)*Aky(i, j)*dvdk(i, j)*fourthird*v(i, j)

            dmufdk_term(4) = Aix(i, j)*Akx(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                             Aix(i, j)*Akx(i, j)*dudi(i, j)*fourthird*u(i, j)+Aix(i, j)*Akx(i, j)*dvdi(i, j)*v(i, j)- &
                             Aix(i, j)*Aky(i, j)*dudi(i, j)*twothird*v(i, j)+Aix(i, j)*Aky(i, j)*dvdi(i, j)*u(i, j)+ &
                             Aiy(i, j)*Akx(i, j)*dudi(i, j)*v(i, j)-Aiy(i, j)*Akx(i, j)*dvdi(i, j)*twothird*u(i, j)+ &
                             Aiy(i, j)*Aky(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+Aiy(i, j)*Aky(i, j)*dudi(i, j)*u(i, j)+ &
                             Aiy(i, j)*Aky(i, j)*dvdi(i, j)*fourthird*v(i, j)+Akx(i, j)**2*Para%Cp*dTdk(i, j)/Para%Pr+ &
                             Akx(i, j)**2*dudk(i, j)*fourthird*u(i, j)+Akx(i, j)**2*dvdk(i, j)*v(i, j)+ &
                             Akx(i, j)*Aky(i, j)*dudk(i, j)*onethird*v(i, j)+Akx(i, j)*Aky(i, j)*dvdk(i, j)*onethird*u(i, j)+ &
                             Aky(i, j)**2*Para%Cp*dTdk(i, j)/Para%Pr+Aky(i, j)**2*dudk(i, j)*u(i, j)+ &
                             Aky(i, j)**2*dvdk(i, j)*fourthird*v(i, j)

            muf_term(4) = Aix(i, j)**2*Para%Cp*d2Tdi2(i, j)/Para%Pr+ &
                          Aix(i, j)**2*d2udi2(i, j)*fourthird*u(i, j)+Aix(i, j)**2*d2vdi2(i, j)*v(i, j)+ &
                          Aix(i, j)**2*dudi(i, j)**2*fourthird+Aix(i, j)**2*dvdi(i, j)**2+ &
                          Aix(i, j)*Aixi(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+Aix(i, j)*Aixi(i, j)*dudi(i, j)*fourthird*u(i, j)+ &
                          Aix(i, j)*Aixi(i, j)*dvdi(i, j)*v(i, j)+Aix(i, j)*Aiy(i, j)*d2udi2(i, j)*onethird*v(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2vdi2(i, j)*onethird*u(i, j)+Aix(i, j)*Aiy(i, j)*dudi(i, j)*dvdi(i, j)*twothird+ &
                          Aix(i, j)*Aiyi(i, j)*dudi(i, j)*v(i, j)-Aix(i, j)*Aiyi(i, j)*dvdi(i, j)*twothird*u(i, j)+ &
                          2.0d0*Aix(i, j)*Akx(i, j)*Para%Cp*d2Tdkdi(i, j)/Para%Pr+ &
                          Aix(i, j)*Akx(i, j)*d2udkdi(i, j)*eighthird*u(i, j)+2.0d0*Aix(i, j)*Akx(i, j)*d2vdkdi(i, j)*v(i, j)+ &
                          Aix(i, j)*Akx(i, j)*dudi(i, j)*dudk(i, j)*eighthird+2.0d0*Aix(i, j)*Akx(i, j)*dvdi(i, j)*dvdk(i, j)+ &
                          Aix(i, j)*Akxi(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+Aix(i, j)*Akxi(i, j)*dudk(i, j)*fourthird*u(i, j)+ &
                          Aix(i, j)*Akxi(i, j)*dvdk(i, j)*v(i, j)+Aix(i, j)*Aky(i, j)*d2udkdi(i, j)*onethird*v(i, j)+ &
                          Aix(i, j)*Aky(i, j)*d2vdkdi(i, j)*onethird*u(i, j)- &
                          Aix(i, j)*Aky(i, j)*dudi(i, j)*dvdk(i, j)*fourthird+2.0d0*Aix(i, j)*Aky(i, j)*dudk(i, j)*dvdi(i, j)+ &
                          Aix(i, j)*Akyi(i, j)*dudk(i, j)*v(i, j)-Aix(i, j)*Akyi(i, j)*dvdk(i, j)*twothird*u(i, j)- &
                          Aixi(i, j)*Aiy(i, j)*dudi(i, j)*twothird*v(i, j)+Aixi(i, j)*Aiy(i, j)*dvdi(i, j)*u(i, j)+ &
                          Aixk(i, j)*Akx(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+Aixk(i, j)*Akx(i, j)*dudi(i, j)*fourthird*u(i, j)+ &
                          Aixk(i, j)*Akx(i, j)*dvdi(i, j)*v(i, j)-Aixk(i, j)*Aky(i, j)*dudi(i, j)*twothird*v(i, j)+ &
                          Aixk(i, j)*Aky(i, j)*dvdi(i, j)*u(i, j)+Aiy(i, j)**2*Para%Cp*d2Tdi2(i, j)/Para%Pr+ &
                          Aiy(i, j)**2*d2udi2(i, j)*u(i, j)+Aiy(i, j)**2*d2vdi2(i, j)*fourthird*v(i, j)+ &
                          Aiy(i, j)**2*dudi(i, j)**2+Aiy(i, j)**2*dvdi(i, j)**2*fourthird+ &
                          Aiy(i, j)*Aiyi(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+Aiy(i, j)*Aiyi(i, j)*dudi(i, j)*u(i, j)+ &
                         Aiy(i, j)*Aiyi(i, j)*dvdi(i, j)*fourthird*v(i, j)+Aiy(i, j)*Akx(i, j)*d2udkdi(i, j)*onethird*v(i, j)+ &
                          Aiy(i, j)*Akx(i, j)*d2vdkdi(i, j)*onethird*u(i, j)+2.0d0*Aiy(i, j)*Akx(i, j)*dudi(i, j)*dvdk(i, j)- &
                         Aiy(i, j)*Akx(i, j)*dudk(i, j)*dvdi(i, j)*fourthird-Aiy(i, j)*Akxi(i, j)*dudk(i, j)*twothird*v(i, j)+ &
                          Aiy(i, j)*Akxi(i, j)*dvdk(i, j)*u(i, j)+2.0d0*Aiy(i, j)*Aky(i, j)*Para%Cp*d2Tdkdi(i, j)/Para%Pr+ &
                          2.0d0*Aiy(i, j)*Aky(i, j)*d2udkdi(i, j)*u(i, j)+Aiy(i, j)*Aky(i, j)*d2vdkdi(i, j)*eighthird*v(i, j)+ &
                          2.0d0*Aiy(i, j)*Aky(i, j)*dudi(i, j)*dudk(i, j)+Aiy(i, j)*Aky(i, j)*dvdi(i, j)*dvdk(i, j)*eighthird+ &
                          Aiy(i, j)*Akyi(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+Aiy(i, j)*Akyi(i, j)*dudk(i, j)*u(i, j)+ &
                          Aiy(i, j)*Akyi(i, j)*dvdk(i, j)*fourthird*v(i, j)+Aiyk(i, j)*Akx(i, j)*dudi(i, j)*v(i, j)- &
                          Aiyk(i, j)*Akx(i, j)*dvdi(i, j)*twothird*u(i, j)+Aiyk(i, j)*Aky(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                          Aiyk(i, j)*Aky(i, j)*dudi(i, j)*u(i, j)+Aiyk(i, j)*Aky(i, j)*dvdi(i, j)*fourthird*v(i, j)+ &
                          Akx(i, j)**2*Para%Cp*d2Tdk2(i, j)/Para%Pr+Akx(i, j)**2*d2udk2(i, j)*fourthird*u(i, j)+ &
                          Akx(i, j)**2*d2vdk2(i, j)*v(i, j)+Akx(i, j)**2*dudk(i, j)**2*fourthird+Akx(i, j)**2*dvdk(i, j)**2+ &
                          Akx(i, j)*Akxk(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+Akx(i, j)*Akxk(i, j)*dudk(i, j)*fourthird*u(i, j)+ &
                          Akx(i, j)*Akxk(i, j)*dvdk(i, j)*v(i, j)+Akx(i, j)*Aky(i, j)*d2udk2(i, j)*onethird*v(i, j)+ &
                         Akx(i, j)*Aky(i, j)*d2vdk2(i, j)*onethird*u(i, j)+Akx(i, j)*Aky(i, j)*dudk(i, j)*dvdk(i, j)*twothird+ &
                          Akx(i, j)*Akyk(i, j)*dudk(i, j)*v(i, j)-Akx(i, j)*Akyk(i, j)*dvdk(i, j)*twothird*u(i, j)- &
                          Akxk(i, j)*Aky(i, j)*dudk(i, j)*twothird*v(i, j)+Akxk(i, j)*Aky(i, j)*dvdk(i, j)*u(i, j)+ &
                          Aky(i, j)**2*Para%Cp*d2Tdk2(i, j)/Para%Pr+Aky(i, j)**2*d2udk2(i, j)*u(i, j)+ &
                          Aky(i, j)**2*d2vdk2(i, j)*fourthird*v(i, j)+Aky(i, j)**2*dudk(i, j)**2+ &
                          Aky(i, j)**2*dvdk(i, j)**2*fourthird+Aky(i, j)*Akyk(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                          Aky(i, j)*Akyk(i, j)*dudk(i, j)*u(i, j)+Aky(i, j)*Akyk(i, j)*dvdk(i, j)*fourthird*v(i, j)

            dmudi_term(4) = Aix(i, j)**2*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                            Aix(i, j)**2*dudi(i, j)*fourthird*u_fluc(i, j)+Aix(i, j)**2*dufdi(i, j)*fourthird*u(i, j)+ &
                            Aix(i, j)**2*dvdi(i, j)*v_fluc(i, j)+Aix(i, j)**2*dvfdi(i, j)*v(i, j)+ &
                            Aix(i, j)*Aiy(i, j)*dudi(i, j)*onethird*v_fluc(i, j)+ &
                            Aix(i, j)*Aiy(i, j)*dufdi(i, j)*onethird*v(i, j)+ &
                            Aix(i, j)*Aiy(i, j)*dvdi(i, j)*onethird*u_fluc(i, j)+ &
                            Aix(i, j)*Aiy(i, j)*dvfdi(i, j)*onethird*u(i, j)+Aix(i, j)*Akx(i, j)*Para%Cp*dTfdk(i, j)/Para%Pr+ &
                            Aix(i, j)*Akx(i, j)*dudk(i, j)*fourthird*u_fluc(i, j)+ &
                            Aix(i, j)*Akx(i, j)*dufdk(i, j)*fourthird*u(i, j)+Aix(i, j)*Akx(i, j)*dvdk(i, j)*v_fluc(i, j)+ &
                            Aix(i, j)*Akx(i, j)*dvfdk(i, j)*v(i, j)+Aix(i, j)*Aky(i, j)*dudk(i, j)*v_fluc(i, j)+ &
                            Aix(i, j)*Aky(i, j)*dufdk(i, j)*v(i, j)-Aix(i, j)*Aky(i, j)*dvdk(i, j)*twothird*u_fluc(i, j)- &
                            Aix(i, j)*Aky(i, j)*dvfdk(i, j)*twothird*u(i, j)+Aiy(i, j)**2*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                            Aiy(i, j)**2*dudi(i, j)*u_fluc(i, j)+Aiy(i, j)**2*dufdi(i, j)*u(i, j)+ &
                            Aiy(i, j)**2*dvdi(i, j)*fourthird*v_fluc(i, j)+Aiy(i, j)**2*dvfdi(i, j)*fourthird*v(i, j)- &
                            Aiy(i, j)*Akx(i, j)*dudk(i, j)*twothird*v_fluc(i, j)- &
                            Aiy(i, j)*Akx(i, j)*dufdk(i, j)*twothird*v(i, j)+Aiy(i, j)*Akx(i, j)*dvdk(i, j)*u_fluc(i, j)+ &
                            Aiy(i, j)*Akx(i, j)*dvfdk(i, j)*u(i, j)+Aiy(i, j)*Aky(i, j)*Para%Cp*dTfdk(i, j)/Para%Pr+ &
                            Aiy(i, j)*Aky(i, j)*dudk(i, j)*u_fluc(i, j)+Aiy(i, j)*Aky(i, j)*dufdk(i, j)*u(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvdk(i, j)*fourthird*v_fluc(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvfdk(i, j)*fourthird*v(i, j)

            dmudk_term(4) = Aix(i, j)*Akx(i, j)*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                            Aix(i, j)*Akx(i, j)*dudi(i, j)*fourthird*u_fluc(i, j)+ &
                            Aix(i, j)*Akx(i, j)*dufdi(i, j)*fourthird*u(i, j)+Aix(i, j)*Akx(i, j)*dvdi(i, j)*v_fluc(i, j)+ &
                            Aix(i, j)*Akx(i, j)*dvfdi(i, j)*v(i, j)-Aix(i, j)*Aky(i, j)*dudi(i, j)*twothird*v_fluc(i, j)- &
                            Aix(i, j)*Aky(i, j)*dufdi(i, j)*twothird*v(i, j)+Aix(i, j)*Aky(i, j)*dvdi(i, j)*u_fluc(i, j)+ &
                            Aix(i, j)*Aky(i, j)*dvfdi(i, j)*u(i, j)+Aiy(i, j)*Akx(i, j)*dudi(i, j)*v_fluc(i, j)+ &
                            Aiy(i, j)*Akx(i, j)*dufdi(i, j)*v(i, j)-Aiy(i, j)*Akx(i, j)*dvdi(i, j)*twothird*u_fluc(i, j)- &
                            Aiy(i, j)*Akx(i, j)*dvfdi(i, j)*twothird*u(i, j)+Aiy(i, j)*Aky(i, j)*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                            Aiy(i, j)*Aky(i, j)*dudi(i, j)*u_fluc(i, j)+Aiy(i, j)*Aky(i, j)*dufdi(i, j)*u(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvdi(i, j)*fourthird*v_fluc(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvfdi(i, j)*fourthird*v(i, j)+Akx(i, j)**2*Para%Cp*dTfdk(i, j)/Para%Pr+ &
                            Akx(i, j)**2*dudk(i, j)*fourthird*u_fluc(i, j)+Akx(i, j)**2*dufdk(i, j)*fourthird*u(i, j)+ &
                            Akx(i, j)**2*dvdk(i, j)*v_fluc(i, j)+Akx(i, j)**2*dvfdk(i, j)*v(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dudk(i, j)*onethird*v_fluc(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dufdk(i, j)*onethird*v(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dvdk(i, j)*onethird*u_fluc(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dvfdk(i, j)*onethird*u(i, j)+Aky(i, j)**2*Para%Cp*dTfdk(i, j)/Para%Pr+ &
                            Aky(i, j)**2*dudk(i, j)*u_fluc(i, j)+Aky(i, j)**2*dufdk(i, j)*u(i, j)+ &
                            Aky(i, j)**2*dvdk(i, j)*fourthird*v_fluc(i, j)+Aky(i, j)**2*dvfdk(i, j)*fourthird*v(i, j)

            mu_term(4) = Aix(i, j)**2*Para%Cp*d2Tfdi2(i, j)/Para%Pr+ &
                         Aix(i, j)**2*d2udi2(i, j)*fourthird*u_fluc(i, j)+Aix(i, j)**2*d2ufdi2(i, j)*fourthird*u(i, j)+ &
                         Aix(i, j)**2*d2vdi2(i, j)*v_fluc(i, j)+Aix(i, j)**2*d2vfdi2(i, j)*v(i, j)+ &
                         Aix(i, j)**2*dudi(i, j)*dufdi(i, j)*eighthird+2.0d0*Aix(i, j)**2*dvdi(i, j)*dvfdi(i, j)+ &
                         Aix(i, j)*Aixi(i, j)*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                         Aix(i, j)*Aixi(i, j)*dudi(i, j)*fourthird*u_fluc(i, j)+ &
                         Aix(i, j)*Aixi(i, j)*dufdi(i, j)*fourthird*u(i, j)+Aix(i, j)*Aixi(i, j)*dvdi(i, j)*v_fluc(i, j)+ &
                         Aix(i, j)*Aixi(i, j)*dvfdi(i, j)*v(i, j)+Aix(i, j)*Aiy(i, j)*d2udi2(i, j)*onethird*v_fluc(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2ufdi2(i, j)*onethird*v(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2vdi2(i, j)*onethird*u_fluc(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2vfdi2(i, j)*onethird*u(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*dudi(i, j)*dvfdi(i, j)*twothird+ &
                         Aix(i, j)*Aiy(i, j)*dufdi(i, j)*dvdi(i, j)*twothird+Aix(i, j)*Aiyi(i, j)*dudi(i, j)*v_fluc(i, j)+ &
                         Aix(i, j)*Aiyi(i, j)*dufdi(i, j)*v(i, j)-Aix(i, j)*Aiyi(i, j)*dvdi(i, j)*twothird*u_fluc(i, j)- &
                         Aix(i, j)*Aiyi(i, j)*dvfdi(i, j)*twothird*u(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*Para%Cp*d2Tfdkdi(i, j)/Para%Pr+ &
                         Aix(i, j)*Akx(i, j)*d2udkdi(i, j)*eighthird*u_fluc(i, j)+ &
                         Aix(i, j)*Akx(i, j)*d2ufdkdi(i, j)*eighthird*u(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*d2vdkdi(i, j)*v_fluc(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*d2vfdkdi(i, j)*v(i, j)+ &
                         Aix(i, j)*Akx(i, j)*dudi(i, j)*dufdk(i, j)*eighthird+ &
                         Aix(i, j)*Akx(i, j)*dudk(i, j)*dufdi(i, j)*eighthird+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*dvdi(i, j)*dvfdk(i, j)+2.0d0*Aix(i, j)*Akx(i, j)*dvdk(i, j)*dvfdi(i, j)+ &
                         Aix(i, j)*Akxi(i, j)*Para%Cp*dTfdk(i, j)/Para%Pr+ &
                         Aix(i, j)*Akxi(i, j)*dudk(i, j)*fourthird*u_fluc(i, j)+ &
                         Aix(i, j)*Akxi(i, j)*dufdk(i, j)*fourthird*u(i, j)+Aix(i, j)*Akxi(i, j)*dvdk(i, j)*v_fluc(i, j)+ &
                         Aix(i, j)*Akxi(i, j)*dvfdk(i, j)*v(i, j)+Aix(i, j)*Aky(i, j)*d2udkdi(i, j)*onethird*v_fluc(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2ufdkdi(i, j)*onethird*v(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2vdkdi(i, j)*onethird*u_fluc(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2vfdkdi(i, j)*onethird*u(i, j)- &
                         Aix(i, j)*Aky(i, j)*dudi(i, j)*dvfdk(i, j)*fourthird+ &
                         2.0d0*Aix(i, j)*Aky(i, j)*dudk(i, j)*dvfdi(i, j)- &
                         Aix(i, j)*Aky(i, j)*dufdi(i, j)*dvdk(i, j)*fourthird+ &
                         2.0d0*Aix(i, j)*Aky(i, j)*dufdk(i, j)*dvdi(i, j)+Aix(i, j)*Akyi(i, j)*dudk(i, j)*v_fluc(i, j)+ &
                         Aix(i, j)*Akyi(i, j)*dufdk(i, j)*v(i, j)-Aix(i, j)*Akyi(i, j)*dvdk(i, j)*twothird*u_fluc(i, j)- &
                         Aix(i, j)*Akyi(i, j)*dvfdk(i, j)*twothird*u(i, j)- &
                         Aixi(i, j)*Aiy(i, j)*dudi(i, j)*twothird*v_fluc(i, j)- &
                         Aixi(i, j)*Aiy(i, j)*dufdi(i, j)*twothird*v(i, j)+Aixi(i, j)*Aiy(i, j)*dvdi(i, j)*u_fluc(i, j)+ &
                         Aixi(i, j)*Aiy(i, j)*dvfdi(i, j)*u(i, j)+Aixk(i, j)*Akx(i, j)*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                         Aixk(i, j)*Akx(i, j)*dudi(i, j)*fourthird*u_fluc(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*dufdi(i, j)*fourthird*u(i, j)+Aixk(i, j)*Akx(i, j)*dvdi(i, j)*v_fluc(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*dvfdi(i, j)*v(i, j)-Aixk(i, j)*Aky(i, j)*dudi(i, j)*twothird*v_fluc(i, j)- &
                         Aixk(i, j)*Aky(i, j)*dufdi(i, j)*twothird*v(i, j)+Aixk(i, j)*Aky(i, j)*dvdi(i, j)*u_fluc(i, j)+ &
                         Aixk(i, j)*Aky(i, j)*dvfdi(i, j)*u(i, j)+Aiy(i, j)**2*Para%Cp*d2Tfdi2(i, j)/Para%Pr+ &
                         Aiy(i, j)**2*d2udi2(i, j)*u_fluc(i, j)+Aiy(i, j)**2*d2ufdi2(i, j)*u(i, j)+ &
                         Aiy(i, j)**2*d2vdi2(i, j)*fourthird*v_fluc(i, j)+Aiy(i, j)**2*d2vfdi2(i, j)*fourthird*v(i, j)+ &
                         2.0d0*Aiy(i, j)**2*dudi(i, j)*dufdi(i, j)+Aiy(i, j)**2*dvdi(i, j)*dvfdi(i, j)*eighthird+ &
                         Aiy(i, j)*Aiyi(i, j)*Para%Cp*dTfdi(i, j)/Para%Pr+Aiy(i, j)*Aiyi(i, j)*dudi(i, j)*u_fluc(i, j)+ &
                         Aiy(i, j)*Aiyi(i, j)*dufdi(i, j)*u(i, j)+Aiy(i, j)*Aiyi(i, j)*dvdi(i, j)*fourthird*v_fluc(i, j)+ &
                         Aiy(i, j)*Aiyi(i, j)*dvfdi(i, j)*fourthird*v(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2udkdi(i, j)*onethird*v_fluc(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2ufdkdi(i, j)*onethird*v(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2vdkdi(i, j)*onethird*u_fluc(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2vfdkdi(i, j)*onethird*u(i, j)+2.0d0*Aiy(i, j)*Akx(i, j)*dudi(i, j)*dvfdk(i, j)- &
                         Aiy(i, j)*Akx(i, j)*dudk(i, j)*dvfdi(i, j)*fourthird+ &
                         2.0d0*Aiy(i, j)*Akx(i, j)*dufdi(i, j)*dvdk(i, j)- &
                         Aiy(i, j)*Akx(i, j)*dufdk(i, j)*dvdi(i, j)*fourthird- &
                         Aiy(i, j)*Akxi(i, j)*dudk(i, j)*twothird*v_fluc(i, j)- &
                         Aiy(i, j)*Akxi(i, j)*dufdk(i, j)*twothird*v(i, j)+Aiy(i, j)*Akxi(i, j)*dvdk(i, j)*u_fluc(i, j)+ &
                         Aiy(i, j)*Akxi(i, j)*dvfdk(i, j)*u(i, j)+2.0d0*Aiy(i, j)*Aky(i, j)*Para%Cp*d2Tfdkdi(i, j)/Para%Pr+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*d2udkdi(i, j)*u_fluc(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*d2ufdkdi(i, j)*u(i, j)+ &
                         Aiy(i, j)*Aky(i, j)*d2vdkdi(i, j)*eighthird*v_fluc(i, j)+ &
                         Aiy(i, j)*Aky(i, j)*d2vfdkdi(i, j)*eighthird*v(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*dudi(i, j)*dufdk(i, j)+2.0d0*Aiy(i, j)*Aky(i, j)*dudk(i, j)*dufdi(i, j)+ &
                         Aiy(i, j)*Aky(i, j)*dvdi(i, j)*dvfdk(i, j)*eighthird+ &
                         Aiy(i, j)*Aky(i, j)*dvdk(i, j)*dvfdi(i, j)*eighthird+ &
                         Aiy(i, j)*Akyi(i, j)*Para%Cp*dTfdk(i, j)/Para%Pr+Aiy(i, j)*Akyi(i, j)*dudk(i, j)*u_fluc(i, j)+ &
                         Aiy(i, j)*Akyi(i, j)*dufdk(i, j)*u(i, j)+Aiy(i, j)*Akyi(i, j)*dvdk(i, j)*fourthird*v_fluc(i, j)+ &
                         Aiy(i, j)*Akyi(i, j)*dvfdk(i, j)*fourthird*v(i, j)+Aiyk(i, j)*Akx(i, j)*dudi(i, j)*v_fluc(i, j)+ &
                         Aiyk(i, j)*Akx(i, j)*dufdi(i, j)*v(i, j)-Aiyk(i, j)*Akx(i, j)*dvdi(i, j)*twothird*u_fluc(i, j)- &
                         Aiyk(i, j)*Akx(i, j)*dvfdi(i, j)*twothird*u(i, j)+Aiyk(i, j)*Aky(i, j)*Para%Cp*dTfdi(i, j)/Para%Pr+ &
                         Aiyk(i, j)*Aky(i, j)*dudi(i, j)*u_fluc(i, j)+Aiyk(i, j)*Aky(i, j)*dufdi(i, j)*u(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*dvdi(i, j)*fourthird*v_fluc(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*dvfdi(i, j)*fourthird*v(i, j)+Akx(i, j)**2*Para%Cp*d2Tfdk2(i, j)/Para%Pr+ &
                         Akx(i, j)**2*d2udk2(i, j)*fourthird*u_fluc(i, j)+Akx(i, j)**2*d2ufdk2(i, j)*fourthird*u(i, j)+ &
                         Akx(i, j)**2*d2vdk2(i, j)*v_fluc(i, j)+Akx(i, j)**2*d2vfdk2(i, j)*v(i, j)+ &
                         Akx(i, j)**2*dudk(i, j)*dufdk(i, j)*eighthird+2.0d0*Akx(i, j)**2*dvdk(i, j)*dvfdk(i, j)+ &
                         Akx(i, j)*Akxk(i, j)*Para%Cp*dTfdk(i, j)/Para%Pr+ &
                         Akx(i, j)*Akxk(i, j)*dudk(i, j)*fourthird*u_fluc(i, j)+ &
                         Akx(i, j)*Akxk(i, j)*dufdk(i, j)*fourthird*u(i, j)+Akx(i, j)*Akxk(i, j)*dvdk(i, j)*v_fluc(i, j)+ &
                         Akx(i, j)*Akxk(i, j)*dvfdk(i, j)*v(i, j)+Akx(i, j)*Aky(i, j)*d2udk2(i, j)*onethird*v_fluc(i, j)+ &
                         Akx(i, j)*Aky(i, j)*d2ufdk2(i, j)*onethird*v(i, j)+ &
                         Akx(i, j)*Aky(i, j)*d2vdk2(i, j)*onethird*u_fluc(i, j)+ &
                         Akx(i, j)*Aky(i, j)*d2vfdk2(i, j)*onethird*u(i, j)+ &
                         Akx(i, j)*Aky(i, j)*dudk(i, j)*dvfdk(i, j)*twothird+ &
                         Akx(i, j)*Aky(i, j)*dufdk(i, j)*dvdk(i, j)*twothird+Akx(i, j)*Akyk(i, j)*dudk(i, j)*v_fluc(i, j)+ &
                         Akx(i, j)*Akyk(i, j)*dufdk(i, j)*v(i, j)-Akx(i, j)*Akyk(i, j)*dvdk(i, j)*twothird*u_fluc(i, j)- &
                         Akx(i, j)*Akyk(i, j)*dvfdk(i, j)*twothird*u(i, j)- &
                         Akxk(i, j)*Aky(i, j)*dudk(i, j)*twothird*v_fluc(i, j)- &
                         Akxk(i, j)*Aky(i, j)*dufdk(i, j)*twothird*v(i, j)+Akxk(i, j)*Aky(i, j)*dvdk(i, j)*u_fluc(i, j)+ &
                         Akxk(i, j)*Aky(i, j)*dvfdk(i, j)*u(i, j)+Aky(i, j)**2*Para%Cp*d2Tfdk2(i, j)/Para%Pr+ &
                         Aky(i, j)**2*d2udk2(i, j)*u_fluc(i, j)+Aky(i, j)**2*d2ufdk2(i, j)*u(i, j)+ &
                         Aky(i, j)**2*d2vdk2(i, j)*fourthird*v_fluc(i, j)+Aky(i, j)**2*d2vfdk2(i, j)*fourthird*v(i, j)+ &
                         2.0d0*Aky(i, j)**2*dudk(i, j)*dufdk(i, j)+Aky(i, j)**2*dvdk(i, j)*dvfdk(i, j)*eighthird+ &
                         Aky(i, j)*Akyk(i, j)*Para%Cp*dTfdk(i, j)/Para%Pr+Aky(i, j)*Akyk(i, j)*dudk(i, j)*u_fluc(i, j)+ &
                         Aky(i, j)*Akyk(i, j)*dufdk(i, j)*u(i, j)+Aky(i, j)*Akyk(i, j)*dvdk(i, j)*fourthird*v_fluc(i, j)+ &
                         Aky(i, j)*Akyk(i, j)*dvfdk(i, j)*fourthird*v(i, j)

            du(i, j, 2) = du(i, j, 2)+ &
                          dmufdi(i, j)*dmufdi_term(1)+ &
                          dmufdk(i, j)*dmufdk_term(1)+ &
                          mu_fluc(i, j)*muf_term(1)+ &
                          dmudi(i, j)*dmudi_term(1)+ &
                          dmudk(i, j)*dmudk_term(1)+ &
                          mu(i, j)*mu_term(1)

            du(i, j, 3) = du(i, j, 3)+ &
                          dmufdi(i, j)*dmufdi_term(2)+ &
                          dmufdk(i, j)*dmufdk_term(2)+ &
                          mu_fluc(i, j)*muf_term(2)+ &
                          dmudi(i, j)*dmudi_term(2)+ &
                          dmudk(i, j)*dmudk_term(2)+ &
                          mu(i, j)*mu_term(2)

            du(i, j, 4) = du(i, j, 4)+ &
                          dmufdi(i, j)*dmufdi_term(3)+ &
                          dmufdk(i, j)*dmufdk_term(3)+ &
                          mu_fluc(i, j)*muf_term(3)+ &
                          dmudi(i, j)*dmudi_term(3)+ &
                          dmudk(i, j)*dmudk_term(3)+ &
                          mu(i, j)*mu_term(3)

            du(i, j, 5) = du(i, j, 5)+ &
                          dmufdi(i, j)*dmufdi_term(4)+ &
                          dmufdk(i, j)*dmufdk_term(4)+ &
                          mu_fluc(i, j)*muf_term(4)+ &
                          dmudi(i, j)*dmudi_term(4)+ &
                          dmudk(i, j)*dmudk_term(4)+ &
                          mu(i, j)*mu_term(4)

            if (abs(Ayy(i, j)) > 1.0d-10) then
                inv_y = 1.0d0/Ayy(i, j)
                dmufdi_term(1) = -Aix(i, j)*inv_y*twothird*v(i, j)

                dmufdk_term(1) = -Akx(i, j)*inv_y*twothird*v(i, j)

                muf_term(1) = Aix(i, j)*dvdi(i, j)*inv_y*onethird+ &
                              Aix(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)+Aiy(i, j)*dudi(i, j)*inv_y+ &
                              Akx(i, j)*dvdk(i, j)*inv_y*onethird+Akx(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)+ &
                              Aky(i, j)*dudk(i, j)*inv_y

                dmudi_term(1) = -Aix(i, j)*inv_y*twothird*v_fluc(i, j)- &
                                Aix(i, j)*inv_y*twothird*w_fluc(i, j)*Para%wavenumber

                dmudk_term(1) = -Akx(i, j)*inv_y*twothird*v_fluc(i, j)- &
                                Akx(i, j)*inv_y*twothird*w_fluc(i, j)*Para%wavenumber

                mu_term(1) = Aix(i, j)*dvfdi(i, j)*inv_y*onethird+ &
                             Aix(i, j)*dwfdi(i, j)*inv_y*onethird*Para%wavenumber+ &
                             Aix(i, j)*dydi(i, j)*inv_y**2*twothird*v_fluc(i, j)+ &
                             Aix(i, j)*dydi(i, j)*inv_y**2*twothird*w_fluc(i, j)*Para%wavenumber+Aiy(i, j)*dufdi(i, j)*inv_y+ &
                             Akx(i, j)*dvfdk(i, j)*inv_y*onethird+Akx(i, j)*dwfdk(i, j)*inv_y*onethird*Para%wavenumber+ &
                             Akx(i, j)*dydk(i, j)*inv_y**2*twothird*v_fluc(i, j)+ &
                             Akx(i, j)*dydk(i, j)*inv_y**2*twothird*w_fluc(i, j)*Para%wavenumber+Aky(i, j)*dufdk(i, j)*inv_y- &
                             inv_y**2*u_fluc(i, j)*Para%wavenumber**2

                dmufdi_term(2) = -Aiy(i, j)*inv_y*twothird*v(i, j)

                dmufdk_term(2) = -Aky(i, j)*inv_y*twothird*v(i, j)

                muf_term(2) = Aiy(i, j)*dvdi(i, j)*fourthird*inv_y+ &
                              Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)+Aky(i, j)*dvdk(i, j)*fourthird*inv_y+ &
                              Aky(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)-2.0d0*inv_y**2*v(i, j)

                dmudi_term(2) = -Aiy(i, j)*inv_y*twothird*v_fluc(i, j)- &
                                Aiy(i, j)*inv_y*twothird*w_fluc(i, j)*Para%wavenumber

                dmudk_term(2) = -Aky(i, j)*inv_y*twothird*v_fluc(i, j)- &
                                Aky(i, j)*inv_y*twothird*w_fluc(i, j)*Para%wavenumber

                mu_term(2) = Aiy(i, j)*dvfdi(i, j)*fourthird*inv_y+ &
                             Aiy(i, j)*dwfdi(i, j)*inv_y*onethird*Para%wavenumber+ &
                             Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*v_fluc(i, j)+ &
                             Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*w_fluc(i, j)*Para%wavenumber+ &
                             Aky(i, j)*dvfdk(i, j)*fourthird*inv_y+Aky(i, j)*dwfdk(i, j)*inv_y*onethird*Para%wavenumber+ &
                             Aky(i, j)*dydk(i, j)*inv_y**2*twothird*v_fluc(i, j)+ &
                             Aky(i, j)*dydk(i, j)*inv_y**2*twothird*w_fluc(i, j)*Para%wavenumber- &
                             inv_y**2*v_fluc(i, j)*Para%wavenumber**2-2.0d0*inv_y**2*v_fluc(i, j)- &
                             3.0d0*inv_y**2*w_fluc(i, j)*Para%wavenumber

                dmufdi_term(3) = 0.0d0

                dmufdk_term(3) = 0.0d0

                muf_term(3) = Aix(i, j)*dudi(i, j)*inv_y*twothird*Para%wavenumber+ &
                              Aiy(i, j)*dvdi(i, j)*inv_y*twothird*Para%wavenumber+ &
                              Akx(i, j)*dudk(i, j)*inv_y*twothird*Para%wavenumber+ &
                              Aky(i, j)*dvdk(i, j)*inv_y*twothird*Para%wavenumber-fourthird*inv_y**2*v(i, j)*Para%wavenumber

                dmudi_term(3) = -Aix(i, j)*inv_y*u_fluc(i, j)*Para%wavenumber- &
                                Aiy(i, j)*inv_y*v_fluc(i, j)*Para%wavenumber-Aiy(i, j)*inv_y*w_fluc(i, j)

                dmudk_term(3) = -Akx(i, j)*inv_y*u_fluc(i, j)*Para%wavenumber- &
                                Aky(i, j)*inv_y*v_fluc(i, j)*Para%wavenumber-Aky(i, j)*inv_y*w_fluc(i, j)

                mu_term(3) = -Aix(i, j)*dufdi(i, j)*inv_y*onethird*Para%wavenumber+ &
                             Aix(i, j)*dydi(i, j)*inv_y**2*u_fluc(i, j)*Para%wavenumber- &
                             Aiy(i, j)*dvfdi(i, j)*inv_y*onethird*Para%wavenumber+Aiy(i, j)*dwfdi(i, j)*inv_y+ &
                             Aiy(i, j)*dydi(i, j)*inv_y**2*v_fluc(i, j)*Para%wavenumber+ &
                             Aiy(i, j)*dydi(i, j)*inv_y**2*w_fluc(i, j)-Akx(i, j)*dufdk(i, j)*inv_y*onethird*Para%wavenumber+ &
                             Akx(i, j)*dydk(i, j)*inv_y**2*u_fluc(i, j)*Para%wavenumber- &
                             Aky(i, j)*dvfdk(i, j)*inv_y*onethird*Para%wavenumber+Aky(i, j)*dwfdk(i, j)*inv_y+ &
                             Aky(i, j)*dydk(i, j)*inv_y**2*v_fluc(i, j)*Para%wavenumber+ &
                             Aky(i, j)*dydk(i, j)*inv_y**2*w_fluc(i, j)-fourthird*inv_y**2*w_fluc(i, j)*Para%wavenumber**2- &
                             inv_y**2*tenthird*v_fluc(i, j)*Para%wavenumber-2.0d0*inv_y**2*w_fluc(i, j)

                dmufdi_term(4) = -Aix(i, j)*inv_y*twothird*u(i, j)*v(i, j)- &
                                 Aiy(i, j)*inv_y*twothird*v(i, j)**2

                dmufdk_term(4) = -Akx(i, j)*inv_y*twothird*u(i, j)*v(i, j)- &
                                 Aky(i, j)*inv_y*twothird*v(i, j)**2

                muf_term(4) = -Aix(i, j)*dudi(i, j)*fourthird*inv_y*v(i, j)+ &
                           Aix(i, j)*dvdi(i, j)*inv_y*onethird*u(i, j)+Aix(i, j)*dydi(i, j)*inv_y**2*twothird*u(i, j)*v(i, j)+ &
                              Aiy(i, j)*Para%Cp*dTdi(i, j)*inv_y/Para%Pr+Aiy(i, j)*dudi(i, j)*inv_y*u(i, j)+ &
                              Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)**2-Akx(i, j)*dudk(i, j)*fourthird*inv_y*v(i, j)+ &
                           Akx(i, j)*dvdk(i, j)*inv_y*onethird*u(i, j)+Akx(i, j)*dydk(i, j)*inv_y**2*twothird*u(i, j)*v(i, j)+ &
                              Aky(i, j)*Para%Cp*dTdk(i, j)*inv_y/Para%Pr+Aky(i, j)*dudk(i, j)*inv_y*u(i, j)+ &
                              Aky(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)**2-inv_y**2*twothird*v(i, j)**2

                dmudi_term(4) = -Aix(i, j)*inv_y*twothird*u(i, j)*v_fluc(i, j)- &
                                Aix(i, j)*inv_y*twothird*u(i, j)*w_fluc(i, j)*Para%wavenumber- &
                                Aix(i, j)*inv_y*twothird*u_fluc(i, j)*v(i, j)-Aiy(i, j)*fourthird*inv_y*v(i, j)*v_fluc(i, j)- &
                                Aiy(i, j)*inv_y*twothird*v(i, j)*w_fluc(i, j)*Para%wavenumber

                dmudk_term(4) = -Akx(i, j)*inv_y*twothird*u(i, j)*v_fluc(i, j)- &
                                Akx(i, j)*inv_y*twothird*u(i, j)*w_fluc(i, j)*Para%wavenumber- &
                                Akx(i, j)*inv_y*twothird*u_fluc(i, j)*v(i, j)-Aky(i, j)*fourthird*inv_y*v(i, j)*v_fluc(i, j)- &
                                Aky(i, j)*inv_y*twothird*v(i, j)*w_fluc(i, j)*Para%wavenumber

                mu_term(4) = -Aix(i, j)*dudi(i, j)*fourthird*inv_y*v_fluc(i, j)- &
                             Aix(i, j)*dudi(i, j)*fourthird*inv_y*w_fluc(i, j)*Para%wavenumber- &
                             Aix(i, j)*dufdi(i, j)*fourthird*inv_y*v(i, j)+Aix(i, j)*dvdi(i, j)*inv_y*onethird*u_fluc(i, j)+ &
                             Aix(i, j)*dvfdi(i, j)*inv_y*onethird*u(i, j)+ &
                             Aix(i, j)*dwfdi(i, j)*inv_y*onethird*u(i, j)*Para%wavenumber+ &
                             Aix(i, j)*dydi(i, j)*inv_y**2*twothird*u(i, j)*v_fluc(i, j)+ &
                             Aix(i, j)*dydi(i, j)*inv_y**2*twothird*u(i, j)*w_fluc(i, j)*Para%wavenumber+ &
                             Aix(i, j)*dydi(i, j)*inv_y**2*twothird*u_fluc(i, j)*v(i, j)+ &
                             Aiy(i, j)*Para%Cp*dTfdi(i, j)*inv_y/Para%Pr+Aiy(i, j)*dudi(i, j)*inv_y*u_fluc(i, j)+ &
                             Aiy(i, j)*dufdi(i, j)*inv_y*u(i, j)- &
                             Aiy(i, j)*dvdi(i, j)*fourthird*inv_y*w_fluc(i, j)*Para%wavenumber+ &
                             Aiy(i, j)*dwfdi(i, j)*inv_y*onethird*v(i, j)*Para%wavenumber+ &
                             Aiy(i, j)*dydi(i, j)*fourthird*inv_y**2*v(i, j)*v_fluc(i, j)+ &
                             Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)*w_fluc(i, j)*Para%wavenumber- &
                             Akx(i, j)*dudk(i, j)*fourthird*inv_y*v_fluc(i, j)- &
                             Akx(i, j)*dudk(i, j)*fourthird*inv_y*w_fluc(i, j)*Para%wavenumber- &
                             Akx(i, j)*dufdk(i, j)*fourthird*inv_y*v(i, j)+Akx(i, j)*dvdk(i, j)*inv_y*onethird*u_fluc(i, j)+ &
                             Akx(i, j)*dvfdk(i, j)*inv_y*onethird*u(i, j)+ &
                             Akx(i, j)*dwfdk(i, j)*inv_y*onethird*u(i, j)*Para%wavenumber+ &
                             Akx(i, j)*dydk(i, j)*inv_y**2*twothird*u(i, j)*v_fluc(i, j)+ &
                             Akx(i, j)*dydk(i, j)*inv_y**2*twothird*u(i, j)*w_fluc(i, j)*Para%wavenumber+ &
                             Akx(i, j)*dydk(i, j)*inv_y**2*twothird*u_fluc(i, j)*v(i, j)+ &
                             Aky(i, j)*Para%Cp*dTfdk(i, j)*inv_y/Para%Pr+Aky(i, j)*dudk(i, j)*inv_y*u_fluc(i, j)+ &
                             Aky(i, j)*dufdk(i, j)*inv_y*u(i, j)- &
                             Aky(i, j)*dvdk(i, j)*fourthird*inv_y*w_fluc(i, j)*Para%wavenumber+ &
                             Aky(i, j)*dwfdk(i, j)*inv_y*onethird*v(i, j)*Para%wavenumber+ &
                             Aky(i, j)*dydk(i, j)*fourthird*inv_y**2*v(i, j)*v_fluc(i, j)+ &
                             Aky(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)*w_fluc(i, j)*Para%wavenumber- &
                             Para%Cp*T_fluc(i, j)*inv_y**2*Para%wavenumber**2/Para%Pr- &
                             fourthird*inv_y**2*v(i, j)*v_fluc(i, j)-inv_y**2*onethird*v(i, j)*w_fluc(i, j)*Para%wavenumber- &
                             inv_y**2*u(i, j)*u_fluc(i, j)*Para%wavenumber**2- &
                             inv_y**2*v(i, j)*v_fluc(i, j)*Para%wavenumber**2

                Source_Vis(i, j, 2) = dmufdi(i, j)*dmufdi_term(1)+ &
                                      dmufdk(i, j)*dmufdk_term(1)+ &
                                      mu_fluc(i, j)*muf_term(1)+ &
                                      dmudi(i, j)*dmudi_term(1)+ &
                                      dmudk(i, j)*dmudk_term(1)+ &
                                      mu(i, j)*mu_term(1)

                Source_Vis(i, j, 3) = dmufdi(i, j)*dmufdi_term(2)+ &
                                      dmufdk(i, j)*dmufdk_term(2)+ &
                                      mu_fluc(i, j)*muf_term(2)+ &
                                      dmudi(i, j)*dmudi_term(2)+ &
                                      dmudk(i, j)*dmudk_term(2)+ &
                                      mu(i, j)*mu_term(2)

                Source_Vis(i, j, 4) = dmufdi(i, j)*dmufdi_term(3)+ &
                                      dmufdk(i, j)*dmufdk_term(3)+ &
                                      mu_fluc(i, j)*muf_term(3)+ &
                                      dmudi(i, j)*dmudi_term(3)+ &
                                      dmudk(i, j)*dmudk_term(3)+ &
                                      mu(i, j)*mu_term(3)

                Source_Vis(i, j, 5) = dmufdi(i, j)*dmufdi_term(4)+ &
                                      dmufdk(i, j)*dmufdk_term(4)+ &
                                      mu_fluc(i, j)*muf_term(4)+ &
                                      dmudi(i, j)*dmudi_term(4)+ &
                                      dmudk(i, j)*dmudk_term(4)+ &
                                      mu(i, j)*mu_term(4)

            else
                Source_Vis(i, j, 2) = 0.0d0
                Source_Vis(i, j, 3) = 0.0d0
                Source_Vis(i, j, 4) = 0.0d0
                Source_Vis(i, j, 5) = 0.0d0
            end if

        end do
        end do

        ! axis treatment
        if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
            if (mod(m_int, 2) == 0) then
                do j = 1, ny
                    Source_Vis(1, j, 2) = (56.0d0*Source_Vis(2, j, 2)-28.0d0*Source_Vis(3, j, 2)+ &
                                           8.0d0*Source_Vis(4, j, 2)-Source_Vis(5, j, 2))/35.0d0
                    Source_Vis(1, j, 3) = 0.0d0
                    Source_Vis(1, j, 4) = 0.0d0
                    Source_Vis(1, j, 5) = (56.0d0*Source_Vis(2, j, 5)-28.0d0*Source_Vis(3, j, 5)+ &
                                           8.0d0*Source_Vis(4, j, 5)-Source_Vis(5, j, 5))/35.0d0
                end do
            else
                do j = 1, ny
                    Source_Vis(1, j, 2) = 0.0d0
                    Source_Vis(1, j, 3) = (56.0d0*Source_Vis(2, j, 3)-28.0d0*Source_Vis(3, j, 3)+ &
                                           8.0d0*Source_Vis(4, j, 3)-Source_Vis(5, j, 3))/35.0d0
                    Source_Vis(1, j, 4) = (56.0d0*Source_Vis(2, j, 4)-28.0d0*Source_Vis(3, j, 4)+ &
                                           8.0d0*Source_Vis(4, j, 4)-Source_Vis(5, j, 4))/35.0d0
                    Source_Vis(1, j, 5) = 0.0d0
                end do
            end if
        end if

        if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
            if (mod(m_int, 2) == 0) then
                do j = 1, ny
                    Source_Vis(nx, j, 2) = (56.0d0*Source_Vis(nx-1, j, 2)-28.0d0*Source_Vis(nx-2, j, 2)+ &
                                            8.0d0*Source_Vis(nx-3, j, 2)-Source_Vis(nx-4, j, 2))/35.0d0
                    Source_Vis(nx, j, 3) = 0.0d0
                    Source_Vis(nx, j, 4) = 0.0d0
                    Source_Vis(nx, j, 5) = (56.0d0*Source_Vis(nx-1, j, 5)-28.0d0*Source_Vis(nx-2, j, 5)+ &
                                            8.0d0*Source_Vis(nx-3, j, 5)-Source_Vis(nx-4, j, 5))/35.0d0
                end do
            else
                do j = 1, ny
                    Source_Vis(nx, j, 2) = 0.0d0
                    Source_Vis(nx, j, 3) = (56.0d0*Source_Vis(nx-1, j, 3)-28.0d0*Source_Vis(nx-2, j, 3)+ &
                                            8.0d0*Source_Vis(nx-3, j, 3)-Source_Vis(nx-4, j, 3))/35.0d0
                    Source_Vis(nx, j, 4) = (56.0d0*Source_Vis(nx-1, j, 4)-28.0d0*Source_Vis(nx-2, j, 4)+ &
                                            8.0d0*Source_Vis(nx-3, j, 4)-Source_Vis(nx-4, j, 4))/35.0d0
                    Source_Vis(nx, j, 5) = 0.0d0
                end do
            end if
        end if

        if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
            if (mod(m_int, 2) == 0) then
                do i = 1, nx
                    Source_Vis(i, 1, 2) = (56.0d0*Source_Vis(i, 2, 2)-28.0d0*Source_Vis(i, 3, 2)+ &
                                           8.0d0*Source_Vis(i, 4, 2)-Source_Vis(i, 5, 2))/35.0d0
                    Source_Vis(i, 1, 3) = 0.0d0
                    Source_Vis(i, 1, 4) = 0.0d0
                    Source_Vis(i, 1, 5) = (56.0d0*Source_Vis(i, 2, 5)-28.0d0*Source_Vis(i, 3, 5)+ &
                                           8.0d0*Source_Vis(i, 4, 5)-Source_Vis(i, 5, 5))/35.0d0
                end do
            else
                do i = 1, nx
                    Source_Vis(i, 1, 2) = 0.0d0
                    Source_Vis(i, 1, 3) = (56.0d0*Source_Vis(i, 2, 3)-28.0d0*Source_Vis(i, 3, 3)+ &
                                           8.0d0*Source_Vis(i, 4, 3)-Source_Vis(i, 5, 3))/35.0d0
                    Source_Vis(i, 1, 4) = (56.0d0*Source_Vis(i, 2, 4)-28.0d0*Source_Vis(i, 3, 4)+ &
                                           8.0d0*Source_Vis(i, 4, 4)-Source_Vis(i, 5, 4))/35.0d0
                    Source_Vis(i, 1, 5) = 0.0d0
                end do
            end if
        end if

        if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
            if (mod(m_int, 2) == 0) then
                do i = 1, nx
                    Source_Vis(i, ny, 2) = (56.0d0*Source_Vis(i, ny-1, 2)-28.0d0*Source_Vis(i, ny-2, 2)+ &
                                            8.0d0*Source_Vis(i, ny-3, 2)-Source_Vis(i, ny-4, 2))/35.0d0
                    Source_Vis(i, ny, 3) = 0.0d0
                    Source_Vis(i, ny, 4) = 0.0d0
                    Source_Vis(i, ny, 5) = (56.0d0*Source_Vis(i, ny-1, 5)-28.0d0*Source_Vis(i, ny-2, 5)+ &
                                            8.0d0*Source_Vis(i, ny-3, 5)-Source_Vis(i, ny-4, 5))/35.0d0
                end do
            else
                do i = 1, nx
                    Source_Vis(i, ny, 2) = 0.0d0
                    Source_Vis(i, ny, 3) = (56.0d0*Source_Vis(i, ny-1, 3)-28.0d0*Source_Vis(i, ny-2, 3)+ &
                                            8.0d0*Source_Vis(i, ny-3, 3)-Source_Vis(i, ny-4, 3))/35.0d0
                    Source_Vis(i, ny, 4) = (56.0d0*Source_Vis(i, ny-1, 4)-28.0d0*Source_Vis(i, ny-2, 4)+ &
                                            8.0d0*Source_Vis(i, ny-3, 4)-Source_Vis(i, ny-4, 4))/35.0d0
                    Source_Vis(i, ny, 5) = 0.0d0
                end do
            end if
        end if

        do j = 1, ny
        do i = 1, nx
            du(i, j, 2) = du(i, j, 2)+Source_Vis(i, j, 2)
            du(i, j, 3) = du(i, j, 3)+Source_Vis(i, j, 3)
            du(i, j, 4) = du(i, j, 4)+Source_Vis(i, j, 4)
            du(i, j, 5) = du(i, j, 5)+Source_Vis(i, j, 5)
        end do
        end do

    end subroutine flux_viscous2d_cyl_d

end module Flux_Viscous_cyl_d
