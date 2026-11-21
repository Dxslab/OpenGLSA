! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Copyright by Li Xinliang
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module Flux_Viscous
    implicit none

contains

    subroutine flux_viscous2d
        use mesh_data
        use flow_data
        use Centrals
        implicit none

        integer :: i, j
        real(dp) :: inv_y
        real(dp) :: dmudk_term(3), dmudi_term(3), mu_term(3)

        real(dp), parameter :: onethird = 1.0d0/3.0d0
        real(dp), parameter :: twothird = 2.0d0/3.0d0
        real(dp), parameter :: fourthird = 4.0d0/3.0d0
        real(dp), parameter :: eighthird = 8.0d0/3.0d0

        call compute_mu
        !
        if (Para%IBC /= BC_Lid_Driven_Cavity) then
            call OCFD2d_dx0(u, dudk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
            call OCFD2d_dx0(v, dvdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
            call OCFD2d_dx0(T, dTdk(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcI)
            !
            call OCFD2d_dy0(u, dudi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
            call OCFD2d_dy0(v, dvdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
            call OCFD2d_dy0(T, dTdi(1:nx, 1:ny), nx, ny, Scheme%Scheme_Vis, BcJ)
            !
            ! k-derivatives are differentiated only in y; x halos are unused.
            call mpi_exchange_y2d(dudk)
            call mpi_exchange_y2d(dvdk)
            call mpi_exchange_y2d(dTdk)
            !
            call OCFD2d_dy0(dudk, d2udkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
            call OCFD2d_dy0(dvdk, d2vdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
            call OCFD2d_dy0(dTdk, d2Tdkdi, nx, ny, Scheme%Scheme_Vis, BcJ)
            !
            call OCFD2d_dx2(u, d2udk2, nx, ny, Scheme%Scheme_Vis, BcI)
            call OCFD2d_dx2(v, d2vdk2, nx, ny, Scheme%Scheme_Vis, BcI)
            call OCFD2d_dx2(T, d2Tdk2, nx, ny, Scheme%Scheme_Vis, BcI)
            call OCFD2d_dy2(u, d2udi2, nx, ny, Scheme%Scheme_Vis, BcJ)
            call OCFD2d_dy2(v, d2vdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
            call OCFD2d_dy2(T, d2Tdi2, nx, ny, Scheme%Scheme_Vis, BcJ)
        end if
        !
        dmudk = diff_mu(dTdk, T)
        dmudi = diff_mu(dTdi, T)
        do j = 1, ny
        do i = 1, nx
            ! x momentum equation viscous source terms for variable viscosity case
            dmudi_term(1) = Aix(i, j)**2*dudi(i, j)*fourthird+ &
                            Aix(i, j)*Aiy(i, j)*dvdi(i, j)*onethird+ &
                            Aix(i, j)*Akx(i, j)*dudk(i, j)*fourthird- &
                            Aix(i, j)*Aky(i, j)*dvdk(i, j)*twothird+Aiy(i, j)**2*dudi(i, j)+ &
                            Aiy(i, j)*Akx(i, j)*dvdk(i, j)+Aiy(i, j)*Aky(i, j)*dudk(i, j)

            dmudk_term(1) = Aix(i, j)*Akx(i, j)*dudi(i, j)*fourthird+ &
                            Aix(i, j)*Aky(i, j)*dvdi(i, j)-Aiy(i, j)*Akx(i, j)*dvdi(i, j)*twothird+ &
                            Aiy(i, j)*Aky(i, j)*dudi(i, j)+Akx(i, j)**2*dudk(i, j)*fourthird+ &
                            Akx(i, j)*Aky(i, j)*dvdk(i, j)*onethird+Aky(i, j)**2*dudk(i, j)

            mu_term(1) = Aix(i, j)**2*d2udi2(i, j)*fourthird+ &
                         Aix(i, j)*Aixi(i, j)*dudi(i, j)*fourthird+ &
                         Aix(i, j)*Aiy(i, j)*d2vdi2(i, j)*onethird- &
                         Aix(i, j)*Aiyi(i, j)*dvdi(i, j)*twothird+ &
                         Aix(i, j)*Akx(i, j)*d2udkdi(i, j)*eighthird+ &
                         Aix(i, j)*Akxi(i, j)*dudk(i, j)*fourthird+ &
                         Aix(i, j)*Aky(i, j)*d2vdkdi(i, j)*onethird- &
                         Aix(i, j)*Akyi(i, j)*dvdk(i, j)*twothird+Aixi(i, j)*Aiy(i, j)*dvdi(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*dudi(i, j)*fourthird+Aixk(i, j)*Aky(i, j)*dvdi(i, j)+ &
                         Aiy(i, j)**2*d2udi2(i, j)+Aiy(i, j)*Aiyi(i, j)*dudi(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2vdkdi(i, j)*onethird+Aiy(i, j)*Akxi(i, j)*dvdk(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*d2udkdi(i, j)+Aiy(i, j)*Akyi(i, j)*dudk(i, j)- &
                         Aiyk(i, j)*Akx(i, j)*dvdi(i, j)*twothird+Aiyk(i, j)*Aky(i, j)*dudi(i, j)+ &
                         Akx(i, j)**2*d2udk2(i, j)*fourthird+Akx(i, j)*Akxk(i, j)*dudk(i, j)*fourthird+ &
                         Akx(i, j)*Aky(i, j)*d2vdk2(i, j)*onethird- &
                         Akx(i, j)*Akyk(i, j)*dvdk(i, j)*twothird+Akxk(i, j)*Aky(i, j)*dvdk(i, j)+ &
                         Aky(i, j)**2*d2udk2(i, j)+Aky(i, j)*Akyk(i, j)*dudk(i, j)

            ! y momentum equation viscous source terms for variable viscosity case
            dmudi_term(2) = Aix(i, j)**2*dvdi(i, j)+Aix(i, j)*Aiy(i, j)*dudi(i, j)*onethird+ &
                            Aix(i, j)*Akx(i, j)*dvdk(i, j)+Aix(i, j)*Aky(i, j)*dudk(i, j)+ &
                            Aiy(i, j)**2*dvdi(i, j)*fourthird-Aiy(i, j)*Akx(i, j)*dudk(i, j)*twothird+ &
                            Aiy(i, j)*Aky(i, j)*dvdk(i, j)*fourthird

            dmudk_term(2) = Aix(i, j)*Akx(i, j)*dvdi(i, j)- &
                            Aix(i, j)*Aky(i, j)*dudi(i, j)*twothird+Aiy(i, j)*Akx(i, j)*dudi(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvdi(i, j)*fourthird+Akx(i, j)**2*dvdk(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dudk(i, j)*onethird+Aky(i, j)**2*dvdk(i, j)*fourthird

            mu_term(2) = Aix(i, j)**2*d2vdi2(i, j)+Aix(i, j)*Aixi(i, j)*dvdi(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2udi2(i, j)*onethird+Aix(i, j)*Aiyi(i, j)*dudi(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*d2vdkdi(i, j)+Aix(i, j)*Akxi(i, j)*dvdk(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2udkdi(i, j)*onethird+Aix(i, j)*Akyi(i, j)*dudk(i, j)- &
                         Aixi(i, j)*Aiy(i, j)*dudi(i, j)*twothird+Aixk(i, j)*Akx(i, j)*dvdi(i, j)- &
                         Aixk(i, j)*Aky(i, j)*dudi(i, j)*twothird+Aiy(i, j)**2*d2vdi2(i, j)*fourthird+ &
                         Aiy(i, j)*Aiyi(i, j)*dvdi(i, j)*fourthird+ &
                         Aiy(i, j)*Akx(i, j)*d2udkdi(i, j)*onethird- &
                         Aiy(i, j)*Akxi(i, j)*dudk(i, j)*twothird+ &
                         Aiy(i, j)*Aky(i, j)*d2vdkdi(i, j)*eighthird+ &
                         Aiy(i, j)*Akyi(i, j)*dvdk(i, j)*fourthird+Aiyk(i, j)*Akx(i, j)*dudi(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*dvdi(i, j)*fourthird+Akx(i, j)**2*d2vdk2(i, j)+ &
                         Akx(i, j)*Akxk(i, j)*dvdk(i, j)+Akx(i, j)*Aky(i, j)*d2udk2(i, j)*onethird+ &
                         Akx(i, j)*Akyk(i, j)*dudk(i, j)-Akxk(i, j)*Aky(i, j)*dudk(i, j)*twothird+ &
                         Aky(i, j)**2*d2vdk2(i, j)*fourthird+Aky(i, j)*Akyk(i, j)*dvdk(i, j)*fourthird

            ! energy equation viscous source terms for variable viscosity case
            dmudi_term(3) = Aix(i, j)**2*Para%Cp*dTdi(i, j)/Para%Pr+ &
                            Aix(i, j)**2*dudi(i, j)*fourthird*u(i, j)+Aix(i, j)**2*dvdi(i, j)*v(i, j)+ &
                            Aix(i, j)*Aiy(i, j)*dudi(i, j)*onethird*v(i, j)+ &
                            Aix(i, j)*Aiy(i, j)*dvdi(i, j)*onethird*u(i, j)+ &
                            Aix(i, j)*Akx(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                            Aix(i, j)*Akx(i, j)*dudk(i, j)*fourthird*u(i, j)+ &
                            Aix(i, j)*Akx(i, j)*dvdk(i, j)*v(i, j)+ &
                            Aix(i, j)*Aky(i, j)*dudk(i, j)*v(i, j)- &
                            Aix(i, j)*Aky(i, j)*dvdk(i, j)*twothird*u(i, j)+ &
                            Aiy(i, j)**2*Para%Cp*dTdi(i, j)/Para%Pr+Aiy(i, j)**2*dudi(i, j)*u(i, j)+ &
                            Aiy(i, j)**2*dvdi(i, j)*fourthird*v(i, j)- &
                            Aiy(i, j)*Akx(i, j)*dudk(i, j)*twothird*v(i, j)+ &
                            Aiy(i, j)*Akx(i, j)*dvdk(i, j)*u(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                            Aiy(i, j)*Aky(i, j)*dudk(i, j)*u(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvdk(i, j)*fourthird*v(i, j)

            dmudk_term(3) = Aix(i, j)*Akx(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                            Aix(i, j)*Akx(i, j)*dudi(i, j)*fourthird*u(i, j)+ &
                            Aix(i, j)*Akx(i, j)*dvdi(i, j)*v(i, j)- &
                            Aix(i, j)*Aky(i, j)*dudi(i, j)*twothird*v(i, j)+ &
                            Aix(i, j)*Aky(i, j)*dvdi(i, j)*u(i, j)+ &
                            Aiy(i, j)*Akx(i, j)*dudi(i, j)*v(i, j)- &
                            Aiy(i, j)*Akx(i, j)*dvdi(i, j)*twothird*u(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                            Aiy(i, j)*Aky(i, j)*dudi(i, j)*u(i, j)+ &
                            Aiy(i, j)*Aky(i, j)*dvdi(i, j)*fourthird*v(i, j)+ &
                            Akx(i, j)**2*Para%Cp*dTdk(i, j)/Para%Pr+ &
                            Akx(i, j)**2*dudk(i, j)*fourthird*u(i, j)+ &
                            Akx(i, j)**2*dvdk(i, j)*v(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dudk(i, j)*onethird*v(i, j)+ &
                            Akx(i, j)*Aky(i, j)*dvdk(i, j)*onethird*u(i, j)+ &
                            Aky(i, j)**2*Para%Cp*dTdk(i, j)/Para%Pr+Aky(i, j)**2*dudk(i, j)*u(i, j)+ &
                            Aky(i, j)**2*dvdk(i, j)*fourthird*v(i, j)

            mu_term(3) = Aix(i, j)**2*Para%Cp*d2Tdi2(i, j)/Para%Pr+ &
                         Aix(i, j)**2*d2udi2(i, j)*fourthird*u(i, j)+Aix(i, j)**2*d2vdi2(i, j)*v(i, j)+ &
                         Aix(i, j)**2*dudi(i, j)**2*fourthird+Aix(i, j)**2*dvdi(i, j)**2+ &
                         Aix(i, j)*Aixi(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                         Aix(i, j)*Aixi(i, j)*dudi(i, j)*fourthird*u(i, j)+ &
                         Aix(i, j)*Aixi(i, j)*dvdi(i, j)*v(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2udi2(i, j)*onethird*v(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*d2vdi2(i, j)*onethird*u(i, j)+ &
                         Aix(i, j)*Aiy(i, j)*dudi(i, j)*dvdi(i, j)*twothird+ &
                         Aix(i, j)*Aiyi(i, j)*dudi(i, j)*v(i, j)- &
                         Aix(i, j)*Aiyi(i, j)*dvdi(i, j)*twothird*u(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*Para%Cp*d2Tdkdi(i, j)/Para%Pr+ &
                         Aix(i, j)*Akx(i, j)*d2udkdi(i, j)*eighthird*u(i, j)+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*d2vdkdi(i, j)*v(i, j)+ &
                         Aix(i, j)*Akx(i, j)*dudi(i, j)*dudk(i, j)*eighthird+ &
                         2.0d0*Aix(i, j)*Akx(i, j)*dvdi(i, j)*dvdk(i, j)+ &
                         Aix(i, j)*Akxi(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                         Aix(i, j)*Akxi(i, j)*dudk(i, j)*fourthird*u(i, j)+ &
                         Aix(i, j)*Akxi(i, j)*dvdk(i, j)*v(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2udkdi(i, j)*onethird*v(i, j)+ &
                         Aix(i, j)*Aky(i, j)*d2vdkdi(i, j)*onethird*u(i, j)- &
                         Aix(i, j)*Aky(i, j)*dudi(i, j)*dvdk(i, j)*fourthird+ &
                         2.0d0*Aix(i, j)*Aky(i, j)*dudk(i, j)*dvdi(i, j)+ &
                         Aix(i, j)*Akyi(i, j)*dudk(i, j)*v(i, j)- &
                         Aix(i, j)*Akyi(i, j)*dvdk(i, j)*twothird*u(i, j)- &
                         Aixi(i, j)*Aiy(i, j)*dudi(i, j)*twothird*v(i, j)+ &
                         Aixi(i, j)*Aiy(i, j)*dvdi(i, j)*u(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                         Aixk(i, j)*Akx(i, j)*dudi(i, j)*fourthird*u(i, j)+ &
                         Aixk(i, j)*Akx(i, j)*dvdi(i, j)*v(i, j)- &
                         Aixk(i, j)*Aky(i, j)*dudi(i, j)*twothird*v(i, j)+ &
                         Aixk(i, j)*Aky(i, j)*dvdi(i, j)*u(i, j)+ &
                         Aiy(i, j)**2*Para%Cp*d2Tdi2(i, j)/Para%Pr+Aiy(i, j)**2*d2udi2(i, j)*u(i, j)+ &
                         Aiy(i, j)**2*d2vdi2(i, j)*fourthird*v(i, j)+Aiy(i, j)**2*dudi(i, j)**2+ &
                         Aiy(i, j)**2*dvdi(i, j)**2*fourthird+ &
                         Aiy(i, j)*Aiyi(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                         Aiy(i, j)*Aiyi(i, j)*dudi(i, j)*u(i, j)+ &
                         Aiy(i, j)*Aiyi(i, j)*dvdi(i, j)*fourthird*v(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2udkdi(i, j)*onethird*v(i, j)+ &
                         Aiy(i, j)*Akx(i, j)*d2vdkdi(i, j)*onethird*u(i, j)+ &
                         2.0d0*Aiy(i, j)*Akx(i, j)*dudi(i, j)*dvdk(i, j)- &
                         Aiy(i, j)*Akx(i, j)*dudk(i, j)*dvdi(i, j)*fourthird- &
                         Aiy(i, j)*Akxi(i, j)*dudk(i, j)*twothird*v(i, j)+ &
                         Aiy(i, j)*Akxi(i, j)*dvdk(i, j)*u(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*Para%Cp*d2Tdkdi(i, j)/Para%Pr+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*d2udkdi(i, j)*u(i, j)+ &
                         Aiy(i, j)*Aky(i, j)*d2vdkdi(i, j)*eighthird*v(i, j)+ &
                         2.0d0*Aiy(i, j)*Aky(i, j)*dudi(i, j)*dudk(i, j)+ &
                         Aiy(i, j)*Aky(i, j)*dvdi(i, j)*dvdk(i, j)*eighthird+ &
                         Aiy(i, j)*Akyi(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                         Aiy(i, j)*Akyi(i, j)*dudk(i, j)*u(i, j)+ &
                         Aiy(i, j)*Akyi(i, j)*dvdk(i, j)*fourthird*v(i, j)+ &
                         Aiyk(i, j)*Akx(i, j)*dudi(i, j)*v(i, j)- &
                         Aiyk(i, j)*Akx(i, j)*dvdi(i, j)*twothird*u(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*Para%Cp*dTdi(i, j)/Para%Pr+ &
                         Aiyk(i, j)*Aky(i, j)*dudi(i, j)*u(i, j)+ &
                         Aiyk(i, j)*Aky(i, j)*dvdi(i, j)*fourthird*v(i, j)+ &
                         Akx(i, j)**2*Para%Cp*d2Tdk2(i, j)/Para%Pr+ &
                         Akx(i, j)**2*d2udk2(i, j)*fourthird*u(i, j)+Akx(i, j)**2*d2vdk2(i, j)*v(i, j)+ &
                         Akx(i, j)**2*dudk(i, j)**2*fourthird+Akx(i, j)**2*dvdk(i, j)**2+ &
                         Akx(i, j)*Akxk(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                         Akx(i, j)*Akxk(i, j)*dudk(i, j)*fourthird*u(i, j)+ &
                         Akx(i, j)*Akxk(i, j)*dvdk(i, j)*v(i, j)+ &
                         Akx(i, j)*Aky(i, j)*d2udk2(i, j)*onethird*v(i, j)+ &
                         Akx(i, j)*Aky(i, j)*d2vdk2(i, j)*onethird*u(i, j)+ &
                         Akx(i, j)*Aky(i, j)*dudk(i, j)*dvdk(i, j)*twothird+ &
                         Akx(i, j)*Akyk(i, j)*dudk(i, j)*v(i, j)- &
                         Akx(i, j)*Akyk(i, j)*dvdk(i, j)*twothird*u(i, j)- &
                         Akxk(i, j)*Aky(i, j)*dudk(i, j)*twothird*v(i, j)+ &
                         Akxk(i, j)*Aky(i, j)*dvdk(i, j)*u(i, j)+ &
                         Aky(i, j)**2*Para%Cp*d2Tdk2(i, j)/Para%Pr+Aky(i, j)**2*d2udk2(i, j)*u(i, j)+ &
                         Aky(i, j)**2*d2vdk2(i, j)*fourthird*v(i, j)+Aky(i, j)**2*dudk(i, j)**2+ &
                         Aky(i, j)**2*dvdk(i, j)**2*fourthird+ &
                         Aky(i, j)*Akyk(i, j)*Para%Cp*dTdk(i, j)/Para%Pr+ &
                         Aky(i, j)*Akyk(i, j)*dudk(i, j)*u(i, j)+ &
                         Aky(i, j)*Akyk(i, j)*dvdk(i, j)*fourthird*v(i, j)
            ! add viscous source terms to du
            du(i, j, 2) = du(i, j, 2)+ &
                          dmudi(i, j)*dmudi_term(1)+dmudk(i, j)*dmudk_term(1)+mu(i, j)*mu_term(1)
            du(i, j, 3) = du(i, j, 3)+ &
                          dmudi(i, j)*dmudi_term(2)+dmudk(i, j)*dmudk_term(2)+mu(i, j)*mu_term(2)
            du(i, j, 4) = du(i, j, 4)+ &
                          dmudi(i, j)*dmudi_term(3)+dmudk(i, j)*dmudk_term(3)+mu(i, j)*mu_term(3)
        end do
        end do


        ! Cylindrical coordinate treatment
        if (Para%IF_Cylindrical_Coordinate == 1) then
            do j = 1, ny
            do i = 1, nx
                if (abs(Ayy(i, j)) > 1.0d-10) then
                    inv_y = 1.0d0/Ayy(i, j)

                    dmudi_term(1) = -Aix(i, j)*inv_y*twothird*v(i, j)
                    dmudk_term(1) = -Akx(i, j)*inv_y*twothird*v(i, j)
                    mu_term(1) = Aix(i, j)*dvdi(i, j)*inv_y*onethird+ &
                                 Aix(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)+ &
                                 Aiy(i, j)*dudi(i, j)*inv_y+ &
                                 Akx(i, j)*dvdk(i, j)*inv_y*onethird+ &
                                 Akx(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)+ &
                                 Aky(i, j)*dudk(i, j)*inv_y

                    dmudi_term(2) = -Aiy(i, j)*inv_y*twothird*v(i, j)
                    dmudk_term(2) = -Aky(i, j)*inv_y*twothird*v(i, j)
                    mu_term(2) = Aiy(i, j)*dvdi(i, j)*fourthird*inv_y+ &
                                 Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)+ &
                                 Aky(i, j)*dvdk(i, j)*fourthird*inv_y+ &
                                 Aky(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)-2.0d0*inv_y**2*v(i, j)

                    dmudi_term(3) = -Aix(i, j)*inv_y*twothird*u(i, j)*v(i, j)- &
                                    Aiy(i, j)*inv_y*twothird*v(i, j)**2
                    dmudk_term(3) = -Akx(i, j)*inv_y*twothird*u(i, j)*v(i, j)- &
                                    Aky(i, j)*inv_y*twothird*v(i, j)**2
                    mu_term(3) = -Aix(i, j)*dudi(i, j)*fourthird*inv_y*v(i, j)+ &
                                 Aix(i, j)*dvdi(i, j)*inv_y*onethird*u(i, j)+ &
                                 Aix(i, j)*dydi(i, j)*inv_y**2*twothird*u(i, j)*v(i, j)+ &
                                 Aiy(i, j)*Para%Cp*dTdi(i, j)*inv_y/Para%Pr+Aiy(i, j)*dudi(i, j)*inv_y*u(i, j)+ &
                                 Aiy(i, j)*dydi(i, j)*inv_y**2*twothird*v(i, j)**2- &
                                 Akx(i, j)*dudk(i, j)*fourthird*inv_y*v(i, j)+ &
                                 Akx(i, j)*dvdk(i, j)*inv_y*onethird*u(i, j)+ &
                                 Akx(i, j)*dydk(i, j)*inv_y**2*twothird*u(i, j)*v(i, j)+ &
                                 Aky(i, j)*Para%Cp*dTdk(i, j)*inv_y/Para%Pr+Aky(i, j)*dudk(i, j)*inv_y*u(i, j)+ &
                                 Aky(i, j)*dydk(i, j)*inv_y**2*twothird*v(i, j)**2-inv_y**2*twothird*v(i, j)**2

                    Source_Vis(i, j, 2) = dmudi(i, j)*dmudi_term(1)+dmudk(i, j)*dmudk_term(1)+mu(i, j)*mu_term(1)
                    Source_Vis(i, j, 3) = dmudi(i, j)*dmudi_term(2)+dmudk(i, j)*dmudk_term(2)+mu(i, j)*mu_term(2)
                    Source_Vis(i, j, 4) = dmudi(i, j)*dmudi_term(3)+dmudk(i, j)*dmudk_term(3)+mu(i, j)*mu_term(3)

                else
                    Source_Vis(i, j, 2) = 0.0d0
                    Source_Vis(i, j, 3) = 0.0d0
                    Source_Vis(i, j, 4) = 0.0d0
                end if
            end do
            end do

            if (Para%Isymmetric_X(1) == 1 .and. npx == 0) then
                do j = 1, ny
                    Source_Vis(1, j, 2) = (56.0d0*Source_Vis(2, j, 2)-28.0d0*Source_Vis(3, j, 2)+ &
                                           8.0d0*Source_Vis(4, j, 2)-Source_Vis(5, j, 2))/35.0d0
                    Source_Vis(1, j, 3) = 0.0d0
                    Source_Vis(1, j, 4) = (56.0d0*Source_Vis(2, j, 4)-28.0d0*Source_Vis(3, j, 4)+ &
                                           8.0d0*Source_Vis(4, j, 4)-Source_Vis(5, j, 4))/35.0d0
                end do
            end if

            if (Para%Isymmetric_X(2) == 1 .and. npx == npx0-1) then
                do j = 1, ny
                    Source_Vis(nx, j, 2) = (56.0d0*Source_Vis(nx-1, j, 2)-28.0d0*Source_Vis(nx-2, j, 2)+ &
                                            8.0d0*Source_Vis(nx-3, j, 2)-Source_Vis(nx-4, j, 2))/35.0d0
                    Source_Vis(nx, j, 3) = 0.0d0
                    Source_Vis(nx, j, 4) = (56.0d0*Source_Vis(nx-1, j, 4)-28.0d0*Source_Vis(nx-2, j, 4)+ &
                                            8.0d0*Source_Vis(nx-3, j, 4)-Source_Vis(nx-4, j, 4))/35.0d0
                end do
            end if

            if (Para%Isymmetric_Y(1) == 1 .and. npy == 0) then
                do i = 1, nx
                    Source_Vis(i, 1, 2) = (56.0d0*Source_Vis(i, 2, 2)-28.0d0*Source_Vis(i, 3, 2)+ &
                                           8.0d0*Source_Vis(i, 4, 2)-Source_Vis(i, 5, 2))/35.0d0
                    Source_Vis(i, 1, 3) = 0.0d0
                    Source_Vis(i, 1, 4) = (56.0d0*Source_Vis(i, 2, 4)-28.0d0*Source_Vis(i, 3, 4)+ &
                                           8.0d0*Source_Vis(i, 4, 4)-Source_Vis(i, 5, 4))/35.0d0
                end do
            end if

            if (Para%Isymmetric_Y(2) == 1 .and. npy == npy0-1) then
                do i = 1, nx
                    Source_Vis(i, ny, 2) = (56.0d0*Source_Vis(i, ny-1, 2)-28.0d0*Source_Vis(i, ny-2, 2)+ &
                                            8.0d0*Source_Vis(i, ny-3, 2)-Source_Vis(i, ny-4, 2))/35.0d0
                    Source_Vis(i, ny, 3) = 0.0d0
                    Source_Vis(i, ny, 4) = (56.0d0*Source_Vis(i, ny-1, 4)-28.0d0*Source_Vis(i, ny-2, 4)+ &
                                            8.0d0*Source_Vis(i, ny-3, 4)-Source_Vis(i, ny-4, 4))/35.0d0
                end do
            end if

            do j = 1, ny
            do i = 1, nx
                du(i, j, 2) = du(i, j, 2)+Source_Vis(i, j, 2)
                du(i, j, 3) = du(i, j, 3)+Source_Vis(i, j, 3)
                du(i, j, 4) = du(i, j, 4)+Source_Vis(i, j, 4)
            end do
            end do

        end if

    end subroutine flux_viscous2d

    subroutine compute_mu
        use mesh_data
        use flow_data
        implicit none

        real(dp):: Tsb, tmpR
        integer:: i, j

        Tsb = 110.4d0/Para%T_inf
        TmpR = 1.0d0/Para%Re

        do j = 1, ny
        do i = 1, nx
            mu(i, j) = TmpR*(1.0d0+Tsb)*sqrt(T(i, j)**3)/(Tsb+T(i, j))
        end do
        end do

    end subroutine compute_mu

    elemental function diff_mu(dTd, T) result(dmud)
        use flow_para, only: Para, dp
        implicit none

        real(dp) :: dmud
        real(dp), intent(IN)  :: dTd, T
        real(dp) :: Tsb

        Tsb = 110.4d0/Para%T_inf
        dmud = (1.0d0+Tsb)*(1.5d0*sqrt(T)*(Tsb+T)-T*sqrt(T))*dTd/(Para%Re*(Tsb+T)**2)

    end function diff_mu

end module flux_viscous
