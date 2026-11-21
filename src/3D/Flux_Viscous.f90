! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Copyright by Li Xinliang
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module Flux_Viscous
    implicit none

contains

    subroutine flux_viscous3d
        use mesh_data
        use flow_data
        use Centrals
        implicit none

        integer :: i, j, k
        real(dp) :: dmudk_term(4), dmudi_term(4), dmuds_term(4), mu_term(4)

        real(dp), parameter :: onethird = 1.0d0/3.0d0
        real(dp), parameter :: twothird = 2.0d0/3.0d0
        real(dp), parameter :: fourthird = 4.0d0/3.0d0
        real(dp), parameter :: eighthird = 8.0d0/3.0d0
        !
        call compute_mu
        !
        ! For the lid-driven cavity the derivatives were already computed in the wall condition
        if (Para%IBC /= BC_Lid_Driven_Cavity) call compute_derivatives
        !
        dmudk = diff_mu(dTdk, T)
        dmudi = diff_mu(dTdi, T)
        dmuds = diff_mu(dTds, T)
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            ! x momentum equation viscous source terms for variable viscosity case
            dmudk_term(1) = Aix(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird+Aix(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)+ &
                            Aix(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dudi(i, j, k)-Aiz(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dudi(i, j, k)+Akx(i, j, k)**2*dudk(i, j, k)*fourthird+ &
                            Akx(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*onethird+Akx(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird+ &
                            Akx(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird-Akx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird- &
                            Akx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+Aky(i, j, k)**2*dudk(i, j, k)+ &
                            Aky(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*duds(i, j, k)+ &
                            Akz(i, j, k)**2*dudk(i, j, k)+Akz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*duds(i, j, k)

            dmudi_term(1) = Aix(i, j, k)**2*dudi(i, j, k)*fourthird+Aix(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)*onethird+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird+Aix(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*fourthird- &
                            Aix(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird-Aix(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird+ &
                            Aix(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird-Aix(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird- &
                            Aix(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+Aiy(i, j, k)**2*dudi(i, j, k)+ &
                            Aiy(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dudk(i, j, k)+ &
                            Aiy(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*duds(i, j, k)+ &
                            Aiz(i, j, k)**2*dudi(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dudk(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*duds(i, j, k)

            dmuds_term(1) = Aix(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird+Aix(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)+ &
                            Aix(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dudi(i, j, k)-Aiz(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dudi(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird+ &
                            Akx(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)+Akx(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*dudk(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird+Akz(i, j, k)*Asz(i, j, k)*dudk(i, j, k)+ &
                            Asx(i, j, k)**2*duds(i, j, k)*fourthird+Asx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*onethird+ &
                            Asx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird+Asy(i, j, k)**2*duds(i, j, k)+ &
                            Asz(i, j, k)**2*duds(i, j, k)

            mu_term(1) = Aix(i, j, k)**2*d2udi2(i, j, k)*fourthird+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dudi(i, j, k)*fourthird+Aix(i, j, k)*Aiy(i, j, k)*d2vdi2(i, j, k)*onethird- &
                         Aix(i, j, k)*Aiyi(i, j, k)*dvdi(i, j, k)*twothird+Aix(i, j, k)*Aiz(i, j, k)*d2wdi2(i, j, k)*onethird- &
                         Aix(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*twothird+Aix(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*eighthird+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*fourthird+Aix(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*onethird- &
                         Aix(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*twothird+Aix(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*onethird- &
                         Aix(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*twothird+Aix(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*eighthird+ &
                         Aix(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*fourthird+Aix(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*onethird- &
                         Aix(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*twothird+Aix(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*onethird- &
                         Aix(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*twothird+Aixi(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)+ &
                         Aixi(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)+Aixk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird+ &
                         Aixk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)+Aixk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird+Aixs(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)+ &
                         Aixs(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)+Aiy(i, j, k)**2*d2udi2(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dudi(i, j, k)+Aiy(i, j, k)*Akx(i, j, k)*d2vdkdi(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Akxi(i, j, k)*dvdk(i, j, k)+2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2udkdi(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dudk(i, j, k)+Aiy(i, j, k)*Asx(i, j, k)*d2vdids(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Asxi(i, j, k)*dvds(i, j, k)+2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2udids(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*duds(i, j, k)-Aiyk(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dudi(i, j, k)-Aiys(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dudi(i, j, k)+Aiz(i, j, k)**2*d2udi2(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dudi(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*d2wdkdi(i, j, k)*onethird+ &
                         Aiz(i, j, k)*Akxi(i, j, k)*dwdk(i, j, k)+2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2udkdi(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dudk(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*d2wdids(i, j, k)*onethird+ &
                         Aiz(i, j, k)*Asxi(i, j, k)*dwds(i, j, k)+2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2udids(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*duds(i, j, k)-Aizk(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dudi(i, j, k)-Aizs(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dudi(i, j, k)+Akx(i, j, k)**2*d2udk2(i, j, k)*fourthird+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dudk(i, j, k)*fourthird+Akx(i, j, k)*Aky(i, j, k)*d2vdk2(i, j, k)*onethird- &
                         Akx(i, j, k)*Akyk(i, j, k)*dvdk(i, j, k)*twothird+Akx(i, j, k)*Akz(i, j, k)*d2wdk2(i, j, k)*onethird- &
                         Akx(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*twothird+Akx(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*eighthird+ &
                         Akx(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*fourthird+Akx(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*onethird- &
                         Akx(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*twothird+Akx(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*onethird- &
                         Akx(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*twothird+Akxk(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)+ &
                         Akxk(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)+Akxs(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird+ &
                         Akxs(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)+Akxs(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)+ &
                         Aky(i, j, k)**2*d2udk2(i, j, k)+Aky(i, j, k)*Akyk(i, j, k)*dudk(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2vdkds(i, j, k)*onethird+Aky(i, j, k)*Asxk(i, j, k)*dvds(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2udkds(i, j, k)+Aky(i, j, k)*Asyk(i, j, k)*duds(i, j, k)- &
                         Akys(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird+Akys(i, j, k)*Asy(i, j, k)*dudk(i, j, k)+ &
                         Akz(i, j, k)**2*d2udk2(i, j, k)+Akz(i, j, k)*Akzk(i, j, k)*dudk(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2wdkds(i, j, k)*onethird+Akz(i, j, k)*Asxk(i, j, k)*dwds(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2udkds(i, j, k)+Akz(i, j, k)*Aszk(i, j, k)*duds(i, j, k)- &
                         Akzs(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird+Akzs(i, j, k)*Asz(i, j, k)*dudk(i, j, k)+ &
                         Asx(i, j, k)**2*d2uds2(i, j, k)*fourthird+Asx(i, j, k)*Asxs(i, j, k)*duds(i, j, k)*fourthird+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2vds2(i, j, k)*onethird-Asx(i, j, k)*Asys(i, j, k)*dvds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2wds2(i, j, k)*onethird-Asx(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*twothird+ &
                         Asxs(i, j, k)*Asy(i, j, k)*dvds(i, j, k)+Asxs(i, j, k)*Asz(i, j, k)*dwds(i, j, k)+ &
                         Asy(i, j, k)**2*d2uds2(i, j, k)+Asy(i, j, k)*Asys(i, j, k)*duds(i, j, k)+Asz(i, j, k)**2*d2uds2(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*duds(i, j, k)

            ! y momentum equation viscous source terms for variable viscosity case
            dmudk_term(2) = Aix(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)-Aix(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Akx(i, j, k)*dudi(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird+ &
                            Aiy(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)-Aiz(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)+Akx(i, j, k)**2*dvdk(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+ &
                            Akx(i, j, k)*Asy(i, j, k)*duds(i, j, k)+Aky(i, j, k)**2*dvdk(i, j, k)*fourthird+ &
                            Aky(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird-Aky(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird-Aky(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+ &
                            Akz(i, j, k)**2*dvdk(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)

            dmudi_term(2) = Aix(i, j, k)**2*dvdi(i, j, k)+Aix(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*onethird+ &
                            Aix(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)+Aix(i, j, k)*Aky(i, j, k)*dudk(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+Aix(i, j, k)*Asy(i, j, k)*duds(i, j, k)+ &
                            Aiy(i, j, k)**2*dvdi(i, j, k)*fourthird+Aiy(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird- &
                            Aiy(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird+Aiy(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*fourthird- &
                            Aiy(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird-Aiy(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird-Aiy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+ &
                            Aiz(i, j, k)**2*dvdi(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)

            dmuds_term(2) = Aix(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)-Aix(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asx(i, j, k)*dudi(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird+ &
                            Aiy(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)-Aiz(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird+Aky(i, j, k)*Asx(i, j, k)*dudk(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird+Aky(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird+Akz(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)+ &
                            Asx(i, j, k)**2*dvds(i, j, k)+Asx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*onethird+ &
                            Asy(i, j, k)**2*dvds(i, j, k)*fourthird+Asy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird+ &
                            Asz(i, j, k)**2*dvds(i, j, k)

            mu_term(2) = Aix(i, j, k)**2*d2vdi2(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dvdi(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2udi2(i, j, k)*onethird+Aix(i, j, k)*Aiyi(i, j, k)*dudi(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2vdkdi(i, j, k)+Aix(i, j, k)*Akxi(i, j, k)*dvdk(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2udkdi(i, j, k)*onethird+Aix(i, j, k)*Akyi(i, j, k)*dudk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2vdids(i, j, k)+Aix(i, j, k)*Asxi(i, j, k)*dvds(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2udids(i, j, k)*onethird+Aix(i, j, k)*Asyi(i, j, k)*duds(i, j, k)- &
                         Aixi(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*twothird+Aixk(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)- &
                         Aixk(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird+Aixs(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)- &
                         Aixs(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird+Aiy(i, j, k)**2*d2vdi2(i, j, k)*fourthird+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dvdi(i, j, k)*fourthird+Aiy(i, j, k)*Aiz(i, j, k)*d2wdi2(i, j, k)*onethird- &
                         Aiy(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*twothird+Aiy(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*onethird- &
                         Aiy(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*twothird+Aiy(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*eighthird+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*fourthird+Aiy(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*onethird- &
                         Aiy(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*twothird+Aiy(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*onethird- &
                         Aiy(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*twothird+Aiy(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*eighthird+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*fourthird+Aiy(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*onethird- &
                         Aiy(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*twothird+Aiyi(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)+ &
                         Aiyk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)+Aiyk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird+ &
                         Aiyk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)+Aiys(i, j, k)*Asx(i, j, k)*dudi(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird+Aiys(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)+ &
                         Aiz(i, j, k)**2*d2vdi2(i, j, k)+Aiz(i, j, k)*Aizi(i, j, k)*dvdi(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2wdkdi(i, j, k)*onethird+Aiz(i, j, k)*Akyi(i, j, k)*dwdk(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2vdkdi(i, j, k)+Aiz(i, j, k)*Akzi(i, j, k)*dvdk(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2wdids(i, j, k)*onethird+Aiz(i, j, k)*Asyi(i, j, k)*dwds(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2vdids(i, j, k)+Aiz(i, j, k)*Aszi(i, j, k)*dvds(i, j, k)- &
                         Aizk(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird+Aizk(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)- &
                         Aizs(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird+Aizs(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)+ &
                         Akx(i, j, k)**2*d2vdk2(i, j, k)+Akx(i, j, k)*Akxk(i, j, k)*dvdk(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2udk2(i, j, k)*onethird+Akx(i, j, k)*Akyk(i, j, k)*dudk(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2vdkds(i, j, k)+Akx(i, j, k)*Asxk(i, j, k)*dvds(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2udkds(i, j, k)*onethird+Akx(i, j, k)*Asyk(i, j, k)*duds(i, j, k)- &
                         Akxk(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*twothird+Akxs(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)- &
                         Akxs(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird+Aky(i, j, k)**2*d2vdk2(i, j, k)*fourthird+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dvdk(i, j, k)*fourthird+Aky(i, j, k)*Akz(i, j, k)*d2wdk2(i, j, k)*onethird- &
                         Aky(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*twothird+Aky(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*onethird- &
                         Aky(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*eighthird+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*fourthird+Aky(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*onethird- &
                         Aky(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*twothird+Akyk(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)+ &
                         Akys(i, j, k)*Asx(i, j, k)*dudk(i, j, k)+Akys(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird+ &
                         Akys(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)+Akz(i, j, k)**2*d2vdk2(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dvdk(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*d2wdkds(i, j, k)*onethird+ &
                         Akz(i, j, k)*Asyk(i, j, k)*dwds(i, j, k)+2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2vdkds(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dvds(i, j, k)-Akzs(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)+Asx(i, j, k)**2*d2vds2(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dvds(i, j, k)+Asx(i, j, k)*Asy(i, j, k)*d2uds2(i, j, k)*onethird+ &
                         Asx(i, j, k)*Asys(i, j, k)*duds(i, j, k)-Asxs(i, j, k)*Asy(i, j, k)*duds(i, j, k)*twothird+ &
                         Asy(i, j, k)**2*d2vds2(i, j, k)*fourthird+Asy(i, j, k)*Asys(i, j, k)*dvds(i, j, k)*fourthird+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2wds2(i, j, k)*onethird-Asy(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*twothird+ &
                         Asys(i, j, k)*Asz(i, j, k)*dwds(i, j, k)+Asz(i, j, k)**2*d2vds2(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dvds(i, j, k)

            ! z momentum equation viscous source terms for variable viscosity case
            dmudk_term(3) = Aix(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)-Aix(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Akx(i, j, k)*dudi(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird+Akx(i, j, k)**2*dwdk(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+ &
                            Akx(i, j, k)*Asz(i, j, k)*duds(i, j, k)+Aky(i, j, k)**2*dwdk(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*onethird+Aky(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+ &
                            Aky(i, j, k)*Asz(i, j, k)*dvds(i, j, k)+Akz(i, j, k)**2*dwdk(i, j, k)*fourthird- &
                            Akz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird-Akz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird

            dmudi_term(3) = Aix(i, j, k)**2*dwdi(i, j, k)+Aix(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*onethird+ &
                            Aix(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)+Aix(i, j, k)*Akz(i, j, k)*dudk(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+Aix(i, j, k)*Asz(i, j, k)*duds(i, j, k)+ &
                            Aiy(i, j, k)**2*dwdi(i, j, k)+Aiy(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*onethird+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)+ &
                            Aiz(i, j, k)**2*dwdi(i, j, k)*fourthird-Aiz(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird- &
                            Aiz(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird+Aiz(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*fourthird- &
                            Aiz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird-Aiz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird

            dmuds_term(3) = Aix(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)-Aix(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asx(i, j, k)*dudi(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird+Akx(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird+Akz(i, j, k)*Asx(i, j, k)*dudk(i, j, k)+ &
                            Akz(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird+ &
                            Asx(i, j, k)**2*dwds(i, j, k)+Asx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*onethird+ &
                            Asy(i, j, k)**2*dwds(i, j, k)+Asy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*onethird+ &
                            Asz(i, j, k)**2*dwds(i, j, k)*fourthird

            mu_term(3) = Aix(i, j, k)**2*d2wdi2(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dwdi(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2udi2(i, j, k)*onethird+Aix(i, j, k)*Aizi(i, j, k)*dudi(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2wdkdi(i, j, k)+Aix(i, j, k)*Akxi(i, j, k)*dwdk(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2udkdi(i, j, k)*onethird+Aix(i, j, k)*Akzi(i, j, k)*dudk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2wdids(i, j, k)+Aix(i, j, k)*Asxi(i, j, k)*dwds(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2udids(i, j, k)*onethird+Aix(i, j, k)*Aszi(i, j, k)*duds(i, j, k)- &
                         Aixi(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*twothird+Aixk(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)- &
                         Aixk(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird+Aixs(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)- &
                         Aixs(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird+Aiy(i, j, k)**2*d2wdi2(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dwdi(i, j, k)+Aiy(i, j, k)*Aiz(i, j, k)*d2vdi2(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Aizi(i, j, k)*dvdi(i, j, k)+2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2wdkdi(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dwdk(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*d2vdkdi(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Akzi(i, j, k)*dvdk(i, j, k)+2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2wdids(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dwds(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*d2vdids(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Aszi(i, j, k)*dvds(i, j, k)-Aiyi(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)-Aiyk(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)-Aiys(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aiz(i, j, k)**2*d2wdi2(i, j, k)*fourthird+Aiz(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*fourthird+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*onethird-Aiz(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*twothird+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*onethird-Aiz(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*twothird+ &
                         Aiz(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*eighthird+Aiz(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*fourthird+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*onethird-Aiz(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*twothird+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*onethird-Aiz(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*twothird+ &
                         Aiz(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*eighthird+Aiz(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*fourthird+ &
                         Aizk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)+Aizk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird+Aizs(i, j, k)*Asx(i, j, k)*dudi(i, j, k)+ &
                         Aizs(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)+Aizs(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird+ &
                         Akx(i, j, k)**2*d2wdk2(i, j, k)+Akx(i, j, k)*Akxk(i, j, k)*dwdk(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2udk2(i, j, k)*onethird+Akx(i, j, k)*Akzk(i, j, k)*dudk(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2wdkds(i, j, k)+Akx(i, j, k)*Asxk(i, j, k)*dwds(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2udkds(i, j, k)*onethird+Akx(i, j, k)*Aszk(i, j, k)*duds(i, j, k)- &
                         Akxk(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*twothird+Akxs(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)- &
                         Akxs(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird+Aky(i, j, k)**2*d2wdk2(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dwdk(i, j, k)+Aky(i, j, k)*Akz(i, j, k)*d2vdk2(i, j, k)*onethird+ &
                         Aky(i, j, k)*Akzk(i, j, k)*dvdk(i, j, k)+2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2wdkds(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dwds(i, j, k)+Aky(i, j, k)*Asz(i, j, k)*d2vdkds(i, j, k)*onethird+ &
                         Aky(i, j, k)*Aszk(i, j, k)*dvds(i, j, k)-Akyk(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*twothird+ &
                         Akys(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)-Akys(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird+ &
                         Akz(i, j, k)**2*d2wdk2(i, j, k)*fourthird+Akz(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*fourthird+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*onethird-Akz(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*twothird+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*onethird-Akz(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*twothird+ &
                         Akz(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*eighthird+Akz(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*fourthird+ &
                         Akzs(i, j, k)*Asx(i, j, k)*dudk(i, j, k)+Akzs(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird+Asx(i, j, k)**2*d2wds2(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dwds(i, j, k)+Asx(i, j, k)*Asz(i, j, k)*d2uds2(i, j, k)*onethird+ &
                         Asx(i, j, k)*Aszs(i, j, k)*duds(i, j, k)-Asxs(i, j, k)*Asz(i, j, k)*duds(i, j, k)*twothird+ &
                         Asy(i, j, k)**2*d2wds2(i, j, k)+Asy(i, j, k)*Asys(i, j, k)*dwds(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2vds2(i, j, k)*onethird+Asy(i, j, k)*Aszs(i, j, k)*dvds(i, j, k)- &
                         Asys(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*twothird+Asz(i, j, k)**2*d2wds2(i, j, k)*fourthird+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*fourthird

            ! energy equation viscous source terms for variable viscosity case
            dmudk_term(4) = Aix(i, j, k)*Akx(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird*u(i, j, k)+ &
                            Aix(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*v(i, j, k)+Aix(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird*v(i, j, k)+ &
                            Aix(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*u(i, j, k)- &
                            Aix(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird*w(i, j, k)+ &
                            Aix(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*u(i, j, k)+Aiy(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*v(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird*u(i, j, k)+ &
                       Aiy(i, j, k)*Aky(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+Aiy(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird*v(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird*w(i, j, k)+ &
                            Aiy(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*v(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird*u(i, j, k)+ &
                            Aiz(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird*v(i, j, k)+ &
                       Aiz(i, j, k)*Akz(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+Aiz(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird*w(i, j, k)+ &
                            Akx(i, j, k)**2*Para%Cp*dTdk(i, j, k)/Para%Pr+Akx(i, j, k)**2*dudk(i, j, k)*fourthird*u(i, j, k)+ &
                            Akx(i, j, k)**2*dvdk(i, j, k)*v(i, j, k)+Akx(i, j, k)**2*dwdk(i, j, k)*w(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*onethird*v(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*onethird*u(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*onethird*w(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird*u(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                            Akx(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird*u(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*v(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                            Akx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*v(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*u(i, j, k)+ &
                            Akx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*w(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*u(i, j, k)+ &
                            Aky(i, j, k)**2*Para%Cp*dTdk(i, j, k)/Para%Pr+Aky(i, j, k)**2*dudk(i, j, k)*u(i, j, k)+ &
                            Aky(i, j, k)**2*dvdk(i, j, k)*fourthird*v(i, j, k)+Aky(i, j, k)**2*dwdk(i, j, k)*w(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*onethird*w(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird*v(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*v(i, j, k)+ &
                       Aky(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*u(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                            Aky(i, j, k)*Asy(i, j, k)*duds(i, j, k)*u(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird*v(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*w(i, j, k)+Aky(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*w(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*v(i, j, k)+ &
                            Akz(i, j, k)**2*Para%Cp*dTdk(i, j, k)/Para%Pr+Akz(i, j, k)**2*dudk(i, j, k)*u(i, j, k)+ &
                            Akz(i, j, k)**2*dvdk(i, j, k)*v(i, j, k)+Akz(i, j, k)**2*dwdk(i, j, k)*fourthird*w(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*w(i, j, k)+ &
                            Akz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*u(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*w(i, j, k)+ &
                       Akz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*v(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                            Akz(i, j, k)*Asz(i, j, k)*duds(i, j, k)*u(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*v(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird*w(i, j, k)

            dmudi_term(4) = Aix(i, j, k)**2*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)**2*dudi(i, j, k)*fourthird*u(i, j, k)+Aix(i, j, k)**2*dvdi(i, j, k)*v(i, j, k)+ &
                            Aix(i, j, k)**2*dwdi(i, j, k)*w(i, j, k)+Aix(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*onethird*v(i, j, k)+ &
                            Aix(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)*onethird*u(i, j, k)+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*onethird*w(i, j, k)+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird*u(i, j, k)+ &
                            Aix(i, j, k)*Akx(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*fourthird*u(i, j, k)+ &
                            Aix(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)*v(i, j, k)+Aix(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)*w(i, j, k)+ &
                            Aix(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*v(i, j, k)- &
                            Aix(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird*u(i, j, k)+ &
                            Aix(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird*u(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird*u(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*v(i, j, k)+Aix(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                            Aix(i, j, k)*Asy(i, j, k)*duds(i, j, k)*v(i, j, k)- &
                            Aix(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*u(i, j, k)+ &
                            Aix(i, j, k)*Asz(i, j, k)*duds(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*u(i, j, k)+ &
                            Aiy(i, j, k)**2*Para%Cp*dTdi(i, j, k)/Para%Pr+Aiy(i, j, k)**2*dudi(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)**2*dvdi(i, j, k)*fourthird*v(i, j, k)+Aiy(i, j, k)**2*dwdi(i, j, k)*w(i, j, k)+ &
                            Aiy(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*onethird*w(i, j, k)+ &
                            Aiy(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird*v(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird*v(i, j, k)+ &
                       Aiy(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)*u(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*fourthird*v(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)*w(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird*v(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*v(i, j, k)+ &
                       Aiy(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*u(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                            Aiy(i, j, k)*Asy(i, j, k)*duds(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird*v(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*w(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*v(i, j, k)+ &
                            Aiz(i, j, k)**2*Para%Cp*dTdi(i, j, k)/Para%Pr+Aiz(i, j, k)**2*dudi(i, j, k)*u(i, j, k)+ &
                            Aiz(i, j, k)**2*dvdi(i, j, k)*v(i, j, k)+Aiz(i, j, k)**2*dwdi(i, j, k)*fourthird*w(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird*w(i, j, k)+ &
                            Aiz(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)*u(i, j, k)- &
                            Aiz(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird*w(i, j, k)+ &
                       Aiz(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)*v(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*u(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*fourthird*w(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*w(i, j, k)+ &
                            Aiz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*u(i, j, k)- &
                            Aiz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*w(i, j, k)+ &
                       Aiz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*v(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                            Aiz(i, j, k)*Asz(i, j, k)*duds(i, j, k)*u(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird*w(i, j, k)

            dmuds_term(4) = Aix(i, j, k)*Asx(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird*u(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*v(i, j, k)+Aix(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird*v(i, j, k)+ &
                            Aix(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*u(i, j, k)- &
                            Aix(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird*w(i, j, k)+ &
                            Aix(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*u(i, j, k)+Aiy(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*v(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird*u(i, j, k)+ &
                       Aiy(i, j, k)*Asy(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+Aiy(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird*v(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird*w(i, j, k)+ &
                            Aiy(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*v(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird*u(i, j, k)+ &
                            Aiz(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird*v(i, j, k)+ &
                       Aiz(i, j, k)*Asz(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+Aiz(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird*w(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                            Akx(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird*u(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*v(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*w(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird*v(i, j, k)+ &
                            Akx(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*u(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird*w(i, j, k)+ &
                            Akx(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*u(i, j, k)+Aky(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*v(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird*u(i, j, k)+ &
                       Aky(i, j, k)*Asy(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+Aky(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*u(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird*v(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*w(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird*w(i, j, k)+ &
                            Aky(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*v(i, j, k)+Akz(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*w(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird*u(i, j, k)+ &
                            Akz(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*w(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird*v(i, j, k)+ &
                       Akz(i, j, k)*Asz(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+Akz(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*u(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*v(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird*w(i, j, k)+ &
                            Asx(i, j, k)**2*Para%Cp*dTds(i, j, k)/Para%Pr+Asx(i, j, k)**2*duds(i, j, k)*fourthird*u(i, j, k)+ &
                            Asx(i, j, k)**2*dvds(i, j, k)*v(i, j, k)+Asx(i, j, k)**2*dwds(i, j, k)*w(i, j, k)+ &
                            Asx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*onethird*v(i, j, k)+ &
                            Asx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*onethird*u(i, j, k)+ &
                            Asx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*onethird*w(i, j, k)+ &
                            Asx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird*u(i, j, k)+ &
                            Asy(i, j, k)**2*Para%Cp*dTds(i, j, k)/Para%Pr+Asy(i, j, k)**2*duds(i, j, k)*u(i, j, k)+ &
                            Asy(i, j, k)**2*dvds(i, j, k)*fourthird*v(i, j, k)+Asy(i, j, k)**2*dwds(i, j, k)*w(i, j, k)+ &
                            Asy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*onethird*w(i, j, k)+ &
                            Asy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird*v(i, j, k)+ &
                            Asz(i, j, k)**2*Para%Cp*dTds(i, j, k)/Para%Pr+Asz(i, j, k)**2*duds(i, j, k)*u(i, j, k)+ &
                            Asz(i, j, k)**2*dvds(i, j, k)*v(i, j, k)+Asz(i, j, k)**2*dwds(i, j, k)*fourthird*w(i, j, k)

            mu_term(4) = Aix(i, j, k)**2*Para%Cp*d2Tdi2(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)**2*d2udi2(i, j, k)*fourthird*u(i, j, k)+Aix(i, j, k)**2*d2vdi2(i, j, k)*v(i, j, k)+ &
                         Aix(i, j, k)**2*d2wdi2(i, j, k)*w(i, j, k)+Aix(i, j, k)**2*dudi(i, j, k)**2*fourthird+ &
                         Aix(i, j, k)**2*dvdi(i, j, k)**2+Aix(i, j, k)**2*dwdi(i, j, k)**2+ &
                         Aix(i, j, k)*Aixi(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dudi(i, j, k)*fourthird*u(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dvdi(i, j, k)*v(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dwdi(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2udi2(i, j, k)*onethird*v(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2vdi2(i, j, k)*onethird*u(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aix(i, j, k)*Aiyi(i, j, k)*dudi(i, j, k)*v(i, j, k)- &
                         Aix(i, j, k)*Aiyi(i, j, k)*dvdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2udi2(i, j, k)*onethird*w(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2wdi2(i, j, k)*onethird*u(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*dwdi(i, j, k)*twothird+ &
                         Aix(i, j, k)*Aizi(i, j, k)*dudi(i, j, k)*w(i, j, k)- &
                         Aix(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*twothird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*Para%Cp*d2Tdkdi(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*eighthird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2vdkdi(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2wdkdi(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*dudk(i, j, k)*eighthird+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*dvdk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*dwdk(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*fourthird*u(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dvdk(i, j, k)*v(i, j, k)+Aix(i, j, k)*Akxi(i, j, k)*dwdk(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2udkdi(i, j, k)*onethird*v(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*dvdk(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*dvdi(i, j, k)+ &
                         Aix(i, j, k)*Akyi(i, j, k)*dudk(i, j, k)*v(i, j, k)- &
                         Aix(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*twothird*u(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2udkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*dwdk(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*dwdi(i, j, k)+ &
                         Aix(i, j, k)*Akzi(i, j, k)*dudk(i, j, k)*w(i, j, k)- &
                         Aix(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*twothird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*Para%Cp*d2Tdids(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*eighthird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2vdids(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2wdids(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*duds(i, j, k)*eighthird+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*dvds(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*dwds(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*fourthird*u(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dvds(i, j, k)*v(i, j, k)+Aix(i, j, k)*Asxi(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2udids(i, j, k)*onethird*v(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*dvds(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dvdi(i, j, k)+ &
                         Aix(i, j, k)*Asyi(i, j, k)*duds(i, j, k)*v(i, j, k)- &
                         Aix(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*twothird*u(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2udids(i, j, k)*onethird*w(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dwdi(i, j, k)+ &
                         Aix(i, j, k)*Aszi(i, j, k)*duds(i, j, k)*w(i, j, k)- &
                         Aix(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*twothird*u(i, j, k)- &
                         Aixi(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*twothird*v(i, j, k)+ &
                         Aixi(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)*u(i, j, k)- &
                         Aixi(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*twothird*w(i, j, k)+ &
                         Aixi(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*u(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird*u(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*v(i, j, k)+Aixk(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                         Aixk(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird*v(i, j, k)+ &
                         Aixk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*u(i, j, k)- &
                         Aixk(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird*w(i, j, k)+ &
                         Aixk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*u(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird*u(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*v(i, j, k)+Aixs(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                         Aixs(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird*v(i, j, k)+ &
                         Aixs(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*u(i, j, k)- &
                         Aixs(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird*w(i, j, k)+ &
                         Aixs(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*u(i, j, k)+Aiy(i, j, k)**2*Para%Cp*d2Tdi2(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)**2*d2udi2(i, j, k)*u(i, j, k)+Aiy(i, j, k)**2*d2vdi2(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)**2*d2wdi2(i, j, k)*w(i, j, k)+Aiy(i, j, k)**2*dudi(i, j, k)**2+ &
                         Aiy(i, j, k)**2*dvdi(i, j, k)**2*fourthird+Aiy(i, j, k)**2*dwdi(i, j, k)**2+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dvdi(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dwdi(i, j, k)*w(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*d2vdi2(i, j, k)*onethird*w(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*d2wdi2(i, j, k)*onethird*v(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*dwdi(i, j, k)*twothird+ &
                         Aiy(i, j, k)*Aizi(i, j, k)*dvdi(i, j, k)*w(i, j, k)- &
                         Aiy(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*onethird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akx(i, j, k)*d2vdkdi(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*dvdk(i, j, k)- &
                         Aiy(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*dvdi(i, j, k)*fourthird- &
                         Aiy(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akxi(i, j, k)*dvdk(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*Para%Cp*d2Tdkdi(i, j, k)/Para%Pr+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2udkdi(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*eighthird*v(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2wdkdi(i, j, k)*w(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*dudk(i, j, k)+ &
                         Aiy(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*dvdk(i, j, k)*eighthird+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*dwdk(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dudk(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dwdk(i, j, k)*w(i, j, k)+ &
                         Aiy(i, j, k)*Akz(i, j, k)*d2vdkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aiy(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*onethird*v(i, j, k)- &
                         Aiy(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*dwdk(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*dwdi(i, j, k)+ &
                         Aiy(i, j, k)*Akzi(i, j, k)*dvdk(i, j, k)*w(i, j, k)- &
                         Aiy(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*onethird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asx(i, j, k)*d2vdids(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*dvds(i, j, k)- &
                         Aiy(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dvdi(i, j, k)*fourthird- &
                         Aiy(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asxi(i, j, k)*dvds(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*Para%Cp*d2Tdids(i, j, k)/Para%Pr+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2udids(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*eighthird*v(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2wdids(i, j, k)*w(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*duds(i, j, k)+ &
                         Aiy(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*dvds(i, j, k)*eighthird+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*dwds(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*duds(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                         Aiy(i, j, k)*Asz(i, j, k)*d2vdids(i, j, k)*onethird*w(i, j, k)+ &
                         Aiy(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*onethird*v(i, j, k)- &
                         Aiy(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dwdi(i, j, k)+ &
                         Aiy(i, j, k)*Aszi(i, j, k)*dvds(i, j, k)*w(i, j, k)- &
                         Aiy(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*twothird*v(i, j, k)- &
                         Aiyi(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aiyi(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*v(i, j, k)+Aiyk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*v(i, j, k)- &
                         Aiyk(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                         Aiyk(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aiyk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*v(i, j, k)+Aiys(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*v(i, j, k)- &
                         Aiys(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*u(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*w(i, j, k)- &
                         Aiys(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aiys(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*v(i, j, k)+Aiz(i, j, k)**2*Para%Cp*d2Tdi2(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)**2*d2udi2(i, j, k)*u(i, j, k)+Aiz(i, j, k)**2*d2vdi2(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)**2*d2wdi2(i, j, k)*fourthird*w(i, j, k)+Aiz(i, j, k)**2*dudi(i, j, k)**2+ &
                         Aiz(i, j, k)**2*dvdi(i, j, k)**2+Aiz(i, j, k)**2*dwdi(i, j, k)**2*fourthird+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dudi(i, j, k)*u(i, j, k)+Aiz(i, j, k)*Aizi(i, j, k)*dvdi(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*fourthird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2wdkdi(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*dwdk(i, j, k)- &
                         Aiz(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akxi(i, j, k)*dwdk(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2wdkdi(i, j, k)*onethird*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*dwdk(i, j, k)- &
                         Aiz(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akyi(i, j, k)*dwdk(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*Para%Cp*d2Tdkdi(i, j, k)/Para%Pr+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2udkdi(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2vdkdi(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*eighthird*w(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*dudk(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*dvdk(i, j, k)+ &
                         Aiz(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*dwdk(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dudk(i, j, k)*u(i, j, k)+Aiz(i, j, k)*Akzi(i, j, k)*dvdk(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*fourthird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2wdids(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*dwds(i, j, k)- &
                         Aiz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asxi(i, j, k)*dwds(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2wdids(i, j, k)*onethird*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*dwds(i, j, k)- &
                         Aiz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asyi(i, j, k)*dwds(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*Para%Cp*d2Tdids(i, j, k)/Para%Pr+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2udids(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2vdids(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*eighthird*w(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*duds(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*dvds(i, j, k)+ &
                         Aiz(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*dwds(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*duds(i, j, k)*u(i, j, k)+Aiz(i, j, k)*Aszi(i, j, k)*dvds(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*fourthird*w(i, j, k)+ &
                         Aizk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*w(i, j, k)- &
                         Aizk(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aizk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*w(i, j, k)- &
                         Aizk(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*u(i, j, k)+Aizk(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*v(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird*w(i, j, k)+ &
                         Aizs(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*w(i, j, k)- &
                         Aizs(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aizs(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*w(i, j, k)- &
                         Aizs(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*u(i, j, k)+Aizs(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*v(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird*w(i, j, k)+ &
                         Akx(i, j, k)**2*Para%Cp*d2Tdk2(i, j, k)/Para%Pr+Akx(i, j, k)**2*d2udk2(i, j, k)*fourthird*u(i, j, k)+ &
                         Akx(i, j, k)**2*d2vdk2(i, j, k)*v(i, j, k)+Akx(i, j, k)**2*d2wdk2(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)**2*dudk(i, j, k)**2*fourthird+Akx(i, j, k)**2*dvdk(i, j, k)**2+ &
                         Akx(i, j, k)**2*dwdk(i, j, k)**2+Akx(i, j, k)*Akxk(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dudk(i, j, k)*fourthird*u(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dvdk(i, j, k)*v(i, j, k)+Akx(i, j, k)*Akxk(i, j, k)*dwdk(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2udk2(i, j, k)*onethird*v(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2vdk2(i, j, k)*onethird*u(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*dvdk(i, j, k)*twothird+ &
                         Akx(i, j, k)*Akyk(i, j, k)*dudk(i, j, k)*v(i, j, k)- &
                         Akx(i, j, k)*Akyk(i, j, k)*dvdk(i, j, k)*twothird*u(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2udk2(i, j, k)*onethird*w(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2wdk2(i, j, k)*onethird*u(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*dwdk(i, j, k)*twothird+ &
                         Akx(i, j, k)*Akzk(i, j, k)*dudk(i, j, k)*w(i, j, k)- &
                         Akx(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*twothird*u(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*Para%Cp*d2Tdkds(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*eighthird*u(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2vdkds(i, j, k)*v(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2wdkds(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*duds(i, j, k)*eighthird+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*dvds(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*dwds(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*fourthird*u(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dvds(i, j, k)*v(i, j, k)+Akx(i, j, k)*Asxk(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2udkds(i, j, k)*onethird*v(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*onethird*u(i, j, k)- &
                         Akx(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*dvds(i, j, k)*fourthird+ &
                         2.0d0*Akx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dvdk(i, j, k)+ &
                         Akx(i, j, k)*Asyk(i, j, k)*duds(i, j, k)*v(i, j, k)- &
                         Akx(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*twothird*u(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2udkds(i, j, k)*onethird*w(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*onethird*u(i, j, k)- &
                         Akx(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Akx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dwdk(i, j, k)+ &
                         Akx(i, j, k)*Aszk(i, j, k)*duds(i, j, k)*w(i, j, k)- &
                         Akx(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*twothird*u(i, j, k)- &
                         Akxk(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*twothird*v(i, j, k)+ &
                         Akxk(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*u(i, j, k)- &
                         Akxk(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*twothird*w(i, j, k)+ &
                         Akxk(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*u(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird*u(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*v(i, j, k)+Akxs(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*w(i, j, k)- &
                         Akxs(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird*v(i, j, k)+ &
                         Akxs(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*u(i, j, k)- &
                         Akxs(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird*w(i, j, k)+ &
                         Akxs(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*u(i, j, k)+Aky(i, j, k)**2*Para%Cp*d2Tdk2(i, j, k)/Para%Pr+ &
                         Aky(i, j, k)**2*d2udk2(i, j, k)*u(i, j, k)+Aky(i, j, k)**2*d2vdk2(i, j, k)*fourthird*v(i, j, k)+ &
                         Aky(i, j, k)**2*d2wdk2(i, j, k)*w(i, j, k)+Aky(i, j, k)**2*dudk(i, j, k)**2+ &
                         Aky(i, j, k)**2*dvdk(i, j, k)**2*fourthird+Aky(i, j, k)**2*dwdk(i, j, k)**2+ &
                         Aky(i, j, k)*Akyk(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dudk(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dvdk(i, j, k)*fourthird*v(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dwdk(i, j, k)*w(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2vdk2(i, j, k)*onethird*w(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2wdk2(i, j, k)*onethird*v(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*dwdk(i, j, k)*twothird+ &
                         Aky(i, j, k)*Akzk(i, j, k)*dvdk(i, j, k)*w(i, j, k)- &
                         Aky(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*twothird*v(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*onethird*v(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2vdkds(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*dvds(i, j, k)- &
                         Aky(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dvdk(i, j, k)*fourthird- &
                         Aky(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*twothird*v(i, j, k)+ &
                         Aky(i, j, k)*Asxk(i, j, k)*dvds(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*Para%Cp*d2Tdkds(i, j, k)/Para%Pr+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2udkds(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*eighthird*v(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2wdkds(i, j, k)*w(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*duds(i, j, k)+ &
                         Aky(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*dvds(i, j, k)*eighthird+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*dwds(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Aky(i, j, k)*Asyk(i, j, k)*duds(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*fourthird*v(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2vdkds(i, j, k)*onethird*w(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*onethird*v(i, j, k)- &
                         Aky(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Aky(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dwdk(i, j, k)+ &
                         Aky(i, j, k)*Aszk(i, j, k)*dvds(i, j, k)*w(i, j, k)- &
                         Aky(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*twothird*v(i, j, k)- &
                         Akyk(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*twothird*w(i, j, k)+ &
                         Akyk(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*v(i, j, k)+Akys(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*v(i, j, k)- &
                         Akys(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird*u(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Akys(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*u(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird*v(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*w(i, j, k)- &
                         Akys(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird*w(i, j, k)+ &
                         Akys(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*v(i, j, k)+Akz(i, j, k)**2*Para%Cp*d2Tdk2(i, j, k)/Para%Pr+ &
                         Akz(i, j, k)**2*d2udk2(i, j, k)*u(i, j, k)+Akz(i, j, k)**2*d2vdk2(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)**2*d2wdk2(i, j, k)*fourthird*w(i, j, k)+Akz(i, j, k)**2*dudk(i, j, k)**2+ &
                         Akz(i, j, k)**2*dvdk(i, j, k)**2+Akz(i, j, k)**2*dwdk(i, j, k)**2*fourthird+ &
                         Akz(i, j, k)*Akzk(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dudk(i, j, k)*u(i, j, k)+Akz(i, j, k)*Akzk(i, j, k)*dvdk(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*fourthird*w(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*onethird*w(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2wdkds(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*dwds(i, j, k)- &
                         Akz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dwdk(i, j, k)*fourthird- &
                         Akz(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*twothird*w(i, j, k)+ &
                         Akz(i, j, k)*Asxk(i, j, k)*dwds(i, j, k)*u(i, j, k)+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*onethird*w(i, j, k)+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2wdkds(i, j, k)*onethird*v(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*dwds(i, j, k)- &
                         Akz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*dwdk(i, j, k)*fourthird- &
                         Akz(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*twothird*w(i, j, k)+ &
                         Akz(i, j, k)*Asyk(i, j, k)*dwds(i, j, k)*v(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*Para%Cp*d2Tdkds(i, j, k)/Para%Pr+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2udkds(i, j, k)*u(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2vdkds(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*eighthird*w(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*duds(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*dvds(i, j, k)+ &
                         Akz(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*dwds(i, j, k)*eighthird+ &
                         Akz(i, j, k)*Aszk(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Akz(i, j, k)*Aszk(i, j, k)*duds(i, j, k)*u(i, j, k)+Akz(i, j, k)*Aszk(i, j, k)*dvds(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*fourthird*w(i, j, k)+ &
                         Akzs(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*w(i, j, k)- &
                         Akzs(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird*u(i, j, k)+ &
                         Akzs(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*w(i, j, k)- &
                         Akzs(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird*v(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*Para%Cp*dTdk(i, j, k)/Para%Pr+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*u(i, j, k)+Akzs(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*v(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird*w(i, j, k)+ &
                         Asx(i, j, k)**2*Para%Cp*d2Tds2(i, j, k)/Para%Pr+Asx(i, j, k)**2*d2uds2(i, j, k)*fourthird*u(i, j, k)+ &
                         Asx(i, j, k)**2*d2vds2(i, j, k)*v(i, j, k)+Asx(i, j, k)**2*d2wds2(i, j, k)*w(i, j, k)+ &
                         Asx(i, j, k)**2*duds(i, j, k)**2*fourthird+Asx(i, j, k)**2*dvds(i, j, k)**2+ &
                         Asx(i, j, k)**2*dwds(i, j, k)**2+Asx(i, j, k)*Asxs(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Asx(i, j, k)*Asxs(i, j, k)*duds(i, j, k)*fourthird*u(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dvds(i, j, k)*v(i, j, k)+Asx(i, j, k)*Asxs(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2uds2(i, j, k)*onethird*v(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2vds2(i, j, k)*onethird*u(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dvds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Asys(i, j, k)*duds(i, j, k)*v(i, j, k)- &
                         Asx(i, j, k)*Asys(i, j, k)*dvds(i, j, k)*twothird*u(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2uds2(i, j, k)*onethird*w(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2wds2(i, j, k)*onethird*u(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dwds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Aszs(i, j, k)*duds(i, j, k)*w(i, j, k)- &
                         Asx(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*twothird*u(i, j, k)- &
                         Asxs(i, j, k)*Asy(i, j, k)*duds(i, j, k)*twothird*v(i, j, k)+ &
                         Asxs(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*u(i, j, k)- &
                         Asxs(i, j, k)*Asz(i, j, k)*duds(i, j, k)*twothird*w(i, j, k)+ &
                         Asxs(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*u(i, j, k)+Asy(i, j, k)**2*Para%Cp*d2Tds2(i, j, k)/Para%Pr+ &
                         Asy(i, j, k)**2*d2uds2(i, j, k)*u(i, j, k)+Asy(i, j, k)**2*d2vds2(i, j, k)*fourthird*v(i, j, k)+ &
                         Asy(i, j, k)**2*d2wds2(i, j, k)*w(i, j, k)+Asy(i, j, k)**2*duds(i, j, k)**2+ &
                         Asy(i, j, k)**2*dvds(i, j, k)**2*fourthird+Asy(i, j, k)**2*dwds(i, j, k)**2+ &
                         Asy(i, j, k)*Asys(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Asy(i, j, k)*Asys(i, j, k)*duds(i, j, k)*u(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dvds(i, j, k)*fourthird*v(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dwds(i, j, k)*w(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2vds2(i, j, k)*onethird*w(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2wds2(i, j, k)*onethird*v(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dwds(i, j, k)*twothird+ &
                         Asy(i, j, k)*Aszs(i, j, k)*dvds(i, j, k)*w(i, j, k)- &
                         Asy(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*twothird*v(i, j, k)- &
                         Asys(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*twothird*w(i, j, k)+ &
                         Asys(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*v(i, j, k)+Asz(i, j, k)**2*Para%Cp*d2Tds2(i, j, k)/Para%Pr+ &
                         Asz(i, j, k)**2*d2uds2(i, j, k)*u(i, j, k)+Asz(i, j, k)**2*d2vds2(i, j, k)*v(i, j, k)+ &
                         Asz(i, j, k)**2*d2wds2(i, j, k)*fourthird*w(i, j, k)+Asz(i, j, k)**2*duds(i, j, k)**2+ &
                         Asz(i, j, k)**2*dvds(i, j, k)**2+Asz(i, j, k)**2*dwds(i, j, k)**2*fourthird+ &
                         Asz(i, j, k)*Aszs(i, j, k)*Para%Cp*dTds(i, j, k)/Para%Pr+ &
                         Asz(i, j, k)*Aszs(i, j, k)*duds(i, j, k)*u(i, j, k)+Asz(i, j, k)*Aszs(i, j, k)*dvds(i, j, k)*v(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*fourthird*w(i, j, k)

            ! add viscous source terms to du
            du(i, j, k, 2) = du(i, j, k, 2)+ &
                             dmudk(i, j, k)*dmudk_term(1)+ &
                             dmudi(i, j, k)*dmudi_term(1)+ &
                             dmuds(i, j, k)*dmuds_term(1)+ &
                             mu(i, j, k)*mu_term(1)
            du(i, j, k, 3) = du(i, j, k, 3)+ &
                             dmudk(i, j, k)*dmudk_term(2)+ &
                             dmudi(i, j, k)*dmudi_term(2)+ &
                             dmuds(i, j, k)*dmuds_term(2)+ &
                             mu(i, j, k)*mu_term(2)
            du(i, j, k, 4) = du(i, j, k, 4)+ &
                             dmudk(i, j, k)*dmudk_term(3)+ &
                             dmudi(i, j, k)*dmudi_term(3)+ &
                             dmuds(i, j, k)*dmuds_term(3)+ &
                             mu(i, j, k)*mu_term(3)
            du(i, j, k, 5) = du(i, j, k, 5)+ &
                             dmudk(i, j, k)*dmudk_term(4)+ &
                             dmudi(i, j, k)*dmudi_term(4)+ &
                             dmuds(i, j, k)*dmuds_term(4)+ &
                             mu(i, j, k)*mu_term(4)
        end do
        end do
        end do

    end subroutine flux_viscous3d

    subroutine compute_mu
        use mesh_data
        use flow_data
        implicit none

        real(dp):: Tsb, tmpR
        integer:: i, j, k

        Tsb = 110.4d0/Para%T_inf
        TmpR = 1.0d0/Para%Re

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            mu(i, j, k) = TmpR*(1.0d0+Tsb)*sqrt(T(i, j, k)**3)/(Tsb+T(i, j, k))
        end do
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


    subroutine compute_derivatives
        ! Velocity and temperature derivatives used by the viscous terms; also reused by the
        ! lid-driven cavity wall-pressure condition, which computes them right after setting the walls
        use mesh_data
        use flow_data
        use Centrals
        implicit none

        call OCFD_dx0(u, dudk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(v, dvdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(w, dwdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(T, dTdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        !
        call OCFD_dy0(u, dudi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(v, dvdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(w, dwdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(T, dTdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        !
        call OCFD_dz0(u, duds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(v, dvds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(w, dwds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(T, dTds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        !
        ! Mixed derivatives use k-derivatives in y/z and i-derivatives in z.
        call mpi_exchange_yz3d(dudk)
        call mpi_exchange_yz3d(dvdk)
        call mpi_exchange_yz3d(dwdk)
        call mpi_exchange_yz3d(dTdk)
        call mpi_exchange_z3d(dudi)
        call mpi_exchange_z3d(dvdi)
        call mpi_exchange_z3d(dwdi)
        call mpi_exchange_z3d(dTdi)
        !
        call OCFD_dy0(dudk, d2udkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dvdk, d2vdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dwdk, d2wdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dTdk, d2Tdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dz0(dudk, d2udkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dvdk, d2vdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dwdk, d2wdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dTdk, d2Tdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dudi, d2udids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dvdi, d2vdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dwdi, d2wdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dTdi, d2Tdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        !
        call OCFD_dx2(u, d2udk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(v, d2vdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(w, d2wdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(T, d2Tdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dy2(u, d2udi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(v, d2vdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(w, d2wdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(T, d2Tdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dz2(u, d2uds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(v, d2vds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(w, d2wds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(T, d2Tds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)

    end subroutine compute_derivatives

end module flux_viscous
