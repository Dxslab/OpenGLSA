! SPDX-License-Identifier: Apache-2.0
! Copyright 2021 Xinliang Li
! Copyright by Li Xinliang
! Modifications Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
! Portions adapted from OpenCFD-SC and modified for OpenGLSA.
! See NOTICE.

module Flux_Viscous_d
    implicit none

contains

    subroutine flux_viscous3d_d
        use mesh_data
        use flow_data
        use fluc_data
        use Centrals
        implicit none

        integer :: i, j, k
        real(dp) :: dmufdk_term(4), dmufdi_term(4), dmufds_term(4), muf_term(4)
        real(dp) :: dmudk_term(4), dmudi_term(4), dmuds_term(4), mu_term(4)

        real(dp), parameter :: onethird = 1.0d0/3.0d0
        real(dp), parameter :: twothird = 2.0d0/3.0d0
        real(dp), parameter :: fourthird = 4.0d0/3.0d0
        real(dp), parameter :: eighthird = 8.0d0/3.0d0

        ! For the lid-driven cavity the derivatives were already computed in the wall condition
        if (Para%IBC /= BC_Lid_Driven_Cavity) call compute_derivatives_d

        call comput_mu_d
        dmufdk = diff_mu_d(dTfdk, T_fluc, dTdk, T)
        dmufdi = diff_mu_d(dTfdi, T_fluc, dTdi, T)
        dmufds = diff_mu_d(dTfds, T_fluc, dTds, T)

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            dmufdk_term(1) = Aix(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird+Aix(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)+ &
                             Aix(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Aky(i, j, k)*dudi(i, j, k)-Aiz(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Akz(i, j, k)*dudi(i, j, k)+Akx(i, j, k)**2*dudk(i, j, k)*fourthird+ &
                             Akx(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*onethird+Akx(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird+ &
                             Akx(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird-Akx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird- &
                             Akx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+Aky(i, j, k)**2*dudk(i, j, k)+ &
                             Aky(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*duds(i, j, k)+ &
                             Akz(i, j, k)**2*dudk(i, j, k)+Akz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+ &
                             Akz(i, j, k)*Asz(i, j, k)*duds(i, j, k)

            dmufdi_term(1) = Aix(i, j, k)**2*dudi(i, j, k)*fourthird+Aix(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)*onethird+ &
                             Aix(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird+Aix(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*fourthird- &
                             Aix(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird-Aix(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird+ &
                             Aix(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird-Aix(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird- &
                             Aix(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+Aiy(i, j, k)**2*dudi(i, j, k)+ &
                             Aiy(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dudk(i, j, k)+ &
                             Aiy(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*duds(i, j, k)+ &
                             Aiz(i, j, k)**2*dudi(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)+ &
                             Aiz(i, j, k)*Akz(i, j, k)*dudk(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+ &
                             Aiz(i, j, k)*Asz(i, j, k)*duds(i, j, k)

            dmufds_term(1) = Aix(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird+Aix(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)+ &
                             Aix(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Asy(i, j, k)*dudi(i, j, k)-Aiz(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Asz(i, j, k)*dudi(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird+ &
                             Akx(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)+Akx(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)- &
                             Aky(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*dudk(i, j, k)- &
                             Akz(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird+Akz(i, j, k)*Asz(i, j, k)*dudk(i, j, k)+ &
                             Asx(i, j, k)**2*duds(i, j, k)*fourthird+Asx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*onethird+ &
                             Asx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird+Asy(i, j, k)**2*duds(i, j, k)+ &
                             Asz(i, j, k)**2*duds(i, j, k)

            muf_term(1) = Aix(i, j, k)**2*d2udi2(i, j, k)*fourthird+ &
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

            dmudk_term(1) = Aix(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*fourthird+ &
                            Aix(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)+Aix(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)*twothird+Aiy(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)*twothird+Aiz(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)+ &
                            Akx(i, j, k)**2*dufdk(i, j, k)*fourthird+Akx(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*onethird+ &
                            Akx(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*fourthird- &
                            Akx(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird-Akx(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird+ &
                            Aky(i, j, k)**2*dufdk(i, j, k)+Aky(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dufds(i, j, k)+Akz(i, j, k)**2*dufdk(i, j, k)+ &
                            Akz(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dufds(i, j, k)

            dmudi_term(1) = Aix(i, j, k)**2*dufdi(i, j, k)*fourthird+ &
                            Aix(i, j, k)*Aiy(i, j, k)*dvfdi(i, j, k)*onethird+Aix(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)*onethird+ &
                            Aix(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*fourthird-Aix(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*twothird- &
                            Aix(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*twothird+Aix(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*fourthird- &
                            Aix(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird-Aix(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird+ &
                            Aiy(i, j, k)**2*dufdi(i, j, k)+Aiy(i, j, k)*Akx(i, j, k)*dvfdk(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)+Aiy(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dufds(i, j, k)+Aiz(i, j, k)**2*dufdi(i, j, k)+ &
                            Aiz(i, j, k)*Akx(i, j, k)*dwfdk(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)+ &
                            Aiz(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*dufds(i, j, k)

            dmuds_term(1) = Aix(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*fourthird+ &
                            Aix(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)+Aix(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)*twothird+Aiy(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)*twothird+Aiz(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*fourthird+Akx(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)+ &
                            Akx(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)-Aky(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)*twothird+ &
                            Aky(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)-Akz(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)*twothird+ &
                            Akz(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)+Asx(i, j, k)**2*dufds(i, j, k)*fourthird+ &
                            Asx(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*onethird+Asx(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*onethird+ &
                            Asy(i, j, k)**2*dufds(i, j, k)+Asz(i, j, k)**2*dufds(i, j, k)

            mu_term(1) = Aix(i, j, k)**2*d2ufdi2(i, j, k)*fourthird+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dufdi(i, j, k)*fourthird+Aix(i, j, k)*Aiy(i, j, k)*d2vfdi2(i, j, k)*onethird- &
                         Aix(i, j, k)*Aiyi(i, j, k)*dvfdi(i, j, k)*twothird+Aix(i, j, k)*Aiz(i, j, k)*d2wfdi2(i, j, k)*onethird- &
                         Aix(i, j, k)*Aizi(i, j, k)*dwfdi(i, j, k)*twothird+Aix(i, j, k)*Akx(i, j, k)*d2ufdkdi(i, j, k)*eighthird+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dufdk(i, j, k)*fourthird+Aix(i, j, k)*Aky(i, j, k)*d2vfdkdi(i, j, k)*onethird- &
                         Aix(i, j, k)*Akyi(i, j, k)*dvfdk(i, j, k)*twothird+Aix(i, j, k)*Akz(i, j, k)*d2wfdkdi(i, j, k)*onethird- &
                         Aix(i, j, k)*Akzi(i, j, k)*dwfdk(i, j, k)*twothird+Aix(i, j, k)*Asx(i, j, k)*d2ufdids(i, j, k)*eighthird+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dufds(i, j, k)*fourthird+Aix(i, j, k)*Asy(i, j, k)*d2vfdids(i, j, k)*onethird- &
                         Aix(i, j, k)*Asyi(i, j, k)*dvfds(i, j, k)*twothird+Aix(i, j, k)*Asz(i, j, k)*d2wfdids(i, j, k)*onethird- &
                         Aix(i, j, k)*Aszi(i, j, k)*dwfds(i, j, k)*twothird+Aixi(i, j, k)*Aiy(i, j, k)*dvfdi(i, j, k)+ &
                         Aixi(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)+Aixk(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*fourthird+ &
                         Aixk(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)+Aixk(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*fourthird+Aixs(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)+ &
                         Aixs(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)+Aiy(i, j, k)**2*d2ufdi2(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dufdi(i, j, k)+Aiy(i, j, k)*Akx(i, j, k)*d2vfdkdi(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Akxi(i, j, k)*dvfdk(i, j, k)+2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2ufdkdi(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dufdk(i, j, k)+Aiy(i, j, k)*Asx(i, j, k)*d2vfdids(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Asxi(i, j, k)*dvfds(i, j, k)+2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2ufdids(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dufds(i, j, k)-Aiyk(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)*twothird+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)-Aiys(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)*twothird+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)+Aiz(i, j, k)**2*d2ufdi2(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dufdi(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*d2wfdkdi(i, j, k)*onethird+ &
                         Aiz(i, j, k)*Akxi(i, j, k)*dwfdk(i, j, k)+2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2ufdkdi(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dufdk(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*d2wfdids(i, j, k)*onethird+ &
                         Aiz(i, j, k)*Asxi(i, j, k)*dwfds(i, j, k)+2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2ufdids(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dufds(i, j, k)-Aizk(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)*twothird+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)-Aizs(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)*twothird+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)+Akx(i, j, k)**2*d2ufdk2(i, j, k)*fourthird+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dufdk(i, j, k)*fourthird+Akx(i, j, k)*Aky(i, j, k)*d2vfdk2(i, j, k)*onethird- &
                         Akx(i, j, k)*Akyk(i, j, k)*dvfdk(i, j, k)*twothird+Akx(i, j, k)*Akz(i, j, k)*d2wfdk2(i, j, k)*onethird- &
                         Akx(i, j, k)*Akzk(i, j, k)*dwfdk(i, j, k)*twothird+Akx(i, j, k)*Asx(i, j, k)*d2ufdkds(i, j, k)*eighthird+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dufds(i, j, k)*fourthird+Akx(i, j, k)*Asy(i, j, k)*d2vfdkds(i, j, k)*onethird- &
                         Akx(i, j, k)*Asyk(i, j, k)*dvfds(i, j, k)*twothird+Akx(i, j, k)*Asz(i, j, k)*d2wfdkds(i, j, k)*onethird- &
                         Akx(i, j, k)*Aszk(i, j, k)*dwfds(i, j, k)*twothird+Akxk(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)+ &
                         Akxk(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)+Akxs(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*fourthird+ &
                         Akxs(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)+Akxs(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)+ &
                         Aky(i, j, k)**2*d2ufdk2(i, j, k)+Aky(i, j, k)*Akyk(i, j, k)*dufdk(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2vfdkds(i, j, k)*onethird+Aky(i, j, k)*Asxk(i, j, k)*dvfds(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2ufdkds(i, j, k)+Aky(i, j, k)*Asyk(i, j, k)*dufds(i, j, k)- &
                         Akys(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)*twothird+Akys(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)+ &
                         Akz(i, j, k)**2*d2ufdk2(i, j, k)+Akz(i, j, k)*Akzk(i, j, k)*dufdk(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2wfdkds(i, j, k)*onethird+Akz(i, j, k)*Asxk(i, j, k)*dwfds(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2ufdkds(i, j, k)+Akz(i, j, k)*Aszk(i, j, k)*dufds(i, j, k)- &
                         Akzs(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)*twothird+Akzs(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)+ &
                         Asx(i, j, k)**2*d2ufds2(i, j, k)*fourthird+Asx(i, j, k)*Asxs(i, j, k)*dufds(i, j, k)*fourthird+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2vfds2(i, j, k)*onethird-Asx(i, j, k)*Asys(i, j, k)*dvfds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2wfds2(i, j, k)*onethird-Asx(i, j, k)*Aszs(i, j, k)*dwfds(i, j, k)*twothird+ &
                         Asxs(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)+Asxs(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)+ &
                         Asy(i, j, k)**2*d2ufds2(i, j, k)+Asy(i, j, k)*Asys(i, j, k)*dufds(i, j, k)+ &
                         Asz(i, j, k)**2*d2ufds2(i, j, k)+Asz(i, j, k)*Aszs(i, j, k)*dufds(i, j, k)

            dmufdk_term(2) = Aix(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)-Aix(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Akx(i, j, k)*dudi(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird+ &
                             Aiy(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)-Aiz(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)+Akx(i, j, k)**2*dvdk(i, j, k)+ &
                             Akx(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+ &
                             Akx(i, j, k)*Asy(i, j, k)*duds(i, j, k)+Aky(i, j, k)**2*dvdk(i, j, k)*fourthird+ &
                             Aky(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird-Aky(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird+ &
                             Aky(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird-Aky(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+ &
                             Akz(i, j, k)**2*dvdk(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+ &
                             Akz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)

            dmufdi_term(2) = Aix(i, j, k)**2*dvdi(i, j, k)+Aix(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*onethird+ &
                             Aix(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)+Aix(i, j, k)*Aky(i, j, k)*dudk(i, j, k)+ &
                             Aix(i, j, k)*Asx(i, j, k)*dvds(i, j, k)+Aix(i, j, k)*Asy(i, j, k)*duds(i, j, k)+ &
                             Aiy(i, j, k)**2*dvdi(i, j, k)*fourthird+Aiy(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird- &
                             Aiy(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird+Aiy(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*fourthird- &
                             Aiy(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird-Aiy(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird-Aiy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird+ &
                             Aiz(i, j, k)**2*dvdi(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)+ &
                             Aiz(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+ &
                             Aiz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)

            dmufds_term(2) = Aix(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)-Aix(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Asx(i, j, k)*dudi(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird+ &
                             Aiy(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)-Aiz(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)- &
                             Akx(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird+Aky(i, j, k)*Asx(i, j, k)*dudk(i, j, k)+ &
                             Aky(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird+Aky(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)- &
                             Akz(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird+Akz(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)+ &
                             Asx(i, j, k)**2*dvds(i, j, k)+Asx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*onethird+ &
                             Asy(i, j, k)**2*dvds(i, j, k)*fourthird+Asy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird+ &
                             Asz(i, j, k)**2*dvds(i, j, k)

            muf_term(2) = Aix(i, j, k)**2*d2vdi2(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dvdi(i, j, k)+ &
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

            dmudk_term(2) = Aix(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)-Aix(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*fourthird+ &
                            Aiy(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)-Aiz(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)+Akx(i, j, k)**2*dvfdk(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)+ &
                            Akx(i, j, k)*Asy(i, j, k)*dufds(i, j, k)+Aky(i, j, k)**2*dvfdk(i, j, k)*fourthird+ &
                            Aky(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*onethird-Aky(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*fourthird-Aky(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird+ &
                            Akz(i, j, k)**2*dvfdk(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)

            dmudi_term(2) = Aix(i, j, k)**2*dvfdi(i, j, k)+Aix(i, j, k)*Aiy(i, j, k)*dufdi(i, j, k)*onethird+ &
                            Aix(i, j, k)*Akx(i, j, k)*dvfdk(i, j, k)+Aix(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)+Aix(i, j, k)*Asy(i, j, k)*dufds(i, j, k)+ &
                            Aiy(i, j, k)**2*dvfdi(i, j, k)*fourthird+Aiy(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)*onethird- &
                            Aiy(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*twothird+Aiy(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*fourthird- &
                            Aiy(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*twothird-Aiy(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*fourthird-Aiy(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird+ &
                            Aiz(i, j, k)**2*dvfdi(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dwfdk(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)

            dmuds_term(2) = Aix(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)-Aix(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*fourthird+ &
                            Aiy(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)-Aiz(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*twothird+Aky(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*fourthird+Aky(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)*twothird+Akz(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)+ &
                            Asx(i, j, k)**2*dvfds(i, j, k)+Asx(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*onethird+ &
                            Asy(i, j, k)**2*dvfds(i, j, k)*fourthird+Asy(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*onethird+ &
                            Asz(i, j, k)**2*dvfds(i, j, k)

            mu_term(2) = Aix(i, j, k)**2*d2vfdi2(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dvfdi(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2ufdi2(i, j, k)*onethird+Aix(i, j, k)*Aiyi(i, j, k)*dufdi(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2vfdkdi(i, j, k)+Aix(i, j, k)*Akxi(i, j, k)*dvfdk(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2ufdkdi(i, j, k)*onethird+Aix(i, j, k)*Akyi(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2vfdids(i, j, k)+Aix(i, j, k)*Asxi(i, j, k)*dvfds(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2ufdids(i, j, k)*onethird+Aix(i, j, k)*Asyi(i, j, k)*dufds(i, j, k)- &
                         Aixi(i, j, k)*Aiy(i, j, k)*dufdi(i, j, k)*twothird+Aixk(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)- &
                         Aixk(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*twothird+Aixs(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)- &
                         Aixs(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*twothird+Aiy(i, j, k)**2*d2vfdi2(i, j, k)*fourthird+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dvfdi(i, j, k)*fourthird+Aiy(i, j, k)*Aiz(i, j, k)*d2wfdi2(i, j, k)*onethird- &
                         Aiy(i, j, k)*Aizi(i, j, k)*dwfdi(i, j, k)*twothird+Aiy(i, j, k)*Akx(i, j, k)*d2ufdkdi(i, j, k)*onethird- &
                         Aiy(i, j, k)*Akxi(i, j, k)*dufdk(i, j, k)*twothird+Aiy(i, j, k)*Aky(i, j, k)*d2vfdkdi(i, j, k)*eighthird+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dvfdk(i, j, k)*fourthird+Aiy(i, j, k)*Akz(i, j, k)*d2wfdkdi(i, j, k)*onethird- &
                         Aiy(i, j, k)*Akzi(i, j, k)*dwfdk(i, j, k)*twothird+Aiy(i, j, k)*Asx(i, j, k)*d2ufdids(i, j, k)*onethird- &
                         Aiy(i, j, k)*Asxi(i, j, k)*dufds(i, j, k)*twothird+Aiy(i, j, k)*Asy(i, j, k)*d2vfdids(i, j, k)*eighthird+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dvfds(i, j, k)*fourthird+Aiy(i, j, k)*Asz(i, j, k)*d2wfdids(i, j, k)*onethird- &
                         Aiy(i, j, k)*Aszi(i, j, k)*dwfds(i, j, k)*twothird+Aiyi(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)+ &
                         Aiyk(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)+Aiyk(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*fourthird+ &
                         Aiyk(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)+Aiys(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*fourthird+Aiys(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)+ &
                         Aiz(i, j, k)**2*d2vfdi2(i, j, k)+Aiz(i, j, k)*Aizi(i, j, k)*dvfdi(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2wfdkdi(i, j, k)*onethird+Aiz(i, j, k)*Akyi(i, j, k)*dwfdk(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2vfdkdi(i, j, k)+Aiz(i, j, k)*Akzi(i, j, k)*dvfdk(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2wfdids(i, j, k)*onethird+Aiz(i, j, k)*Asyi(i, j, k)*dwfds(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2vfdids(i, j, k)+Aiz(i, j, k)*Aszi(i, j, k)*dvfds(i, j, k)- &
                         Aizk(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)*twothird+Aizk(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)- &
                         Aizs(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)*twothird+Aizs(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)+ &
                         Akx(i, j, k)**2*d2vfdk2(i, j, k)+Akx(i, j, k)*Akxk(i, j, k)*dvfdk(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2ufdk2(i, j, k)*onethird+Akx(i, j, k)*Akyk(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2vfdkds(i, j, k)+Akx(i, j, k)*Asxk(i, j, k)*dvfds(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2ufdkds(i, j, k)*onethird+Akx(i, j, k)*Asyk(i, j, k)*dufds(i, j, k)- &
                         Akxk(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*twothird+Akxs(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)- &
                         Akxs(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*twothird+Aky(i, j, k)**2*d2vfdk2(i, j, k)*fourthird+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dvfdk(i, j, k)*fourthird+Aky(i, j, k)*Akz(i, j, k)*d2wfdk2(i, j, k)*onethird- &
                         Aky(i, j, k)*Akzk(i, j, k)*dwfdk(i, j, k)*twothird+Aky(i, j, k)*Asx(i, j, k)*d2ufdkds(i, j, k)*onethird- &
                         Aky(i, j, k)*Asxk(i, j, k)*dufds(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*d2vfdkds(i, j, k)*eighthird+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dvfds(i, j, k)*fourthird+Aky(i, j, k)*Asz(i, j, k)*d2wfdkds(i, j, k)*onethird- &
                         Aky(i, j, k)*Aszk(i, j, k)*dwfds(i, j, k)*twothird+Akyk(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)+ &
                         Akys(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)+Akys(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*fourthird+ &
                         Akys(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)+Akz(i, j, k)**2*d2vfdk2(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dvfdk(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*d2wfdkds(i, j, k)*onethird+ &
                         Akz(i, j, k)*Asyk(i, j, k)*dwfds(i, j, k)+2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2vfdkds(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dvfds(i, j, k)-Akzs(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)*twothird+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)+Asx(i, j, k)**2*d2vfds2(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dvfds(i, j, k)+Asx(i, j, k)*Asy(i, j, k)*d2ufds2(i, j, k)*onethird+ &
                         Asx(i, j, k)*Asys(i, j, k)*dufds(i, j, k)-Asxs(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*twothird+ &
                         Asy(i, j, k)**2*d2vfds2(i, j, k)*fourthird+Asy(i, j, k)*Asys(i, j, k)*dvfds(i, j, k)*fourthird+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2wfds2(i, j, k)*onethird-Asy(i, j, k)*Aszs(i, j, k)*dwfds(i, j, k)*twothird+ &
                         Asys(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)+Asz(i, j, k)**2*d2vfds2(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dvfds(i, j, k)

            dmufdk_term(3) = Aix(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)-Aix(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Akx(i, j, k)*dudi(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)+ &
                             Aiz(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird+Akx(i, j, k)**2*dwdk(i, j, k)+ &
                             Akx(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+ &
                             Akx(i, j, k)*Asz(i, j, k)*duds(i, j, k)+Aky(i, j, k)**2*dwdk(i, j, k)+ &
                             Aky(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*onethird+Aky(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+ &
                             Aky(i, j, k)*Asz(i, j, k)*dvds(i, j, k)+Akz(i, j, k)**2*dwdk(i, j, k)*fourthird- &
                             Akz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird-Akz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird+ &
                             Akz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird

            dmufdi_term(3) = Aix(i, j, k)**2*dwdi(i, j, k)+Aix(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*onethird+ &
                             Aix(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)+Aix(i, j, k)*Akz(i, j, k)*dudk(i, j, k)+ &
                             Aix(i, j, k)*Asx(i, j, k)*dwds(i, j, k)+Aix(i, j, k)*Asz(i, j, k)*duds(i, j, k)+ &
                             Aiy(i, j, k)**2*dwdi(i, j, k)+Aiy(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*onethird+ &
                             Aiy(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)+ &
                             Aiy(i, j, k)*Asy(i, j, k)*dwds(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)+ &
                             Aiz(i, j, k)**2*dwdi(i, j, k)*fourthird-Aiz(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird- &
                             Aiz(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird+Aiz(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*fourthird- &
                             Aiz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird-Aiz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird

            dmufds_term(3) = Aix(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)-Aix(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird+ &
                             Aiy(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)-Aiy(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird+ &
                             Aiz(i, j, k)*Asx(i, j, k)*dudi(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)+ &
                             Aiz(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird+Akx(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)- &
                             Akx(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)- &
                             Aky(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird+Akz(i, j, k)*Asx(i, j, k)*dudk(i, j, k)+ &
                             Akz(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird+ &
                             Asx(i, j, k)**2*dwds(i, j, k)+Asx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*onethird+ &
                             Asy(i, j, k)**2*dwds(i, j, k)+Asy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*onethird+ &
                             Asz(i, j, k)**2*dwds(i, j, k)*fourthird

            muf_term(3) = Aix(i, j, k)**2*d2wdi2(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dwdi(i, j, k)+ &
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

            dmudk_term(3) = Aix(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)-Aix(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)-Aiy(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*fourthird+Akx(i, j, k)**2*dwfdk(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*onethird+Akx(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)+ &
                            Akx(i, j, k)*Asz(i, j, k)*dufds(i, j, k)+Aky(i, j, k)**2*dwfdk(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*onethird+Aky(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)+ &
                            Aky(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)+Akz(i, j, k)**2*dwfdk(i, j, k)*fourthird- &
                            Akz(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird-Akz(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*fourthird

            dmudi_term(3) = Aix(i, j, k)**2*dwfdi(i, j, k)+Aix(i, j, k)*Aiz(i, j, k)*dufdi(i, j, k)*onethird+ &
                            Aix(i, j, k)*Akx(i, j, k)*dwfdk(i, j, k)+Aix(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)+Aix(i, j, k)*Asz(i, j, k)*dufds(i, j, k)+ &
                            Aiy(i, j, k)**2*dwfdi(i, j, k)+Aiy(i, j, k)*Aiz(i, j, k)*dvfdi(i, j, k)*onethird+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dwfdk(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)+ &
                            Aiz(i, j, k)**2*dwfdi(i, j, k)*fourthird-Aiz(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*twothird- &
                            Aiz(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*twothird+Aiz(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*fourthird- &
                            Aiz(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird-Aiz(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*fourthird

            dmuds_term(3) = Aix(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)-Aix(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*twothird+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)-Aiy(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*twothird+ &
                            Aiz(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*fourthird+Akx(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*twothird+Aky(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*twothird+Akz(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)+ &
                            Akz(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*fourthird+ &
                            Asx(i, j, k)**2*dwfds(i, j, k)+Asx(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*onethird+ &
                            Asy(i, j, k)**2*dwfds(i, j, k)+Asy(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*onethird+ &
                            Asz(i, j, k)**2*dwfds(i, j, k)*fourthird

            mu_term(3) = Aix(i, j, k)**2*d2wfdi2(i, j, k)+Aix(i, j, k)*Aixi(i, j, k)*dwfdi(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2ufdi2(i, j, k)*onethird+Aix(i, j, k)*Aizi(i, j, k)*dufdi(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2wfdkdi(i, j, k)+Aix(i, j, k)*Akxi(i, j, k)*dwfdk(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2ufdkdi(i, j, k)*onethird+Aix(i, j, k)*Akzi(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2wfdids(i, j, k)+Aix(i, j, k)*Asxi(i, j, k)*dwfds(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2ufdids(i, j, k)*onethird+Aix(i, j, k)*Aszi(i, j, k)*dufds(i, j, k)- &
                         Aixi(i, j, k)*Aiz(i, j, k)*dufdi(i, j, k)*twothird+Aixk(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)- &
                         Aixk(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*twothird+Aixs(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)- &
                         Aixs(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*twothird+Aiy(i, j, k)**2*d2wfdi2(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dwfdi(i, j, k)+Aiy(i, j, k)*Aiz(i, j, k)*d2vfdi2(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Aizi(i, j, k)*dvfdi(i, j, k)+2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2wfdkdi(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dwfdk(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*d2vfdkdi(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Akzi(i, j, k)*dvfdk(i, j, k)+2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2wfdids(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dwfds(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*d2vfdids(i, j, k)*onethird+ &
                         Aiy(i, j, k)*Aszi(i, j, k)*dvfds(i, j, k)-Aiyi(i, j, k)*Aiz(i, j, k)*dvfdi(i, j, k)*twothird+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)-Aiyk(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*twothird+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)-Aiys(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*twothird+ &
                         Aiz(i, j, k)**2*d2wfdi2(i, j, k)*fourthird+Aiz(i, j, k)*Aizi(i, j, k)*dwfdi(i, j, k)*fourthird+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2ufdkdi(i, j, k)*onethird-Aiz(i, j, k)*Akxi(i, j, k)*dufdk(i, j, k)*twothird+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2vfdkdi(i, j, k)*onethird-Aiz(i, j, k)*Akyi(i, j, k)*dvfdk(i, j, k)*twothird+ &
                         Aiz(i, j, k)*Akz(i, j, k)*d2wfdkdi(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dwfdk(i, j, k)*fourthird+Aiz(i, j, k)*Asx(i, j, k)*d2ufdids(i, j, k)*onethird- &
                         Aiz(i, j, k)*Asxi(i, j, k)*dufds(i, j, k)*twothird+Aiz(i, j, k)*Asy(i, j, k)*d2vfdids(i, j, k)*onethird- &
                         Aiz(i, j, k)*Asyi(i, j, k)*dvfds(i, j, k)*twothird+Aiz(i, j, k)*Asz(i, j, k)*d2wfdids(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dwfds(i, j, k)*fourthird+Aizk(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)+ &
                         Aizk(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)+Aizk(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*fourthird+ &
                         Aizs(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)+Aizs(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*fourthird+Akx(i, j, k)**2*d2wfdk2(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dwfdk(i, j, k)+Akx(i, j, k)*Akz(i, j, k)*d2ufdk2(i, j, k)*onethird+ &
                         Akx(i, j, k)*Akzk(i, j, k)*dufdk(i, j, k)+2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2wfdkds(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dwfds(i, j, k)+Akx(i, j, k)*Asz(i, j, k)*d2ufdkds(i, j, k)*onethird+ &
                         Akx(i, j, k)*Aszk(i, j, k)*dufds(i, j, k)-Akxk(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*twothird+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)-Akxs(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*twothird+ &
                         Aky(i, j, k)**2*d2wfdk2(i, j, k)+Aky(i, j, k)*Akyk(i, j, k)*dwfdk(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2vfdk2(i, j, k)*onethird+Aky(i, j, k)*Akzk(i, j, k)*dvfdk(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2wfdkds(i, j, k)+Aky(i, j, k)*Asyk(i, j, k)*dwfds(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2vfdkds(i, j, k)*onethird+Aky(i, j, k)*Aszk(i, j, k)*dvfds(i, j, k)- &
                         Akyk(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*twothird+Akys(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)- &
                         Akys(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*twothird+Akz(i, j, k)**2*d2wfdk2(i, j, k)*fourthird+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dwfdk(i, j, k)*fourthird+Akz(i, j, k)*Asx(i, j, k)*d2ufdkds(i, j, k)*onethird- &
                         Akz(i, j, k)*Asxk(i, j, k)*dufds(i, j, k)*twothird+Akz(i, j, k)*Asy(i, j, k)*d2vfdkds(i, j, k)*onethird- &
                         Akz(i, j, k)*Asyk(i, j, k)*dvfds(i, j, k)*twothird+Akz(i, j, k)*Asz(i, j, k)*d2wfdkds(i, j, k)*eighthird+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dwfds(i, j, k)*fourthird+Akzs(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)+ &
                         Akzs(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)+Akzs(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*fourthird+ &
                         Asx(i, j, k)**2*d2wfds2(i, j, k)+Asx(i, j, k)*Asxs(i, j, k)*dwfds(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2ufds2(i, j, k)*onethird+Asx(i, j, k)*Aszs(i, j, k)*dufds(i, j, k)- &
                         Asxs(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*twothird+Asy(i, j, k)**2*d2wfds2(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dwfds(i, j, k)+Asy(i, j, k)*Asz(i, j, k)*d2vfds2(i, j, k)*onethird+ &
                         Asy(i, j, k)*Aszs(i, j, k)*dvfds(i, j, k)-Asys(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*twothird+ &
                         Asz(i, j, k)**2*d2wfds2(i, j, k)*fourthird+Asz(i, j, k)*Aszs(i, j, k)*dwfds(i, j, k)*fourthird

            dmufdk_term(4) = Aix(i, j, k)*Akx(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
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

            dmufdi_term(4) = Aix(i, j, k)**2*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
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

            dmufds_term(4) = Aix(i, j, k)*Asx(i, j, k)*Para%Cp*dTdi(i, j, k)/Para%Pr+ &
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

            muf_term(4) = Aix(i, j, k)**2*Para%Cp*d2Tdi2(i, j, k)/Para%Pr+ &
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

            dmudk_term(4) = Aix(i, j, k)*Akx(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*fourthird*u(i, j, k)+ &
                      Aix(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+Aix(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                      Aix(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+Aix(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aix(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*twothird*v(i, j, k)+ &
                      Aix(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*u_fluc(i, j, k)+Aix(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*u(i, j, k)- &
                            Aix(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aix(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*twothird*w(i, j, k)+ &
                      Aix(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*u_fluc(i, j, k)+Aix(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*u(i, j, k)+ &
                      Aiy(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*v_fluc(i, j, k)+Aiy(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*v(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)*twothird*u(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                      Aiy(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*fourthird*v(i, j, k)+ &
                      Aiy(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aiy(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*twothird*w(i, j, k)+ &
                      Aiy(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*v_fluc(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*v(i, j, k)+ &
                      Aiz(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*w_fluc(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)*twothird*u(i, j, k)+ &
                      Aiz(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*w_fluc(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aiz(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)*twothird*v(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                      Aiz(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                      Aiz(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*fourthird*w(i, j, k)+ &
                           Akx(i, j, k)**2*Para%Cp*dTfdk(i, j, k)/Para%Pr+Akx(i, j, k)**2*dudk(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Akx(i, j, k)**2*dufdk(i, j, k)*fourthird*u(i, j, k)+Akx(i, j, k)**2*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                            Akx(i, j, k)**2*dvfdk(i, j, k)*v(i, j, k)+Akx(i, j, k)**2*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                       Akx(i, j, k)**2*dwfdk(i, j, k)*w(i, j, k)+Akx(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*onethird*v_fluc(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*onethird*v(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*onethird*u_fluc(i, j, k)+ &
                            Akx(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*onethird*u(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*onethird*w_fluc(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*onethird*w(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird*u_fluc(i, j, k)+ &
                            Akx(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*onethird*u(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                            Akx(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*fourthird*u(i, j, k)+ &
                      Akx(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                      Akx(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                      Akx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*v_fluc(i, j, k)+Akx(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*v(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird*u(i, j, k)+ &
                      Akx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*w_fluc(i, j, k)+Akx(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*w(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird*u(i, j, k)+ &
                            Aky(i, j, k)**2*Para%Cp*dTfdk(i, j, k)/Para%Pr+Aky(i, j, k)**2*dudk(i, j, k)*u_fluc(i, j, k)+ &
                            Aky(i, j, k)**2*dufdk(i, j, k)*u(i, j, k)+Aky(i, j, k)**2*dvdk(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aky(i, j, k)**2*dvfdk(i, j, k)*fourthird*v(i, j, k)+Aky(i, j, k)**2*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                       Aky(i, j, k)**2*dwfdk(i, j, k)*w(i, j, k)+Aky(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*onethird*w_fluc(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*onethird*w(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*onethird*v_fluc(i, j, k)+ &
                            Aky(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*onethird*v(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird*v(i, j, k)+ &
                      Aky(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*u_fluc(i, j, k)+Aky(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)*u(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                      Aky(i, j, k)*Asy(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*fourthird*v(i, j, k)+ &
                      Aky(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                      Aky(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*w_fluc(i, j, k)+Aky(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*w(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird*v(i, j, k)+ &
                            Akz(i, j, k)**2*Para%Cp*dTfdk(i, j, k)/Para%Pr+Akz(i, j, k)**2*dudk(i, j, k)*u_fluc(i, j, k)+ &
                            Akz(i, j, k)**2*dufdk(i, j, k)*u(i, j, k)+Akz(i, j, k)**2*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                            Akz(i, j, k)**2*dvfdk(i, j, k)*v(i, j, k)+Akz(i, j, k)**2*dwdk(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Akz(i, j, k)**2*dwfdk(i, j, k)*fourthird*w(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird*w(i, j, k)+ &
                      Akz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*u_fluc(i, j, k)+Akz(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)*u(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird*w(i, j, k)+ &
                      Akz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*v_fluc(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)*v(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                      Akz(i, j, k)*Asz(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                      Akz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*fourthird*w(i, j, k)

            dmudi_term(4) = Aix(i, j, k)**2*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)**2*dudi(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)**2*dufdi(i, j, k)*fourthird*u(i, j, k)+Aix(i, j, k)**2*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                            Aix(i, j, k)**2*dvfdi(i, j, k)*v(i, j, k)+Aix(i, j, k)**2*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                       Aix(i, j, k)**2*dwfdi(i, j, k)*w(i, j, k)+Aix(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*onethird*v_fluc(i, j, k)+ &
                            Aix(i, j, k)*Aiy(i, j, k)*dufdi(i, j, k)*onethird*v(i, j, k)+ &
                            Aix(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)*onethird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)*Aiy(i, j, k)*dvfdi(i, j, k)*onethird*u(i, j, k)+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*onethird*w_fluc(i, j, k)+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dufdi(i, j, k)*onethird*w(i, j, k)+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)*onethird*u(i, j, k)+ &
                            Aix(i, j, k)*Akx(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*fourthird*u(i, j, k)+ &
                      Aix(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+Aix(i, j, k)*Akx(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                      Aix(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+Aix(i, j, k)*Akx(i, j, k)*dwfdk(i, j, k)*w(i, j, k)+ &
                      Aix(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*v_fluc(i, j, k)+Aix(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*v(i, j, k)- &
                            Aix(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aix(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*twothird*u(i, j, k)+ &
                      Aix(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*w_fluc(i, j, k)+Aix(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aix(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*twothird*u(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Asx(i, j, k)*duds(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*fourthird*u(i, j, k)+ &
                      Aix(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+Aix(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                      Aix(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+Aix(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                      Aix(i, j, k)*Asy(i, j, k)*duds(i, j, k)*v_fluc(i, j, k)+Aix(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*v(i, j, k)- &
                            Aix(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aix(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird*u(i, j, k)+ &
                      Aix(i, j, k)*Asz(i, j, k)*duds(i, j, k)*w_fluc(i, j, k)+Aix(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aix(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird*u(i, j, k)+ &
                            Aiy(i, j, k)**2*Para%Cp*dTfdi(i, j, k)/Para%Pr+Aiy(i, j, k)**2*dudi(i, j, k)*u_fluc(i, j, k)+ &
                            Aiy(i, j, k)**2*dufdi(i, j, k)*u(i, j, k)+Aiy(i, j, k)**2*dvdi(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aiy(i, j, k)**2*dvfdi(i, j, k)*fourthird*v(i, j, k)+Aiy(i, j, k)**2*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                       Aiy(i, j, k)**2*dwfdi(i, j, k)*w(i, j, k)+Aiy(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*onethird*w_fluc(i, j, k)+ &
                            Aiy(i, j, k)*Aiz(i, j, k)*dvfdi(i, j, k)*onethird*w(i, j, k)+ &
                            Aiy(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*onethird*v_fluc(i, j, k)+ &
                            Aiy(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)*onethird*v(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aiy(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*twothird*v(i, j, k)+ &
                      Aiy(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)*Akx(i, j, k)*dvfdk(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                      Aiy(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aiy(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*fourthird*v(i, j, k)+ &
                      Aiy(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+Aiy(i, j, k)*Aky(i, j, k)*dwfdk(i, j, k)*w(i, j, k)+ &
                      Aiy(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*w_fluc(i, j, k)+Aiy(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aiy(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*twothird*v(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird*v(i, j, k)+ &
                      Aiy(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)*Asx(i, j, k)*dvfds(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                      Aiy(i, j, k)*Asy(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*fourthird*v(i, j, k)+ &
                      Aiy(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                      Aiy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*w_fluc(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aiy(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*twothird*v(i, j, k)+ &
                            Aiz(i, j, k)**2*Para%Cp*dTfdi(i, j, k)/Para%Pr+Aiz(i, j, k)**2*dudi(i, j, k)*u_fluc(i, j, k)+ &
                            Aiz(i, j, k)**2*dufdi(i, j, k)*u(i, j, k)+Aiz(i, j, k)**2*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                            Aiz(i, j, k)**2*dvfdi(i, j, k)*v(i, j, k)+Aiz(i, j, k)**2*dwdi(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Aiz(i, j, k)**2*dwfdi(i, j, k)*fourthird*w(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aiz(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*twothird*w(i, j, k)+ &
                      Aiz(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)*Akx(i, j, k)*dwfdk(i, j, k)*u(i, j, k)- &
                            Aiz(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aiz(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*twothird*w(i, j, k)+ &
                      Aiz(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)*Aky(i, j, k)*dwfdk(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                      Aiz(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                      Aiz(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Aiz(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*fourthird*w(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*twothird*w(i, j, k)+ &
                      Aiz(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*dwfds(i, j, k)*u(i, j, k)- &
                            Aiz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aiz(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*twothird*w(i, j, k)+ &
                      Aiz(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dwfds(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                      Aiz(i, j, k)*Asz(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                      Aiz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*fourthird*w(i, j, k)

            dmuds_term(4) = Aix(i, j, k)*Asx(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                            Aix(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Aix(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*fourthird*u(i, j, k)+ &
                      Aix(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+Aix(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                      Aix(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+Aix(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                            Aix(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aix(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*twothird*v(i, j, k)+ &
                      Aix(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*u_fluc(i, j, k)+Aix(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*u(i, j, k)- &
                            Aix(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aix(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*twothird*w(i, j, k)+ &
                      Aix(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*u_fluc(i, j, k)+Aix(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*u(i, j, k)+ &
                      Aiy(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*v_fluc(i, j, k)+Aiy(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*v(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aiy(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)*twothird*u(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                      Aiy(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aiy(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*fourthird*v(i, j, k)+ &
                      Aiy(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+Aiy(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                            Aiy(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aiy(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*twothird*w(i, j, k)+ &
                      Aiy(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*v_fluc(i, j, k)+Aiy(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*v(i, j, k)+ &
                      Aiz(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*w_fluc(i, j, k)+Aiz(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aiz(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)*twothird*u(i, j, k)+ &
                      Aiz(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*w_fluc(i, j, k)+Aiz(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*w(i, j, k)- &
                            Aiz(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Aiz(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)*twothird*v(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                      Aiz(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                      Aiz(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Aiz(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*fourthird*w(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                            Akx(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Akx(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*fourthird*u(i, j, k)+ &
                      Akx(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                      Akx(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+Akx(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)*w(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Akx(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*twothird*v(i, j, k)+ &
                      Akx(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*u_fluc(i, j, k)+Akx(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*u(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Akx(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*twothird*w(i, j, k)+ &
                      Akx(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*u_fluc(i, j, k)+Akx(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*u(i, j, k)+ &
                      Aky(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*v_fluc(i, j, k)+Aky(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*v(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Aky(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)*twothird*u(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                      Aky(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Aky(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*fourthird*v(i, j, k)+ &
                      Aky(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+Aky(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)*w(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird*w_fluc(i, j, k)- &
                            Aky(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*twothird*w(i, j, k)+ &
                      Aky(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*v_fluc(i, j, k)+Aky(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*v(i, j, k)+ &
                      Akz(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*w_fluc(i, j, k)+Akz(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*w(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                            Akz(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)*twothird*u(i, j, k)+ &
                      Akz(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*w_fluc(i, j, k)+Akz(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*w(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird*v_fluc(i, j, k)- &
                            Akz(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)*twothird*v(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                      Akz(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                      Akz(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+Akz(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Akz(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*fourthird*w(i, j, k)+ &
                           Asx(i, j, k)**2*Para%Cp*dTfds(i, j, k)/Para%Pr+Asx(i, j, k)**2*duds(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                            Asx(i, j, k)**2*dufds(i, j, k)*fourthird*u(i, j, k)+Asx(i, j, k)**2*dvds(i, j, k)*v_fluc(i, j, k)+ &
                            Asx(i, j, k)**2*dvfds(i, j, k)*v(i, j, k)+Asx(i, j, k)**2*dwds(i, j, k)*w_fluc(i, j, k)+ &
                       Asx(i, j, k)**2*dwfds(i, j, k)*w(i, j, k)+Asx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*onethird*v_fluc(i, j, k)+ &
                            Asx(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*onethird*v(i, j, k)+ &
                            Asx(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*onethird*u_fluc(i, j, k)+ &
                            Asx(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*onethird*u(i, j, k)+ &
                            Asx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*onethird*w_fluc(i, j, k)+ &
                            Asx(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*onethird*w(i, j, k)+ &
                            Asx(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird*u_fluc(i, j, k)+ &
                            Asx(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*onethird*u(i, j, k)+ &
                            Asy(i, j, k)**2*Para%Cp*dTfds(i, j, k)/Para%Pr+Asy(i, j, k)**2*duds(i, j, k)*u_fluc(i, j, k)+ &
                            Asy(i, j, k)**2*dufds(i, j, k)*u(i, j, k)+Asy(i, j, k)**2*dvds(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                            Asy(i, j, k)**2*dvfds(i, j, k)*fourthird*v(i, j, k)+Asy(i, j, k)**2*dwds(i, j, k)*w_fluc(i, j, k)+ &
                       Asy(i, j, k)**2*dwfds(i, j, k)*w(i, j, k)+Asy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*onethird*w_fluc(i, j, k)+ &
                            Asy(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*onethird*w(i, j, k)+ &
                            Asy(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*onethird*v_fluc(i, j, k)+ &
                            Asy(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*onethird*v(i, j, k)+ &
                            Asz(i, j, k)**2*Para%Cp*dTfds(i, j, k)/Para%Pr+Asz(i, j, k)**2*duds(i, j, k)*u_fluc(i, j, k)+ &
                            Asz(i, j, k)**2*dufds(i, j, k)*u(i, j, k)+Asz(i, j, k)**2*dvds(i, j, k)*v_fluc(i, j, k)+ &
                            Asz(i, j, k)**2*dvfds(i, j, k)*v(i, j, k)+Asz(i, j, k)**2*dwds(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                            Asz(i, j, k)**2*dwfds(i, j, k)*fourthird*w(i, j, k)

            mu_term(4) = Aix(i, j, k)**2*Para%Cp*d2Tfdi2(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)**2*d2udi2(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)**2*d2ufdi2(i, j, k)*fourthird*u(i, j, k)+Aix(i, j, k)**2*d2vdi2(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)**2*d2vfdi2(i, j, k)*v(i, j, k)+Aix(i, j, k)**2*d2wdi2(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)**2*d2wfdi2(i, j, k)*w(i, j, k)+Aix(i, j, k)**2*dudi(i, j, k)*dufdi(i, j, k)*eighthird+ &
                         2.0d0*Aix(i, j, k)**2*dvdi(i, j, k)*dvfdi(i, j, k)+2.0d0*Aix(i, j, k)**2*dwdi(i, j, k)*dwfdi(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dudi(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dufdi(i, j, k)*fourthird*u(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aixi(i, j, k)*dwfdi(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2udi2(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2ufdi2(i, j, k)*onethird*v(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2vdi2(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*d2vfdi2(i, j, k)*onethird*u(i, j, k)+ &
                         Aix(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*dvfdi(i, j, k)*twothird+ &
                         Aix(i, j, k)*Aiy(i, j, k)*dufdi(i, j, k)*dvdi(i, j, k)*twothird+ &
                         Aix(i, j, k)*Aiyi(i, j, k)*dudi(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aiyi(i, j, k)*dufdi(i, j, k)*v(i, j, k)- &
                         Aix(i, j, k)*Aiyi(i, j, k)*dvdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aix(i, j, k)*Aiyi(i, j, k)*dvfdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2udi2(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2ufdi2(i, j, k)*onethird*w(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2wdi2(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*d2wfdi2(i, j, k)*onethird*u(i, j, k)+ &
                         Aix(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*dwfdi(i, j, k)*twothird+ &
                         Aix(i, j, k)*Aiz(i, j, k)*dufdi(i, j, k)*dwdi(i, j, k)*twothird+ &
                         Aix(i, j, k)*Aizi(i, j, k)*dudi(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aizi(i, j, k)*dufdi(i, j, k)*w(i, j, k)- &
                         Aix(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aix(i, j, k)*Aizi(i, j, k)*dwfdi(i, j, k)*twothird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*Para%Cp*d2Tfdkdi(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*eighthird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akx(i, j, k)*d2ufdkdi(i, j, k)*eighthird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2vdkdi(i, j, k)*v_fluc(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2vfdkdi(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2wdkdi(i, j, k)*w_fluc(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*d2wfdkdi(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*dufdk(i, j, k)*eighthird+ &
                         Aix(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*dufdi(i, j, k)*eighthird+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*dvfdk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*dvdk(i, j, k)*dvfdi(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*dwfdk(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Akx(i, j, k)*dwdk(i, j, k)*dwfdi(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dufdk(i, j, k)*fourthird*u(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akxi(i, j, k)*dwfdk(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2udkdi(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2ufdkdi(i, j, k)*onethird*v(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aky(i, j, k)*d2vfdkdi(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*dvfdk(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*dvfdi(i, j, k)- &
                         Aix(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*dvdk(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*dvdi(i, j, k)+ &
                         Aix(i, j, k)*Akyi(i, j, k)*dudk(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akyi(i, j, k)*dufdk(i, j, k)*v(i, j, k)- &
                         Aix(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aix(i, j, k)*Akyi(i, j, k)*dvfdk(i, j, k)*twothird*u(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2udkdi(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2ufdkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akz(i, j, k)*d2wfdkdi(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*dwfdk(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*dwfdi(i, j, k)- &
                         Aix(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*dwdk(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*dwdi(i, j, k)+ &
                         Aix(i, j, k)*Akzi(i, j, k)*dudk(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Akzi(i, j, k)*dufdk(i, j, k)*w(i, j, k)- &
                         Aix(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aix(i, j, k)*Akzi(i, j, k)*dwfdk(i, j, k)*twothird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*Para%Cp*d2Tfdids(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*eighthird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asx(i, j, k)*d2ufdids(i, j, k)*eighthird*u(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2vdids(i, j, k)*v_fluc(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2vfdids(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2wdids(i, j, k)*w_fluc(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*d2wfdids(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*dufds(i, j, k)*eighthird+ &
                         Aix(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dufdi(i, j, k)*eighthird+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*dvfds(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*dvfdi(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*dwfds(i, j, k)+ &
                         2.0d0*Aix(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*dwfdi(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Aix(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dufds(i, j, k)*fourthird*u(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asxi(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2udids(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2ufdids(i, j, k)*onethird*v(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asy(i, j, k)*d2vfdids(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*dvfds(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dvfdi(i, j, k)- &
                         Aix(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*dvds(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*dvdi(i, j, k)+ &
                         Aix(i, j, k)*Asyi(i, j, k)*duds(i, j, k)*v_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asyi(i, j, k)*dufds(i, j, k)*v(i, j, k)- &
                         Aix(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aix(i, j, k)*Asyi(i, j, k)*dvfds(i, j, k)*twothird*u(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2udids(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2ufdids(i, j, k)*onethird*w(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aix(i, j, k)*Asz(i, j, k)*d2wfdids(i, j, k)*onethird*u(i, j, k)- &
                         Aix(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*dwfds(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dwfdi(i, j, k)- &
                         Aix(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Aix(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*dwdi(i, j, k)+ &
                         Aix(i, j, k)*Aszi(i, j, k)*duds(i, j, k)*w_fluc(i, j, k)+ &
                         Aix(i, j, k)*Aszi(i, j, k)*dufds(i, j, k)*w(i, j, k)- &
                         Aix(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aix(i, j, k)*Aszi(i, j, k)*dwfds(i, j, k)*twothird*u(i, j, k)- &
                         Aixi(i, j, k)*Aiy(i, j, k)*dudi(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aixi(i, j, k)*Aiy(i, j, k)*dufdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aixi(i, j, k)*Aiy(i, j, k)*dvdi(i, j, k)*u_fluc(i, j, k)+ &
                         Aixi(i, j, k)*Aiy(i, j, k)*dvfdi(i, j, k)*u(i, j, k)- &
                         Aixi(i, j, k)*Aiz(i, j, k)*dudi(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aixi(i, j, k)*Aiz(i, j, k)*dufdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aixi(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*u_fluc(i, j, k)+ &
                         Aixi(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)*u(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*fourthird*u(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aixk(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                         Aixk(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aixk(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aixk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*u_fluc(i, j, k)+ &
                         Aixk(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*u(i, j, k)- &
                         Aixk(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aixk(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aixk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*u_fluc(i, j, k)+ &
                         Aixk(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*u(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*fourthird*u(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aixs(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                         Aixs(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aixs(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aixs(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*u_fluc(i, j, k)+ &
                         Aixs(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*u(i, j, k)- &
                         Aixs(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aixs(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aixs(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*u_fluc(i, j, k)+ &
                         Aixs(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*u(i, j, k)+Aiy(i, j, k)**2*Para%Cp*d2Tfdi2(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)**2*d2udi2(i, j, k)*u_fluc(i, j, k)+Aiy(i, j, k)**2*d2ufdi2(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)**2*d2vdi2(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)**2*d2vfdi2(i, j, k)*fourthird*v(i, j, k)+Aiy(i, j, k)**2*d2wdi2(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)**2*d2wfdi2(i, j, k)*w(i, j, k)+2.0d0*Aiy(i, j, k)**2*dudi(i, j, k)*dufdi(i, j, k)+ &
                         Aiy(i, j, k)**2*dvdi(i, j, k)*dvfdi(i, j, k)*eighthird+ &
                         2.0d0*Aiy(i, j, k)**2*dwdi(i, j, k)*dwfdi(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dvdi(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dvfdi(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aiyi(i, j, k)*dwfdi(i, j, k)*w(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*d2vdi2(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*d2vfdi2(i, j, k)*onethird*w(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*d2wdi2(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*d2wfdi2(i, j, k)*onethird*v(i, j, k)+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*dwfdi(i, j, k)*twothird+ &
                         Aiy(i, j, k)*Aiz(i, j, k)*dvfdi(i, j, k)*dwdi(i, j, k)*twothird+ &
                         Aiy(i, j, k)*Aizi(i, j, k)*dvdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aizi(i, j, k)*dvfdi(i, j, k)*w(i, j, k)- &
                         Aiy(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aiy(i, j, k)*Aizi(i, j, k)*dwfdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akx(i, j, k)*d2ufdkdi(i, j, k)*onethird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akx(i, j, k)*d2vdkdi(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akx(i, j, k)*d2vfdkdi(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*dvfdk(i, j, k)- &
                         Aiy(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*dvfdi(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*dvdk(i, j, k)- &
                         Aiy(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*dvdi(i, j, k)*fourthird- &
                         Aiy(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aiy(i, j, k)*Akxi(i, j, k)*dufdk(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akxi(i, j, k)*dvdk(i, j, k)*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akxi(i, j, k)*dvfdk(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*Para%Cp*d2Tfdkdi(i, j, k)/Para%Pr+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2udkdi(i, j, k)*u_fluc(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2ufdkdi(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*eighthird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aky(i, j, k)*d2vfdkdi(i, j, k)*eighthird*v(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2wdkdi(i, j, k)*w_fluc(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*d2wfdkdi(i, j, k)*w(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*dufdi(i, j, k)+ &
                         Aiy(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*dvfdk(i, j, k)*eighthird+ &
                         Aiy(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*dvfdi(i, j, k)*eighthird+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*dwfdk(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Aky(i, j, k)*dwdk(i, j, k)*dwfdi(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dvfdk(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akyi(i, j, k)*dwfdk(i, j, k)*w(i, j, k)+ &
                         Aiy(i, j, k)*Akz(i, j, k)*d2vdkdi(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akz(i, j, k)*d2vfdkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aiy(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akz(i, j, k)*d2wfdkdi(i, j, k)*onethird*v(i, j, k)- &
                         Aiy(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*dwfdk(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*dwfdi(i, j, k)- &
                         Aiy(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*dwdk(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*dwdi(i, j, k)+ &
                         Aiy(i, j, k)*Akzi(i, j, k)*dvdk(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Akzi(i, j, k)*dvfdk(i, j, k)*w(i, j, k)- &
                         Aiy(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aiy(i, j, k)*Akzi(i, j, k)*dwfdk(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asx(i, j, k)*d2ufdids(i, j, k)*onethird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asx(i, j, k)*d2vdids(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asx(i, j, k)*d2vfdids(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*dvfds(i, j, k)- &
                         Aiy(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dvfdi(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*dvds(i, j, k)- &
                         Aiy(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*dvdi(i, j, k)*fourthird- &
                         Aiy(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aiy(i, j, k)*Asxi(i, j, k)*dufds(i, j, k)*twothird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asxi(i, j, k)*dvds(i, j, k)*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asxi(i, j, k)*dvfds(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*Para%Cp*d2Tfdids(i, j, k)/Para%Pr+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2udids(i, j, k)*u_fluc(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2ufdids(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*eighthird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asy(i, j, k)*d2vfdids(i, j, k)*eighthird*v(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2wdids(i, j, k)*w_fluc(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*d2wfdids(i, j, k)*w(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*dufds(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dufdi(i, j, k)+ &
                         Aiy(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*dvfds(i, j, k)*eighthird+ &
                         Aiy(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*dvfdi(i, j, k)*eighthird+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*dwfds(i, j, k)+ &
                         2.0d0*Aiy(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*dwfdi(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dvfds(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asyi(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                         Aiy(i, j, k)*Asz(i, j, k)*d2vdids(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asz(i, j, k)*d2vfdids(i, j, k)*onethird*w(i, j, k)+ &
                         Aiy(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Asz(i, j, k)*d2wfdids(i, j, k)*onethird*v(i, j, k)- &
                         Aiy(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*dwfds(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dwfdi(i, j, k)- &
                         Aiy(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Aiy(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*dwdi(i, j, k)+ &
                         Aiy(i, j, k)*Aszi(i, j, k)*dvds(i, j, k)*w_fluc(i, j, k)+ &
                         Aiy(i, j, k)*Aszi(i, j, k)*dvfds(i, j, k)*w(i, j, k)- &
                         Aiy(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aiy(i, j, k)*Aszi(i, j, k)*dwfds(i, j, k)*twothird*v(i, j, k)- &
                         Aiyi(i, j, k)*Aiz(i, j, k)*dvdi(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiyi(i, j, k)*Aiz(i, j, k)*dvfdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aiyi(i, j, k)*Aiz(i, j, k)*dwdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aiyi(i, j, k)*Aiz(i, j, k)*dwfdi(i, j, k)*v(i, j, k)+ &
                         Aiyk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*v_fluc(i, j, k)+ &
                         Aiyk(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*v(i, j, k)- &
                         Aiyk(i, j, k)*Akx(i, j, k)*dvdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aiyk(i, j, k)*Akx(i, j, k)*dvfdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aiyk(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                         Aiyk(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiyk(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aiyk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aiyk(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*v(i, j, k)+ &
                         Aiys(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*v_fluc(i, j, k)+ &
                         Aiys(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*v(i, j, k)- &
                         Aiys(i, j, k)*Asx(i, j, k)*dvdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aiys(i, j, k)*Asx(i, j, k)*dvfdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*fourthird*v(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aiys(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)*w(i, j, k)- &
                         Aiys(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiys(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*twothird*w(i, j, k)+ &
                         Aiys(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aiys(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*v(i, j, k)+Aiz(i, j, k)**2*Para%Cp*d2Tfdi2(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)**2*d2udi2(i, j, k)*u_fluc(i, j, k)+Aiz(i, j, k)**2*d2ufdi2(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)**2*d2vdi2(i, j, k)*v_fluc(i, j, k)+Aiz(i, j, k)**2*d2vfdi2(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)**2*d2wdi2(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)**2*d2wfdi2(i, j, k)*fourthird*w(i, j, k)+2.0d0*Aiz(i, j, k)**2*dudi(i, j, k)*dufdi(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)**2*dvdi(i, j, k)*dvfdi(i, j, k)+ &
                         Aiz(i, j, k)**2*dwdi(i, j, k)*dwfdi(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dwdi(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aizi(i, j, k)*dwfdi(i, j, k)*fourthird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2udkdi(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2ufdkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2wdkdi(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akx(i, j, k)*d2wfdkdi(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*dwfdk(i, j, k)- &
                         Aiz(i, j, k)*Akx(i, j, k)*dudk(i, j, k)*dwfdi(i, j, k)*fourthird+ &
                         2.0d0*Aiz(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*dwdk(i, j, k)- &
                         Aiz(i, j, k)*Akx(i, j, k)*dufdk(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Akxi(i, j, k)*dudk(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiz(i, j, k)*Akxi(i, j, k)*dufdk(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akxi(i, j, k)*dwdk(i, j, k)*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akxi(i, j, k)*dwfdk(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2vdkdi(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2vfdkdi(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2wdkdi(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aky(i, j, k)*d2wfdkdi(i, j, k)*onethird*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*dwfdk(i, j, k)- &
                         Aiz(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*dwfdi(i, j, k)*fourthird+ &
                         2.0d0*Aiz(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*dwdk(i, j, k)- &
                         Aiz(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Akyi(i, j, k)*dvdk(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiz(i, j, k)*Akyi(i, j, k)*dvfdk(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Akyi(i, j, k)*dwdk(i, j, k)*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akyi(i, j, k)*dwfdk(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*Para%Cp*d2Tfdkdi(i, j, k)/Para%Pr+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2udkdi(i, j, k)*u_fluc(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2ufdkdi(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2vdkdi(i, j, k)*v_fluc(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*d2vfdkdi(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Akz(i, j, k)*d2wdkdi(i, j, k)*eighthird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akz(i, j, k)*d2wfdkdi(i, j, k)*eighthird*w(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*dufdi(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*dvfdk(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*dvfdi(i, j, k)+ &
                         Aiz(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*dwfdk(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*dwfdi(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dwdk(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Akzi(i, j, k)*dwfdk(i, j, k)*fourthird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2udids(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2ufdids(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2wdids(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asx(i, j, k)*d2wfdids(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*dwfds(i, j, k)- &
                         Aiz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dwfdi(i, j, k)*fourthird+ &
                         2.0d0*Aiz(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*dwds(i, j, k)- &
                         Aiz(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Asxi(i, j, k)*duds(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiz(i, j, k)*Asxi(i, j, k)*dufds(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asxi(i, j, k)*dwds(i, j, k)*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asxi(i, j, k)*dwfds(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2vdids(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2vfdids(i, j, k)*onethird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2wdids(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asy(i, j, k)*d2wfdids(i, j, k)*onethird*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*dwfds(i, j, k)- &
                         Aiz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*dwfdi(i, j, k)*fourthird+ &
                         2.0d0*Aiz(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*dwds(i, j, k)- &
                         Aiz(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*dwdi(i, j, k)*fourthird- &
                         Aiz(i, j, k)*Asyi(i, j, k)*dvds(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Aiz(i, j, k)*Asyi(i, j, k)*dvfds(i, j, k)*twothird*w(i, j, k)+ &
                         Aiz(i, j, k)*Asyi(i, j, k)*dwds(i, j, k)*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asyi(i, j, k)*dwfds(i, j, k)*v(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*Para%Cp*d2Tfdids(i, j, k)/Para%Pr+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2udids(i, j, k)*u_fluc(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2ufdids(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2vdids(i, j, k)*v_fluc(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*d2vfdids(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Asz(i, j, k)*d2wdids(i, j, k)*eighthird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Asz(i, j, k)*d2wfdids(i, j, k)*eighthird*w(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*dufds(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dufdi(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*dvfds(i, j, k)+ &
                         2.0d0*Aiz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dvfdi(i, j, k)+ &
                         Aiz(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*dwfds(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*dwfdi(i, j, k)*eighthird+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dwds(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Aiz(i, j, k)*Aszi(i, j, k)*dwfds(i, j, k)*fourthird*w(i, j, k)+ &
                         Aizk(i, j, k)*Akx(i, j, k)*dudi(i, j, k)*w_fluc(i, j, k)+ &
                         Aizk(i, j, k)*Akx(i, j, k)*dufdi(i, j, k)*w(i, j, k)- &
                         Aizk(i, j, k)*Akx(i, j, k)*dwdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aizk(i, j, k)*Akx(i, j, k)*dwfdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aizk(i, j, k)*Aky(i, j, k)*dvdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aizk(i, j, k)*Aky(i, j, k)*dvfdi(i, j, k)*w(i, j, k)- &
                         Aizk(i, j, k)*Aky(i, j, k)*dwdi(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aizk(i, j, k)*Aky(i, j, k)*dwfdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dwdi(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Aizk(i, j, k)*Akz(i, j, k)*dwfdi(i, j, k)*fourthird*w(i, j, k)+ &
                         Aizs(i, j, k)*Asx(i, j, k)*dudi(i, j, k)*w_fluc(i, j, k)+ &
                         Aizs(i, j, k)*Asx(i, j, k)*dufdi(i, j, k)*w(i, j, k)- &
                         Aizs(i, j, k)*Asx(i, j, k)*dwdi(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Aizs(i, j, k)*Asx(i, j, k)*dwfdi(i, j, k)*twothird*u(i, j, k)+ &
                         Aizs(i, j, k)*Asy(i, j, k)*dvdi(i, j, k)*w_fluc(i, j, k)+ &
                         Aizs(i, j, k)*Asy(i, j, k)*dvfdi(i, j, k)*w(i, j, k)- &
                         Aizs(i, j, k)*Asy(i, j, k)*dwdi(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aizs(i, j, k)*Asy(i, j, k)*dwfdi(i, j, k)*twothird*v(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*Para%Cp*dTfdi(i, j, k)/Para%Pr+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dudi(i, j, k)*u_fluc(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dufdi(i, j, k)*u(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dvdi(i, j, k)*v_fluc(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dvfdi(i, j, k)*v(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dwdi(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Aizs(i, j, k)*Asz(i, j, k)*dwfdi(i, j, k)*fourthird*w(i, j, k)+ &
                         Akx(i, j, k)**2*Para%Cp*d2Tfdk2(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)**2*d2udk2(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)**2*d2ufdk2(i, j, k)*fourthird*u(i, j, k)+Akx(i, j, k)**2*d2vdk2(i, j, k)*v_fluc(i, j, k)+ &
                         Akx(i, j, k)**2*d2vfdk2(i, j, k)*v(i, j, k)+Akx(i, j, k)**2*d2wdk2(i, j, k)*w_fluc(i, j, k)+ &
                         Akx(i, j, k)**2*d2wfdk2(i, j, k)*w(i, j, k)+Akx(i, j, k)**2*dudk(i, j, k)*dufdk(i, j, k)*eighthird+ &
                         2.0d0*Akx(i, j, k)**2*dvdk(i, j, k)*dvfdk(i, j, k)+2.0d0*Akx(i, j, k)**2*dwdk(i, j, k)*dwfdk(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dudk(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dufdk(i, j, k)*fourthird*u(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akxk(i, j, k)*dwfdk(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2udk2(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2ufdk2(i, j, k)*onethird*v(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2vdk2(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*d2vfdk2(i, j, k)*onethird*u(i, j, k)+ &
                         Akx(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*dvfdk(i, j, k)*twothird+ &
                         Akx(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*dvdk(i, j, k)*twothird+ &
                         Akx(i, j, k)*Akyk(i, j, k)*dudk(i, j, k)*v_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akyk(i, j, k)*dufdk(i, j, k)*v(i, j, k)- &
                         Akx(i, j, k)*Akyk(i, j, k)*dvdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Akx(i, j, k)*Akyk(i, j, k)*dvfdk(i, j, k)*twothird*u(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2udk2(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2ufdk2(i, j, k)*onethird*w(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2wdk2(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*d2wfdk2(i, j, k)*onethird*u(i, j, k)+ &
                         Akx(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*dwfdk(i, j, k)*twothird+ &
                         Akx(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*dwdk(i, j, k)*twothird+ &
                         Akx(i, j, k)*Akzk(i, j, k)*dudk(i, j, k)*w_fluc(i, j, k)+ &
                         Akx(i, j, k)*Akzk(i, j, k)*dufdk(i, j, k)*w(i, j, k)- &
                         Akx(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Akx(i, j, k)*Akzk(i, j, k)*dwfdk(i, j, k)*twothird*u(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*Para%Cp*d2Tfdkds(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*eighthird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asx(i, j, k)*d2ufdkds(i, j, k)*eighthird*u(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2vdkds(i, j, k)*v_fluc(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2vfdkds(i, j, k)*v(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2wdkds(i, j, k)*w_fluc(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*d2wfdkds(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*dufds(i, j, k)*eighthird+ &
                         Akx(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dufdk(i, j, k)*eighthird+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*dvfds(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*dvds(i, j, k)*dvfdk(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*dwfds(i, j, k)+ &
                         2.0d0*Akx(i, j, k)*Asx(i, j, k)*dwds(i, j, k)*dwfdk(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Akx(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dufds(i, j, k)*fourthird*u(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asxk(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2udkds(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2ufdkds(i, j, k)*onethird*v(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asy(i, j, k)*d2vfdkds(i, j, k)*onethird*u(i, j, k)- &
                         Akx(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*dvfds(i, j, k)*fourthird+ &
                         2.0d0*Akx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dvfdk(i, j, k)- &
                         Akx(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*dvds(i, j, k)*fourthird+ &
                         2.0d0*Akx(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*dvdk(i, j, k)+ &
                         Akx(i, j, k)*Asyk(i, j, k)*duds(i, j, k)*v_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asyk(i, j, k)*dufds(i, j, k)*v(i, j, k)- &
                         Akx(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Akx(i, j, k)*Asyk(i, j, k)*dvfds(i, j, k)*twothird*u(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2udkds(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2ufdkds(i, j, k)*onethird*w(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Akx(i, j, k)*Asz(i, j, k)*d2wfdkds(i, j, k)*onethird*u(i, j, k)- &
                         Akx(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*dwfds(i, j, k)*fourthird+ &
                         2.0d0*Akx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dwfdk(i, j, k)- &
                         Akx(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Akx(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*dwdk(i, j, k)+ &
                         Akx(i, j, k)*Aszk(i, j, k)*duds(i, j, k)*w_fluc(i, j, k)+ &
                         Akx(i, j, k)*Aszk(i, j, k)*dufds(i, j, k)*w(i, j, k)- &
                         Akx(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Akx(i, j, k)*Aszk(i, j, k)*dwfds(i, j, k)*twothird*u(i, j, k)- &
                         Akxk(i, j, k)*Aky(i, j, k)*dudk(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Akxk(i, j, k)*Aky(i, j, k)*dufdk(i, j, k)*twothird*v(i, j, k)+ &
                         Akxk(i, j, k)*Aky(i, j, k)*dvdk(i, j, k)*u_fluc(i, j, k)+ &
                         Akxk(i, j, k)*Aky(i, j, k)*dvfdk(i, j, k)*u(i, j, k)- &
                         Akxk(i, j, k)*Akz(i, j, k)*dudk(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Akxk(i, j, k)*Akz(i, j, k)*dufdk(i, j, k)*twothird*w(i, j, k)+ &
                         Akxk(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*u_fluc(i, j, k)+ &
                         Akxk(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*u(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*fourthird*u(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                         Akxs(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)*w(i, j, k)- &
                         Akxs(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Akxs(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*twothird*v(i, j, k)+ &
                         Akxs(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*u_fluc(i, j, k)+ &
                         Akxs(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*u(i, j, k)- &
                         Akxs(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Akxs(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*twothird*w(i, j, k)+ &
                         Akxs(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*u_fluc(i, j, k)+ &
                         Akxs(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*u(i, j, k)+Aky(i, j, k)**2*Para%Cp*d2Tfdk2(i, j, k)/Para%Pr+ &
                         Aky(i, j, k)**2*d2udk2(i, j, k)*u_fluc(i, j, k)+Aky(i, j, k)**2*d2ufdk2(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)**2*d2vdk2(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)**2*d2vfdk2(i, j, k)*fourthird*v(i, j, k)+Aky(i, j, k)**2*d2wdk2(i, j, k)*w_fluc(i, j, k)+ &
                         Aky(i, j, k)**2*d2wfdk2(i, j, k)*w(i, j, k)+2.0d0*Aky(i, j, k)**2*dudk(i, j, k)*dufdk(i, j, k)+ &
                         Aky(i, j, k)**2*dvdk(i, j, k)*dvfdk(i, j, k)*eighthird+ &
                         2.0d0*Aky(i, j, k)**2*dwdk(i, j, k)*dwfdk(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dvdk(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dvfdk(i, j, k)*fourthird*v(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                         Aky(i, j, k)*Akyk(i, j, k)*dwfdk(i, j, k)*w(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2vdk2(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2vfdk2(i, j, k)*onethird*w(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2wdk2(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*d2wfdk2(i, j, k)*onethird*v(i, j, k)+ &
                         Aky(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*dwfdk(i, j, k)*twothird+ &
                         Aky(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*dwdk(i, j, k)*twothird+ &
                         Aky(i, j, k)*Akzk(i, j, k)*dvdk(i, j, k)*w_fluc(i, j, k)+ &
                         Aky(i, j, k)*Akzk(i, j, k)*dvfdk(i, j, k)*w(i, j, k)- &
                         Aky(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aky(i, j, k)*Akzk(i, j, k)*dwfdk(i, j, k)*twothird*v(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2ufdkds(i, j, k)*onethird*v(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2vdkds(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asx(i, j, k)*d2vfdkds(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*dvfds(i, j, k)- &
                         Aky(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dvfdk(i, j, k)*fourthird+ &
                         2.0d0*Aky(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*dvds(i, j, k)- &
                         Aky(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*dvdk(i, j, k)*fourthird- &
                         Aky(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aky(i, j, k)*Asxk(i, j, k)*dufds(i, j, k)*twothird*v(i, j, k)+ &
                         Aky(i, j, k)*Asxk(i, j, k)*dvds(i, j, k)*u_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asxk(i, j, k)*dvfds(i, j, k)*u(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*Para%Cp*d2Tfdkds(i, j, k)/Para%Pr+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2udkds(i, j, k)*u_fluc(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2ufdkds(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*eighthird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asy(i, j, k)*d2vfdkds(i, j, k)*eighthird*v(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2wdkds(i, j, k)*w_fluc(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*d2wfdkds(i, j, k)*w(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*dufds(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dufdk(i, j, k)+ &
                         Aky(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*dvfds(i, j, k)*eighthird+ &
                         Aky(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*dvfdk(i, j, k)*eighthird+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*dwfds(i, j, k)+ &
                         2.0d0*Aky(i, j, k)*Asy(i, j, k)*dwds(i, j, k)*dwfdk(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Aky(i, j, k)*Asyk(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dvfds(i, j, k)*fourthird*v(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asyk(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2vdkds(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2vfdkds(i, j, k)*onethird*w(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Aky(i, j, k)*Asz(i, j, k)*d2wfdkds(i, j, k)*onethird*v(i, j, k)- &
                         Aky(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*dwfds(i, j, k)*fourthird+ &
                         2.0d0*Aky(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dwfdk(i, j, k)- &
                         Aky(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*dwds(i, j, k)*fourthird+ &
                         2.0d0*Aky(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*dwdk(i, j, k)+ &
                         Aky(i, j, k)*Aszk(i, j, k)*dvds(i, j, k)*w_fluc(i, j, k)+ &
                         Aky(i, j, k)*Aszk(i, j, k)*dvfds(i, j, k)*w(i, j, k)- &
                         Aky(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Aky(i, j, k)*Aszk(i, j, k)*dwfds(i, j, k)*twothird*v(i, j, k)- &
                         Akyk(i, j, k)*Akz(i, j, k)*dvdk(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Akyk(i, j, k)*Akz(i, j, k)*dvfdk(i, j, k)*twothird*w(i, j, k)+ &
                         Akyk(i, j, k)*Akz(i, j, k)*dwdk(i, j, k)*v_fluc(i, j, k)+ &
                         Akyk(i, j, k)*Akz(i, j, k)*dwfdk(i, j, k)*v(i, j, k)+ &
                         Akys(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*v_fluc(i, j, k)+ &
                         Akys(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*v(i, j, k)- &
                         Akys(i, j, k)*Asx(i, j, k)*dvdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Akys(i, j, k)*Asx(i, j, k)*dvfdk(i, j, k)*twothird*u(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Akys(i, j, k)*Asy(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*fourthird*v(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*w_fluc(i, j, k)+ &
                         Akys(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)*w(i, j, k)- &
                         Akys(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Akys(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*twothird*w(i, j, k)+ &
                         Akys(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*v_fluc(i, j, k)+ &
                         Akys(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*v(i, j, k)+Akz(i, j, k)**2*Para%Cp*d2Tfdk2(i, j, k)/Para%Pr+ &
                         Akz(i, j, k)**2*d2udk2(i, j, k)*u_fluc(i, j, k)+Akz(i, j, k)**2*d2ufdk2(i, j, k)*u(i, j, k)+ &
                         Akz(i, j, k)**2*d2vdk2(i, j, k)*v_fluc(i, j, k)+Akz(i, j, k)**2*d2vfdk2(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)**2*d2wdk2(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Akz(i, j, k)**2*d2wfdk2(i, j, k)*fourthird*w(i, j, k)+2.0d0*Akz(i, j, k)**2*dudk(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Akz(i, j, k)**2*dvdk(i, j, k)*dvfdk(i, j, k)+ &
                         Akz(i, j, k)**2*dwdk(i, j, k)*dwfdk(i, j, k)*eighthird+ &
                         Akz(i, j, k)*Akzk(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dwdk(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Akz(i, j, k)*Akzk(i, j, k)*dwfdk(i, j, k)*fourthird*w(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2udkds(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2ufdkds(i, j, k)*onethird*w(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2wdkds(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asx(i, j, k)*d2wfdkds(i, j, k)*onethird*u(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*dwfds(i, j, k)- &
                         Akz(i, j, k)*Asx(i, j, k)*duds(i, j, k)*dwfdk(i, j, k)*fourthird+ &
                         2.0d0*Akz(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*dwds(i, j, k)- &
                         Akz(i, j, k)*Asx(i, j, k)*dufds(i, j, k)*dwdk(i, j, k)*fourthird- &
                         Akz(i, j, k)*Asxk(i, j, k)*duds(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Akz(i, j, k)*Asxk(i, j, k)*dufds(i, j, k)*twothird*w(i, j, k)+ &
                         Akz(i, j, k)*Asxk(i, j, k)*dwds(i, j, k)*u_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asxk(i, j, k)*dwfds(i, j, k)*u(i, j, k)+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2vdkds(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2vfdkds(i, j, k)*onethird*w(i, j, k)+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2wdkds(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asy(i, j, k)*d2wfdkds(i, j, k)*onethird*v(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*dwfds(i, j, k)- &
                         Akz(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*dwfdk(i, j, k)*fourthird+ &
                         2.0d0*Akz(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*dwds(i, j, k)- &
                         Akz(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*dwdk(i, j, k)*fourthird- &
                         Akz(i, j, k)*Asyk(i, j, k)*dvds(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Akz(i, j, k)*Asyk(i, j, k)*dvfds(i, j, k)*twothird*w(i, j, k)+ &
                         Akz(i, j, k)*Asyk(i, j, k)*dwds(i, j, k)*v_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asyk(i, j, k)*dwfds(i, j, k)*v(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*Para%Cp*d2Tfdkds(i, j, k)/Para%Pr+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2udkds(i, j, k)*u_fluc(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2ufdkds(i, j, k)*u(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2vdkds(i, j, k)*v_fluc(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*d2vfdkds(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)*Asz(i, j, k)*d2wdkds(i, j, k)*eighthird*w_fluc(i, j, k)+ &
                         Akz(i, j, k)*Asz(i, j, k)*d2wfdkds(i, j, k)*eighthird*w(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*dufds(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dufdk(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*dvfds(i, j, k)+ &
                         2.0d0*Akz(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dvfdk(i, j, k)+ &
                         Akz(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*dwfds(i, j, k)*eighthird+ &
                         Akz(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*dwfdk(i, j, k)*eighthird+ &
                         Akz(i, j, k)*Aszk(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Akz(i, j, k)*Aszk(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dwds(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Akz(i, j, k)*Aszk(i, j, k)*dwfds(i, j, k)*fourthird*w(i, j, k)+ &
                         Akzs(i, j, k)*Asx(i, j, k)*dudk(i, j, k)*w_fluc(i, j, k)+ &
                         Akzs(i, j, k)*Asx(i, j, k)*dufdk(i, j, k)*w(i, j, k)- &
                         Akzs(i, j, k)*Asx(i, j, k)*dwdk(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Akzs(i, j, k)*Asx(i, j, k)*dwfdk(i, j, k)*twothird*u(i, j, k)+ &
                         Akzs(i, j, k)*Asy(i, j, k)*dvdk(i, j, k)*w_fluc(i, j, k)+ &
                         Akzs(i, j, k)*Asy(i, j, k)*dvfdk(i, j, k)*w(i, j, k)- &
                         Akzs(i, j, k)*Asy(i, j, k)*dwdk(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Akzs(i, j, k)*Asy(i, j, k)*dwfdk(i, j, k)*twothird*v(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*Para%Cp*dTfdk(i, j, k)/Para%Pr+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dudk(i, j, k)*u_fluc(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dufdk(i, j, k)*u(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dvdk(i, j, k)*v_fluc(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dvfdk(i, j, k)*v(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dwdk(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Akzs(i, j, k)*Asz(i, j, k)*dwfdk(i, j, k)*fourthird*w(i, j, k)+ &
                         Asx(i, j, k)**2*Para%Cp*d2Tfds2(i, j, k)/Para%Pr+ &
                         Asx(i, j, k)**2*d2uds2(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Asx(i, j, k)**2*d2ufds2(i, j, k)*fourthird*u(i, j, k)+Asx(i, j, k)**2*d2vds2(i, j, k)*v_fluc(i, j, k)+ &
                         Asx(i, j, k)**2*d2vfds2(i, j, k)*v(i, j, k)+Asx(i, j, k)**2*d2wds2(i, j, k)*w_fluc(i, j, k)+ &
                         Asx(i, j, k)**2*d2wfds2(i, j, k)*w(i, j, k)+Asx(i, j, k)**2*duds(i, j, k)*dufds(i, j, k)*eighthird+ &
                         2.0d0*Asx(i, j, k)**2*dvds(i, j, k)*dvfds(i, j, k)+2.0d0*Asx(i, j, k)**2*dwds(i, j, k)*dwfds(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Asx(i, j, k)*Asxs(i, j, k)*duds(i, j, k)*fourthird*u_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dufds(i, j, k)*fourthird*u(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asxs(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2uds2(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2ufds2(i, j, k)*onethird*v(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2vds2(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*d2vfds2(i, j, k)*onethird*u(i, j, k)+ &
                         Asx(i, j, k)*Asy(i, j, k)*duds(i, j, k)*dvfds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*dvds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Asys(i, j, k)*duds(i, j, k)*v_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asys(i, j, k)*dufds(i, j, k)*v(i, j, k)- &
                         Asx(i, j, k)*Asys(i, j, k)*dvds(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Asx(i, j, k)*Asys(i, j, k)*dvfds(i, j, k)*twothird*u(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2uds2(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2ufds2(i, j, k)*onethird*w(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2wds2(i, j, k)*onethird*u_fluc(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*d2wfds2(i, j, k)*onethird*u(i, j, k)+ &
                         Asx(i, j, k)*Asz(i, j, k)*duds(i, j, k)*dwfds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*dwds(i, j, k)*twothird+ &
                         Asx(i, j, k)*Aszs(i, j, k)*duds(i, j, k)*w_fluc(i, j, k)+ &
                         Asx(i, j, k)*Aszs(i, j, k)*dufds(i, j, k)*w(i, j, k)- &
                         Asx(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*twothird*u_fluc(i, j, k)- &
                         Asx(i, j, k)*Aszs(i, j, k)*dwfds(i, j, k)*twothird*u(i, j, k)- &
                         Asxs(i, j, k)*Asy(i, j, k)*duds(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Asxs(i, j, k)*Asy(i, j, k)*dufds(i, j, k)*twothird*v(i, j, k)+ &
                         Asxs(i, j, k)*Asy(i, j, k)*dvds(i, j, k)*u_fluc(i, j, k)+ &
                         Asxs(i, j, k)*Asy(i, j, k)*dvfds(i, j, k)*u(i, j, k)- &
                         Asxs(i, j, k)*Asz(i, j, k)*duds(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Asxs(i, j, k)*Asz(i, j, k)*dufds(i, j, k)*twothird*w(i, j, k)+ &
                         Asxs(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*u_fluc(i, j, k)+ &
                         Asxs(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*u(i, j, k)+Asy(i, j, k)**2*Para%Cp*d2Tfds2(i, j, k)/Para%Pr+ &
                         Asy(i, j, k)**2*d2uds2(i, j, k)*u_fluc(i, j, k)+Asy(i, j, k)**2*d2ufds2(i, j, k)*u(i, j, k)+ &
                         Asy(i, j, k)**2*d2vds2(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Asy(i, j, k)**2*d2vfds2(i, j, k)*fourthird*v(i, j, k)+Asy(i, j, k)**2*d2wds2(i, j, k)*w_fluc(i, j, k)+ &
                         Asy(i, j, k)**2*d2wfds2(i, j, k)*w(i, j, k)+2.0d0*Asy(i, j, k)**2*duds(i, j, k)*dufds(i, j, k)+ &
                         Asy(i, j, k)**2*dvds(i, j, k)*dvfds(i, j, k)*eighthird+ &
                         2.0d0*Asy(i, j, k)**2*dwds(i, j, k)*dwfds(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Asy(i, j, k)*Asys(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dvds(i, j, k)*fourthird*v_fluc(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dvfds(i, j, k)*fourthird*v(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dwds(i, j, k)*w_fluc(i, j, k)+ &
                         Asy(i, j, k)*Asys(i, j, k)*dwfds(i, j, k)*w(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2vds2(i, j, k)*onethird*w_fluc(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2vfds2(i, j, k)*onethird*w(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2wds2(i, j, k)*onethird*v_fluc(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*d2wfds2(i, j, k)*onethird*v(i, j, k)+ &
                         Asy(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*dwfds(i, j, k)*twothird+ &
                         Asy(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*dwds(i, j, k)*twothird+ &
                         Asy(i, j, k)*Aszs(i, j, k)*dvds(i, j, k)*w_fluc(i, j, k)+ &
                         Asy(i, j, k)*Aszs(i, j, k)*dvfds(i, j, k)*w(i, j, k)- &
                         Asy(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*twothird*v_fluc(i, j, k)- &
                         Asy(i, j, k)*Aszs(i, j, k)*dwfds(i, j, k)*twothird*v(i, j, k)- &
                         Asys(i, j, k)*Asz(i, j, k)*dvds(i, j, k)*twothird*w_fluc(i, j, k)- &
                         Asys(i, j, k)*Asz(i, j, k)*dvfds(i, j, k)*twothird*w(i, j, k)+ &
                         Asys(i, j, k)*Asz(i, j, k)*dwds(i, j, k)*v_fluc(i, j, k)+ &
                         Asys(i, j, k)*Asz(i, j, k)*dwfds(i, j, k)*v(i, j, k)+Asz(i, j, k)**2*Para%Cp*d2Tfds2(i, j, k)/Para%Pr+ &
                         Asz(i, j, k)**2*d2uds2(i, j, k)*u_fluc(i, j, k)+Asz(i, j, k)**2*d2ufds2(i, j, k)*u(i, j, k)+ &
                         Asz(i, j, k)**2*d2vds2(i, j, k)*v_fluc(i, j, k)+Asz(i, j, k)**2*d2vfds2(i, j, k)*v(i, j, k)+ &
                         Asz(i, j, k)**2*d2wds2(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Asz(i, j, k)**2*d2wfds2(i, j, k)*fourthird*w(i, j, k)+2.0d0*Asz(i, j, k)**2*duds(i, j, k)*dufds(i, j, k)+ &
                         2.0d0*Asz(i, j, k)**2*dvds(i, j, k)*dvfds(i, j, k)+ &
                         Asz(i, j, k)**2*dwds(i, j, k)*dwfds(i, j, k)*eighthird+ &
                         Asz(i, j, k)*Aszs(i, j, k)*Para%Cp*dTfds(i, j, k)/Para%Pr+ &
                         Asz(i, j, k)*Aszs(i, j, k)*duds(i, j, k)*u_fluc(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dufds(i, j, k)*u(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dvds(i, j, k)*v_fluc(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dvfds(i, j, k)*v(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dwds(i, j, k)*fourthird*w_fluc(i, j, k)+ &
                         Asz(i, j, k)*Aszs(i, j, k)*dwfds(i, j, k)*fourthird*w(i, j, k)

            du(i, j, k, 2) = du(i, j, k, 2)+ &
                             dmufdk(i, j, k)*dmufdk_term(1)+ &
                             dmufdi(i, j, k)*dmufdi_term(1)+ &
                             dmufds(i, j, k)*dmufds_term(1)+ &
                             mu_fluc(i, j, k)*muf_term(1)+ &
                             dmudk(i, j, k)*dmudk_term(1)+ &
                             dmudi(i, j, k)*dmudi_term(1)+ &
                             dmuds(i, j, k)*dmuds_term(1)+ &
                             mu(i, j, k)*mu_term(1)

            du(i, j, k, 3) = du(i, j, k, 3)+ &
                             dmufdk(i, j, k)*dmufdk_term(2)+ &
                             dmufdi(i, j, k)*dmufdi_term(2)+ &
                             dmufds(i, j, k)*dmufds_term(2)+ &
                             mu_fluc(i, j, k)*muf_term(2)+ &
                             dmudk(i, j, k)*dmudk_term(2)+ &
                             dmudi(i, j, k)*dmudi_term(2)+ &
                             dmuds(i, j, k)*dmuds_term(2)+ &
                             mu(i, j, k)*mu_term(2)

            du(i, j, k, 4) = du(i, j, k, 4)+ &
                             dmufdk(i, j, k)*dmufdk_term(3)+ &
                             dmufdi(i, j, k)*dmufdi_term(3)+ &
                             dmufds(i, j, k)*dmufds_term(3)+ &
                             mu_fluc(i, j, k)*muf_term(3)+ &
                             dmudk(i, j, k)*dmudk_term(3)+ &
                             dmudi(i, j, k)*dmudi_term(3)+ &
                             dmuds(i, j, k)*dmuds_term(3)+ &
                             mu(i, j, k)*mu_term(3)

            du(i, j, k, 5) = du(i, j, k, 5)+ &
                             dmufdk(i, j, k)*dmufdk_term(4)+ &
                             dmufdi(i, j, k)*dmufdi_term(4)+ &
                             dmufds(i, j, k)*dmufds_term(4)+ &
                             mu_fluc(i, j, k)*muf_term(4)+ &
                             dmudk(i, j, k)*dmudk_term(4)+ &
                             dmudi(i, j, k)*dmudi_term(4)+ &
                             dmuds(i, j, k)*dmuds_term(4)+ &
                             mu(i, j, k)*mu_term(4)

        end do
        end do
        end do

    end subroutine flux_viscous3d_d

    subroutine comput_mu_d
        use precision, only: dp
        use mesh_data, only: nx, ny, nz
        use flow_para, only: Para
        use flow_data, only: T
        use fluc_data, only: mu_fluc, T_fluc
        implicit none

        real(dp):: Tsb, tmpR, C_const
        integer:: i, j, k

        Tsb = 110.4d0/Para%T_inf
        TmpR = 1.0d0/Para%Re
        C_const = TmpR*(1.0d0+Tsb)/2.0d0

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            mu_fluc(i, j, k) = C_const*sqrt(T(i, j, k))*(T(i, j, k)+3.0d0*Tsb)*T_fluc(i, j, k)/((T(i, j, k)+Tsb)**2)
        end do
        end do
        end do

    end subroutine comput_mu_d

    elemental function diff_mu_d(dTfd, Tf, dTd, T) result(dmufd)
        use flow_para, only: Para, dp
        implicit none

        real(dp) :: dmufd
        real(dp), intent(IN)  :: dTfd, Tf, dTd, T
        real(dp) :: Tsb, TmpR, Term1, Term2

        Tsb = 110.4d0/Para%T_inf
        TmpR = 1.0d0/Para%Re

        Term1 = Tf*dTd*(-T**2-6.0d0*T*Tsb+3.0d0*Tsb**2)
        Term2 = 2.0d0*dTfd*T*(T+Tsb)*(T+3.0d0*Tsb)

        dmufd = TmpR*(Tsb+1.0d0)*(Term1+Term2)/(4.0d0*sqrt(T)*(T+Tsb)**3)

    end function diff_mu_d


    subroutine compute_derivatives_d
        ! Velocity and temperature derivatives used by the viscous terms; also reused by the
        ! lid-driven cavity wall-pressure condition, which computes them right after setting the walls
        use mesh_data
        use flow_data
        use fluc_data
        use Centrals
        implicit none

        call OCFD_dx0(u_fluc, dufdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(v_fluc, dvfdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(w_fluc, dwfdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx0(T_fluc, dTfdk(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcI)

        call OCFD_dy0(u_fluc, dufdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(v_fluc, dvfdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(w_fluc, dwfdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(T_fluc, dTfdi(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcJ)

        call OCFD_dz0(u_fluc, dufds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(v_fluc, dvfds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(w_fluc, dwfds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(T_fluc, dTfds(1:nx, 1:ny, 1:nz), nx, ny, nz, Scheme%Scheme_Vis, BcK)

        ! Mixed derivatives use k-derivatives in y/z and i-derivatives in z.
        call mpi_exchange_yz3d(dufdk)
        call mpi_exchange_yz3d(dvfdk)
        call mpi_exchange_yz3d(dwfdk)
        call mpi_exchange_yz3d(dTfdk)
        call mpi_exchange_z3d(dufdi)
        call mpi_exchange_z3d(dvfdi)
        call mpi_exchange_z3d(dwfdi)
        call mpi_exchange_z3d(dTfdi)

        call OCFD_dy0(dufdk, d2ufdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dvfdk, d2vfdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dwfdk, d2wfdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy0(dTfdk, d2Tfdkdi, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dz0(dufdk, d2ufdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dvfdk, d2vfdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dwfdk, d2wfdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dTfdk, d2Tfdkds, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dufdi, d2ufdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dvfdi, d2vfdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dwfdi, d2wfdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz0(dTfdi, d2Tfdids, nx, ny, nz, Scheme%Scheme_Vis, BcK)

        call OCFD_dx2(u_fluc, d2ufdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(v_fluc, d2vfdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(w_fluc, d2wfdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dx2(T_fluc, d2Tfdk2, nx, ny, nz, Scheme%Scheme_Vis, BcI)
        call OCFD_dy2(u_fluc, d2ufdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(v_fluc, d2vfdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(w_fluc, d2wfdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dy2(T_fluc, d2Tfdi2, nx, ny, nz, Scheme%Scheme_Vis, BcJ)
        call OCFD_dz2(u_fluc, d2ufds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(v_fluc, d2vfds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(w_fluc, d2wfds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)
        call OCFD_dz2(T_fluc, d2Tfds2, nx, ny, nz, Scheme%Scheme_Vis, BcK)

    end subroutine compute_derivatives_d

end module Flux_Viscous_d
