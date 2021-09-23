

# Module for Algal Growth, Grazing, and Sinking
# Adapted from Cerco and Coel 1993; Clark et al. 2020
# JB Clark NASA GSFC Aug 2021

import numpy as np
import Nutrients as nt
import DefVars as dv
import functions as fct
import pandas as pd

Algae1 = pd.DataFrame()
Algae2 = pd.DataFrame()


def DTB1(B1, NH4, NO3, T, PAR, FRate, RiverAlg, Q, V):

    # calculate temperaturee scaled production rate
    P1 = dv.PM1*fct.get_bell(dv.kta1a, dv.kta1b, dv.TrefB1, T)

    # calculate the temperature scaled metabolic raaate
    BM1 = dv.BMR1*fct.get_ft(dv.ktb1, dv.Tb1, T)

    # calculate the predation
    PR1 = dv.BPR1*B1*B1*fct.get_ft(dv.ktpr1, dv.Tpr1, T)

    # calculate the nutrient preference for NH4 or NO3
    PN1 = nt.Npref(NH4, NO3, dv.khn1)
    NL1 = Nlim(NH4, NO3, PN1, dv.khn1)

    # calculate light limitation
    FI1 = LightLim(PAR, P1, dv.alpha1)

    MyLim = np.minimum(NL1, FI1)

    # now calculate the scaled production rate
    # from nutrient and light limitation
    P1 = P1*MyLim/dv.cchla1

    # calculate NPP substracting out photorespiration and Metabolism
    NPP1 = (P1*(1-dv.PRSP1)-BM1)*B1
    # calculatee grosss primary production
    GPP1 = P1*B1

    # calculate the loss of phytoplankton due to sinking
    SinkingRate = dv.WB1/(dv.Z*0.5)

    # now calculate the total change in B1 over time, g m^-3 s^-1
    DTB1out = (NPP1-PR1-SinkingRate*B1)/86400-B1*FRate+RiverAlg*Q/V

    # calculate ammonium and nitrate uptake  of algae
    NH4A = AlgNH4(PN1, P1, B1, PR1, dv.FNP1, dv.ANC1)
    NO3A = AlgNO3(PN1, P1, B1, dv.ANC1)

    # now fill the data structure with variables that can be used elsewhere

    Algae1.NPP = NPP1 * dv.Z * 0.5
    Algae1.GPP = GPP1
    Algae1.DTB1 = DTB1out
    Algae1.NH4A = NH4A
    Algae1.NO3A = NO3A
    Algae1.PN = PN1
    Algae1.NL = NL1
    Algae1.PR1 = PR1
    Algae1.P1 = P1
    Algae1.FI = FI1

    return Algae1


def DTB2(B2, NH4, NO3, T, PAR, FRate, RiverAlg, Q, V):

    # calculate temperaturee scaled production rate
    P2 = dv.PM2*fct.get_bell(dv.kta2a, dv.kta2b, dv.TrefB2, T)

    # calculate the temperature scaled metabolic raaate
    BM2 = dv.BMR2*fct.get_ft(dv.ktb2, dv.Tb2, T)

    # calculate the predation
    PR2 = dv.BPR2*B2*B2*fct.get_ft(dv.ktpr2, dv.Tpr2, T)

    # calculate the nutrient preference for NH4 or NO3
    PN2 = nt.Npref(NH4, NO3, dv.khn2)
    NL2 = Nlim(NH4, NO3, PN2, dv.khn2)

    # calculate light limitation
    FI2 = LightLim(PAR, P2, dv.alpha2)

    MyLim = np.minimum(NL2, FI2)

# now calculate the scaled production rate from nutrient and light limitation
    P2 = P2*MyLim/dv.cchla2

    # calculate NPP substracting out photorespiration and Metabolism
    NPP2 = (P2*(1-dv.PRSP2)-BM2)*B2
    # calculatee grosss primary production
    GPP2 = P2*B2

    # calculate the loss of phytoplankton due to sinking
    SinkingRate = dv.WB2/(dv.Z*0.5)

    # now calculate the total change in B1 over time
    DTB2out = (NPP2-PR2-SinkingRate*B2)/86400-B2*FRate+RiverAlg*Q/V

    # calculate ammonium balance of algae
    NH4A = AlgNH4(PN2, P2, B2, PR2, dv.FNP2, dv.ANC2)

    NO3A = AlgNO3(PN2, P2, B2, dv.ANC2)

    # now fill the data structure with variables that can be used elsewhere
    Algae2.NPP = NPP2 * dv.Z * 0.5
    Algae2.GPP = GPP2
    Algae2.DTB2 = DTB2out
    Algae2.NH4A = NH4A
    Algae2.NO3A = NO3A
    Algae2.PN = PN2
    Algae2.NL = NL2
    Algae2.PR2 = PR2
    Algae2.P2 = P2
    Algae2.FI = FI2

    return Algae2


def AlgNH4(PN, P, B, PR, FNP, ANC):
    # net change in NH4 due to algal uptake and algal predation
    # PN is defined above,
    # P is temperature scaled algal production rate
    # B is the algal biomass, in Carbon
    # PR is the Algal predation rate
    # FNP is the fractino of algal predation that is lost as ammonium (defined in DefVars)
    # ANC is the algal nitrogen to carbon ration
    NH4A = -PN*P*B*ANC+FNP*PR*ANC
    return NH4A


def AlgNO3(PN, P, B, ANC):
    # change in NO3 due to algal uptake
    # PN is defined above,
    # P is temperature scaled algal production rate
    # B is the algal biomass, in Carbon
    NO3A = (PN-1.0)*P*B*ANC
    return NO3A


def Nlim(NH4, NO3, PN, KHN):

    topside = 2*NH4+NO3-PN*(NO3+NH4)
    bottomside = KHN+2*NH4+NO3-PN*(NO3+NH4)

    NL = topside/bottomside
    return NL


def LightLim(PAR, P, alpha):

    IK = P/alpha

    FI = PAR/np.sqrt(IK*IK+PAR*PAR)

    return FI
