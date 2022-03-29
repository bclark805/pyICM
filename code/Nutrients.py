#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:45:19 2021

@author: jbclark8
"""


# module to calculate reactions related to Nutrientss (currently just inorganic nitrogen)
# adapted from Clark et al. 2020
# JB Clark, NASA GSFC Aug 2021

import numpy as np
import functions as fct
import DefVars as dv
# import Algae
import pandas as pd


def DTNO3(NH4, NO3in, T, DO2, NO3A, Q, V, FR, RivNO3):
    deltaNO3 = pd.DataFrame
    # function to calculate the change in NO3 over time due to biological and chemical changes
    # will deal with inputs and outputs in the main program
    DTNO3out = 0.

    # calculate the input from the river/upstream
    DTNO3out = RivNO3*Q/V
    # calculate the loss out downstream
    DTNO3out = DTNO3out-NO3in*FR

    NT = Nitrification(NH4, T, DO2)
    # DNT=Denitrification(NO3,T,DO2)
    DNT = 0.

    DTNO3out = (NT-DNT+NO3A)/86400+DTNO3out

    deltaNO3.DTNO3 = DTNO3out
    deltaNO3.NT = NT

    return deltaNO3


def DTNH4(NH4in, NO3in, T, DO2, MNLDON, NH4A, Q, V, FR, RivNH4):
    deltaNH4 = pd.DataFrame
    # function to calculate the change in ammonium due to biological
    # and chemical properties
    # NH4 is ammonium concentration coming in
    # MNLDON is the remineralization of DON to ammonium

    # calculate the input from the river/upstream g m^-3 s^-1
    DTNH4out = RivNH4*Q/V

    # calculate the loss out downstream,  g m^-3 s^-1
    DTNH4out = DTNH4out-NH4in*FR

    # call nitrification
    NT = Nitrification(NH4in, T, DO2)

    DTNH4out = (-NT+MNLDON+NH4A)/86400+DTNH4out

    deltaNH4.DTNH4 = DTNH4out
    deltaNH4.NT = -NT

    return deltaNH4

# nitrification


def Nitrification(NH4in, T, DO2):

    # temperature effect on nitrification
    FTNT = fct.get_bell(dv.ktnt1, dv.ktnt2, dv.TMNT, T)
    # oxygen effect on nitrification
    O2lim = DO2/(dv.KHONT+DO2)
    # NH4 limitation on nitrification
    NH4lim = NH4in/(dv.KHNNT+NH4in)
    # nitrifcation in units of g N m^-3 d^-1
    NT = O2lim*NH4lim*FTNT*dv.NTM*NH4in

    return NT


def Npref(NH4in, NO3in, KHN):
    # function to calculate the algal preference for nitrogen
    # brings in NH4, NO3 and half saturation for algal nitrogen uptake, KHN
    # is used heree and in algae growth
    # total dissolved inorganic nitrogen
    DIN = NH4in+NO3in
    # Left side of Function (Clark et al. 2020)
    LeftSide = NO3in/((KHN+NH4in)*(KHN+NO3in))
    # Right side of Function (Clark et al. 2020)
    RightSide = KHN/(DIN*(KHN+NO3in))
    # Nitrogen preference equation
    PN = NH4in*(LeftSide+RightSide)

    return PN
