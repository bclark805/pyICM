#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:45:19 2021

@author: jbclark8
"""


# module to calculate the interactions between
# ISS and the water column

# JB Clark, NASA GSFC Aug 2021


import DefVars as dv


def deltaISS(Q, Volume, FlushingRate, Input_ISS, ISS):
    # import flushing rate, the outside concentration of ISS,
    #  and the current concentration of ISS
    # now calculate the change in ISS due to inputs, sinking, and outputs
    DTISS = 0.0
    # change in ISS due to inputs into the box, g m^-3 s^-1
    DTISS = Input_ISS*Q/Volume
    # change in ISS due to flushing of the box
    DTISS = DTISS-ISS*FlushingRate
    # change in ISS due to sinking
    SinkingRate = dv.WISS/(dv.Z*0.5)
    DTISS = DTISS-DTISS*SinkingRate/86400

    return DTISS


###########
