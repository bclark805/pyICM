#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:45:19 2021

@author: jbclark8
"""


import numpy as np

# functions utilized throughout the ICM code
# will be called from other modules

# JB Clark, NASA GSFC Aug 2021


def get_bell(A1, A2, Ref, Val):
    # this is a function to create a bell shaped function, utilized
    # initially in phytoplankton temperature response to growth
    # but also more recently in the resposne of sediment sinking
    # velocity to concentration
    # can also be be utilized for the flocculation of DOM
    # (or anything else the user can think of)
    #######
    # A1=sub-optimal shape parameter
    # A2=super-optimal shape parameter
    # Ref=reference value where shape changes
    # Val=actual value of state variable

    if Val < Ref:
        theta = np.exp(-A1*(Val-Ref)**2)
    else:
        theta = np.exp(-A2*(Ref-Val)**2)
    return theta


def get_ft(kt, Tref, T):
    # exponentially based temperature function used in many
    # biogeochemical processes
    # kt is the shape parameter for the temperature function
    # Tref is the reference temperature
    # T is the actual temperature
    theta = np.exp(kt*(T-Tref))
    return theta


def CalcVolume(X, Y, Z):
    # function to calculate the volume of the box
    Vol = X*Y*Z
    return Vol


def CalcFlushingRate(Q, Vol):
    # calculate the fluxing rate in days
    # Q is the flow, in m^3/s
    # Vol is the volume of the box

    Qseconds = Q  # convert discharge from seconds to days
    FlushingRate = 1/(Vol/Qseconds)
    return FlushingRate
