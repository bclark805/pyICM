# module to calculate reactions related to DOM (C and N)
# adapted from Clark et al. 2020
# JB Clark, NASA GSFC Aug 2021
import pandas as pd
import numpy as np
import functions as fct
import DefVars as dv
import POM
import Nutrients as nt
import Light as lt
def DTDOC(CDOC,NCDOC,CP1,CP2,PR1,PR2,Ed,S,T):
    # DTCDOC isthe change in colored DOC concentration over time
    # DTNCDOC is the change in non-colored DOC concentration over time
    # CP is the algal photorespiration and metabolism for algae 1 and 2
    # PR is the algal predation for algae 1 and 2
    # Ed is the spectral downwelling irradiance
    # S is salinity
    # T is temperature
    # LPOC is the labile POC concentration
    # RPOC is the refractory POC concentration
    C1=dv.FCD1*CP1 + dv.FCDP*PR1
    C2=dv.FCD2*CP2 + dv.FCDP*PR2
    AlgDOC=C1+C2

    MNLCDOC1=CDOC[1]
    MNLCDOC2=CDOC[2]
    MNLCDOC3=CDOC[3]

    MNLNCDOC1=NCDOC[1]
    MNLNCDOC2=NCDOC[2]
    MNLNCDOC3=NCDOC[3]
    
    COAGC=1.0
    DTPDOC=1.0

    return DTCDOC, DTNCDOC, COAGC, DTPDOC

#######