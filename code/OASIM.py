# OASIM module (Greg and Casey, 2009)
# to calculate water leaving radiance
# from the fully spectral downwelling radiance
# total absorption and total scattering

import numpy as np

def OASIM(lam, T, S, Ed, aT, bbT):
    # main function
    rho = 0.021 # surface reflectance
    rn = 1.341 # index of refraction
    rn2 = rn*rn
    Q = np.pi
    bbrw = 0.5 # backscattering to total scattering ratio of water

    bsw = zhangbsw(lam, T, S)

    bb = bbrw * bsw + bbT
    sfceu = eu0n(Ed, aT, bb)
    rlwn = (1.0 - rho) * sfceu / (rn2 * Q)
    Rrs=rlwn / Ed

    return rlwn, Rrs

def zhangbsw(lam, Tc, S):
    # functions to calculate scattering of sea water with varying salinity
    # from
    # Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
    # seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710
    # values of the constants
    Na = 6.0221417930E23 # Avogadro's constant
    Kbz = 1.3806503E-23  # Boltzmann constant
    Tk = Tc + 273.15       #  Absolute tempearture
    rM0 = 18.0E-3        #  Molecular weight of water in kg/mol
    pi = np.arccos(-1.0)
    rlambda = lam

    delta = 0.039 # . Farinato and Roswell (1976)

    nsw,dnds=getnswdnds(rlambda,Tc,S)
    
    IsoComp=getIsoComp(Tc,S)

    density_sw=get_density_sw(Tc,S)
    
    dlnawds=getdlnawds(Tc,S)
    
    DFRI=getDFRI(nsw)

    # volume scattering at 90 degree due to the density fluctuation
    # beta_df = pi*pi/2.0*((lambda*1e-9)**(-4))*Kbz*Tk*IsoComp*DFRI**2*(6+6*delta)/(6-7*delta)
    beta_df = pi * pi / 2.0 * ((rlambda * 1.0E-9)**(-4)) * Kbz \
    * Tk * IsoComp * DFRI**2 * (6.0 + 6.0 * delta) / (6.0-7.0 * delta)

    # volume scattering at 90 degree due to the concentration fluctuation
    # flu_con = S*M0*dnds.^2/density_sw/(-dlnawds)/Na;
    flu_con = S * rM0 * dnds**2 / density_sw / (-dlnawds) / Na
    #beta_cf = 2*pi*pi*((lambda*1e-9).^(-4)).*nsw.^2.*(flu_con)*(6+6*delta)/(6-7*delta);
    beta_cf = 2.0 * pi * pi * ((rlambda * 1.0E-9)**(-4)) * nsw**2 * (flu_con) \
    * (6.0 + 6.0 * delta) / (6.0 - 7.0 * delta);
    # total volume scattering at 90 degree
    beta90sw = beta_df + beta_cf;
    # bsw=8*pi/3*beta90sw*(2+delta)/(1+delta);
    bsw = 8.0 * pi / 3.0 * beta90sw * (2.0 + delta) / (1.0 + delta);

    return bsw

def getnswdnds(rlam,Tc,S):
    # refractive index of air is from Ciddor (1996,Applied Optics)

    rn_air = 1.0 + (5792105.0 / (238.0185 - 1.0 / (rlam / 1.0E3)**2) 
    +167917.0 / (57.362 - 1.0 / (rlam /1.0E3)**2)) / 1.0E8

    # refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    rn0 = 1.31405
    rn1 = 1.779E-4 
    rn2 = -1.05E-6 
    rn3 = 1.6E-8 
    rn4 = -2.02E-6 
    rn5 = 15.868
    rn6 = 0.01155
    rn7 = -0.00423
    rn8 = -4382.0 
    rn9 = 1.1455E6

    # pure seawater
    nsw = rn0 + (rn1 + rn2 * Tc + rn3 * Tc**2) * S + rn4 * Tc**2
    + (rn5 + rn6 * S + rn7 * Tc) / rlam + rn8 / rlam**2 + rn9 / rlam**3 

    nsw = nsw * rn_air

    dnswds = (rn1 + rn2 * Tc + rn3 * Tc**2 + rn6 / rlam) * rn_air

    return nsw, dnswds

def getIsoComp(Tc,S):
    # pure water secant bulk Millero (1980, Deep-sea Research)

    rkw = 19652.21 + 148.4206 * Tc - 2.327105 * Tc**2 + 1.360477E-2 * Tc**3
    -5.155288E-5 * Tc**4  

    # isothermal compressibility from Kell sound measurement in pure water
    # Btw = (50.88630+0.717582.*Tc+0.7819867E-3.*Tc.^2+31.62214E-6.*Tc.^3       ...
    #     -0.1323594E-6.*Tc.^4+0.634575E-9.*Tc.^5)./(1.0+21.65928E-3.*Tc).*1.0E-6;

    # seawater secant bulk
    a0 = 54.6746 - 0.603459 * Tc + 1.09987E-2 * Tc**2 - 6.167E-5 * Tc**3
    b0 = 7.944E-2 + 1.6483E-2 * Tc - 5.3009E-4 * Tc**2

    rKs = rkw + a0 * S + b0 * S**1.5

    # calculate seawater isothermal compressibility from the secant bulk
    IsoComp = 1.0 / rKs * 1.0E-5; # unit is pa
    return IsoComp

def get_density_sw(Tc,S):
    # density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981

    a0 = 8.24493E-1;
    a1 = -4.0899E-3;
    a2 = 7.6438E-5;
    a3 = -8.2467E-7;
    a4 = 5.3875E-9;
    a5 = -5.72466E-3;
    a6 = 1.0227E-4;
    a7 = -1.6546E-6;
    a8 = 4.8314E-4;
    b0 = 999.842594;
    b1 = 6.793952E-2;
    b2 = -9.09529E-3;
    b3 = 1.001685E-4;
    b4 = -1.120083E-6;
    b5 = 6.536332E-9;

    # density for pure water
    density_w = b0 + b1 * Tc + b2 * Tc**2 + b3 * Tc**3 + b4 * Tc**4 + b5 * Tc**5;
    
    # density for pure seawater
    density_sw = density_w \
    +((a0 + a1 * Tc + a2 * Tc**2 + a3 * Tc**3 + a4 * Tc**4) * S \
    +(a5 + a6 * Tc + a7 * Tc**2) * S**1.5 + a8 * S**2)

    return density_sw

def getdlnawds(Tc,S):
    # water activity data of seawater is from Millero and Leung (1976,American
    # Journal of Science,276,1035-1077). Table 19 was reproduced using
    # Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
    # dlnawds is partial derivative of natural logarithm of water activity
    # w.r.t.salinity

    dlnawds = (-5.58651E-4 + 2.40452E-7 * Tc - 3.12165E-9 * Tc**2 + 2.40808E-11 * Tc**3) \
    +1.5 * (1.79613E-5 - 9.9422E-8 * Tc + 2.08919E-9 * Tc**2 - 1.39872E-11 * Tc**3) * S**0.5 \
    +2.0 * (-2.31065E-6 - 1.37674E-9 * Tc - 1.93316E-11 * Tc**2) * S
    
    return dlnawds

def getDFRI(n_wat):

    n_wat2 = n_wat**2
    n_density_derivative = (n_wat2 - 1.0) * (1.0 + 2.0 / 3.0 * (n_wat2 + 2.0) * (n_wat / 3.0 - 1.0 / 3.0 / n_wat)**2)

    return n_density_derivative 

def eu0n(edtop,a,bb):
    # Computes surface normalized upwelling irradiance.  
    # A full treatment of diffuse and direct path lengths is required.  

    # Constants
      rmud = 1.0 ;                # avg cosine direct down
      rmuu = 1.0 / 0.4 ;            # avg cosine diffuse up
      bbw = 0.5 ;   # backscattering to forward scattering ratio
      rd = 1.5 ;  # these are taken from Ackleson, et al. 1994 (JGR)
      ru = 3.0;

    # Compute Eu
    # Eu from Ed
      ad = a * rmud
      bd = rd * bb * rmud
      au = a * rmuu
      bu = ru * bb * rmuu
      cd = ad + bd
      cu = au + bu
      bquad = cd - cu
      cquad = bd * bu - cd * cu
      sqarg = bquad * bquad - 4.0 * cquad
      a2 = 0.5 * ( -bquad - np.sqrt(sqarg))
      sfceun = (a2 + cd) / bu

      sfceu = edtop * sfceun;
      
      return sfceu



