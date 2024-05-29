import numpy as np

temp_offset = 21

"""========== General Equations =========="""
def Density(h, rho0=1.225,  constlambda=-0.0065, g0=9.80665, T0=288.15, R=287):

    return rho0*(1+(constlambda*h/T0))**(-(g0/(R*constlambda)+1))

def Power(Pto, h, rho0=1.225, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h) # temperature correction is for temperature offset
    Correction = Temp / (Temp + temp_offset)
    return Pto*(Correction*Density(h)/rho0)**(3/4)

def Cd(cdo, A, e, Cl):
    return cdo+Cl**2/(np.pi*A*e)


"""========== Flight Operations =========="""
def Stallspeedx(h, Vs, Clmax, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h) # temperature correction is for temperature offset
    Correction = Temp / (Temp + temp_offset)
    return 0.5*Correction*Density(h)*Vs**2*Clmax

def Takeoff(TOP, WS, h, ClmaxTO, rho0=1.225, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h)
    Correction = Temp / (Temp + temp_offset)
    constsigma=Correction*Density(h)/rho0
    Clto=ClmaxTO/1.21 # roskam pt1 p95
    return TOP*Clto*constsigma/WS

def Landingx(Clmax_land, h, sland, f, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h)
    Correction = Temp / (Temp + temp_offset)
    return (Clmax_land*Correction*Density(h)*sland/0.305199384478051392)/(2*f) # was 0.5915 changed to 0.305199384478051392

def Cruise(etap, h, Cd0, Vcruise, WS, A, e, cruisepowersetting, cruiseweightwrtMTOW, rho0=1.225, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h)
    Correction = Temp / (Temp + temp_offset)
    a=(cruisepowersetting/cruiseweightwrtMTOW)*etap*(Correction*Density(h)/rho0)**0.75
    b=((Cd0*0.5*Correction*Density(h)*Vcruise**3)/(cruiseweightwrtMTOW*WS)+(cruiseweightwrtMTOW*WS)/(np.pi*A*e*0.5*Correction*Density(h)*Vcruise))**(-1)
    return a*b

def Climbrate(etap, A, e, h, Cdo, climbrate, WS, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h)
    Correction = Temp / (Temp + temp_offset)
    a=climbrate+(np.sqrt(WS)*np.sqrt(2/(Correction*Density(h))))/(1.345*(A*e)**0.75/(Cdo**0.25))
    return etap/a

def Climbgradient(etap, WS, climbrate, V, A, e, Clmax, Clsafetyfactor, Cdo, h, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h)
    Correction = Temp / (Temp + temp_offset)
    return etap/(np.sqrt(WS)*(climbrate/V+Cd(Cdo, A, e, Clmax/Clsafetyfactor)/(Clmax/Clsafetyfactor))*np.sqrt(2/(Correction*Density(h)*(Clmax/Clsafetyfactor))))

def Manouvering(Cdo, h, V, WS, nmax, A, e, etap, T0=288.15, constlambda=-0.0065):
    Temp = (T0 + constlambda * h)
    Correction = Temp / (Temp + temp_offset)
    return ((Cdo*0.5*Correction*Density(h)*V**3/WS+WS*nmax**2/(np.pi*A*e*0.5*Correction*Density(h)*V))/etap)**(-1)