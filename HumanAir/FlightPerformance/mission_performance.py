# from engine specs:
# T/O fuel consumption 0.17 kg/hp/hr -> 0.34532749982929480304 N/W/s
# averaged non-TO fuel consumption 0.16 kg/hp -> 0.3250141174863951048 N/W/s

def fuel_rate(acf, P_shaft=None, TO=False):
    """
    Calculates the fuel consumed for given power required in N/s.
    
    When 'TO' is set to 'True' this function disregards 'P_shaft' and assumes
    max takeoff power.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    P_shaft : float
        Shaft power required [W].
    TO  : boolean
        Whether to use T/O power.

    Returns
    -------
    float
        Fuel consumption rate [N/s].

    """
    if TO == False:
        return P_shaft / acf.eff_powertrain * acf.fuel_cons_flight
    else:
        return acf.electric_takeoff_power * acf.fuel_cons_TO
    
    
def bat_cap_rate(acf, P_shaft=None, TO=False):
    """
    Calculates the rate at which battery capacity is depleted in Wh/s.
    
    When 'TO' is set to 'True' this function disregards P_shaft and assumes
    max takeoff power.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    P_shaft : float
        Shaft power required [W].
    TO  : boolean
        Whether to use T/O power.

    Returns
    -------
    float
        Fuel consumption rate [N/s].

    """
    if TO == False:
        return P_shaft / acf.eff_powertrain / acf.eff_electric_motor / acf.eff_battery
    else:
        return acf.electric_takeoff_power/ acf.eff_powertrain / acf.eff_electric_motor / acf.eff_battery
    

def take_off(acf, h=None, dT=None, surface="grass", electric=False):
    if h == None and dT == None:
        h = acf.h_TO
        dT = acf.dT_default
        
    dt, acf.TAS = acf.takeoff_ground_run(acf.W_current, h, dT, 0, surface, electric=electric, calc_time=True)
    
    if electric:
        acf.bat_cap -= acf.bat_cap_rate(TO=True) * dt
    else:
        acf.fuel -= acf.fuel_rate(TO=True) * dt
        
    acf.alt = 0 # neglect what little alt has been gained
    
    
def climb(acf, h=None, dT=None, rate=2.5):
    ground_distance = 0
    return ground_distance
    
def cruise(acf, ground_distance, h=None, dT=None):
    # TODO: check if flying at max CL/CD is possible
    foo = 1
    
def loiter(acf, t, h=None, dT=None):
    foo = 2
    
def descent(acf, h=None, dT=None, rate=2.5):
    ground_distance = 0
    return ground_distance
    
def land(acf, h=None, dT=None, surface="grass"):
    

