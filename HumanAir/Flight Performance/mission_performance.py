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
        return acf.takeoff_power_sealevel * acf.fuel_cons_TO
    
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

def take_off(acf, h=None, dT=None, surface="grass"):
    ground_distance = 0
    return ground_distance    
    
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
    ground_distance = 0
    return ground_distance
    

