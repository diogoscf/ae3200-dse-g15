import numpy as np
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import unit_conversions as conv
import aircraft
from helper import density

def P_shaft(acf, h, dT, use_takeoff_power=False):
    """
    Returns the maximum available shaft power (therefore excluding propeller
    efficiency) for either continuous operation or short operation.
    
    It is assumed that shaft power decreases with air density but is constant
    with velocity. The decrease with altitude is computed using the Gagg and
    Ferrar model (sourced from Gudmundseo).
    
    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    h : float
        The geopotential altitude [m].
    dT : float
        ISA temperature offset [deg C].
    use_takeoff_power : boolean, optional
        Whether to use takeoff power. The default is False, which will result
        in the maximum continuous power to be returned.

    Returns
    -------
    float
        The maximum available shaft power [W].

    """
    
    alt_correction = 1.132 * (density(h, dT)/1.225) - 0.132 # Gagg and Ferrar model, Gudmundsen eq 7-16
   
    if use_takeoff_power:
        return acf.takeoff_power_sealevel * acf.eff_powertrain * alt_correction
    else:
        return acf.max_cont_power_sealevel * acf.eff_powertrain * alt_correction

def P_a(acf, h, dT, use_takeoff_power=False, V=None):
    """
    Max continuous power available for given conditions. If no airspeed is
    given then the maximum propeller efficiency is used instead of the speed
    dependent one.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    h : float
        The geopotential altitude [m].
    dT : float
        The ISA temperature offset.
    use_takeoff_power : boolean, optional
        Wether to use takeoff power instead of max continuous power. The
        default is False.
    V : float, optional
        Airspeed. The default is None, which results in max prop efficiency
        being used.

    Returns
    -------
    float
        Power available [W].

    """

    if V == None:
        return acf.max_eff_prop * P_shaft(acf, h, dT, use_takeoff_power=use_takeoff_power)
    else:
        return acf.prop_eff(V, h, dT) * P_shaft(acf, h, dT, use_takeoff_power=use_takeoff_power)

    """
    TODO: RC_max() assumes constant P_a with velocity, if this changes here
    then the implementation of that function must be changed too """

 
def T(acf, V_ms, h, dT, use_takeoff_power=False):
    """
    Approximates the thrust for given conditions. It uses the cubic spline
    interpolation as derived in Gudmundsun section 14.4.2 (method #3).

    https://www.sciencedirect.com/science/article/pii/B9780123973085000143#s0235
                                                           
    The main purpose of this function is for calculating thrust at low
    airspeeds; for cruise speed and beyond one can simply use T=P/V*eff_prop.
                                                           
    Basically it interpolates between the static thrust at zero airspeed,
    P/V*Max_prop_efficiency at cruise and P/V*Max_prop_efficiency at maximum
    airspeed. This is done since the efficiency for constant-speed propellers
    can be assumed to be constant for high air speeds. The limitation is that
    we do not know when the maximum efficiency is reached: that specific speed
    may be lower than the cruise speed, giving an underestimated thrust for
    V < V_cruise.
    
    TODO: check with prop specs when max efficiency is reached.
    
    This function will throw an exception for V>V_max as beyond the maximum
    speed the interpolation is no longer valid (see Gudmundson fig 14-42).
    
    For calculating V_max it uses the OEW.
    TODO: think over this decision again

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    V : float
        True airspeed [m/s].
    h : float
        Geopotential altitude [m].
    dT : float
        ISA temperature offset.
    use_takeoff_power : boolean, optional
        Whether to use takeoff power instead of maximum continuous power. The
        default is False.

    Returns
    -------
    Thrust [N].

    """
    
    rho = density(h, dT)
    
    # max prop efficiency
    eff_prop_max = acf.max_eff_prop
    
    # get shaft power
    P_sh_W = P_shaft(acf, h, dT, use_takeoff_power)
    P_sh_hp = conv.W_to_hp(P_sh_W)
    
    # calculate static thrust
    A_prop = np.pi * (acf.propeller_diameter/2)**2
    A_spinner = np.pi * (acf.spinner_diameter/2)**2
    T_static_N = 0.85 * P_sh_W**(2/3) * (2 * rho * A_prop)**(1/3) * (1 - A_spinner/A_prop) # Gudmundson eq 14-64
    T_static_lb = conv.N_to_lbs(T_static_N)
    
    # cruise speed
    V_C_ms = acf.V_cruise
    V_C_kt = conv.m_s_to_kt(V_C_ms)
    
    # max speed
    V_H_ms = acf.V_max(acf.W_OE, h, dT)
    V_H_kt = conv.m_s_to_kt(V_H_ms)
    
    if V_ms > V_H_ms:
        raise Exception("Velocity given to thrust function exceeds max speed, \
                        this is out of bounds.")
    
    # cruise thrust
    T_cruise_N = P_sh_W * eff_prop_max / V_C_ms
    T_cruise_lb = conv.N_to_lbs(T_cruise_N)

    # max speed thrust
    T_Vmax_N = P_sh_W * eff_prop_max / V_H_ms
    T_Vmax_lb = conv.N_to_lbs(T_Vmax_N)
    
    # speed in knots for which thrust need to be calculated
    V_kt = conv.m_s_to_kt(V_ms)
    
    
    #
    # now performing the actual calculations according to Gudmundson eq 14-41 and eq 14-42
    #
    
    # matrix
    MTX = [[0,           0,         0,         1],
           [V_C_kt**3,   V_C_kt**2, V_C_kt,    1],
           [3*V_C_kt**2, 2*V_C_kt,  1,         0],
           [V_H_kt**3,   V_H_kt**2, V_H_kt**2, 1]]
    
    # vecotr
    VCT = [T_static_lb, T_cruise_lb, -eff_prop_max*325.8*P_sh_hp/V_C_kt**2, T_Vmax_lb]
    
    # first multiplication
    result = np.dot(np.linalg.inv(MTX), VCT)
    
    # velocity vector
    VCT_2 = [V_kt**3, V_kt**2, V_kt, 1]
    
    print(f"Static thrust: {T_static_N:.2f}\n \
            Cruise thrust: {T_cruise_N:.2f}\n \
            V_max thrust:  {T_Vmax_N:.2f}")

    return conv.lbs_to_N(np.dot(result, VCT_2))

def prop_eff(acf, V, h, dT, use_takeoff_power=False):
    return T(acf, V, h, dT, use_takeoff_power=use_takeoff_power)*V/P_shaft(acf, h, dT, use_takeoff_power=use_takeoff_power)
    