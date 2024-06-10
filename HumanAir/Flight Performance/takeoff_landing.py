import numpy as np
from scipy.integrate import quad

from helper import density

def takeoff_ground_run(acf, W, h, dT, slope, surface):
    """
    Calculates the take-off ground run for given conditions using the Ruijgrok
    method (p372-376).
    
    Does not take rotation phase into account, simply calculates distance from
    V=0 to V=lift off speed.
    
    It is inaccurate for high postive slopes and will throw an exception when
    the integration error estimate exceeds 0.001 m.
        
    Parameters
    ----------
    W : float
        Gross aircraft weight at takeoff [N].
    h : float
        Elevation of airstrip [m].
    dT : float
        ISA temperature offset at airstrip [deg C].
    slope : float
        Slope of runway [%] (vertical/horizontal*100), 
        positive = sloped upwards at takeoff.
    surface : string
        Either "grass" or "paved"

    Returns
    -------
    s_run : float
        Ground run [m]
    accuracy : float
        Integration absolute error estimate [m]

    """
    
    # TODO: multiply by safety factor to ensure takeoff not at very end of runway?
    
    # TODO: ground effect currently only accounted for by reduced induced drag
    # there is a good discussion on ground effect in  
    # Fundamentals of Aircraft and Airship Design: Volume 1 by nicolai
    # pages 257-260

    rho = density(h, dT)
    g = 9.80665
    
    CL_ground = acf.CL_ground_TO # C_L during ground run
    CD_ground = acf.CD(CL_ground, gear="down", flaps="TO", ground_effect_h=20) # C_D during ground run
    
    
    # Torenbeek estimates V_LOF to be 1.2*stall speed in T/O 
    # configuration. This is a conservative estimate since CS23 requires the
    # rotation speed to be >=V_S1 or >=1.05*V_MC and the speed at 15m >=1.2*V_S1
    V_S_TO = np.sqrt(W / (0.5 * rho * acf.S * acf.CLmax_TO))
    V_LOF = max(1.2 * V_S_TO, 0)#35) # 35 at 750+18 was chosen by control department as minimum
    # TODO: check if 35 is really required
    xi = np.arctan(slope/100) # angle of slope
    
    if surface == "paved":
        mu_r = 0.02 # concrete, ruijgrok p372, see also gudmundson table 17-3
    elif surface == "grass":
        mu_r = 0.05 # short cut grass, ruijgrok p372
    
    def f(V):
        # eq 16.2-6 with acceleration from eq 16.2-9
        T = acf.T(V, h, dT)
        D = 0.5 * rho * V**2 * acf.S * CD_ground
        L = 0.5 * rho * V**2 * acf.S * CL_ground
        return V / (g/W * (T - D - mu_r*(W*np.cos(xi) - L) - W*np.sin(xi)))
    
    s_run, accuracy = quad(f, 0, V_LOF) # integrate from V=0 to V=V_LOF
    
    if accuracy > 10**-3:
        raise Exception(f"Low accuracy in integrating takeoff distance: \
                        s_run = {s_run}, accuracy = {accuracy}. \
                        This is likely due to a high positive slope. Slope={slope}%")
    #
    # Torenbeek method, torenbeek p167-168
    # Assumes constant speed propeller for Torenbeek method.
    # Note that Torenbeek uses kg for weights, not newton
    # however it really is not accurate enough, it was found that
    # T_bar underestimates avg thrust (for unknown reason, could be a mistake
    # in units on my part) causing much greater landing distances. Also the
    # drag approximation mu_quote is not as accurate as calculating drag force
    # manually.
    #
    
    # P_to = acf.P_a(h, dT, use_takeoff_power=True)
    # sigma = rho/1.225
    # CD0 = acf.CD(0, gear="down", flaps="TO")
    # P_to_kg = P_to / g # formula requires kgm/s
    # W_kg = W / g # torenbeek uses kg
    
    # # adjusted friction coefficient
    # mu_quote = 0.05 + .72 * CD0/CL_max # 0.04-0.05 for short grass and 0.02 for concrete, if a fixed pitch prop is used use 0.576
    
    # # avg thrust in kg
    # Tbar_kg = acf.T(0.707 * V_LOF, h, dT, use_takeoff_power=True) / g #0.321 * P_to_kg * (sigma * acf.number_of_engines * acf.propeller_diameter**2 / P_to_kg)**(1/3)
    
    # # now calculate ground run
    # s_run2 = V_LOF**2 / (2*g) / (Tbar_kg/W_kg - mu_quote)
        
    #
    # Nicolai method, calculate acceleration at average speed 0.707V_LOF
    # and assume this remains constant.
    # from Fundamental of Aircraft and Airship design: volume 1 pages 257-260
    # it was found that this method is very close to the ruijgrok integration
    # when the speed-adjusted thrust is used. When T = T_to/V is used it
    # overestimates avg thrust by quite a bit.
    #
    
    # V = 0.707 * V_LOF # 0.707 to get the average velocity (when assuming const acceleration)
    # T = acf.T(V, h, dT, use_takeoff_power=True) # P_to / V
    # D = 0.5 * rho * V**2 * acf.S * CD_ground
    # L = 0.5 * rho * V**2 * acf.S * CL_ground
    # a = g/W * (T - D - 0.05*(W - L)) # avg acceleration - nicolai has no friction coefficient for dry grass so using ruijgrok's value
    # s_run3 = 0.5 * V_LOF**2 / a
    
    
    return s_run, accuracy

    
def landing_ground_distance(acf, W, h, dT, slope, surface, reversible_pitch=False):
    """
    Estimates landing ground run distance using gudmundsen method (chapter 22).
    
    Assumes one second of free roll after touchdown before braking. 
    
    Raises an exception when estimateed absolute error of integration drops 
    below 0.001.

    Parameters
    ----------
    acf : Aircraft
        The aircraft object.
    W : float
        Aircraft gross weight [N].
    h : float
        Geopotential altitude [m].
    dT : float
        ISA temperature offset.
    slope : float
        Slope of runway, positive = upwards, in %.
    surface : string
        "grass" or "paved".

    Returns
    -------
    float
        Ground run distance [m].

    """
            
    rho = density(h, dT)
    g = 9.80665
    
    V_T = 1.1 * np.sqrt(W/( 0.5 * rho * acf.S * acf.CLmax_land)) # touchdown velocity 

    CL = acf.CL_ground_land # TODO: take into account ground effect
    CD = acf.CD(CL, gear="down", flaps="land", ground_effect_h=0)
    
    xi = np.arctan(slope/100) # angle of slope

    if surface == "paved":
        mu = 0.4
    elif surface == "grass":
        mu = 0.2 # wet grass, gudmundsen provides no value for dry grass, = conservative
        
    if reversible_pitch:
        T = -0.4 * acf.T(0, h, dT) # -40% of static thrust according to gudmundsen
    else:
        T = 0.07 * acf.T(0, h, dT) # -+7% of static thrust according to gudmundsen
    
    V_avg = V_T / np.sqrt(2)   
    L = 0.5 * rho * V_avg**2 * acf.S * CL
    D = 0.5 * rho * V_avg**2 * acf.S * CD
    
    D_g = mu*(W-L)*acf.weight_on_MLG # assumes zero pitching moment

    s_ground = - V_T**2 * W / (2 * g * (T - D - D_g - W * np.sin(xi))) # gudmundsen eq 22-2 and 22-13
    
    s_ground += 3* V_T # add 2 seconds of overflight, 1 second of free roll at touchdown speed
    
   #  def f(V):
   #      L = 0.5 * rho * V**2 * acf.S * CL
   #      D = 0.5 * rho * V**2 * acf.S * CD
   #      D_g = mu*(W-L)*acf.weight_on_MLG # assumes zero pitching moment
   #      return (W/g * V) / (T - D - D_g - W*np.sin(xi))
   #  s_ground2, accuracy = quad(f, V_T, 0)
   # # s_ground2 += V_T
    
   #  if accuracy > 10**-3:
   #      raise Exception(f"Low accuracy in integrating landing distance: \
   #                      s_ground = {s_ground}, accuracy = {accuracy}.")
        
    # TODO: verwijder "MAC" van design.json als deze niet nodig is want deze is dubbel
    
    return s_ground#, s_ground2
    