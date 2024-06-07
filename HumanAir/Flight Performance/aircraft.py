import json
import numpy as np
from scipy.optimize import root_scalar

from helper import density
import thrust_power
import takeoff_landing

class Aircraft:
    """
    Aircraft object for evaluating performance.
    
    IMPORTANT: the thrust function stores intermediate values for faster
    calculations. If aircraft parameters change it should be ensured that
    these intermediate calculations are ereased.
    """
    
    def __init__(self, FILE="design.json"):
        with open("../Configurations/" + FILE,'r') as f:
            dat = json.load(f)
            
        # assumptions:
        # - variable pitch propeller (in takeoff run calculation)
            
        # settings
        self.CL_climb_safetyfactor = 1.2 # from ADSEE I lecture 3
        
        # performance
        self.V_cruise       = dat["Performance"]["Vc_m/s"]
        self.h_cruise       = dat["Performance"]["Altitude_Cruise_m"]
        self.h_cruise_max   = dat["Performance"]["Altitude_Max_Cruise_m"]
        self.dT_cruise      = dat["Performance"]["Temp_offset_TO_Land_cruise"]

        # aerodynamic data
        self.CD0_clean      = dat["Aero"]["CD0"]
        self.AR             = dat["Aero"]["AR"]
        self.e_clean        = dat["Aero"]["e"]
        self.S              = dat["Aero"]["S_Wing"]
        self.b              = dat["Aero"]["b_Wing"]
        
        self.CLmax_clean    = dat["Aero"]["CLmax_clean"]
        self.CLmax_TO       = dat["Aero"]["CLmax_TO"]
        self.CLmax_land     = dat["Aero"]["CLmax_TO"] # TODO: check if still correct
        
        # TODO: fill in
        self.CL_ground_TO   = .5 # C_L while on level runway with flaps in TO position 
        self.CL_ground_land = .7 # C_L while on level runway with flaps in land position
        self.number_of_engines = 1 # TODO: use design.json
        self.propeller_diameter = 2 # TODO: use deisng.json
        self.spinner_diameter = 0.3 # TODO: use deisng.json

        self.wing_h_above_ground = dat["Geometry"]["fus_height_m"] + \
            dat["Landing_gear"]["Hs_m"] + dat["Landing_gear"]["Dwm_m"]/2
        
        self.retractable_gear = dat["Landing_gear"]["Retractable"]

        # weights
        self.W_OE           = dat["Weights"]["OEW_N"]
        self.W_F            = dat["Weights"]["WF_N"]
        self.W_bat          = dat["Weights"]["Wbat_N"]
        self.W_pl_des       = dat["Weights"]["Wpl_des_kg"] * 9.80655
        self.W_pl_max       = dat["Weights"]["Wpl_max_kg"] * 9.80655
        self.W_MTO          = dat["Weights"]["MTOW_N"]
        
        x_cg = dat["Stability"]["Cg_Front"]*dat["Aero"]["MAC_wing"] + dat["Geometry"]["XLEMAC_m"]
        x_main, x_nose = dat["Landing_gear"]["Xmw_m"], dat["Landing_gear"]["Xnw_m"]
        self.weight_on_MLG  = (x_cg - x_nose) / (x_main - x_nose) # fraction of weight on MLG with front most possible CG position
        print(self.weight_on_MLG)
        # propulsion
        # TODO: correct
        self.max_cont_power_sealevel  = 250500 # engine power without losses 
        # TODO: correct
        self.takeoff_power_sealevel   = 250500 # engine power without losses
        # TODO: check definition
        self.eff_powertrain = dat["Power_prop"]["eta_powertrain"]
        self.max_eff_prop       = dat["Power_prop"]["eta_p"] # max eff prop
        
        # calculated data
        self.LDmax = 0.5 * np.sqrt(np.pi*self.AR*self.e_clean/self.CD0_clean) # ruijgrok p107
        self.CD_LDmax = 2 * self.CD0_clean # ruijgrok p106
        self.CL_LDmax = self.CD_LDmax * self.LDmax
        
        self.L3D2max = 3 * np.sqrt(3) / 16 * np.pi * self.AR * self.e_clean \
            * np.sqrt(np.pi * self.AR * self.e_clean / self.CD0_clean) # for P_r_min for max climbrate; ruijgrok p107
        self.CD_L3D2max = 4 * self.CD0_clean # ruijgrok p107
        self.CL_L3D2max = (self.L3D2max * self.CD_L3D2max**2)**(1/3)
                
    def CD(self, CL, gear="up", flaps="up", ground_effect_h=None): # this parabolic approximation is about accurate until CL=1 (ruijgrok p226)
        """ Parabolic estimation of C_D for given C_L for given gear/flaps conditions (default=clean) """
        
        # TODO: replace by more accurate representation for high-alpha CD estimation
        
        CD0 = self.CD0_clean
        e = self.e_clean
        
        if gear == "down" and self.retractable_gear:
            CD0 += 0.011 # p264 of Nicolai: Fundamentals of Aircraft and Airship Design: Volume 1
            
        if flaps == "TO":
            CD0 += 0.01 # roskam pt1 p127 and nicolai p242
            e -= 0.05 # roskam pt1 p127
            
        elif flaps == "land":
            # TODO: improve
            CD0 += 0.03 # guess based on nicolai p242
            
        # ground effect reduces induced drag, using nicolai section 10.2
        # (credit to Wieselberger for the equation), according to gudmundson
        # a good approximation for 0.033 < h_b < 0.25
        sigma = 0
        if ground_effect_h is not None:
            h_b = (self.wing_h_above_ground + ground_effect_h) / self.b
            sigma = (1 - 1.32 * h_b)/(1.05 + 7.4*h_b)
            sigma = max(sigma,0) # for around h_b > .75 sigma will become negative
            
        return CD0 + CL**2 / (np.pi * self.AR * e) * (1-sigma)
    
    def P_shaft(self, h, dT, use_takeoff_power=False):
        return thrust_power.P_shaft(self, h, dT, use_takeoff_power=use_takeoff_power)
    
    def P_a(self, h, dT, use_takeoff_power=False, V=None):
        return thrust_power.P_a(self, h, dT, use_takeoff_power=use_takeoff_power, V=V)
        
    def T(self, V_ms, h, dT, use_takeoff_power=False):
        return thrust_power.T(self, V_ms, h, dT, use_takeoff_power=use_takeoff_power)
    
    def prop_eff(acf, V, h, dT, use_takeoff_power=False):
        return thrust_power.prop_eff(acf, V, h, dT, use_takeoff_power=use_takeoff_power)
        
    def V_Dmin(self, W, h, dT):
        """ Calculates velocity corresponding with minimum drag in clean configuration """
        a = 1 / (self.CD0_clean * np.pi * self.AR * self.e_clean)
        return np.sqrt( W/self.S * 2/density(h, dT) * a )
    
    def D_min(self, W):
        """ Minimum drag (in Newtons) for given weight in clean configuration. Corresponds to V_Dmin """
        return W / self.LDmax
    

        
    def P_r_min(self, W, h, dT):
        """ Minimum possible value for  power required for given conditions in clean configuration,
            note that this assumes a small flight path angle. """
        return W * np.sqrt(W/self.S * 2/density(h, dT) * 1/self.L3D2max) # ruijgrok p221
        
    def RC_max(self, W, h, dT):
        """ Maximum possible climb rate for given conditions in clean configuration"""
        # TODO: take into account effect of low airspeed on thrust or P_a ?
        return (self.P_a(h, dT) - self.P_r_min(W, h, dT)) / W
    
    def climb_angle_max(self, W, h, dT, gear="up", flaps="up"):
        """ Maximum possible climb slope (RC/V) for given conditions in given config (default clean) """
        # TODO: take into account effect of low airspeed on thrust or P_a ?
        # adsee I q3 lecture 3 slide 93
        if flaps == "up":
            CL = self.CLmax_clean / self.CL_climb_safetyfactor
        elif flaps == "TO":
            CL = self.CLmax_TO / self.CL_climb_safetyfactor
        elif flaps == "land":
            CL = self.CLmax_land / self.CL_climb_safetyfactor
        return self.P_a(h, dT)/W * 1/np.sqrt(W/self.S * 2/density(h, dT) * 1/CL) - self.CD(CL, gear=gear, flaps=flaps)/CL
    
    def stall_speed(self, W, h, dT, flaps="up"):
        """ Calculates stall speed for given conditions in given config (default clean) """
        if flaps == "up":
            CL = self.CLmax_clean
        elif flaps == "TO":
            CL = self.CLmax_TO
        elif flaps == "land":
            CL = self.CLmax_land
        return np.sqrt(W / (0.5 * density(h, dT) * self.S * CL))
    
    def V_max(self, W, h, dT):
        """
        Calculates the maximum true airspeed in clean configuration.
        
        It does so by finding the point where P_a() equals the drag force times
        velocity.
                
        Note that the maximum possible P_a() is used. this means that the
        calculated value is inaccurate for high altitudes.
        
        TODO: fix this
        
        This is because the calculation of power available for low speeds
        needs V_max, causing an infinite loop.
        

        Parameters
        ----------
        W : float
            Aircraft weight [N].
        h : float
            Geopotential altitude [m].
        dT : float
            ISA temperature offset [deg C].

        Returns
        -------
        float
            The maximum true airspeed in m/s.
        """

        rho = density(h, dT)
        def diff(V): # returns difference between power available and power required
            CL = W / (0.5 * rho * self.S * V**2)
            CD = self.CD(CL)
            D = 0.5 * rho * self.S * CD * V**2
            return self.P_a(h, dT) - D*V
        return root_scalar(diff, x0=75).root
    
    def power_setting(self, W, V, h, dT):
        """ Gives power setting required to maintain cruise (0-1, but also goes beyond 1) """
        CL = W / (0.5 * density(h, dT) * self.S * V**2)
        D = 0.5 * density(h, dT) * self.S * self.CD(CL) * V**2
        if CL > 1: return 0
        return V * D / self.P_a(h, dT)
    
    
    #
    # Take-off and landing
    #

    def takeoff_ground_run(self, W, h, dT, slope, surface):
        return takeoff_landing.takeoff_ground_run(self, W, h, dT, slope, surface)
       
    def landing_ground_distance(self, W, h, dT, slope, surface):
        return takeoff_landing.landing_ground_distance(self, W, h, dT, slope, surface)