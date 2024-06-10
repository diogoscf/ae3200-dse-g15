import numpy as np
from scipy.optimize import root_scalar, minimize_scalar

from helper import density
import thrust_power
import takeoff_landing 

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import aircraft_data

# old span: 36.27269326905173
#        "Fuel_engine_P_TO_W": 368000,
#        "Fuel_engine_P_max_cont_W": 338000


class Aircraft:
    """
    Aircraft object for evaluating performance.
    
    Currently it only uses only the fuel engine.
    
    IMPORTANT: the thrust function stores intermediate values for faster
    calculations. If aircraft parameters change it should be ensured that
    these intermediate calculations are ereased.
    """
    
    def __init__(self, FILE="design.json"):
        dat = aircraft_data.aircraft_data
            
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
        self.CLmax_land     = dat["Aero"]["CLmax_Land"]
        
        self.CL_ground_TO   = dat["Flaps"]["CL_AoA0_landing"] 
        self.CL_ground_land = dat["Flaps"]["CL_AoA0_takeoff"]
        self.number_of_engines = 1
        self.propeller_diameter = dat["Power_prop"]["Dp_m"]
        self.spinner_diameter = 0.3

        self.wing_h_above_ground = dat["Geometry"]["fus_height_m"] + \
            dat["Landing_gear"]["Hs_m"] + dat["Landing_gear"]["Dwm_m"]/2
        
        self.retractable_gear = dat["Landing_gear"]["Retractable"]

        # weights
        self.W_OE           = dat["CL2Weight"]["OEW"]
        self.W_MF           = dat["CL2Weight"]["Wfuel_N"] # max fuel weight
        self.W_bat          = dat["CL2Weight"]["Wbat_N"]
        self.W_pl_no_pilot  = dat["CL2Weight"]["Wpl_w/o_pilot"] # class ii
        self.W_pl_des       = dat["Weights"]["Wpl_des_kg"] * 9.80655 # class i
        self.W_pl_max       = dat["Weights"]["Wpl_max_kg"] * 9.80655 # class i
        self.W_MTO          = dat["CL2Weight"]["MTOW_N"] #- self.W_pl_max/2 - self.W_MF / 2
        
        x_cg = dat["Stability"]["Cg_Front"]*dat["Aero"]["MAC_wing"] + dat["Geometry"]["XLEMAC_m"]
        x_main, x_nose = dat["Landing_gear"]["Xmw_m"], dat["Landing_gear"]["Xnw_m"]
        self.weight_on_MLG  = (x_cg - x_nose) / (x_main - x_nose) # fraction of weight on MLG with front most possible CG position
        
        # propulsion
        self.max_cont_power_sealevel  = dat["Power_prop"]["Fuel_engine_P_max_cont_W"]
        self.takeoff_power_sealevel   = dat["Power_prop"]["Fuel_engine_P_TO_W"]
        self.eff_powertrain     = dat["Power_prop"]["eta_powertrain"]
        self.max_eff_prop       = dat["Power_prop"]["eta_p"]
        
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
        
        # TODO: replace by more accurate representation for high-alpha CD and flap estimation
        
        CD0 = self.CD0_clean
        e = self.e_clean
        
        if gear == "down" and self.retractable_gear:
            CD0 += 0.011 # p264 of Nicolai: Fundamentals of Aircraft and Airship Design: Volume 1
            
        if flaps == "TO":
            CD0 += 0.003 # guess based on figure 9.25 nicolai p242
            
        elif flaps == "land":
            CD0 += 0.0175 # guess based on figure 9.25 nicolai p242
            
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
        a = 1 / np.sqrt(self.CD0_clean * np.pi * self.AR * self.e_clean)
        b = np.sqrt( W/self.S * 2/density(h, dT) * a ) # ruijgrok p224, same as V at (L/D)max
        c = np.sqrt( W / (0.5 * density(h, dT) * self.S * self.CL_LDmax))
        if not np.isclose(b,c):
            print(f"V_D min error: {b:.2f} - {c:.2f}")
        return b
    
    def D_min(self, W):
        """ Minimum drag (in Newtons) for given weight in clean configuration. Corresponds to V_Dmin """
        return W / self.LDmax # ruijgrok p220
        
    def P_r_min(self, W, h, dT):
        """ Minimum possible value for  power required for given conditions in clean configuration,
            note that this assumes a small flight path angle. """
        return W * np.sqrt(W/self.S * 2/density(h, dT) * 1/self.L3D2max) # ruijgrok p221
        
    def RC_max(self, W, h, dT):
        """
        Calculates maximum possible climb rate for clean configuration in given
        condition. Assumes constant V.
        
        Uses scipy.optimize.minimize_scalar to find the velocity with maximum
        climb rate (to take into account non-constant power available.)

        Parameters
        ----------
        W : float
            Aircraft gross weight [N].
        h : float
            Geopotential altitude [m].
        dT : TYPE
            ISA temperature offset [deg C].

        Returns
        -------
        float
            Maximum climb rate [m/s].

        """
        rho = density(h, dT)
        
        def RC(V): # calculates climbrate *-1
            CL = W / (0.5 * rho * V**2 * self.S)
            CD = self.CD(CL)
            P_r = V * (0.5 * rho * V**2 * self.S * CD)
            P_a = self.P_a(h, dT, V=V)
            return -(P_a - P_r) / W # minus so that minimize() will maximize
        
        bounds = [self.stall_speed(W, h, dT), self.V_max(W, h, dT)]
        if bounds[0] > bounds[1]:
            return 0
        
        res = minimize_scalar(RC, bounds=bounds, method="bounded") # minus to convert back to positive
        rc_max = -res.fun
        
        #
        # alternate method, does not take into account varying prop efficiency,
        # result is almost the same
        #
        #rc2 = (self.P_a(h, dT) - self.P_r_min(W, h, dT)) / W
        #print(f"RC {res.x:.2f} - {rc_max:.2f} - {rc2:.2f}")
        
        return rc_max
    
    def climb_slope_max(self, W, h, dT, gear="up", flaps="up"):
        """
        Calculates maximum possuble climb slope for given conditions, assuming
        constant V.
        
        Takes into account lower power available at low speeds, and since the
        max climb slope is achieved at high CL (=low speeds) this method gives
        signifiantly lower results than when constant P_a is assumed.

        Parameters
        ----------
        W : float
            Aircraft gross weight.
        h : float
            Altitude.
        dT : float
            ISA temperature offset.
        gear : string, optional
            "up" or "down". The default is "up".
        flaps : float, optional
            "up", "TO", or "landing. The default is "up".

        Returns
        -------
        float
            Maximum climb slope ranging 0-100% [no unit].

        """
        rho = density(h, dT)
        
        def slope(V): # calculates slope *-1
            CL = W / (0.5 * rho * V**2 * self.S)
            CD = self.CD(CL, gear=gear, flaps=flaps)
            D = 0.5 * rho * V**2 * self.S * CD
            T = self.T(V, h, dT)
            return -(T - D) / W # minus so that minimize() will maximize
        
        bounds = [self.stall_speed(W, h, dT), self.V_max(W, h, dT)]
        if bounds[0] > bounds[1]:
            return 0
        
        res = minimize_scalar(slope, bounds=bounds, method="bounded") # minus to convert back to positive
        slope_max = np.tan(np.arcsin(-res.fun)) *100 # result is sin(gamma), we need tan(gamma)

        #
        # alternate method for testing (does not take into account varying prop eff)
        #
        
        # if flaps == "up":
        #     CL = self.CLmax_clean / self.CL_climb_safetyfactor
        # elif flaps == "TO":
        #     CL = self.CLmax_TO / self.CL_climb_safetyfactor
        # elif flaps == "land":
        #     CL = self.CLmax_land / self.CL_climb_safetyfactor
        # slope_max2 = 100* self.P_a(h, dT)/W * 1/np.sqrt(W/self.S * 2/density(h, dT) * 1/CL) - self.CD(CL, gear=gear, flaps=flaps)/CL
        # print(f"SL {res.x:.2f} - {slope_max:.2f} - {slope_max2:.2f}")
        
        return slope_max
    
    def stall_speed(self, W, h, dT, flaps="up"):
        """
        Calcualtes stall speed for given conditions in given config (default
        is clean).

        Parameters
        ----------
        W : float
            Aircraft gross weight [N].
        h : float
            Geopotential altitude [m].
        dT : float
            ISA temperature offset [deg C].
        flaps : string, optional
            "up", "TO", or "landing". The default is "up".

        Returns
        -------
        float
            Stall speed [m/s].

        """
        if flaps == "up":
            CL = self.CLmax_clean
        elif flaps == "TO":
            CL = self.CLmax_TO
        elif flaps == "land":
            CL = self.CLmax_land
        
        return np.sqrt(W / (0.5 * density(h, dT) * self.S * CL))
    
    def V_max(self, W, h, dT, use_max_prop_eff=False):
        """
        Calculates the maximum true airspeed in clean configuration.
        
        It does so by finding the point where P_a() equals the drag force times
        velocity.
                
        When use_max_prop_eff is set to true it will not take into account
        low-speed losses. This is to prevent an infinite loop when used by
        the thrust function.   
        
        TODO: root_scalar method employed is not very reliable at high altitudes

        Parameters
        ----------
        W : float
            Aircraft weight [N].
        h : float
            Geopotential altitude [m].
        dT : float
            ISA temperature offset [deg C].
        use_max_prop_efficiency : boolean
            Whether to assume max prop efficiency

        Returns
        -------
        float
            The maximum true airspeed in m/s.
        """
        rho = density(h, dT)
        
        def diff(V): # returns difference between power available and power required
            if type(V) is not np.float64:
                V = V[0]
            CL = W / (0.5 * rho * self.S * V**2)
            CD = self.CD(CL)
            D = 0.5 * rho * self.S * CD * V**2
            if use_max_prop_eff:
                return self.P_a(h, dT) - D*V
            else:
                return self.P_a(h, dT, V=V) - D*V
            
        return root_scalar(diff, x0=150).root
    
    def power_setting(self, V, W, h, dT):
        """
        Returns power setting (0-1 and >1 if impossible speed) to maintain 
        steady horizontal flight at given speed in clean configuration.
        
        Returns 0 if C_L > 1 because with current drag estimation CD is
        inaccurate for C_L>1.
        
        TODO: change this when CD() has been updated.

        Parameters
        ----------
        V : float
            True airpseed [m/s].
        W : float
            Aircraft gross weight [N].
        h : float
            Altitude [m].
        dT : float
            ISA temperature offset [deg C].

        Returns
        -------
        float
            Power setting [0-1 and >1 if impossible].

        """
        CL = W / (0.5 * density(h, dT) * V**2 * self.S)
        D = 0.5 * density(h, dT) * V**2 * self.S * self.CD(CL)
        
        if CL > 1: return 0
        
        return V * D / self.P_a(h, dT, V=V)
    
    
    #
    # Take-off and landing
    #

    def takeoff_ground_run(self, W, h, dT, slope, surface):
        return takeoff_landing.takeoff_ground_run(self, W, h, dT, slope, surface)
       
    def landing_ground_distance(self, W, h, dT, slope, surface, reversible_pitch=False):
        return takeoff_landing.landing_ground_distance(self, W, h, dT, slope, surface, reversible_pitch=reversible_pitch)