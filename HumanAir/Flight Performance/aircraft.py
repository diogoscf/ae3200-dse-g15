from helper import density
import json
import numpy as np
from scipy.optimize import root_scalar
#from scipy.integrate import quad

class Aircraft:
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
        
        self.CLmax_clean    = dat["Aero"]["CLmax_clean"]
        self.CLmax_TO       = dat["Aero"]["CLmax_TO"]
        self.CLmax_land     = dat["Aero"]["CLmax_TO"] # TODO: check if still correct
        
        #self.CL_ground      = .5 # C_L while on level runway TODO: fill in
        self.number_of_engines = 1 # TODO: use design.json
        self.propeller_diameter = 2 # TODO: use deisng.json

        # weights
        self.W_OE           = dat["Weights"]["OEW_N"]
        self.W_F            = dat["Weights"]["WF_N"]
        self.W_bat          = dat["Weights"]["Wbat_N"]
        self.W_pl_des       = dat["Weights"]["Wpl_des_kg"] * 9.80655
        self.W_pl_max       = dat["Weights"]["Wpl_max_kg"] * 9.80655
        self.W_MTO          = dat["Weights"]["MTOW_N"]

        # propulsion
        self.max_cont_power_sealevel  = 250500 # engine power without losses TODO: correct
        self.takeoff_power_sealevel   = 250500 # engine power without losses TODO: correct
        self.eff_powertrain = dat["Power_prop"]["eta_powertrain"] # TODO: check definition
        self.eff_prop       = dat["Power_prop"]["eta_p"]
        
        # calculated data
        self.LDmax = 0.5 * np.sqrt(np.pi*self.AR*self.e_clean/self.CD0_clean) # ruijgrok p107
        self.CD_LDmax = 2 * self.CD0_clean # ruijgrok p106
        self.CL_LDmax = self.CD_LDmax * self.LDmax
        
        self.L3D2max = 3 * np.sqrt(3) / 16 * np.pi * self.AR * self.e_clean \
            * np.sqrt(np.pi * self.AR * self.e_clean / self.CD0_clean) # for P_r_min for max climbrate; ruijgrok p107
        self.CD_L3D2max = 4 * self.CD0_clean # ruijgrok p107
        self.CL_L3D2max = (self.L3D2max * self.CD_L3D2max**2)**(1/3)
                
    def CD(self, CL, gear="up", flaps="up"): # this parabolic approximation is about accurate until CL=1 (ruijgrok p226)
        """ Parabolic estimation of C_D for given C_L for given gear/flaps conditions (default=clean) """
        
        # TODO: replace by more accurate representation for high-alpha CD estimation
        
        CD0 = self.CD0_clean
        e = self.e_clean
        
        if gear == "down":
            CD0 += 0.005 # TODO: estimate
            
        if flaps == "TO":
            CD0 += 0.005 # TODO: estimate
            e -= 0.05 # TODO: estimate
        elif flaps == "land":
            CD0 += 0.01 # TODO: estimate
            e -= 0.1 # TODO: estimate
            
        return CD0 + CL**2 / (np.pi * self.AR * e)
    
    def V_Dmin(self, W, h, dT):
        """ Calculates velocity corresponding with minimum drag in clean configuration """
        a = 1 / (self.CD0_clean * np.pi * self.AR * self.e_clean)
        return np.sqrt( W/self.S * 2/density(h, dT) * a )
    
    def D_min(self, W):
        """ Minimum drag (in Newtons) for given weight in clean configuration. Corresponds to V_Dmin """
        return W / self.LDmax
    
    def P_a(self, h, dT):
        """ Max continuous power for given conditions
        TODO: temporary placeholder implementation
        NOTE: RC_max() assumes constant P_a with velocity, if this changes here
        then the implementation of that function must be changed too """
        alt_correction = (density(h, dT)/1.225)**0.7
        return alt_correction * self.eff_powertrain * self.eff_prop * self.max_cont_power_sealevel
    
    def P_to(self, h, dT):
        """ Available takeoff power for given conditions
        TODO: temporary placeholder implementation
        NOTE: RC_max() assumes constant P_a with velocity, if this changes here
        then the implementation of that function must be changed too """
        alt_correction = (density(h, dT)/1.225)**0.7
        return alt_correction * self.eff_powertrain * self.eff_prop * self.takeoff_power_sealevel

        
    def P_r_min(self, W, h, dT):
        """ Minimum possible value for  power required for given conditions in clean configuration,
            note that this assumes a small flight path angle. """
        return W * np.sqrt(W/self.S * 2/density(h, dT) * 1/self.L3D2max) # ruijgrok p221
        
    def RC_max(self, W, h, dT):
        """ Maximum possible climb rate for given conditions in clean configuration"""
        return (self.P_a(h, dT) - self.P_r_min(W, h, dT)) / W
    
    def climb_angle_max(self, W, h, dT, gear="up", flaps="up"):
        """ Maximum possible climb slope (RC/V) for given conditions in given config (default clean) """
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
        """ Calculates the maximum true airspeed for given conditions in clean config """
        P_a = self.P_a(h, dT)
        rho = density(h, dT)
        def diff(V):
            CL = W / (0.5 * rho * self.S * V**2)
            CD = self.CD(CL)
            D = 0.5 * rho * self.S * CD * V**2
            return P_a - D*V
        return root_scalar(diff, x0=75).root
    
    def power_setting(self, W, V, h, dT):
        """ Gives power setting required to maintain cruise (0-1, but also goes beyond 1) """
        CL = W / (0.5 * density(h, dT) * self.S * V**2)
        D = 0.5 * density(h, dT) * self.S * self.CD(CL) * V**2
        if CL > 1: return 0
        return V * D / self.P_a(h, dT)

    def takeoff_ground_run(self, W, h, dT, slope):

        # this uses the torenbeek method, torenbeek p167-168
        # however it really is not accurate enough
        
        # also ground effect should be taken into account
        
        rho = density(h, dT)
        sigma = rho/1.225
        g = 9.80665
        
        CLmax = self.CLmax_TO
        V_S_TO = np.sqrt(W / (0.5 * rho * self.S * CLmax)) # TODO: ground effect should be taken into account here
        V_LOF = 1.2 * V_S_TO # lift off speed approximation according to torenbeek
        #CL_LOF = W / (0.5 * rho * V_LOF**2 * self.S)
        CD0 = self.CD(0, gear="down", flaps="TO")
        
        P_to = self.P_to(h, dT) / g # in kgm/s
        
        mu_quote = 0.05 + .72 * CD0/CLmax # 0.04-0.05 for short grass and 0.02 for concrete, if a fixed pitch prop is used use 0.576
        
        # avg thrust in kgm/s
        Tbar = .321 * P_to * (sigma / self.number_of_engines * self.propeller_diameter**2 / P_to)**(1/3)
        Tbar *= g # convert to W
        
        s_run = V_LOF**2 / (2*g) / (Tbar/W - mu_quote)
        
        return s_run
    
        #CL_G = self.CL_ground
        #CD_G = self.CD(CL_G)
        
        #def s_g(V):
           # return V / (g * (P_a(h, dT)))
        
        
        #mu_r = 0.05 # coefficient of rolling friction for short-cut grass, ruijgrok p372 and roskam pt7 p40
        
        
        
        