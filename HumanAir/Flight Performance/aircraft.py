from helper import *
import json
import numpy as np

class Aircraft:
    def __init__(self, FILE="design.json"):
        with open("../Configurations/" + FILE,'r') as f:
            dat = json.load(f)
            
        # aerodynamic data
        self.CD0            = dat["Aero"]["CD0"]
        self.AR             = dat["Aero"]["AR"]
        self.e              = dat["Aero"]["e"]
        self.S              = 19.32 # GA8
        
        self.CLmax_clean    = dat["Aero"]["CLmax_clean"]
        self.CLmax_TO       = dat["Aero"]["CLmax_TO"]
        self.CLmax_land     = dat["Aero"]["CLmax_land"]

        # weights
        self.W_OE           = dat["Weights"]["OEW_N"]
        self.W_F            = dat["Weights"]["WF_N"]
        self.W_bat          = dat["Weights"]["Wbat_N"]
        self.W_pl_des       = dat["Weights"]["Wpl_des_kg"] * 9.80655
        self.W_pl_max       = dat["Weights"]["Wpl_max_kg"] * 9.80655
        self.W_MTO          = dat["Weights"]["MTOW_N"]

        # propulsion
        self.sealevelpower  = 220000 # GA8
        self.eff_powertrain = dat["Power_prop"]["eta_powertrain"] # TODO: check definition
        self.eff_prop       = dat["Power_prop"]["eta_p"]
        
        # calculated data
        self.LDmax = 0.5 * np.sqrt(np.pi*self.AR*self.e/self.CD0) # ruijgrok p107
        self.CD_LDmax = 2 * self.CD0 # ruijgrok p106
        self.CL_LDmax = self.CD_LDmax * self.LDmax
        
        self.L3D2max = 3 * np.sqrt(3) / 16 * np.pi * self.AR * self.e \
            * np.sqrt(np.pi * self.AR * self.e / self.CD0) # for P_r_min for max climbrate; ruijgrok p107
        
    def CD(self, CL): # this parabolic approximation is about accurate until CL=1 (ruijgrok p226)
        """ Parabolic estimation of C_D for given C_L """
        return self.CD0 + self.CL**2 / (np.pi * self.AR * self.e)
    
    def V_Dmin(self, W, h, dT):
        """ Calculates velocity corresponding with minimum drag """
        a = 1 / (self.CD0 * np.pi * self.AR * self.e)
        return np.sqrt( W/self.S * 2/density(h, dT) * a )
    
    def D_min(self, W):
        """ Minimum drag (in Newtons) for given weight. Corresponds to V_Dmin """
        return 2 * W * np.sqrt(self.CD0 / (np.pi * self.AR * self.e))
    
    def P_a(self, h, dT):
        """ Available power for given conditions """
        # TODO: temporary placeholder implementation
        # NOTE: RC_max() assumes constant P_a with velocity, if this changes here
        # then the implementation of that function must be changed too
        alt_correction = (density(h, dT)/1.225)**0.7
        return alt_correction * self.eff_powertrain * self.eff_prop * self.sealevelpower
        
    def P_r_min(self, W, h, dT):
        """ Minimum possible power required for given conditions,
            assumes small flight path angle gamme. """
        return W * np.sqrt(W/self.S * 2/density(h, dT) * 1/self.L3D2max) # ruijgrok p221
        
    def RC_max(self, W, h, dT):
        """ Maximum possible climb rate for given conditions """
        return (self.P_a(h, dT) - self.P_r_min(W, h, dT)) / W