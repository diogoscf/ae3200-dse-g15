"""
Input parameters:
---
W_TO [lbs]
S_Wing [ft^2]
n_ult [-]
AR_Wing [-]
Sweep_c/4_Wing [deg]
Sweep_c/2_Wing [deg]
Taper_Wing [-]
(t/c)_max_Wing [-]
V_H (Maximum level speed at sea level) [kts]
b_wing [ft]
t_r_wing (Maximum thickness of wing root chord) [ft]

"""
import sys
import os

import numpy as np


sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data
from unit_conversions import m_to_ft, N_to_lbs, m_squared_to_ft_squared, m_s_to_kt, W_to_hp
class Class_II_Weight:
    def __init__(self, aircraft_data):
        self.W_TO=N_to_lbs(aircraft_data["Weights"]["MTOW_N"])
        self.W_L=N_to_lbs(aircraft_data["Weights"]["W_L"]) #Not in dict yet
        self.W_F=N_to_lbs(aircraft_data["Weights"]["WF_N"])

        self.S_Wing = m_squared_to_ft_squared(aircraft_data["Aero"]["S_Wing"]) #Not in dict yet
        self.S_h = m_squared_to_ft_squared(aircraft_data["Aero"]["S_h"]) #Not in dict yet
        self.S_v = m_squared_to_ft_squared(aircraft_data["Aero"]["S_v"]) #Not in dict yet

        self.n_ult = aircraft_data["Performance"]["n_ult"]  #Not in dict yet
        self.n_ult_l = 5.7

        self.AR_Wing = aircraft_data["Aero"]["AR"]
        self.AR_h = aircraft_data["Aero"]["AR_HS"]
        self.AR_v = aircraft_data["Aero"]["AR_v"] #Not in dict yet

        self.QuarterChordSweep_Wing = np.deg2rad(aircraft_data["Aero"]["QuarterChordSweep_Wing_deg"])
        self.HalfChordSweep_Wing = np.deg2rad(aircraft_data["Aero"]["HalfChordSweep_Wing_deg"]) #Not in dict yet
        self.QuarterChordSweep_v = np.deg2rad(aircraft_data["Aero"]["QuarterChordSweep_v_deg"]) #Not in dict yet

        self.Taper_Wing = aircraft_data["Aero"]["Taper_Wing"]
        self.tc_m_Wing = m_to_ft(aircraft_data["Aero"]["tc_m_Wing"]) #Not in dict yet

        self.b_Wing = m_to_ft(aircraft_data["Aero"]["b_Wing"]) #Not in dict yet
        self.b_h = m_to_ft(aircraft_data["Aero"]["b_h"]) #Not in dict yet
        self.b_v = m_to_ft(aircraft_data["Aero"]["b_v"]) #Not in dict yet

        self.t_root_max_Wing = m_to_ft(aircraft_data["Aero"]["t_root_max_Wing"]) #Not in dict yet
        self.t_root_max_h = m_to_ft(aircraft_data["Aero"]["t_root_max_h"]) #Not in dict yet
        self.t_root_max_v = m_to_ft(aircraft_data["Aero"]["t_root_max_v"]) #Not in dict yet

        self.V_H = m_s_to_kt(aircraft_data["Performance"]["V_H"]) #Not in dict yet
        self.V_c = m_s_to_kt(aircraft_data["Performance"]["Vc_m/s"])
        self.QCW_to_QCh = m_to_ft(aircraft_data["Stability"]["QCW_to_QCh"]) #Not in dict yet
        self.l_f_nonosecone = m_to_ft(aircraft_data["Geometry"]["l_f_nonosecone"]) #Not sure if this is the same as fuselage length
        self.p_max = m_to_ft(aircraft_data["Geometry"]["fuselage_max_perimeter"])
        self.N_pax = aircraft_data["Performance"]["N_pax"] #Includes pilots and is not in dict yet
        self.l_f = m_to_ft(aircraft_data["Geometry"]["fus_length_m"])
        self.w_f = m_to_ft(aircraft_data["Geometry"]["fus_width_m"]) #Not in dict yet
        self.h_f = m_to_ft(aircraft_data["Geometry"]["fus_height_m"]) #Not in dict yet

        self.P_TO = W_to_hp(aircraft_data["Power_prop"]["P_TO"]) #Not in dict yet
        self.K_n = aircraft_data["Power_prop"]["K_n"] #Not in dict yet, 0.37 for radial, 0.24 for horizontally opposed
        self.K_p = 1.1
        self.K_pg = 1.16
        self.K_fsp = 6.55

        self.int = aircraft_data["Power_prop"]["int"] #Fraction of integral fuel tanks, Not in dict yet

        self.l_s_m = aircraft_data["Landing_gear"]["l_s_m"] #Not sure if in dict
        self.l_s_n = aircraft_data["Landing_gear"]["l_s_n"] #Not sure if in dict

        self.N_e = aircraft_data["Power_prop"]["N_e"] #Not in dict yet
        self.N_t = aircraft_data["Power_prop"]["N_t"] #Not in dict yet





        """========== Structure Weight =========="""
    def WingWeight(self):
        results = {}

        results["Cessna"] = 0.002933*self.S_Wing**1.018*self.AR_Wing**2.473*self.n_ult**0.611
        results["USAF"] = 96.948*((self.W_TO*self.n_ult*10**(-5))**0.65*(self.AR_Wing/np.cos(self.QuarterChordSweep_Wing))**0.57*(self.S_Wing/100)**0.61*((1+self.Taper_Wing)/2*self.tc_m_Wing)**0.36*(1+self.V_H/500)**0.5)**0.993
        results["Torenbeek"] = 0.00125*self.W_TO*(self.b_Wing/np.cos(self.HalfChordSweep_Wing))**0.75*(1+(6.3*np.cos(self.HalfChordSweep_Wing)/self.b_Wing)**0.5)*self.n_ult**0.55*(self.b_Wing*self.S_Wing/self.t_root_max_Wing*self.W_TO*np.cos(self.HalfChordSweep_Wing))**0.30

        results["Average"] = np.average([results["Cessna"], results["USAF"], results["Torenbeek"]])
        return results

    def EmpennageWeight(self):
        results = {}
        """ Cessna """
        results["Cessna"]["W_h"] = (3.184*self.W_TO**0.887*self.S_h**0.101*self.AR_h**0.138)/(174.04*self.t_root_max_h**0.223)
        results["Cessna"]["W_v"] = (1.68*self.W_TO**0.567*self.S_v**0.1249*self.AR_v**0.482)/(639.95*self.t_root_max_v**0.747*(np.cos(self.QuarterChordSweep_v))**0.882)
        results["Cessna"]["W_c"] = 0
        results["Cessna"]["Total"] = results["Cessna"]["W_h"]+results["Cessna"]["W_v"]+results["Cessna"]["W_c"]

        """ USAF """
        results["USAF"]["W_h"] = 127*((self.W_TO*self.n_ult*10**(-5))**0.87*(self.S_h/100)**1.2*0.289*(self.QCW_to_QCh/10)**0.483*(self.b_h/self.t_root_max_h)**0.5)**0.458
        results["USAF"]["W_v"] = 98.5*((self.W_TO*self.n_ult*10**(-5))**0.87*(self.S_v/100)**1.2*0.289*(self.b_v/self.t_root_max_v)**0.5)**0.458
        results["USAF"]["W_c"] = 0
        results["USAF"]["Total"] = results["USAF"]["W_h"]+results["USAF"]["W_v"]+results["USAF"]["W_c"]

        """ Torenbeek """
        results["Torenbeek"] = 0.04*(self.n_ult*(self.S_v+self.S_h)**2)**0.75

        results["Average"]=np.average([results["Cessna"]["Total"], results["USAF"]["Total"], results["Torenbeek"]])

        return results

    def FuselageWeight(self):
        results = {}

        results["Cessna"] = 14.86*self.W_TO**0.144*(self.l_f_nonosecone/self.p_max)**0.778*self.l_f_nonosecone**0.383*self.N_pax**0.455
        results["USAF"] = 200*((self.W_TO*self.n_ult*10**(-5))**0.286*(self.l_f/10)**0.857*((self.w_f+self.h_f)/10)*(self.V_c/100)**0.338)**1.1

        results["Average"] = np.average([results["Cessna"], results["USAF"]])
        return results

    def NacelleWeight(self):
        results = {}

        results["Cessna"]=self.P_TO*self.K_n
        results["Torenbeek"] = 2.5*self.P_TO**0.5

        results["Average"] = np.average([results["Cessna"], results["Torenbeek"]])
        return results

    def LandingGearWeight(self, retractable):
        results = {}

        results["Cessna"] = 0.013*self.W_TO+0.362*self.W_L**0.417*self.n_ult_l**0.950*self.l_s_m**0.183*6.2+0.0013*self.W_TO+0.007157*self.W_L**0.749*self.n_ult_l*self.l_s_n**0.788
        if retractable:
            results["Cessna"]+=0.014*self.W_TO

        results["USAF"] = 0.054*self.l_s_m**0.501*self.W_L*self.n_ult_l**0.684

        results["Average"] = np.average([results["Cessna"], results["USAF"]])
        return results


    def StructureWeight_Total(self):
        return self.WingWeight()["Average"]+self.EmpennageWeight()["Average"]+self.FuselageWeight()["Average"]+self.NacelleWeight()["Average"]+self.LandingGearWeight()["Average"]


    """========== Powerplant Weight =========="""
    def PowerplantWeight_Total(self):
        results = {}
        results["USAF"]["WeWaiWpropWp"] = 2.575*(self.K_p*self.P_TO)**0.922*self.N_e
        results["USAF"]["Wfs"] = 2.49*((self.W_F/self.K_fsp)**0.6*(1/(1+self.int))**0.3*self.N_t**0.20*self.N_e**0.13)**1.21
        results["USAF"]["Total"]=results["USAF"]["WeWaiWpropWp"]+results["USAF"]["Wfs"]-self.NacelleWeight()["Average"]

        results["Torenbeek"]["Total"] = self.K_pg*(self.K_p*self.P_TO+0.24*self.P_TO)
        return np.average([results["USAF"]["Total"],results["Torenbeek"]["Total"]])

    """========== Fixed Equipment Weight =========="""


    def FixedEquipmentWeight_Total(self):

        return 1