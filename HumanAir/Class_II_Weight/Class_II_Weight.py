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
        self.W_E=N_to_lbs(aircraft_data["Weights"]["W_E"]) #Not sure if in dict, as only OEW is given

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
        self.M_D = m_s_to_kt(aircraft_data["Performance"]["M_D"]) #Not in dict yet
        self.QCW_to_QCh = m_to_ft(aircraft_data["Stability"]["QCW_to_QCh"]) #Not in dict yet
        self.l_f_nonosecone = m_to_ft(aircraft_data["Geometry"]["l_f_nonosecone"]) #Not sure if this is the same as fuselage length
        self.paint = aircraft_data["General"]["Paint"]
        self.p_max = m_to_ft(aircraft_data["Geometry"]["fuselage_max_perimeter"])
        self.N_pax = aircraft_data["Performance"]["N_pax"] #Includes pilots and is not in dict yet
        self.N_row = aircraft_data["Geometry"]["N_row"] #Not in dict yet
        self.l_f = m_to_ft(aircraft_data["Geometry"]["fus_length_m"])
        self.w_f = m_to_ft(aircraft_data["Geometry"]["fus_width_m"]) #Not in dict yet
        self.h_f = m_to_ft(aircraft_data["Geometry"]["fus_height_m"]) #Not in dict yet

        self.P_TO = W_to_hp(aircraft_data["Power_prop"]["P_TO"]) #Not in dict yet
        self.K_n = aircraft_data["Power_prop"]["K_n"] #Not in dict yet, 0.37 for radial, 0.24 for horizontally opposed
        self.K_p = 1.1
        self.K_pg = 1.16
        self.K_fsp = 6.55

        self.int = aircraft_data["Power_prop"]["int"] #Fraction of integral fuel tanks, Not in dict yet

        self.l_s_m = m_to_ft(aircraft_data["Landing_gear"]["l_s_m"]) #Not sure if in dict
        self.l_s_n = m_to_ft(aircraft_data["Landing_gear"]["l_s_n"]) #Not sure if in dict
        self.retractable = aircraft_data["Landing_gear"]["Retractable"]

        self.N_e = aircraft_data["Power_prop"]["N_e"] #Not in dict yet
        self.N_t = aircraft_data["Power_prop"]["N_t"] #Not in dict yet

        self.PoweredFlightControls = aircraft_data["General"]["PoweredFlightControls"] #Not in dict yet
        self.DuplicatedFlightControls = aircraft_data["General"]["DuplicatedFlightControls"] #Not in dict yet





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

    def LandingGearWeight(self):
        results = {}

        results["Cessna"] = 0.013*self.W_TO+0.362*self.W_L**0.417*self.n_ult_l**0.950*self.l_s_m**0.183*6.2+0.0013*self.W_TO+0.007157*self.W_L**0.749*self.n_ult_l*self.l_s_n**0.788
        if self.retractable:
            results["Cessna"]+=0.014*self.W_TO

        results["USAF"] = 0.054*self.l_s_m**0.501*self.W_L*self.n_ult_l**0.684

        results["Average"] = np.average([results["Cessna"], results["USAF"]])
        return results


    def StructureWeight_Total(self):
        return self.WingWeight()["Average"]+self.EmpennageWeight()["Average"]+self.FuselageWeight()["Average"]+self.NacelleWeight()["Average"]+self.LandingGearWeight()["Average"]


    """========== Powerplant Weight =========="""
    def FuelSystemWeight(self):
        results={}

        results["Cessna"] = 0.40*self.W_F/self.K_fsp
        results["USAF"] = 2.49*((self.W_F/self.K_fsp)**0.6*(1/(1+self.int))**0.3*self.N_t**0.20*self.N_e**0.13)**1.21
        results["Torenbeek"] = 2*(self.W_F/5.87)**0.667

        results["Average"] = np.average([results["Cessna"], results["USAF"], results["Torenbeek"]])

        return results

    def PowerplantWeight_Total(self):
        results = {}
        results["USAF"]["WeWaiWpropWp"] = 2.575*(self.K_p*self.P_TO)**0.922*self.N_e
        results["USAF"]["Total"]=results["USAF"]["WeWaiWpropWp"]+self.FuelSystemWeight()["USAF"]-self.NacelleWeight()["Average"]

        results["Torenbeek"]["Total"] = self.K_pg*(self.K_p*self.P_TO+0.24*self.P_TO)
        return np.average([results["USAF"]["Total"],results["Torenbeek"]["Total"]])

    """========== Fixed Equipment Weight =========="""
    def FlightControlSystem(self):
        results={}

        results["Cessna"]=0.016*self.W_TO
        if not self.PoweredFlightControls:
            results["USAF"] = 1.066*self.W_TO**0.626
            if not self.DuplicatedFlightControls:
                results["Torenbeek"]=0.33*self.W_TO**(2/3)
        else:
            results["USAF"] = 1.08*self.W_TO**0.7
            results["Torenbeek"]=False

        finalresults=[]
        for i in results.values():
            if i!=False:
                finalresults.append(i)
        results["Average"]=np.average(finalresults)
        return results

    def HydraulicsPneumatics(self):
        return {"Average": 0.008*self.W_TO}
    def InstrumentsAvionicsElectronics(self):
        return {"Average": 33*self.N_pax}


    def ElectricalSystemWeight(self):
        results={}
        results["Cessna"] = 0.0268*self.W_TO
        results["USAF"] = 426*((self.FuelSystemWeight()["Average"]+self.InstrumentsAvionicsElectronics()["Average"])/1000)**0.51

        results["Average"] = np.average([results["Cessna"],results["USAF"]])
        return results

    def AirconPressurizationAntiDeicingWeight(self):
        results={}

        results["USAF"] = 0.265*self.W_TO**0.52*self.N_pax**0.68*self.InstrumentsAvionicsElectronics()["Average"]**0.17*self.M_D**0.08
        if self.N_e>1:
            results["Torenbeek"]=0.018*self.W_E
        else:
            results["Torenbeek"] = 2.5*self.N_pax

        results["Average"]= np.average([results["USAF"], results["Torenbeek"]])
        return results

    def OxygenSystem(self):
        return {"Average": 20+0.5*self.N_pax}

    def APU(self):
        results={}
        if self.APU:
            results["Average"] = 0.0085*self.W_TO
        else:
            results["Average"] = 0
        return results

    def Furnishings(self):
        results = {}
        results["Cessna"] = 0.412*self.N_pax**1.145*self.W_TO**0.489
        results["Torenbeek"] = 5+13*self.N_pax+25*self.N_row

        results["Average"] = np.average([results["Cessna"], results["Torenbeek"]])
        return results

    def AuxiliaryGear(self):
        return {"Average": 0.01*self.W_E}

    def Paint(self):
        if self.paint:
            results={"Average": 0.0045*self.W_TO}
        else:
            results={"Average": 0}
        return results



    def FixedEquipmentWeight_Total(self):
        return self.FlightControlSystem()["Average"]+self.HydraulicsPneumatics()["Average"]+self.InstrumentsAvionicsElectronics()["Average"]+self.ElectricalSystemWeight()["Average"]+self.AirconPressurizationAntiDeicingWeight()["Average"]+self.OxygenSystem()["Average"]+self.APU()["Average"]+self.Furnishings()["Average"]+self.AuxiliaryGear()["Average"]+self.Paint()["Average"]

    def NewEmptyWeight(self):
        return self.PowerplantWeight_Total()+self.StructureWeight_Total()+self.FixedEquipmentWeight_Total()


def RunClassII(aircraft_data, check):
    p=Class_II_Weight(aircraft_data)

    if check:
        print("========== Structures Weight ==========")
        print('\nWing Weight = ', p.WingWeight()["Average"])
        print('Empennage Weight =', p.EmpennageWeight()['Average'])
        print('Fuselage Weight = ', p.FuselageWeight()["Average"])
        print('Nacelle Weight = ', p.NacelleWeight()["Average"])
        print('Landing Gear Weight = ', p.LandingGearWeight()["Average"])
        print('\nTotal Structures Weight = ', p.StructureWeight_Total())
        print('\n\n ========== Powerplant Weight ==========')
        print('\n Fuel System Weight = ', p.FuelSystemWeight()["Average"])
        print('\n Total Powerplant Weight = ', p.PowerplantWeight_Total())
        print('\n\n ========== Fixed Equipment Weight ==========')
        print('\n Flight Control Systems Weight = ', p.FlightControlSystem()["Average"])
        print("Hydraulics and/or Pneumatics Weight = ", p.HydraulicsPneumatics()["Average"])
        print("Instruments, Avionics and Electronics Weight = ", p.InstrumentsAvionicsElectronics()["Average"])
        print("Electrical System Weight = ", p.ElectricalSystemWeight()["Average"])
        print("Airconditioning, Pressurization and Anti or Deicing Weight = ", p.AirconPressurizationAntiDeicingWeight()["Average"])
        print("Oxygen System Weight = ", p.OxygenSystem()["Average"])
        print("APU Weight = ", p.APU()["Average"])
        print("Furnishings Weight = ", p.Furnishings()["Average"])
        print("Auxiliary Gear Weight = ", p.AuxiliaryGear()["Average"])
        print("Paint Weight = ", p.Paint()["Average"])
        print("\nTotal Fixed Equipment Weight = ", p.FixedEquipmentWeight_Total())
        print("\n\nTotal Empty Weight = ", p.NewEmptyWeight())

    return p.NewEmptyWeight()