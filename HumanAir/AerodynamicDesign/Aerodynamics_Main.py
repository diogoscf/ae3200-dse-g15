import numpy as np
import json
import os
import sys
import matplotlib.pyplot as plt

from AerodynamicDesign.ISA import ISA
from AerodynamicDesign.PlanformDesign import Planform
from AerodynamicDesign.FlowParameters import Flow
from AerodynamicDesign.LongitudinalStability import LongitudinalStability

sys.path.append(os.path.dirname(os.path.dirname(__file__)))



from aircraft_data import aircraft_data

def aerodynamic_design(ac_data, checkwingplanform=False, checkflowparameters=False, checkstability=False, checkhsplanform=False):
    MTOW=ac_data["Weights"]['MTOW_N']
    WS=ac_data["Performance"]["W/S_N/m2"]
    AR_Wing = ac_data["Aero"]["AR"]
    AR_HS = ac_data["Aero"]["AR_HS"]
    Taper_Wing = ac_data["Aero"]["Taper_Wing"]
    Taper_HS = ac_data["Aero"]["Taper_HS"]
    QuarterChordSweep_Wing = ac_data["Aero"]["QuarterChordSweep_Wing_deg"]
    QuarterChordSweep_HS = ac_data["Aero"]["QuarterChordSweep_HS_deg"]
    CruiseHeight=ac_data["Performance"]["Altitude_Cruise_m"]
    Temp_offset=ac_data["Performance"]["Temp_offset_TO_Land_cruise"]
    TemperatureGradient=-0.0065
    CruiseVelocity=ac_data["Performance"]["Vc_m/s"]
    CLh=ac_data["Stability"]["C_L_h"]
    CLah=ac_data["Stability"]["C_L_AH"]
    Xcgh=ac_data["Stability"]["X_cg_HS"]
    XLEMAC=ac_data["Stability"]["XLEMAC_m"]
    CgAft=ac_data["Stability"]["Cg_Aft"]
    CgFwd=ac_data["Stability"]["Cg_Front"]
    SM=ac_data["Stability"]["Stability_Margin"]
    deda=ac_data["Aero"]["deda"]
    VhV=ac_data["Aero"]["VhV"]
    FuselageLength=ac_data["Geometry"]["fus_length_m"]

    WingPlanform = Planform(AR_Wing, Taper_Wing, QuarterChordSweep_Wing,MTOW=MTOW, WS=WS)

    ISACruise = ISA(CruiseHeight, Temp_offset, TemperatureGradient)

    if checkwingplanform:
        WingPlanform.PlotWingPlanform()

    Flowvariables=Flow(CruiseVelocity, ISACruise, WingPlanform)

    if checkflowparameters:
        print("MACH = ", Flowvariables.Mach())
        print("Reynolds = ", Flowvariables.Reynolds())
        print("Beta =", Flowvariables.Beta())


    """========== Stability Analysis =========="""
    # replace with 'FX_63-137.json' and for pycharm 'NACA0012.json'
    # replace with 'c:\\Users\\nicho\\Documents\\GitHub\\ae3200-dse-g15\\HumanAir\\Configurations\\FX_63-137.json' and same for the other airfoil for vscode (change nicho to your username)
    Stab=LongitudinalStability(CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, WingPlanform, 
                               '../AerodynamicDesign/Airfoils/FX_63-137.json',
                               '../AerodynamicDesign/Airfoils/NACA0012.json')

    if checkstability:
        Stab.Plotting()

    HSPlanform = Planform(AR_HS, Taper_HS, QuarterChordSweep_HS, S=WingPlanform.WingSurfaceArea()*Stab.ShS())
    if checkhsplanform:
        HSPlanform.PlotWingPlanform()

    mac = WingPlanform.MAC()

    return mac

if __name__ == "__main__":
    # replace with 'design.json' for pycharm
    # replace with 'c:\\Users\\nicho\\Documents\\GitHub\\ae3200-dse-g15\\HumanAir\\Configurations\\design.json' for vscode (change nicho to your username)
    # with open("design.json",'r') as f:
    #     dict = json.load(f)
    aerodynamic_design(aircraft_data,checkwingplanform=True, checkflowparameters=True, checkstability=True, checkhsplanform=True)
