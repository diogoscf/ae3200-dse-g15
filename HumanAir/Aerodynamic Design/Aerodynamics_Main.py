from ISA import ISA
import matplotlib.pyplot as plt
from PlanformDesign import Planform
from FlowParameters import Flow
from LongitudinalStability import LongitudinalStability
import numpy as np
import json
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data

MTOW=aircraft_data["Weights"]['MTOW_N']
WS=aircraft_data["Performance"]["W/S_N/m2"]
AR_Wing = aircraft_data["Aero"]["AR_Wing"]
AR_HS = aircraft_data["Aero"]["AR_HS"]
Taper_Wing = aircraft_data["Aero"]["Taper_Wing"]
Taper_HS = aircraft_data["Aero"]["Taper_HS"]
QuarterChordSweep_Wing = aircraft_data["Aero"]["QuarterChordSweep_Wing_deg"]
QuarterChordSweep_HS = aircraft_data["Aero"]["QuarterChordSweep_HS_deg"]
CruiseHeight=aircraft_data["Performance"]["Altitude_Cruise_m"]
TemperatureGradient=-0.0065
CruiseVelocity=aircraft_data["Performance"]["Vc_m/s"]
CLh=aircraft_data["Stability"]["C_L_h"]
CLah=aircraft_data["Stability"]["C_L_AH"]
Xcgh=aircraft_data["Stability"]["X_cg_HS"]
XLEMAC=aircraft_data["Stability"]["XLEMAC_m"]
CgAft=aircraft_data["Stability"]["Cg_Aft"]
CgFwd=aircraft_data["Stability"]["Cg_Front"]
SM=aircraft_data["Stability"]["Stability_Margin"]
deda=aircraft_data["Aero"]["deda"]
VhV=aircraft_data["Aero"]["VhV"]
FuselageLength=aircraft_data["Geometry"]["fus_length_m"]

checkwingplanform=True
checkflowparameters=True
checkstability=True
checkhsplanform=True

WingPlanform = Planform(AR_Wing, Taper_Wing, QuarterChordSweep_Wing,MTOW=MTOW, WS=WS)
ISACruise = ISA(CruiseHeight, TemperatureGradient)

if checkwingplanform:
    WingPlanform.PlotWingPlanform()

Flowvariables=Flow(CruiseVelocity, ISACruise, WingPlanform)

if checkflowparameters:
    print("MACH = ", Flowvariables.Mach())
    print("Reynolds = ", Flowvariables.Reynolds())
    print("Beta =", Flowvariables.Beta())


"""========== Stability Analysis =========="""
Stab=LongitudinalStability(CLh, CLah, Xcgh, XLEMAC, CgAft, CgFwd, SM, deda, VhV, FuselageLength, WingPlanform, 'Airfoils/FX_63-137.json', 'Airfoils/NACA0012.json')

if checkstability:
    Stab.Plotting()

HSPlanform = Planform(AR_HS, Taper_HS, QuarterChordSweep_HS, S=WingPlanform.WingSurfaceArea()*Stab.ShS())
if checkhsplanform:
    HSPlanform.PlotWingPlanform()

