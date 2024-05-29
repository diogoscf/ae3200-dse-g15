from ISA import ISA
import matplotlib.pyplot as plt
from PlanformDesign import Planform
from FlowParameters import Flow
from LongitudinalStability import LongitudinalStability
import numpy as np
import json

MTOW=1946.34
WS=618
AR=9
Taper = 0.4
QuarterChordSweep = 0
CruiseHeight=3000
TemperatureGradient=-0.0065
CruiseVelocity=60
CLh=-0.5
CLah=1.72
Xcgh=0.95
XLEMAC=2.69
CgAft=0.75
CgFwd=0.25
SM=0.05
deda=0
VhV=1
FuselageLength=10.09

checkwingplanform=True
checkflowparameters=True
checkstability=True

WingPlanform = Planform(AR, Taper, QuarterChordSweep,MTOW, WS)
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

HSPlanform = Planform(5, 0.4, 30, S=WingPlanform.WingSurfaceArea()*Stab.ShS())


