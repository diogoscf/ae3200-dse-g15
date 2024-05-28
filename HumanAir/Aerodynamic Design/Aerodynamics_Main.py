from ISA import ISA
import matplotlib.pyplot as plt
from PlanformDesign import Planform
from FlowParameters import Flow
from LongitudinalStability import LongitudinalStability
import numpy as np
import json

MTOW=1946.34
WS=618
AR=7.44
Taper = 0.4
QuarterChordSweep = 0
CruiseHeight=3000
TemperatureGradient=-0.0065
CruiseVelocity=60

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
Stab=LongitudinalStability(CLh, CLah, Xcgh, XLEMAC, VhV, WingPlanform, )

if checkstability:
    Stab.Plotting()

HSPlanform = Planform(5, 0.4, 30, S=WingPlanform.WingSurfaceArea()*Stab.ShS())

# with open('FX_63-137_AirfoilData.json') as airfoildata:
#     p = json.load(airfoildata)
#     print(p)


