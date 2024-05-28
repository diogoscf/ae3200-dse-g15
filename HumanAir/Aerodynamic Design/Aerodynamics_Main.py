from ISA import ISA
import matplotlib.pyplot as plt
from PlanformDesign import Planform
from FlowParameters import Mach, Reynolds
from LongitudinalStability import LongitudinalStability


MTOW=1946.34
WS=618
AR=9.35
Taper = 0.4
QuarterChordSweep = 0
CruiseHeight=3000
TemperatureGradient=-0.0065
CruiseVelocity=60

checkwingplanform=False
checkflowparameters=False

WingPlanform = Planform(MTOW, WS, AR, Taper, QuarterChordSweep)
ISACruise = ISA(CruiseHeight, TemperatureGradient)

if checkwingplanform:
    WingPlanform.PlotWingPlanform()

MACH=Mach(CruiseVelocity, ISACruise.SpeedOfSound())
Reynolds=Reynolds(CruiseVelocity,ISACruise.DynamicViscosity(), WingPlanform.MAC(), ISACruise.Density())

if checkflowparameters:
    print("MACH = ", MACH)
    print("Reynolds = ", Reynolds)

Stab=LongitudinalStability(-0.5, 1.72, 0.343, 0.95, 1, WingPlanform.MAC(), -0.199498857, 0.287701353, 0.05, 5.93792624, 6.237426906, 0, 10.09, 0.75, 0.25)
print(Stab.ShS())
Stab.Plotting()









