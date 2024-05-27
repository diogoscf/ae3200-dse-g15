from ISA import ISA
from PlanformDesign import Planform
from FlowParameters import Mach, Reynolds

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









