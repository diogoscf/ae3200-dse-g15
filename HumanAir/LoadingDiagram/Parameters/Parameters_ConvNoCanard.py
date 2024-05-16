import numpy as np

"""========== Aircraft Design Parameters =========="""
name="Conventional Aircraft, No Canard"
A=9.38 #Cessna 206: 9.38, B2: 5.75
e=0.82
TOP=66.91 # was 165
eta_p=0.85

Clmax_clean=1.6
Clmax_TO=1.6
Clmax_Land=2.5


Cdo=0.024 #Cessna :0.028 B2: 0.0065, 0.0165

"""========== Mission Parameters =========="""
h_TO=1800
h_Cruise=3000
h_Land=1800

V_stall=31.38 #CS23 - at take off weight
V_climb=35
V_cruise=60 #Cessna 206

s_land=500
f= 0.997#W_to/W_land
climbrate=5
print(climbrate)
Cl_SafetyFactor=1.2
nmax=4.5

CruisePower=0.8
CruiseWeight=0.9

