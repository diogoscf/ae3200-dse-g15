import numpy as np

"""========== Aircraft Design Parameters =========="""
name="Conventional Aircraft, No Canard"
A=7.44 #Cessna 206: 9.38, B2: 5.75
e=0.82
TOP=60.44
eta_p=0.85

Clmax_clean=2
Clmax_TO=2.3
Clmax_Land=2.5


Cdo=0.025 #Cessna :0.028 B2: 0.0065, 0.0165

"""========== Mission Parameters =========="""
h_TO=1800
h_Cruise=3000
h_Land=1800

V_stall=31.38 #Cessna 206
V_climb=1.2*V_stall
V_cruise=60 #Cessna 206

s_land=500
f= 0.95#W_to/W_land
climbrate=4.5
Cl_SafetyFactor=1.2
nmax=4.5

CruisePower=0.8
CruiseWeight=0.97

