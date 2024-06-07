import pandas as pd
import numpy as np
import sys
from math import tan, sqrt, pi
import time
import os
import json
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data
from HumanAir.Weights_and_CG.weight_fractions_C2 import calculate_lh, optimised_xlemac_landing_gears
from isa import isa

def Geometry(acd = aircraft_data):
    taper = acd["Aero"]["Taper_Wing"]
    # Root chord
    Cr = 3/2*acd["Aero"]["MAC_wing"]*((1+taper)/(1+taper+taper**2))

    # Net wing area, outside of fuselage
    Cfus = Cr*(1 - (1 - taper)*acd["Geometry"]["fus_width_m"]/acd["Aero"]["b_Wing"])
    Snet = acd["Aero"]["S_Wing"] - (Cr + Cfus)*acd["Geometry"]["fus_width_m"]/2

    # Length to net root chord, outside of fuselage
    Ymgc = acd["Aero"]["b_Wing"]/(2*(1-taper))*(1 - 2/3*((1+taper+taper**2)/(1+taper)))
    tanLE = -4/acd["Aero"]["AR"]*(-0.25*(1 - taper)/(1 + taper))
    l_fn = acd["Geometry"]["XLEMAC_m"] - tanLE*(Ymgc-acd["Geometry"]["fus_width_m"]/2)

    # Stabiliser affected by propwash
    taper_h = acd["Aero"]["Taper_HS"]
    Cr_h = 3/2*acd["Aero"]["MAC_HS"]*((1+taper_h)/(1+taper_h*taper_h**2))
    Cpw = Cr_h * (1 - (1 - taper_h)*acd["Power_prop"]["Dp_m"]/acd["Aero"]["b_h"])
    Shslip = (Cr_h + Cpw)*acd["Power_prop"]["Dp_m"]/2

    # Distance between quarter chord MAC of wing and stabiliser

    return Cr, Snet, l_fn, Shslip

def TailAero(l_H, acd = aircraft_data):
    # Get and convert values to imp*rial
    Shslip = Geometry(acd)[3] * 10.7639104
    Sh = acd["Aero"]["S_h"] * 10.7639104
    U1 = acd["Performance"]["Vc_m/s"] * 3.2808399
    Pav = 0.98 * 400 * acd["Power_prop"]["eta_p"]
    Dp = acd["Power_prop"]["Dp_m"] / 0.3048
    qbar = (0.5 * isa(acd["Performance"]["Altitude_Cruise_m"], acd["Performance"]["Temp_offset_TO_Land_cruise"])[2]*acd["Performance"]["Vc_m/s"]**2) * 0.020885

    # Magic Roskam equation, only god knows how this works
    eta_H = 1 + (Shslip/Sh) * (2200*Pav / (qbar*U1*pi*Dp**2))
    VhVcorr =acd["Aero"]["VhV"] * eta_H

    # Downwash Gradient supposedly from Slingerland whoever that may be
    r = l_H*2/acd["Aero"]["b_Wing"]
    mtv = 0*2/acd["Aero"]["b_Wing"]
    deda = ((r/(r**2 + mtv**2)) * 0.4876/np.sqrt(r**2 + 0.6319 + mtv**2) + (1 + (r**2/(r**2 + 0.7915 + 5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*acd["Aero"]["CLalpha"]/(pi*acd["Aero"]["AR"])

    return VhVcorr, deda

def Liftrate(l_H, acd = aircraft_data):
    # Tail lift rate
    SweepHS_05 = -4/acd["Aero"]["AR_HS"]*(0.25*(1-acd["Aero"]["Taper_HS"])/(1+acd["Aero"]["Taper_HS"]))
    ClaH = 2*pi*acd["Aero"]["AR_HS"] / (2 + sqrt(4 + (acd["Aero"]["AR_HS"]/0.95)**2 * (1 + tan(SweepHS_05)**2)))

    # Aircraft less tail lift rate
    Snet = Geometry(acd)[1]
    CLaAH = acd["Aero"]["CLalpha"]*(1+2.15*acd["Geometry"]["fus_width_m"]/acd["Aero"]["b_Wing"])*Snet/acd["Aero"]["S_Wing"] + pi/2*acd["Geometry"]["fus_width_m"]**2/acd["Aero"]["S_Wing"]

    # Controllability tail lift
    CLH = -0.35*acd["Aero"]["AR_HS"]**(1/3)

    # Controllability aircraft-less-tail lift
    CLAH = acd["Aero"]["CLmax_Land"] - CLH*TailAero(l_H, acd)[0]**2 *acd["Aero"]["S_h"]/acd["Aero"]["S_Wing"]

    return ClaH, CLaAH, CLH, CLAH

def Xacplane(l_H, acd = aircraft_data):
    # From ADSEE data for no wing sweep
    Xac_wing = 0.205 + 0.005*acd["Aero"]["AR"]

    # Fuselage contribution
    CLaAH = Liftrate(l_H, acd)[1]
    l_fn = Geometry(acd)[2]
    Xac_fus = - (1.8/CLaAH)*(acd["Geometry"]["fus_width_m"]*acd["Geometry"]["fus_height_m"]*l_fn)/(acd["Aero"]["S_Wing"]*acd["Aero"]["MAC_wing"])

    # Total Xac
    Xac = Xac_fus + Xac_wing

    return Xac

def MomentCoefficient(l_H, acd = aircraft_data):
    # Wing contribution
    with open("Airfoils/FX_63-137.json") as WingAirfoil:
        data = json.load(WingAirfoil)
        Cm_0_Wing = data['Cm_0']

    Cmacw = Cm_0_Wing*(acd["Aero"]["AR"]*np.cos(acd["Aero"]["QuarterChordSweep_Wing_deg"])**2/(acd["Aero"]["AR"]+2*np.cos(acd["Aero"]["QuarterChordSweep_Wing_deg"])))

    # Flap Contribution
    ClaH, CLaAH, CLH, CLAH = Liftrate(l_H, acd)
    Xac = Xacplane(l_H, acd)
    Cmacflap = acd["Flaps"]["mu2"]*(-acd["Flaps"]["mu1"]*acd["Flaps"]["DeltaClmax"]*acd["Flaps"]["cprime_c_landing"]-(CLAH+acd["Flaps"]["DeltaClmax"]*(1-acd["Flaps"]["Swf"]/acd["Aero"]["S_Wing"]))*1/8*acd["Flaps"]["cprime_c_landing"]*(acd["Flaps"]["cprime_c_landing"]-1))-acd["Aero"]["CLmax_Land"]*(0.25-Xac/acd["Aero"]["MAC_wing"])

    # Fuselage contribution
    Cmacfus = -1.8*(1-2.5*acd["Geometry"]["fus_width_m"]/acd["Geometry"]["fus_length_m"])*(pi*acd["Geometry"]["fus_width_m"]*acd["Geometry"]["fus_height_m"]*acd["Geometry"]["fus_length_m"])/(4*acd["Aero"]["S_Wing"]*acd["Aero"]["MAC_wing"])*acd["Flaps"]["CL0"]/CLaAH

    return Cmacw + Cmacflap + Cmacfus

def StabControl(acd = aircraft_data):
    # Initiate xcg/MAC
    Xcg = np.arange(-0.20, 1.20, 0.01)

    # l_H iteration/calculation, XLEMAC placement
    # TODO: add this
    # TODO: done
    calculate_lh(ac_data = aircraft_data, hinge_chord_percentage = 3/4)
    l_H = acd["Stability"]["QCW_to_QCh"]

    # Get necessary values from functions
    Xac = Xacplane(l_H, acd)
    ClaH, CLaAH, CLH, CLAH = Liftrate(l_H, acd)
    VhVcorr, deda = TailAero(l_H, acd)
    Cmac = MomentCoefficient(l_H, acd)
    # TODO: Fix Cmac, now positive

    # Stability Sh/S
    StabSM = (Xcg - Xac + 0.05)/((ClaH/CLaAH)*(1-deda)*l_H/acd["Aero"]["MAC_wing"]*VhVcorr**2)
    StabNeutral = (Xcg - Xac)/((ClaH/CLaAH)*(1-deda)*l_H/acd["Aero"]["MAC_wing"]*VhVcorr**2)

    # Controllability Sh/S
    Control = (Xcg - Xac + Cmac/CLAH)/(CLH/CLAH*l_H/acd["Aero"]["MAC_wing"]*VhVcorr**2)

    return StabSM, StabNeutral, Control, Xcg

def TailSizing(acd = aircraft_data, begin_value = -1, end_value = 8, step = 1):

    for x_percentage in range(begin_value, end_value, step):

        optimised_xlemac_landing_gears(ac_data = acd, percentage = x_percentage / 10, bat_xcg_init = 0.1)
        StabSM, _, Control, Xcg = StabControl(acd = acd)
        Xcg_excursion = np.array([0,0,0]) # TODO:potatoplot()

        optimisation = False
        S_h = acd["Aero"]["S_h"]
        S = acd["Aero"]["S_Wing"]

        while not optimisation:

            Sh_S = S_h / S

            if round(Sh_S, 4) in [round(i, 4) for i in StabSM] and round(Sh_S, 4) in [round(i, 4) for i in Control]:

                end_idx = np.where(round(Sh_S,4) == np.round(Control, 4))[0][0]
                begin_idx = np.where(round(Sh_S,4) == np.round(StabSM, 4))[0][0]

                if end_idx - begin_idx == np.shape(Xcg_excursion):
                    optimisation = True

            else:
                S_h -= 0.01

        Plotting(acd = acd, show = False)
        Sh_S_list = np.ones(np.shape(Xcg_excursion)[0]) * Sh_S

        plt.plot(Xcg_excursion, Sh_S_list, label = "Optimised Xcg for Landing Gear")

        answer = input("Do you want to continue? [Y/N]")

        if answer.lower() == "n":

            # saving all of the updated horizontal tail values
            acd["Aero"]["S_h"] = S_h
            break

    # TODO: Make this

def Plotting(acd = aircraft_data, show = True):
    # Get data to plot from previous functions
    StabSM, StabNeutral, Control, Xcg = StabControl(acd)
    TailSizing(acd = aircraft_data)
    # Actual plotting
    plt.plot(Xcg, StabSM, label="Stability with Safety Margin", color="limegreen", linewidth=2.2)
    plt.plot(Xcg, StabNeutral, label="Neutral Stability", linestyle="--", color="red")
    plt.plot(Xcg, Control, label="Controllability", color="dodgerblue", linewidth=2.2)
    plt.fill_between(Xcg, 0, StabSM, color='crimson', alpha=.2)
    plt.fill_between(Xcg, 0, Control, color='crimson', alpha=.2)
    plt.xlim(-0.2, 1.2)
    plt.ylim(0, 0.5)
    plt.legend()

    if show:
        plt.show()

if __name__ == "__main__":
    #Plotting()
    TailSizing(acd=aircraft_data)
