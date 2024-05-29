import pandas as pd
import numpy as np
import sys
from math import tan, sqrt
import time
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data

def find_lg(ac_datafile = aircraft_data):
    # Import tyre database
    tyres = pd.read_csv('tiredata.csv', index_col=0).to_numpy()

    # Choose smallest available tyre
    Pmw = 0.90 * ac_datafile["Weights"]["MTOW_N"] / (2 * 9.81)
    Pnw = 0.10 * ac_datafile["Weights"]["MTOW_N"] / 9.81
    for tyre in range(len(tyres[:, 0])):
        Wt_m = tyres[tyre, 0]
        Dw_m = tyres[tyre, 1]
        if tyres[tyre, 2] >= Pmw:
            break

    for tyre in range(len(tyres[:, 0])):
        Wt_n = tyres[tyre, 0]
        Dw_n = tyres[tyre, 1]
        if tyres[tyre, 2] >= Pnw:
            break

    if Pmw > tyres[-1, 2]:
        print("WARNING: NO TYRE AVAILBLE")
        con = input("Continue? (y/n): ")
        if con == "n":
            sys.exit("Too heavy for landing gear")

    # Calculate landing gear geometry
    Hcg = 0.5 * ac_datafile["Geometry"]["fus_height_m"]
    H_s = 1.5 * Dw_m
    l_m = tan(np.radians(15))*(Hcg+H_s+0.5*Dw_m)
    l_n = l_m * Pmw/Pnw
    ymin = (l_m + l_n)/(sqrt(l_n**2 * tan(np.radians(55))**2 / (Hcg+H_s+0.5*Dw_m)**2 - 1))

    # Write values to dict
    ac_datafile["Landing_gear"]["lm_m"] = l_m
    ac_datafile["Landing_gear"]["ln_m"] = l_n
    ac_datafile["Landing_gear"]["ymin_m"] = ymin
    ac_datafile["Landing_gear"]["Hs_m"] = H_s
    ac_datafile["Landing_gear"]["Dwm_m"] = Dw_m
    ac_datafile["Landing_gear"]["Dwn_m"] = Dw_n
    ac_datafile["Landing_gear"]["Wtm_m"] = Wt_m
    ac_datafile["Landing_gear"]["Wtn_m"] = Wt_n

    return l_m, l_n, ymin, H_s, Dw_m, Wt_m, Dw_n, Wt_n

def component_mass(ac_datafile = aircraft_data):
    # Import statistical weight fraction data
    fracs = pd.read_csv('fraction-database.csv', index_col=0).to_numpy()

    # Convert weights to kg and with contingency
    # TODO: Change to actual datafile
    MTOW_cont = ac_datafile["Weights"]["MTOW_N"]/9.81*ac_datafile["Contingency"]
    OEW_cont = ac_datafile["Weights"]["OEW_N"]/9.81*ac_datafile["Contingency"]
    Wbat_cont = ac_datafile["Weights"]["Wbat_N"]/9.81*ac_datafile["Contingency"]
    Ww_cont = ac_datafile["Weights"]["Ww_N"]/9.81*ac_datafile["Contingency"]

    # Set up known fractions
    OEW_frac = OEW_cont/MTOW_cont
    OEW_avg = np.average(fracs[:, -1])
    Ww_frac = Ww_cont/MTOW_cont
    Wbat_frac = Wbat_cont/MTOW_cont

    # Get adjusted component fractions: Wing, MLG, pwtr, NLG, fus, emp, FE, bat, EW
    Ww_diff = Ww_frac - (np.average(fracs[:, 0])*(OEW_frac / OEW_avg) - (Wbat_frac * np.average(fracs[:, 0]) / OEW_avg))
    wcg = np.zeros((3, 9))
    wcg[0, 0] = Ww_frac
    wcg[0, -2] = Wbat_frac
    wcg[0, -1] = OEW_frac
    for i in range(1, 7, 1):
        wcg[0, i] = np.average(fracs[:, i])*(OEW_frac / OEW_avg) - (Wbat_frac * np.average(fracs[:, i]) / OEW_avg) - (Ww_diff * np.average(fracs[:, i]) / (OEW_avg - np.average(fracs[:, 0])))
    for i in range(len(wcg[0, :])):
        wcg[1, i] = wcg[0, i] * MTOW_cont
    fracsum = np.sum(wcg[0, 0:-1])

    # Check whether fractions make sense
    if abs(fracsum/OEW_frac - 1) > 0.01:
        print("WARNING: WEIGHT FRACTIONS DIVERGE")
        con = input("Continue? (y/n): ")
        if con == "n":
            sys.exit("Weight fractions diverged")

    # Return fractions and masses of each component: Wing, MLG, pwtr, NLG, fus, emp, FE, bat, EW
    return wcg

def iterate_cg_lg(ac_datafile = aircraft_data):
    # Get fractions, weights, cg
    wcg = component_mass(ac_datafile)
    WF_cont = ac_datafile["Weights"]["WF_N"] / 9.81 * ac_datafile["Contingency"]
    WPL_cont = ac_datafile["Weights"]["Wpl_des_kg"] * ac_datafile["Contingency"]

    # Get preliminary moving CG locations from the nose
    Xcg_pld = 0.5 * ac_datafile["Geometry"]["fus_length_m"]
    Xcg_f = ac_datafile["Geometry"]["XLEMAC_m"] + 0.4 * ac_datafile["Geometry"]["MGC_m"]

    # Get preliminary component CG locations
    CGw_MAC = 0.4 * ac_datafile["Geometry"]["MGC_m"]
    wcg[2, 0] = ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC
    wcg[2, 2] = 0.05 * ac_datafile["Geometry"]["fus_length_m"]
    wcg[2, 4] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]
    wcg[2, 5] = 0.9 * ac_datafile["Geometry"]["fus_length_m"]
    wcg[2, 6] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]
    wcg[2, 7] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]

    Xcg_OEW = ((ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC)*wcg[1, 0] + np.average(wcg[2, [2, 4, 5, 6, 7]], weights=wcg[1, [2, 4, 5, 6, 7]]) * np.sum(wcg[1, [2, 4, 5, 6, 7]])) / (wcg[1, 0] + np.sum(wcg[1, [2, 4, 5, 6, 7]]))
    CGlist = [Xcg_OEW, (Xcg_OEW*wcg[1, -1] + Xcg_pld*WPL_cont)/(wcg[1, -1]+WPL_cont), (Xcg_OEW*wcg[1, -1] + Xcg_f*WF_cont)/(wcg[1, -1]+WF_cont), (Xcg_OEW*wcg[1, -1] + Xcg_pld*WPL_cont + Xcg_f*WF_cont)/(wcg[1, -1]+WPL_cont+WF_cont)]
    aftcg = np.max(CGlist)

    l_m, l_n = find_lg(ac_datafile)[0:2]
    wcg[2, 1] = aftcg + l_m
    wcg[2, 3] = aftcg - l_n
    wcg[2, -1] = Xcg_OEW

    # Iterate on CG and LEMAC positions
    iter = 1.0
    xlemac = ac_datafile["Geometry"]["XLEMAC_m"]
    wcg[2, 0] = CGw_MAC + xlemac
    while iter > 0.0001:
        xlemacold = xlemac
        Xcg_OEW = np.average(wcg[2, 0:8], weights=wcg[1, 0:8])
        Xcg_f = xlemac + 0.4 * ac_datafile["Geometry"]["MGC_m"]
        CGlist = [Xcg_OEW, (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont) / (wcg[1, -1] + WPL_cont), (Xcg_OEW * wcg[1, -1] + Xcg_f * WF_cont) / (wcg[1, -1] + WF_cont), (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont + Xcg_f * WF_cont) / (wcg[1, -1] + WPL_cont + WF_cont)]
        aftcg = np.max(CGlist)
        l_m, l_n = find_lg(ac_datafile)[0:2]
        wcg[2, 1] = aftcg + l_m
        wcg[2, 3] = aftcg - l_n
        wcg[2, -1] = Xcg_OEW
        cgwg = np.average(wcg[2, 0:2] - xlemac, weights=wcg[1, 0:2])
        xlemac = np.average(wcg[2, 2:8], weights=wcg[1, 2:8])+ac_datafile["Geometry"]["MGC_m"]*((cgwg/ac_datafile["Geometry"]["MGC_m"])*np.sum(wcg[1, 0:2])/np.sum(wcg[1, 2:8])-0.2*(1+np.sum(wcg[1, 0:2])/np.sum(wcg[1, 2:8])))
        wcg[2, 0] = CGw_MAC + xlemac
        iter = abs(xlemacold/xlemac - 1)

    ac_datafile["Geometry"]["XLEMAC_m"] = xlemac
    return wcg, CGlist, xlemac

if __name__ == "__main__":
    init = time.process_time()
    print(iterate_cg_lg(aircraft_data)[2])
    total = time.process_time() - init
    print(total)