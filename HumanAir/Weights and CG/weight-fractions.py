import pandas as pd
import numpy as np
import sys
from math import tan, sqrt
import time
#from aircraft_data import aircraft_data

test_dict = {"cont": 1.2, "MTOW": 25753, "OEW": 17651, "WF": 1922, "WPL": 6180, "Wbat": 5094, "Ww": 2477, "MAC": 1.93, "l_f": 10.09, "h_f": 1.8, "XLEMAC": 3.2}

def find_lg(ac_datafile):
    # Import tyre database
    tires = pd.read_csv('tiredata.csv', index_col=0).to_numpy()

    # Choose smallest available tyre
    Pmw = 0.90 * ac_datafile["MTOW"] / (2 * 9.81)
    Pnw = 0.10 * ac_datafile["MTOW"] / 9.81
    for tire in range(len(tires[:, 0])):
        Wt_m = tires[tire, 0]
        Dw_m = tires[tire, 1]
        if tires[tire, 2] >= Pmw:
            break

    for tire in range(len(tires[:, 0])):
        Wt_n = tires[tire, 0]
        Dw_n = tires[tire, 1]
        if tires[tire, 2] >= Pnw:
            break

    if Pmw > tires[-1, 2]:
        print("WARNING: NO TYRE AVAILBLE")
        con = input("Continue? (y/n): ")
        if con == "n":
            sys.exit("Too heavy for landing gear")

    # Calculate landing gear geometry
    Hcg = 0.5 * ac_datafile["h_f"]
    H_s = 1.5 * Dw_m
    l_m = tan(np.radians(15))*(Hcg+H_s+0.5*Dw_m)
    l_n = l_m * Pmw/Pnw
    ymin = (l_m + l_n)/(sqrt(l_n**2 * tan(np.radians(55))**2 / (Hcg+H_s+0.5*Dw_m)**2 - 1))

    return l_m, l_n, ymin, H_s, Dw_m, Wt_m, Dw_n, Wt_n

def component_mass(ac_datafile):
    # Import statistical weight fraction data
    fracs = pd.read_csv('fraction-database.csv', index_col=0).to_numpy()

    # Convert weights to kg and with contingency
    # TODO: Change to actual datafile
    MTOW_cont = ac_datafile["MTOW"]/9.81*ac_datafile["cont"]
    OEW_cont = ac_datafile["OEW"]/9.81*ac_datafile["cont"]
    Wbat_cont = ac_datafile["Wbat"]/9.81*ac_datafile["cont"]
    Ww_cont = ac_datafile["Ww"]/9.81*ac_datafile["cont"]

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

def iterate_cg_lg(ac_datafile):
    # Get fractions, weights, cg
    wcg = component_mass(ac_datafile)
    WF_cont = ac_datafile["WF"] / 9.81 * ac_datafile["cont"]
    WPL_cont = ac_datafile["WPL"] / 9.81 * ac_datafile["cont"]

    # Get preliminary moving CG locations from the nose
    Xcg_pld = 0.5 * ac_datafile["l_f"]
    Xcg_f = ac_datafile["XLEMAC"] + 0.4 * ac_datafile["MAC"]

    # Get preliminary component CG locations
    CGw_MAC = 0.4 * ac_datafile["MAC"]
    wcg[2, 0] = ac_datafile["XLEMAC"] + CGw_MAC
    wcg[2, 2] = 0.05 * ac_datafile["l_f"]
    wcg[2, 4] = 0.4 * ac_datafile["l_f"]
    wcg[2, 5] = 0.9 * ac_datafile["l_f"]
    wcg[2, 6] = 0.4 * ac_datafile["l_f"]
    wcg[2, 7] = 0.4 * ac_datafile["l_f"]

    Xcg_OEW = ((ac_datafile["XLEMAC"] + CGw_MAC)*wcg[1, 0] + np.average(wcg[2, [2, 4, 5, 6, 7]], weights=wcg[1, [2, 4, 5, 6, 7]]) * np.sum(wcg[1, [2, 4, 5, 6, 7]])) / (wcg[1, 0] + np.sum(wcg[1, [2, 4, 5, 6, 7]]))
    CGlist = [Xcg_OEW, (Xcg_OEW*wcg[1, -1] + Xcg_pld*WPL_cont)/(wcg[1, -1]+WPL_cont), (Xcg_OEW*wcg[1, -1] + Xcg_f*WF_cont)/(wcg[1, -1]+WF_cont), (Xcg_OEW*wcg[1, -1] + Xcg_pld*WPL_cont + Xcg_f*WF_cont)/(wcg[1, -1]+WPL_cont+WF_cont)]
    aftcg = np.max(CGlist)

    l_m, l_n = find_lg(ac_datafile)[0:2]
    wcg[2, 1] = aftcg + l_m
    wcg[2, 3] = aftcg - l_n
    wcg[2, -1] = Xcg_OEW

    # Iterate on CG and LEMAC positions
    iter = 1.0
    xlemac = ac_datafile["XLEMAC"]
    wcg[2, 0] = CGw_MAC + xlemac
    while iter > 0.0001:
        xlemacold = xlemac
        Xcg_OEW = np.average(wcg[2, 0:8], weights=wcg[1, 0:8])
        Xcg_f = xlemac + 0.4 * ac_datafile["MAC"]
        CGlist = [Xcg_OEW, (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont) / (wcg[1, -1] + WPL_cont), (Xcg_OEW * wcg[1, -1] + Xcg_f * WF_cont) / (wcg[1, -1] + WF_cont), (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont + Xcg_f * WF_cont) / (wcg[1, -1] + WPL_cont + WF_cont)]
        aftcg = np.max(CGlist)
        l_m, l_n = find_lg(ac_datafile)[0:2]
        wcg[2, 1] = aftcg + l_m
        wcg[2, 3] = aftcg - l_n
        wcg[2, -1] = Xcg_OEW
        cgwg = np.average(wcg[2, 0:2] - xlemac, weights=wcg[1, 0:2])
        xlemac = np.average(wcg[2, 2:8], weights=wcg[1, 2:8])+ac_datafile["MAC"]*((cgwg/ac_datafile["MAC"])*np.sum(wcg[1, 0:2])/np.sum(wcg[1, 2:8])-0.2*(1+np.sum(wcg[1, 0:2])/np.sum(wcg[1, 2:8])))
        wcg[2, 0] = CGw_MAC + xlemac
        iter = abs(xlemacold/xlemac - 1)
    return wcg, CGlist, xlemac

if __name__ == "__main__":
    init = time.process_time()
    print(iterate_cg_lg(test_dict)[2])
    total = time.process_time() - init
    print(total)