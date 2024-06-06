import pandas as pd
import numpy as np
import sys
from math import tan, sqrt
import time
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data

def find_lg(nose_loading, aftcg, ac_datafile=aircraft_data):
    # Import tyre database
    tyre_file = os.path.join(os.path.dirname(__file__), "tiredata.csv")
    tyres = pd.read_csv(tyre_file, index_col=0).to_numpy()

    # Choose smallest available tyre
    nose_loading = nose_loading
    Pmw = (1 - nose_loading) * ac_datafile["CL2Weight"]["MTOW_N"] / (2 * 9.81)
    Pnw = nose_loading * ac_datafile["CL2Weight"]["MTOW_N"] / 9.81  # accounts for additional load from front CG
    Pmg = (1 - nose_loading) * ac_datafile["CL2Weight"]["MTOW_N"] / (9.81)
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
        print("WARNING: NO TYRE AVAILABLE")
        con = input("Continue? (y/n): ")
        if con == "n":
            sys.exit("Too heavy for landing gear")

    # Calculate landing gear geometry
    Hcg = 0.5 * ac_datafile["Geometry"]["fus_height_m"]
    H_s = 1.5 * Dw_m  # initial guess
    l_m = tan(np.radians(16)) * (Hcg + H_s + 0.5 * Dw_m)
    H_strike = (
        ac_datafile["Geometry"]["fus_length_m"] - ac_datafile["Geometry"]["tail_length_m"] - (aftcg + l_m)
    ) * np.tan(np.radians(18))
    iter = 1.0
    while iter > 0.0001:
        H_s = H_strike
        l_m = tan(np.radians(16)) * (Hcg + H_s + 0.5 * Dw_m)
        H_strike = (
            ac_datafile["Geometry"]["fus_length_m"] - ac_datafile["Geometry"]["tail_length_m"] - (aftcg + l_m)
        ) * np.tan(np.radians(18))
        iter = abs(H_s / H_strike - 1)
        if H_strike < 0.6 * Dw_m:
            H_strike = 0.6 * Dw_m
            iter = 0

    H_s = H_strike
    l_n = l_m * Pmg / Pnw
    ymin = (l_m + l_n) / (sqrt(l_n**2 * tan(np.radians(55)) ** 2 / (Hcg + H_s + 0.5 * Dw_m) ** 2 - 1))

    # Write values to dict
    ac_datafile["Landing_gear"]["lm_m"] = l_m
    ac_datafile["Landing_gear"]["ln_m"] = l_n
    ac_datafile["Landing_gear"]["ymin_m"] = ymin
    ac_datafile["Landing_gear"]["Hs_m"] = H_s
    ac_datafile["Landing_gear"]["Dwm_m"] = Dw_m
    ac_datafile["Landing_gear"]["Dwn_m"] = Dw_n
    ac_datafile["Landing_gear"]["Wtm_m"] = Wt_m
    ac_datafile["Landing_gear"]["Wtn_m"] = Wt_n
    ac_datafile["Landing_gear"]["PMW_N"] = Pmw * 9.81
    ac_datafile["Landing_gear"]["PNW_N"] = Pnw * 9.81

    return l_m, l_n, Pmg, Pnw, H_s, Pmw, ymin, Dw_m, Wt_m, Dw_n, Wt_n


def component_mass(ac_datafile=aircraft_data):
    # Import statistical weight fraction data
    fracs_file = os.path.join(os.path.dirname(__file__), "fraction-database.csv")
    fracs = pd.read_csv(fracs_file, index_col=0).to_numpy()

    # Convert weights to kg and with contingency, put them in a list
    wcg = np.zeros((3, 9))

    MTOW_cont = ac_datafile["CL2Weight"]["MTOW_N"] / 9.81
    OEW_cont = wcg[1, -1] = ac_datafile["CL2Weight"]["OEW"] / 9.81
    Ww_cont = wcg[1, 0] = ac_datafile["CL2Weight"]["Wing Weight"] / 9.81
    Wmg = wcg[1, 1] = ac_datafile["CL2Weight"]["Landing Gear Weight"] / 9.81 * 0.772371311
    Wpwtr = wcg[1, 2] = ac_datafile["CL2Weight"]["Total Powerplant Weight"] / 9.81
    Wnw = wcg[1, 3] = ac_datafile["CL2Weight"]["Landing Gear Weight"] / 9.81 * 0.227628689
    Wfus = wcg[1, 4] = ac_datafile["CL2Weight"]["Fuselage Weight"] / 9.81
    Wemp = wcg[1, 5] = ac_datafile["CL2Weight"]["Empennage Weight"] / 9.81
    Wfe = wcg[1, 6] = ac_datafile["CL2Weight"]["Total Fixed Equipment Weight"] / 9.81
    Wbat_cont = wcg[1, 7] = ac_datafile["CL2Weight"]["Wbat_N"] / 9.81

    # Set up weight fractions
    for weight in range(len(wcg[0, :])):
        wcg[0, weight] = wcg[1, weight] / MTOW_cont

    # Check whether fractions make sense
    fracsum = np.sum(wcg[0, 0:-1])
    if abs(fracsum / wcg[0, -1] - 1) > 0.05:
        print("WARNING: WEIGHT FRACTIONS DIFFER MORE THAN 5%")
        print("Expected OEW/MTOW:", wcg[0, -1])
        print("Summed OEW/MTOW:", fracsum)
        con = input("Continue? (y/n): ")
        if con == "n":
            sys.exit("Weight fractions do not add up")

    # Return fractions and masses of each component: Wing, MLG, pwtr, NLG, fus, emp, FE, bat, EW
    return wcg

def cg_excursion(ftb_list, btf_list, Xcg_OEW, ac_datafile=aircraft_data, plot = False):
    masslist = [ac_datafile["CL2Weight"]["OEW"] / 9.81]
    CGlist = [Xcg_OEW]
    for i in ftb_list:
        CGlist.append((masslist[i]*CGlist[i]+ftb_list[0, i]*ftb_list[1, i])/(masslist[i]+ftb_list[0, i]))

def potato_diagrams(Xcg_OEW, ac_datafile=aircraft_data, plot=False):
    # TODO: Update Xcg with fuselage sizing values
    # Define passenger Xcg positions, 2 passengers per row
    Xcg_row1 = 4
    Xcg_row2 = 5.1
    Xcg_row3 = 6.2

    # Define pilot/cargo/fuel Xcg
    Xcg_pilot = 2.8
    Xcg_luggage = 7
    Xcg_cargomax = 6.5 # Maximum aft position of full cargo load -> dangerous goods
    Xcg_fuel = ac_datafile["Geometry"]["XLEMAC_m"] + 0.4 * ac_datafile["Aero"]["MAC_wing"]

    # Define loading dict
    m_pax = ac_datafile["CL2Weight"]["Passenger Mass"]
    m_lug = ac_datafile["CL2Weight"]["Luggage Mass"]
    m_pil = ac_datafile["CL2Weight"]["Pilot Weight"] / 9.81
    m_pl = ac_datafile["CL2Weight"]["Wpl"] - ac_datafile["CL2Weight"]["Pilot Weight"] / 9.81
    m_f = ac_datafile["CL2Weight"]["Wfuel_N"] / 9.81
    Pax_ftb = np.array([[Xcg_row1, Xcg_row1, Xcg_row2, Xcg_row2, Xcg_row3, Xcg_row3], [m_pax, m_pax, m_pax, m_pax, m_pax, m_pax]])
    Pax_btf = np.array([[Xcg_row3, Xcg_row3, Xcg_row2, Xcg_row2, Xcg_row1, Xcg_row1], [m_pax, m_pax, m_pax, m_pax, m_pax, m_pax]])
    pil = np.array([[Xcg_pilot], [m_pil]])
    lug = np.vstack((Xcg_luggage*np.ones(6), m_lug*np.ones(6)))
    pl = np.array([[Xcg_cargomax], [m_pl]])
    fl = np.array([[Xcg_fuel], [m_f]])

    # Define loading situations
    fcp_1 = np.hstack((fl, lug, pil, Pax_ftb))
    fcp_2 = np.hstack((fl, lug, Pax_btf, pil))

    fpc_1 = np.hstack((fl, pil, Pax_ftb, lug))
    fpc_2 = np.hstack((fl, Pax_btf, pil, lug))

    cpf_1 = np.hstack((lug, pil, Pax_ftb, fl))
    cpf_2 = np.hstack((lug, Pax_btf, pil, fl))

    cfp_1 = np.hstack((lug, fl, pil, Pax_ftb))
    cfp_2 = np.hstack((lug, fl, Pax_btf, pil))

    pcf_1 = np.hstack((pil, Pax_ftb, lug, fl))
    pcf_2 = np.hstack((Pax_btf, pil, lug, fl))

    pfc_1 = np.hstack((pil, Pax_ftb, fl, lug))
    pfc_2 = np.hstack((Pax_btf, pil, fl, lug))

    pmpl_1 = np.hstack((pil, pl))
    pmpl_2 = np.hstack((pl, pil))

    # Get CG excursion


def iterate_cg_lg(ac_datafile=aircraft_data, PERCENTAGE=0.2):
    # Set distance of nosewheel from nose [m]
    nose_distance = 0.4
    nose_loading = 0.08
    MAClimit = 3.2

    # Get fractions, weights, cg
    wcg = component_mass(ac_datafile)
    WF_cont = ac_datafile["CL2Weight"]["Wfuel_N"] / 9.81
    WPL_cont = ac_datafile["CL2Weight"]["Wpl"]

    # Get preliminary component CG locations
    CGw_MAC = 0.4 * ac_datafile["Aero"]["MAC_wing"]
    wcg[2, 0] = ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC  # Wing

    # Powertrain
    Weng = ac_datafile["CL2Weight"]["Engine Weight"]
    Wprop = ac_datafile["CL2Weight"]["Propeller Weight"]
    Wmotor = ac_datafile["CL2Weight"]["Electric Motor Weight"]

    wcg[2, 2] = (
        Weng * ac_datafile["Stability"]["Xcg_engine_m"]
        + Wprop * ac_datafile["Stability"]["Xcg_prop_m"]
        + Wmotor * ac_datafile["Stability"]["Xcg_motor_m"]
    ) / (Weng + Wprop + Wmotor)

    wcg[2, 4] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]  # Fuselage
    wcg[2, 5] = 0.95 * ac_datafile["Geometry"]["fus_length_m"]  # Empennage
    wcg[2, 6] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]  # Fixed equipment
    wcg[2, 7] = ac_datafile["Stability"]["Xcg_battery_m"]  # Battery

    Xcg_OEW = (
        (ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC) * wcg[1, 0]
        + np.average(wcg[2, [2, 4, 5, 6, 7]], weights=wcg[1, [2, 4, 5, 6, 7]]) * np.sum(wcg[1, [2, 4, 5, 6, 7]])
    ) / (wcg[1, 0] + np.sum(wcg[1, [2, 4, 5, 6, 7]]))

    # Get preliminary moving CG locations from the nose
    # TODO: updated potato diagram for better CG list
    Xcg_pld = 0.5 * ac_datafile["Geometry"]["fus_length_m"]
    Xcg_f = ac_datafile["Geometry"]["XLEMAC_m"] + 0.4 * ac_datafile["Aero"]["MAC_wing"]

    CGlist = [
        Xcg_OEW,
        (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont) / (wcg[1, -1] + WPL_cont),
        (Xcg_OEW * wcg[1, -1] + Xcg_f * WF_cont) / (wcg[1, -1] + WF_cont),
        (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont + Xcg_f * WF_cont) / (wcg[1, -1] + WPL_cont + WF_cont),
    ]
    aftcg = np.max(CGlist)

    l_m, l_n, Pmg, Pnw, H_s = find_lg(nose_loading, aftcg, ac_datafile)[0:5]  # Nose loading of 8% initially
    wcg[2, 1] = aftcg + l_m
    wcg[2, 3] = aftcg - l_n
    wcg[2, -1] = Xcg_OEW

    # Iterate on CG and LEMAC positions
    iter = 1.0
    xlemac = ac_datafile["Geometry"]["XLEMAC_m"]
    wcg[2, 0] = CGw_MAC + xlemac

    while iter > 0.0001:
        # Get CG excursion positions
        xlemacold = xlemac
        Xcg_OEW = np.average(wcg[2, 0:8], weights=wcg[1, 0:8])
        Xcg_f = xlemac + 0.4 * ac_datafile["Aero"]["MAC_wing"]

        # TODO: updated potato diagram for better CG list

        CGlist = [
            Xcg_OEW,
            (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont) / (wcg[1, -1] + WPL_cont),
            (Xcg_OEW * wcg[1, -1] + Xcg_f * WF_cont) / (wcg[1, -1] + WF_cont),
            (Xcg_OEW * wcg[1, -1] + Xcg_pld * WPL_cont + Xcg_f * WF_cont) / (wcg[1, -1] + WPL_cont + WF_cont),
        ]
        aftcg = np.max(CGlist)
        # Revise nosewheel loading in case wheel is too far forward
        if wcg[2, 3] < nose_distance:
            wcg[2, 3] = nose_distance
            l_n = aftcg - wcg[2, 3]
            nose_loading = 1 / (l_n / l_m + 1)
            if nose_loading > 0.15:
                print("WARNING: TOO MUCH LOAD ON NOSE WHEEL")
                con = input("Continue? (y/n): ")
                if con == "n":
                    sys.exit("Nose wheel could not be placed")

        # Place nosewheel
        l_m, l_n, Pmg, Pnw, H_s = find_lg(nose_loading, aftcg, ac_datafile)[0:5]
        wcg[2, 1] = aftcg + l_m
        wcg[2, 3] = aftcg - l_n
        if wcg[2, 3] > nose_distance:
            wcg[2, 3] = nose_distance
            l_n = aftcg - wcg[2, 3]
            l_m = l_n * Pnw / Pmg
            wcg[2, 1] = aftcg + l_m

        # Update X LEMAC
        wcg[2, -1] = Xcg_OEW
        cgwg = np.average(wcg[2, 0:2] - xlemac, weights=wcg[1, 0:2])
        xlemac = np.average(wcg[2, 2:8], weights=wcg[1, 2:8]) + ac_datafile["Aero"]["MAC_wing"] * (
            (cgwg / ac_datafile["Aero"]["MAC_wing"]) * np.sum(wcg[1, 0:2]) / np.sum(wcg[1, 2:8])
            - PERCENTAGE * (1 + np.sum(wcg[1, 0:2]) / np.sum(wcg[1, 2:8]))
        )
        wcg[2, 0] = CGw_MAC + xlemac
        iter = abs(xlemacold / xlemac - 1)

    ac_datafile["Geometry"]["XLEMAC_m"] = xlemac
    ac_datafile["Landing_gear"]["Xmw_m"] = wcg[2, 1]
    ac_datafile["Landing_gear"]["Xnw_m"] = wcg[2, 3]
    return wcg, CGlist, xlemac

if __name__ == "__main__":
    init = time.process_time()
    print(iterate_cg_lg(aircraft_data, PERCENTAGE=0.2))
    total = time.process_time() - init
    print(total)
    potato_diagrams(4)
