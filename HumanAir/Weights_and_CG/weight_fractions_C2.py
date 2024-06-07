import pandas as pd
import numpy as np
import sys
from math import tan, sqrt, floor, ceil
import time
import os
import matplotlib.pyplot as plt

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
    OEW_cont = wcg[1, -1] = ac_datafile["CL2Weight"]["OEW"] / 9.81 - ac_datafile["CL2Weight"]["W_pilot"] / 9.81
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

        if con.lower() == "n":
            sys.exit("Weight fractions do not add up")

    # Return fractions and masses of each component: Wing, MLG, pwtr, NLG, fus, emp, FE, bat, EW
    return wcg

def potato_diagrams(ftb_list, btf_list, Xcg_OEW, xlemac, name, ac_datafile=aircraft_data, plot = False):
    masslist_1 = [ac_datafile["CL2Weight"]["OEW"] / 9.81 - ac_datafile["CL2Weight"]["W_pilot"] / 9.81]
    CGlist_1 = [Xcg_OEW]
    for i in range(ftb_list.shape[1]):
        CGlist_1.append((masslist_1[i]*CGlist_1[i]+ftb_list[0, i]*ftb_list[1, i])/(masslist_1[i]+ftb_list[1, i]))
        masslist_1.append(masslist_1[i]+ftb_list[1, i])

    masslist_2 = [ac_datafile["CL2Weight"]["OEW"] / 9.81 - ac_datafile["CL2Weight"]["W_pilot"] / 9.81]
    CGlist_2 = [Xcg_OEW]
    for i in range(btf_list.shape[1]):
        CGlist_2.append((masslist_2[i] * CGlist_2[i] + btf_list[0, i] * btf_list[1, i]) / (masslist_2[i] + btf_list[1, i]))
        masslist_2.append(masslist_2[i] + btf_list[1, i])

    CGfwd = min(CGlist_1 + CGlist_2) - 0.02*ac_datafile["Aero"]["MAC_wing"]
    CGaft = max(CGlist_1 + CGlist_2) + 0.02*ac_datafile["Aero"]["MAC_wing"]
    CGfwdMAC = (CGfwd - xlemac)/ac_datafile["Aero"]["MAC_wing"]
    CGaftMAC = (CGaft - xlemac)/ac_datafile["Aero"]["MAC_wing"]

    if plot == True:
        plt.plot(((np.array(CGlist_1) - xlemac)/ac_datafile["Aero"]["MAC_wing"]), masslist_1, marker='o', color='red', markersize=4)
        plt.plot(((np.array(CGlist_2) - xlemac)/ac_datafile["Aero"]["MAC_wing"]), masslist_2, marker='o', color='blue', markersize=4)
        plt.plot([CGfwdMAC, CGfwdMAC], [0, max(masslist_1)], color='dimgrey', linewidth=1, linestyle='dashdot')
        plt.plot([CGaftMAC, CGaftMAC], [0, max(masslist_1)], color='dimgrey', linewidth=1, linestyle='dashdot')
        plt.ylabel(r"Mass [kg]")
        plt.xlabel(r"$X_{cg}$ [MAC]")
        plt.title(name)
        plt.xlim(floor(CGfwdMAC*20)/20, ceil(CGaftMAC*20)/20)
        plt.ylim(masslist_1[0], 1.1*(max(masslist_1 + masslist_2) - masslist_1[0]) + masslist_1[0])
        plt.tight_layout()
        plt.show()

    return [CGfwd, CGaft]

def cg_excursion(Xcg_OEW, xlemac, ac_datafile=aircraft_data, plot=False):
    # TODO: Update Xcg with fuselage sizing values
    # Define passenger Xcg positions, 2 passengers per row
    Xcg_row1 = 4
    Xcg_row2 = 5.1
    Xcg_row3 = 6.2

    # Define pilot/cargo/fuel Xcg
    Xcg_pilot = 2.8
    Xcg_luggage = 7
    Xcg_cargomax = 6.5 # Maximum aft position of full cargo load -> dangerous goods
    Xcg_fuel = xlemac + 0.4 * ac_datafile["Aero"]["MAC_wing"]

    # Define loading dict
    m_pax = ac_datafile["CL2Weight"]["Passenger Mass"]
    m_lug = ac_datafile["CL2Weight"]["Luggage Mass"]
    m_pil = ac_datafile["CL2Weight"]["W_pilot"] / 9.81
    m_pl = ac_datafile["CL2Weight"]["Wpl_w/o_pilot"] / 9.81
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
    CGlist_fcp = potato_diagrams(fcp_1, fcp_2, Xcg_OEW, xlemac, "Fuel/Luggage/Passengers", ac_datafile, plot)
    CGlist_fpc = potato_diagrams(fpc_1, fpc_2, Xcg_OEW, xlemac, "Fuel/Passengers/Luggage", ac_datafile, plot)
    CGlist_cpf = potato_diagrams(cpf_1, cpf_2, Xcg_OEW, xlemac, "Luggage/Passengers/Fuel", ac_datafile, plot)
    CGlist_cfp = potato_diagrams(cfp_1, cfp_2, Xcg_OEW, xlemac, "Luggage/Fuel/Passengers", ac_datafile, plot)
    CGlist_pcf = potato_diagrams(pcf_1, pcf_2, Xcg_OEW, xlemac, "Passengers/Luggage/Fuel", ac_datafile, plot)
    CGlist_pfc = potato_diagrams(pfc_1, pfc_2, Xcg_OEW, xlemac, "Passengers/Fuel/Luggage", ac_datafile, plot)
    CGlist_pmpl = potato_diagrams(pmpl_1, pmpl_2, Xcg_OEW, xlemac, "Full Cargo Load", ac_datafile, plot)

    return (CGlist_fcp + CGlist_fpc + CGlist_cpf + CGlist_cfp + CGlist_pcf + CGlist_pfc + CGlist_pmpl)

def iterate_cg_lg(ac_datafile=aircraft_data, PERCENTAGE=0.2, bat_xcg=0.5, plot=False):
    # Set distance of nosewheel from nose [m]
    nose_distance = 0.4
    nose_loading = 0.08
    MAClimit = 3.2

    # Get fractions, weights, cg
    wcg = component_mass(ac_datafile)

    # Get preliminary component CG locations
    CGw_MAC = 0.4 * ac_datafile["Aero"]["MAC_wing"]

    wcg[2, 0] = ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC  # Wing distance from tip

    # Powertrain
    Weng = ac_datafile["CL2Weight"]["Engine Weight"]
    Wprop = ac_datafile["CL2Weight"]["Propeller Weight"]
    Wmotor = ac_datafile["CL2Weight"]["Electric Motor Weight"]

    wcg[2, 2] = (
        Weng * ac_datafile["Stability"]["Xcg_engine_m"]
        + Wprop * ac_datafile["Stability"]["Xcg_prop_m"]
        + Wmotor * ac_datafile["Stability"]["Xcg_motor_m"]
    ) / (Weng + Wprop + Wmotor)

    wcg[2, 4] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]  # Fuselage distance from tip
    wcg[2, 5] = 0.95 * ac_datafile["Geometry"]["fus_length_m"]  # Empennage distance from tip
    wcg[2, 6] = 0.4 * ac_datafile["Geometry"]["fus_length_m"]  # Fixed equipment distance from tip
    wcg[2, 7] = bat_xcg * ac_datafile["Geometry"]["fus_length_m"]  # Battery distance from tip

    Xcg_OEW = (
        (ac_datafile["Geometry"]["XLEMAC_m"] + CGw_MAC) * wcg[1, 0]
        + np.average(wcg[2, [2, 4, 5, 6, 7]], weights=wcg[1, [2, 4, 5, 6, 7]]) * np.sum(wcg[1, [2, 4, 5, 6, 7]])
    ) / (wcg[1, 0] + np.sum(wcg[1, [2, 4, 5, 6, 7]]))

    # Get preliminary moving CG locations from the nose
    xlemac = ac_datafile["Geometry"]["XLEMAC_m"]
    aftcg = max(cg_excursion(Xcg_OEW, xlemac))

    l_m, l_n, Pmg, Pnw, H_s = find_lg(nose_loading, aftcg, ac_datafile)[0:5]  # Nose loading of 8% initially
    wcg[2, 1] = aftcg + l_m
    wcg[2, 3] = aftcg - l_n
    wcg[2, -1] = Xcg_OEW

    # Iterate on CG and LEMAC positions
    iter = 1.0
    wcg[2, 0] = CGw_MAC + xlemac

    while iter > 0.0001:  # convergence criterion
        # Get CG excursion positions
        xlemacold = xlemac
        Xcg_OEW = np.average(wcg[2, 0:8], weights=wcg[1, 0:8])

        aftcg = max(cg_excursion(Xcg_OEW, xlemac))

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
        cgwg = np.average(wcg[2, 0:2] - xlemac, weights=wcg[1, 0:2])  # wing group cg location
        xlemac = np.average(wcg[2, 2:8], weights=wcg[1, 2:8]) + ac_datafile["Aero"]["MAC_wing"] * (
            (cgwg / ac_datafile["Aero"]["MAC_wing"]) * np.sum(wcg[1, 0:2]) / np.sum(wcg[1, 2:8])
            - PERCENTAGE * (1 + np.sum(wcg[1, 0:2]) / np.sum(wcg[1, 2:8]))
        )

        wcg[2, 0] = CGw_MAC + xlemac
        iter = abs(xlemacold / xlemac - 1)

    ac_datafile["Geometry"]["XLEMAC_m"] = xlemac
    ac_datafile["Landing_gear"]["Xmw_m"] = wcg[2, 1]
    ac_datafile["Landing_gear"]["Xnw_m"] = wcg[2, 3]

    CGlist = cg_excursion(Xcg_OEW, xlemac, ac_datafile, plot=plot)

    return wcg, CGlist, xlemac


# iterate though lemac such that the x_lemac is larger than 3.2m and that the landing gear can fit with the batteries
def optimised_xlemac_landing_gears(ac_data=aircraft_data, percentage=0.5, bat_xcg_init=0.1):

    # initialise the sizing parameter
    sizing = False
    bat_xcg = bat_xcg_init

    # loop to get the optimised lemac and landing gear data with the batteries in the right position and the lemac larger than 3.2m
    while not sizing:

        # assume that the landing gear retracts forward
        forward_retractable = True

        # initialise the check if the sizing is possible
        check_sizing = False

        _, _, xlemac = iterate_cg_lg(ac_datafile=ac_data, PERCENTAGE=percentage, bat_xcg=bat_xcg)

        if xlemac > 3.2:

            # TODO: call the function from fuselage sizing and check if the boxes overlap and returns the lemac

            #if forward_retractable:
                # check if the batteries are in the right position

                    # check_sizing = True
                    # forward_retractable = False

            #if not forward_retractable:

                    # check if the batteries are in the right position
                    # check_sizing = True

            if check_sizing:
                sizing = True
                optimised_xlemac = xlemac
                ac_data["Geometry"]["Xcg_battery_m"] = bat_xcg * ac_data["Geometry"]["fus_length_m"]

            # check if the batteries are too far forward
            if bat_xcg > 0.9:
                break
        else:
            bat_xcg += 0.01
        sizing = True
        optimised_xlemac = xlemac
        ac_data["Geometry"]["Xcg_battery_m"] = bat_xcg * ac_data["Geometry"]["fus_length_m"]

    ac_data['Geometry']["XLEMAC_m"] = optimised_xlemac

    return sizing


def calculate_lh(ac_data = aircraft_data, hinge_chord_percentage = 3/4):
    # lh is defined as the distance from quarter chord location of the wing to the quarter chord location of the horizontal tail
    QCW_mac = ac_data["Geometry"]["XLEMAC_m"] + 0.25 * ac_data["Aero"]["MAC_wing"]

    # get the horizontal stabiliser data from the aircraft data
    AR_h = ac_data["Aero"]["AR_HS"]
    taper_h = ac_data["Aero"]["Taper_HS"]
    c_root_h = ac_data["Aero"]["c_root_HS"]
    b_h = ac_data["Aero"]["b_h"]

    # calculate the leading edge angle of the horizontal stabiliser and the x lemac
    tan_LE_sweep = tan(0) - 4 / AR_h * ((- hinge_chord_percentage * c_root_h) * (1 - taper_h) / (1 + taper_h))

    # calculate where the mac of the horizontal stabiliser wrt the leading edge
    y_mac_h = b_h / 6 * (1 + 2 * taper_h) /(1 + taper_h)
    x_mac_h = y_mac_h * tan_LE_sweep

    # get the quarter chord location of the horizontal stabiliser
    QCH_mac = x_mac_h + 0.25 * ac_data["Aero"]["MAC_HS"] + ac_data['Geometry']["fus_length_m"] - ac_data["Aero"]["c_root_HS"]

    # update the aircraft data with the new lh
    ac_data['Stability']["QCW_to_QCh"] = QCH_mac - QCW_mac


if __name__ == "__main__":
    init = time.process_time()
    print(iterate_cg_lg(aircraft_data, PERCENTAGE=0.2, plot=True))
    total = time.process_time() - init
    calculate_lh(ac_data=aircraft_data, hinge_chord_percentage=3/4)
    print(total)
