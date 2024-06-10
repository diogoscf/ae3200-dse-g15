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
from isa import isa

def cy_y(acd=aircraft_data, y=1):
    return 2 * acd["Aero"]["S_Wing"] / (1 + acd["Aero"]["Taper_Wing"]) / acd["Aero"]["b_Wing"] * (y**2 / 2 - (1 - acd["Aero"]["Taper_Wing"]) / acd["Aero"]["b_Wing"] * 2 * y**3 / 3)

def cy_y2(acd=aircraft_data, y=1):
    return 2 * acd["Aero"]["S_Wing"] / (1 + acd["Aero"]["Taper_Wing"]) / acd["Aero"]["b_Wing"] * (y**3 / 3 - (1 - acd["Aero"]["Taper_Wing"]) / acd["Aero"]["b_Wing"] * 2 * y**4 / 4)


def AileronSizing(acd = aircraft_data):

    # define slope of aileron effectiveness graph taken from ADSEE slides
    start_x = 0.2
    start_y = 0.4
    end_x = 0.4
    end_y = 0.6
    slope = (end_y - start_y) / (end_x - start_x)

    # positon of aileron hinge can be varied but always needs to be placed behind 
    pos_lst = np.arange(0.25,0.36,0.05)

    design_found = False 

    for pos in pos_lst:

        tau = slope * (pos - start_x) + start_y

        CL_delta_a = 2 * acd["Aero"]["cl_alpha_airfoil_deg"] / (180/3.14) * tau / acd["Aero"]["S_Wing"] / acd["Aero"]["b_Wing"] * (cy_y(acd, 0.95 * acd["Aero"]["b_Wing"] / 2) - cy_y(acd, (acd["Flaps"]["flap_end"] + 0.05) * acd["Aero"]["b_Wing"] / 2)) 

        CL_P = -4 * (acd["Aero"]["cl_alpha_airfoil_deg"] / (180/3.14) + acd["Aero"]["cd_0_airfoil"]) / acd["Aero"]["S_Wing"] / acd["Aero"]["b_Wing"] * (cy_y2(acd, acd["Aero"]["b_Wing"]/2) - cy_y2(acd, 0))

        P_rad = - CL_delta_a / CL_P * 14 * (2* acd["Performance"]["Vc_m/s"]) / acd["Aero"]["b_Wing"]
        
        P_deg = P_rad * (180/3.14) 

        turn_time = round(60 / P_deg,2)

        if turn_time < (acd["CL2Weight"]["MTOW_N"] / 9.81 + 200) / 590:
            design_found = True
            acd["Aileron"]["start"] = acd["Flaps"]["flap_end"] + 0.05
            acd["Aileron"]["end"] = 0.95
            acd["Aileron"]["roll_rate_deg"] = P_deg
            acd["Aileron"]["roll_rate_rad"] = P_rad
            acd["Aileron"]["CL_delta_a"] = CL_delta_a
            acd["Aileron"]["CL_P"] = CL_P
            acd["Aileron"]["turn_time"] = turn_time
            acd['Aileron']['hinge_position'] = pos
            break
    
    if design_found:
        return acd
    else:
        raise Exception("No suitable design found.")

if __name__ == "__main__":
    AileronSizing()
