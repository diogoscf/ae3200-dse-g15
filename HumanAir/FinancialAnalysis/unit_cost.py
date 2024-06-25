import numpy as np
import os
import pickle
import sys
import matplotlib.pyplot as plt


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, c206_data
from HumanAir.CO2_Calculator.co2v2 import calculate_co2_reduction_flightdist


def unit_cost(ref_data=c206_data, ac_data=aircraft_data, per_co2_before=0.37, per_co2_after=0.7):

    distance_lst = np.arange(100, 610, 100)
    frq_lst = [70.07, 18.47, 5.19, 3.13, 1.53, 0.61]
    # frq_lst = [50.07, 38.47, 5.19, 3.13, 1.53, 0.61]

    # frq_lst = [100, 0 , 0, 0, 0, 0]
    # frq_lst = [0, 0, 0, 0, 0, 100]

    idx = 0
    new_unit_cost = 0
    ref_co2 = 0
    ac_co2 = 0
    ac_co2_new_bat = 0
    new_unit_cost_new_bat = 0

    for distance in distance_lst:
        distance_m = distance / 100 * 185200

        ac_flight_time_h = distance_m / ac_data["Performance"]["Vc_m/s"] / 3600
        ref_flight_time_h = distance_m / ref_data["Vc_m/s"] / 3600

        new_unit_cost += (frq_lst[idx] / 100) * ((ac_data["Power_prop"]["E_bat_Wh"] / 1000 * 0.29 / ac_flight_time_h))

        new_unit_cost_new_bat += (frq_lst[idx] / 100) * (
            (ac_data["Power_prop"]["E_bat_Wh"] * 685 / 350 / 1000 * 0.29 / ac_flight_time_h)
        )

        idx += 1

    co2_ratio = calculate_co2_reduction_flightdist(ac_data=ac_data)
    ac_data["Power_prop"]["E_bat_Wh"] = ac_data["Power_prop"]["E_bat_Wh"] * 685 / 350
    co2_ratio_new_bat = calculate_co2_reduction_flightdist(ac_data=ac_data)

    print("Reduction CO2 w.r.t reference a/c with 350 Wh/kg batteries in 2030: ", round(co2_ratio * 100), "%")
    print("Hourly cost a/c with 350 Wh/kg batteries in 2030: ", 210 + round(1.07 * new_unit_cost, 0))
    print("Reduction CO2 w.r.t reference a/c with 685 Wh/kg batteries in 2035: ", round(co2_ratio_new_bat * 100), "%")
    print("Hourly cost a/c with 685 Wh/kg batteries in 2035: ", 210 + round(1.07 * new_unit_cost_new_bat, 0))


if __name__ == "__main__":
    unit_cost(per_co2_before=1.0, per_co2_after=1.0)
