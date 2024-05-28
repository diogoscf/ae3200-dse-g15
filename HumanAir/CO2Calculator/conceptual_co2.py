import numpy as np
import csv
import os
import json
import sys
import time

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data, c206_data
from unit_conversions import m_s_to_kt

MAINTENANCE_CO2_OVERHAUL = 0.2 # 20% of maintenance CO2 is overhaul

# Black Magic Maintenance Costs
a_factor = 68.37332366 # [$/h] at 1 flight hour
exponent = -0.423

maintenance_cost_per_hour = lambda FT: a_factor * (FT ** exponent)
maintenance_cost_per_hour = np.vectorize(maintenance_cost_per_hour)

# Allows caching of relevant data for quick iteration
relevant_c206_data = None
calculated_mission_frequencies = None


def calculate_mission_freqs(mission_file):
    """Calculate the frequency of each mission in the mission file.

    Parameters
    ----------
    mission_file : str
        The name of the mission file to read.
    Returns
    -------
    freqs : np.ndarray
        An array of the frequencies of each mission.
        Each row is formatted as [mission length, number, relative frequency].
    """

    with open(os.path.join(os.path.dirname(__file__), mission_file), "r") as f:
        reader = csv.reader(f)
        data = list(reader)

    data = [k for k in data if len(k) > 0]
    known = []
    for row in data:
        if row[1]:
            known.append((float(row[2]) / float(row[3])))

    height_avg = np.mean(known)
    freqs = [[float(r[0]), np.rint(float(r[3]) * height_avg)] for r in data]
    freqs = np.array([[int(f[0]), int(f[1]), f[1] / np.sum(np.array(freqs)[:, 1])] for f in freqs])
    return freqs

def calculate_standard_co2(mission_freqs, ac_data = c206_data, design_range = 600):
    """
    Calculate the CO2 emissions (per hour) of a standard aircraft based on the mission profile.

    Parameters
    ----------
    mission_freqs : np.ndarray
        An array of the frequencies of each mission.
        Each row is formatted as [mission length, number, relative frequency].
    ac_data : dict
        The aircraft data to use for the calculation. Will use C206 data by default.
    design_range : float
        The design range of the (new) aircraft (in nautical miles). Default is 600.
    Returns
    -------
    co2 : float
        The CO2 emissions of the Cessna 206 based on the mission profile (in kgCO2)
    maintenance_co2 : float
        The maintenance CO2 emissions of the Cessna 206
    """

    vc_kts = m_s_to_kt(ac_data["Vc_m/s"])
    flight_time_h = mission_freqs[:,0] * 2 / vc_kts
    fuel_emissions_return = flight_time_h * ac_data["fuel_burn_L/h"] * ac_data["fuel_density_kg/L"] * ac_data["CO2_emissions_kg/kg"]
    relevant_idx = np.where(mission_freqs[:,0]*2 <= design_range)
    weighted_fuel_co2 = np.sum(fuel_emissions_return[relevant_idx] * mission_freqs[relevant_idx, 2])/np.sum(mission_freqs[relevant_idx, 2])

    co2 = weighted_fuel_co2 / (ac_data["co2_fuel_%"] / 100)
    maintenance_co2 = co2 - weighted_fuel_co2
    return co2, maintenance_co2

def calculate_new_co2(mission_freqs, ac_data = aircraft_data, maintenance_standard_co2 = None, V_standard_kts = None, standard_ac_data = c206_data):
    """
    Calculate the CO2 emissions (per hour) of a new aircraft design based on the mission profile.

    Parameters
    ----------
    mission_freqs : np.ndarray
        An array of the frequencies of each mission.
        Each row is formatted as [mission length, number, relative frequency].
    ac_data : dict
        The aircraft data to use for the calculation. Will use main new design data by default.
    maintenance_standard_co2 : float or None
        The maintenance CO2 emissions of the standard aircraft. If None, will calculate it from `standard_ac_data`.
    V_standard_kts : float or None
        The cruise speed of the standard aircraft. If None, will fetch it from `standard_ac_data`.
    standard_ac_data : dict
        The standard aircraft data to use for the calculation. Will use C206 data by default.
    Returns
    -------
    co2 : float
        The CO2 emissions of the new aircraft design based on the mission profile (in kgCO2)
    """

    vc_kts = m_s_to_kt(ac_data["Vc_m/s"])
    flight_time_h = mission_freqs[:,0] * 2 / vc_kts
    battery_usage_ratio = ac_data["E_bat_Wh"]/(flight_time_h * ac_data["P_req_cruise_W"]) # TODO: Check if E_bat_Wh is before or after efficiency

    # print(flight_time_h)
    # print(battery_usage_ratio)
    # print(" ")
    battery_usage_ratio[battery_usage_ratio > 1] = 1
    fuel_required_return_kg = flight_time_h * ac_data["P_req_cruise_W"] * (1 - battery_usage_ratio) / (ac_data["E_fuel_Wh/kg"] * ac_data["n_generator"]) # TODO: Check if power train efficiency matters
    fuel_emissions_return = fuel_required_return_kg * ac_data["CO2_emissions_kg/kg"]

    if maintenance_standard_co2 is None:
        _, maintenance_standard_co2 = calculate_standard_co2(mission_freqs, standard_ac_data, ac_data["range_nm"])
    if V_standard_kts is None:
        V_standard_kts = m_s_to_kt(standard_ac_data["Vc_m/s"])
    
    FT_ratio = (V_standard_kts / vc_kts) # new_ac_FT / standard_ac_FT
    maintenance_cost_return_standard = maintenance_cost_per_hour(flight_time_h/FT_ratio/2) # V&V Note: the weighted average of this should be equal to $96 for a C206
    maintenance_cost_return_new = maintenance_cost_per_hour(flight_time_h/2)

    overhaul_co2 = MAINTENANCE_CO2_OVERHAUL * maintenance_standard_co2 * FT_ratio

    relevant_idx = np.where(mission_freqs[:,0]*2 <= ac_data["range_nm"])
    # Note: the following 2 don't use `relevant_idx` bc the maintenance cost for the 206 only adds up to 96 if it is calculated for the full range
    # In any case, the influence on the results is minimal
    maintenance_cost_return_standard_avg = np.sum(maintenance_cost_return_standard[:] * mission_freqs[:, 2])/np.sum(mission_freqs[:, 2])
    maintenance_cost_return_new_avg = np.sum(maintenance_cost_return_new[:] * mission_freqs[:, 2])/np.sum(mission_freqs[:, 2])
    cost_ratio = maintenance_cost_return_new_avg / maintenance_cost_return_standard_avg 
    maintenance_co2 = (1 - MAINTENANCE_CO2_OVERHAUL) * maintenance_standard_co2 * cost_ratio * FT_ratio

    weighted_fuel_co2 = np.sum(fuel_emissions_return[relevant_idx] * mission_freqs[relevant_idx, 2])/np.sum(mission_freqs[relevant_idx, 2])
    co2 = weighted_fuel_co2 + overhaul_co2 + maintenance_co2
    return co2

def calculate_co2_reduction(mission_file = "maf_mission_graph.csv", ac_data = aircraft_data, standard_ac_data = c206_data):
    """Calculate the CO2 reduction of the aircraft based on the mission profile.

    Parameters
    ----------
    mission_file : str
        The name of the mission file to read. Default is "maf_mission_graph.csv".
    ac_data : dict
        The aircraft data to use for the calculation. Will use main new design data by default.
    standard_ac_data : dict
        The standard aircraft data to use for the calculation. Will use C206 data by default.
    Returns
    -------
    co2_reduction : float
        The CO2 reduction (i.e. ratio) of the aircraft based on the mission profile.
    """
    global relevant_c206_data, calculated_mission_frequencies

    if calculated_mission_frequencies is None or calculated_mission_frequencies[1] != mission_file:
        mission_freqs = calculate_mission_freqs(mission_file)
        calculated_mission_frequencies = mission_freqs, mission_file
    else:
        mission_freqs = calculated_mission_frequencies[0]

    if relevant_c206_data is None or standard_ac_data != c206_data:
        c206_co2, c206_maintenance_co2 = calculate_standard_co2(mission_freqs, c206_data, ac_data["range_nm"])
        relevant_c206_data = (c206_co2, c206_maintenance_co2)
    else:
        c206_co2, c206_maintenance_co2 = relevant_c206_data

    new_co2 = calculate_new_co2(mission_freqs, ac_data, c206_maintenance_co2, m_s_to_kt(standard_ac_data["Vc_m/s"]), standard_ac_data)
    co2_ratio = 1 - (new_co2 / c206_co2)
    return co2_ratio


if __name__ == "__main__":
    # t1 = time.process_time()
    # n = 10000
    # i = 0
    # while i < n:
    #     co2_ratio = calculate_co2_reduction("maf_mission_graph.csv")
    #     i += 1
    # tot1 = time.process_time() - t1

    # i = 0
    # t2 = time.process_time()

    # while i < n:
    #     co2_ratio = calculate_co2_reduction("maf_mission_graph.csv")
    #     i += 1
    
    # tot2 = time.process_time() - t2
    # print(f"Time 1: {tot1:.8f}s, Time 2: {tot2:.8f}s")

    co2_ratio = calculate_co2_reduction("maf_mission_graph.csv")
    print(f"CO2 reduction: {co2_ratio*100:.2f}%")