# NOT WORKING YET!!!!!!!!!!!!!!!!!

import numpy as np
import csv
import os
import json
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data, c206_data

# Black Magic Maintenance Costs
a_factor = 68.37332366 # [$/h] at 1 flight hour
exponent = -0.423

maintenance_cost_per_hour = lambda FT: a_factor * (FT ** exponent)
maintenance_cost_per_hour = np.vectorize(maintenance_cost_per_hour)


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
    """Calculate the CO2 emissions of the a standard aircraft based on the mission profile.

    Parameters
    ----------
    mission_freqs : np.ndarray
        An array of the frequencies of each mission.
        Each row is formatted as [mission length, number, relative frequency].
    ac_data : dict
        The aircraft data to use for the calculation. Will use C206 data by default.
    design_range : float
        The design range of the aircraft (in nautical miles). Default is 600.
    Returns
    -------
    co2 : float
        The CO2 emissions of the Cessna 206 based on the mission profile (in kgCO2)
    """

    vc_kts = ac_data["Vc_m/s"]
    flight_time_h = np.array([mission_freqs[0] * 2 / vc_kts for mission in mission_freqs])
    fuel_emissions_return = flight_time_h * ac_data["fuel_burn_L/h"] * ac_data["fuel_density_kg/L"] * ac_data["CO2_emissions_kg/kg"]
    maintenance_cost_return = maintenance_cost_per_hour(flight_time_h) # V&V Note: the weighted average of this should be equal to $96 for a C206
    relevant_idx = np.where(mission_freqs[0]*2 <= design_range)
    weighted_fuel_co2 = np.sum(fuel_emissions_return * mission_freqs[:, 2])

    return co2

def calculate_co2_reduction(mission_file):
    """Calculate the CO2 reduction of the aircraft based on the mission profile.

    Parameters
    ----------
    mission_file : str
        The name of the mission file to read.
    Returns
    -------
    co2_reduction : float
        The CO2 reduction of the aircraft based on the mission profile (in %)
    """

    mission_freqs = calculate_mission_freqs(mission_file)
    c206_co2 = calculate_standard_co2(mission_freqs, c206_data)


if __name__ == "__main__":
    print(calculate_co2_reduction("maf_mission_graph.csv"))