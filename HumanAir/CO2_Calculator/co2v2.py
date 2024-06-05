import numpy as np
import os
import pickle
import sys

# import time
# import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, c206_data
from HumanAir.unit_conversions import m_s_to_kt
from HumanAir.CO2_Calculator.conceptual_co2 import maintenance_cost_per_hour, calculate_mission_freqs

MAINTENANCE_CO2_OVERHAUL = 0.2  # 20% of maintenance CO2 is overhaul

# Caching for SPEEEEEEEEEEEED
relevant_c206_data = None
maintenance_overhaul_co2 = None
vc_check = None  # Check if the cruise speed has changed
mission_freqs_global = None
flight_distribution = None


def get_flight_distribution(dist_file="flight_dist.pickle", size=10000):
    with open(os.path.join(os.path.dirname(__file__), dist_file), "rb") as f:
        dist = pickle.load(f)
    return dist.rvs(size=size)


def calculate_standard_co2_flightdist(flightdist, ac_data=c206_data, design_range=600):
    """
    Calculate the CO2 emissions (per hour) of a standard aircraft based on the mission profile,
    as given by the frequency of flight lengths.

    Parameters
    ----------
    flightdist : np.ndarray
        Simulated lengths of flights from a log-normal distribution
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
    flight_time_h = flightdist[np.where(flightdist <= design_range)] / vc_kts
    fuel_emissions_return = (
        flight_time_h * ac_data["fuel_burn_L/h"] * ac_data["fuel_density_kg/L"] * ac_data["CO2_emissions_kg/kg"]
    )
    fuel_co2 = np.mean(fuel_emissions_return)

    co2 = fuel_co2 / (ac_data["co2_fuel_%"] / 100)
    maintenance_co2 = co2 - fuel_co2

    return co2, maintenance_co2


def calculate_new_co2_flightdist(
    flightdist,
    mission_freqs,
    ac_data=aircraft_data,
    maintenance_standard_co2=None,
    V_standard_kts=None,
    standard_ac_data=c206_data,
    mission_file="maf_mission_graph.csv",
    maintenance_overhaul_co2=None,
):
    """
    Calculate the CO2 emissions (per hour) of a new aircraft design based on the mission profile,
    as given by the frequency of flight lengths.

    Parameters
    ----------
    flightdist : np.ndarray
        Simulated lengths of flights from a log-normal distribution
    mission_freqs : np.ndarray
        The frequency of flights of different lengths. Will calculate it if not provided.
    ac_data : dict
        The aircraft data to use for the calculation. Will use main new design data by default.
    maintenance_standard_co2 : float or None
        The maintenance CO2 emissions of the standard aircraft.
        If None, will calculate it from `standard_ac_data`.
    V_standard_kts : float or None
        The cruise speed of the standard aircraft. If None, will fetch it from `standard_ac_data`.
    standard_ac_data : dict
        The standard aircraft data to use for the calculation. Will use C206 data by default.
    mission_file: str, default "maf_mission_graph.csv"
        The file containing the mission profile data. Default is "maf_mission_graph.csv".
        Only relevant if `mission_freqs` is None.
    maintenance_overhaul_co2 : float or None
        The CO2 emissions of the overhaul and maintenance. If None, will calculate it.

    Returns
    -------
    co2 : float
        The CO2 emissions of the new aircraft design based on the mission profile (in kgCO2).
    maintenance_overhaul_co2 : float
        The maintenance and overhaul CO2 emissions of the new aircraft design (in kgCO2).
    """

    vc_kts = m_s_to_kt(ac_data["Performance"]["Vc_m/s"])
    flight_time_h = flightdist[np.where(flightdist <= ac_data["Performance"]["range_nm"])] / vc_kts
    battery_usage_ratio = ac_data["Power_prop"]["E_bat_Wh"] / (
        flight_time_h * ac_data["Power_prop"]["P_req_cruise_W"]
    )  # TODO: Check if E_bat_Wh is before or after efficiency
    battery_usage_ratio[battery_usage_ratio > 1] = 1
    fuel_required_return_kg = (
        flight_time_h
        * ac_data["Power_prop"]["P_req_cruise_W"]
        * (1 - battery_usage_ratio)
        / (ac_data["Power_prop"]["E_fuel_Wh/kg"] * ac_data["Power_prop"]["eta_generator"])
    )  # TODO: Check if power train efficiency matters
    fuel_emissions_return = fuel_required_return_kg * ac_data["Performance"]["CO2_emissions_kg/kg"]

    if maintenance_standard_co2 is None:
        _, maintenance_standard_co2 = calculate_standard_co2_flightdist(
            flightdist, standard_ac_data, ac_data["Performance"]["range_nm"]
        )
    if V_standard_kts is None:
        V_standard_kts = m_s_to_kt(standard_ac_data["Vc_m/s"])

    if maintenance_overhaul_co2 is None:
        FT_ratio = V_standard_kts / vc_kts  # new_ac_FT / standard_ac_FT
        FT_maintenance_h = mission_freqs[:, 0] / vc_kts
        maintenance_cost_return_standard = maintenance_cost_per_hour(
            FT_maintenance_h / FT_ratio
        )  # V&V Note: the weighted average of this should be equal to $96 for a C206
        maintenance_cost_return_new = maintenance_cost_per_hour(FT_maintenance_h)

        overhaul_co2 = MAINTENANCE_CO2_OVERHAUL * maintenance_standard_co2 * FT_ratio

        maintenance_cost_return_standard_avg = np.sum(maintenance_cost_return_standard * mission_freqs[:, 2]) / np.sum(
            mission_freqs[:, 2]
        )
        maintenance_cost_return_new_avg = np.sum(maintenance_cost_return_new * mission_freqs[:, 2]) / np.sum(
            mission_freqs[:, 2]
        )
        cost_ratio = maintenance_cost_return_new_avg / maintenance_cost_return_standard_avg
        maintenance_co2 = (1 - MAINTENANCE_CO2_OVERHAUL) * maintenance_standard_co2 * cost_ratio * FT_ratio
        maintenance_overhaul_co2 = overhaul_co2 + maintenance_co2

    fuel_co2 = np.mean(fuel_emissions_return)
    co2 = fuel_co2 + maintenance_overhaul_co2
    return co2, maintenance_overhaul_co2


def calculate_co2_reduction_flightdist(
    ac_data=aircraft_data,
    standard_ac_data=c206_data,
    flightdist=None,
    flightdist_file="flightdist.pickle",
):
    """Calculate the CO2 reduction of the aircraft based on average flight length.

    Parameters
    ----------
    ac_data : dict, optional
        The aircraft data to use for the calculation. Will use main new design data by default.
    standard_ac_data : dict, optional
        The standard aircraft data to use for the calculation. Will use C206 data by default.
    flightdist : np.ndarray or None, optional
        Distribution of flight lengths. If None, will use the default distribution,
        or calculate it if not yet available.
    flightdist_file : str, optional
        The file containing the flight length distribution. Default is "flightdist.pickle".
        Only relevant if `flightdist` is None.

    Returns
    -------
    co2_reduction : float
        The CO2 reduction (i.e. ratio) of the aircraft based on the mission profile.
    """
    global relevant_c206_data, maintenance_overhaul_co2, mission_freqs_global, vc_check, flight_distribution

    if flightdist is None:
        if flight_distribution is None or flightdist_file != "flightdist.pickle":
            flight_distribution = get_flight_distribution(flightdist_file, 10000)

        flightdist = flight_distribution

    if relevant_c206_data is None or standard_ac_data != c206_data:
        c206_co2, c206_maintenance_co2 = calculate_standard_co2_flightdist(flightdist, c206_data)
        relevant_c206_data = (c206_co2, c206_maintenance_co2)
    else:
        c206_co2, c206_maintenance_co2 = relevant_c206_data

    if vc_check != m_s_to_kt(ac_data["Performance"]["Vc_m/s"]):
        passed_overhaul_co2 = None
    else:
        passed_overhaul_co2 = maintenance_overhaul_co2

    if mission_freqs_global is None:
        mission_freqs_global = calculate_mission_freqs("maf_mission_graph.csv")

    new_co2, mo_co2 = calculate_new_co2_flightdist(
        flightdist,
        mission_freqs_global,
        ac_data=ac_data,
        maintenance_standard_co2=c206_maintenance_co2,
        V_standard_kts=m_s_to_kt(standard_ac_data["Vc_m/s"]),
        standard_ac_data=c206_data,
        mission_file="maf_mission_graph.csv",
        maintenance_overhaul_co2=passed_overhaul_co2,
    )

    maintenance_overhaul_co2 = mo_co2
    vc_check = m_s_to_kt(ac_data["Performance"]["Vc_m/s"])

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

    co2_ratio = calculate_co2_reduction_flightdist(ac_data=aircraft_data, standard_ac_data=c206_data)
    print(f"CO2 reduction: {co2_ratio*100:.2f}%")
