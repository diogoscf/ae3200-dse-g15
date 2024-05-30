import numpy as np
import sys
import os
import math

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from aircraft_data import aircraft_data, c206_data
from CO2_Calculator.conceptual_co2 import calculate_new_co2, calculate_mission_freqs
from unit_conversions import nm_to_m

vol_jet_a1_price = 0.8 #US$/L
jet_a1_dens = 0.8025 #kg/m3

def hourly_operating_cost(mission_file, standard_aircraft_data = c206_data, aircraft_data = aircraft_data):
    """
    Calculate the hourly operating cost of a new aircraft design based on the mission range. Overhaul taken 
    from MAF values, so not yet dependent on the mission profile legs. Should be improved later. 
    Maintenance _is_ dependent on flight time already though. Fuel cost based on inaccurate estimate of 
    fuel burn (assuming cruise speed for all 600 nm), improved fuel burn estimate shall be used.  

    Parameters
    ----------
    mission_file : str
        The name of the mission file to read.
    standard_ac_data : dict
        The standard aircraft data to use for the calculation. Will use C206 data by default.
    ac_data : dict
        The aircraft data to use for the calculation. Will use main new design data by default.

    Returns
    -------
    total_hourly_cost : float
        The total new aircraft operating cost per hour of flight in (US$)
    """
    overhaul_cost = standard_aircraft_data["overhaul_per_hour"] * standard_aircraft_data["Vc_m/s"] / aircraft_data["Performance"]["Vc_m/s"]

    mission_freqs = calculate_mission_freqs(mission_file)
    _, maintenance_cost = calculate_new_co2(mission_freqs, ac_data = aircraft_data, maintenance_standard_co2 = None, V_standard_kts = None, standard_ac_data = standard_aircraft_data)

    endurance = nm_to_m(aircraft_data["Performance"]["range_nm"])/aircraft_data["Performance"]["Vc_m/s"]/3600

    fuel_burn = aircraft_data["Weights"]["MFW_N"]/9.80665/endurance
    fuel_cost = fuel_burn*vol_jet_a1_price/jet_a1_dens

    return overhaul_cost + maintenance_cost + fuel_cost


if __name__ == "__main__":
    print(hourly_operating_cost(f"maf_mission_graph.csv", c206_data, aircraft_data)) 
    


    # V&V
    input = {
        "mission_file": f"maf_mission_graph.csv", 
        "standard_aircraft_data": {'name': 'c206', 'pretty_name': 'Cessna 206', 'Vc_m/s': 73, 'fuel_burn_L/h': 62.5, 'fuel_density_kg/L': 0.717, 'CO2_emissions_kg/kg': 3.05, 'co2_fuel_%': 89, 'fuel_per_hour': 144, 'maintenance_per_hour': 129, 'overhaul_per_hour': 55},
        "aircraft_data": {'name': 'final_design', 'pretty_name': 'Final Design', 'MTOW_N': 15269, 'OEW_N': 8816, 'MFW_N': 1390, 'Wpl_des_kg': 540, 'Wpl_max_kg': 630, 'CLmax_clean': 1.6, 'CLmax_land': 2.5, 'W/S_N/m2': 618, 'W/P_N/W': 0.118, 'Vc_m/s': 60, 'MGC_m': 1.93, 'CLalpha': 6.24, 'E_bat_Wh': 186451, 'E_fuel_Wh/kg': 11972, 'η_bat': 0.85, 'η_generator': 0.45, 'DoD_bat': 0.8, 'η_electricmotor': 0.925, 'η_powertrain': 0.9216, 'CD0': 0.02, 'AR': 9.38, 'e': 0.82, 'range_nm': 600, 'P_req_cruise_W': 217309, 'CO2_emissions_kg/kg': 3.16}
    }

    assert(
            math.isclose(182.7,
                        hourly_operating_cost(input["mission_file"], input["standard_aircraft_data"], input["aircraft_data"]),
                        abs_tol=0.05)
    )
