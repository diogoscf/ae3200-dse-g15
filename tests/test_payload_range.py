import pytest
import json
import os
import sys

# import matplotlib.pyplot as plt  # Importing plt to mock the plot setup

# Append parent directory to sys.path to import payload_range_points
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from HumanAir.PayloadRange_Diagrams.payload_range import payload_range_points

# Mock aircraft data for testing
mock_aircraft_data = {
    "MTOW_N": 1000000.0,
    "OEW_N": 500000.0,
    "Wpl_des_kg": 10000.0,
    "Wpl_max_kg": 20000.0,
    "MFW_N": 450000.0,
    "range_nm": 1500.0,
    "E_bat_Wh": 500000.0,
    "η_bat": 0.9,
    "W/P_N/W": 0.2,
    "Vc_m/s": 200.0,
    "E_fuel_Wh/kg": 12000.0,
    "η_generator": 0.85,
    "pretty_name": "Test Aircraft",
}

# Directory where the mock file will be saved
mock_config_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "HumanAir", "Configurations")
mock_config_file = os.path.join(mock_config_dir, "mock_aircraft_data.json")

# Save the mock aircraft data to a JSON file
if not os.path.exists(mock_config_dir):
    os.makedirs(mock_config_dir)
with open(mock_config_file, "w", encoding="utf-8") as f:
    json.dump(mock_aircraft_data, f)


# Define the test function
def test_payload_range_points():
    # Calculate payload range points using the mock data
    results = payload_range_points("mock_aircraft_data.json")

    # Unpack results
    (
        zero_fuel_range_km,
        maxpayload_maxrange_km,
        design_range_km,
        max_fuel_range_km,
        ferry_range_km,
        Wpl_max_N,
        Wpl_des_N,
        Wpl_maxfuel_N,
        pretty_name,
    ) = results

    # Define expected values with a tolerance (relative error)
    expected_zero_fuel_range_km = 6300.0
    expected_maxpayload_maxrange_km = 2100.0
    expected_design_range_km = 2778.0
    expected_max_fuel_range_km = 3450.0
    expected_ferry_range_km = 4800.0
    expected_Wpl_max_N = 20000.0 * 9.80665
    expected_Wpl_des_N = 10000.0 * 9.80665
    expected_Wpl_maxfuel_N = 1000000.0 - 500000.0 - 300000.0
    expected_pretty_name = "Test Aircraft"

    # Test assertions with approximate equality
    assert zero_fuel_range_km == pytest.approx(expected_zero_fuel_range_km, rel=1e-2)
    # assert maxpayload_maxrange_km == pytest.approx(expected_maxpayload_maxrange_km, rel=1e-2)
    assert design_range_km == pytest.approx(expected_design_range_km, rel=1e-2)
    # assert max_fuel_range_km == pytest.approx(expected_max_fuel_range_km, rel=1e-2)
    # assert ferry_range_km == pytest.approx(expected_ferry_range_km, rel=1e-2)
    assert Wpl_max_N == pytest.approx(expected_Wpl_max_N, rel=1e-2)
    assert Wpl_des_N == pytest.approx(expected_Wpl_des_N, rel=1e-2)
    # assert Wpl_maxfuel_N == pytest.approx(expected_Wpl_maxfuel_N, rel=1e-2)
    assert pretty_name == expected_pretty_name

    # Clean up: remove the mock file after testing
    os.remove(mock_config_file)
