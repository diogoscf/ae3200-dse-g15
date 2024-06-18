import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from HumanAir.StructuralAnalysis.LoadDistributions import interpolate_Cl_Cd_Cm


def test_get_deflection():
    assert False


def test_get_twist():
    assert False


def test_force_distribution():
    assert False


def test_moment_distribution():
    assert False


def test_weight_distribution():
    assert False


def test_axial_distribution_ground():
    assert False


def test_axial_distribution_flight():
    assert False


def test_read_points_from_load_dist():
    assert False


def test_strut_error_calculation():
    assert False


def test_strut_force():
    assert False


def test_strut_error_calc_w_args():
    assert False


def test_internal_loads():
    assert False


def test_interpolate_cl_cd_cm():
    Cl_data = {0: {"y_span": [0.1, 0.2, 0.3, 0.4, 0.5], "coefficient": [0.11, 0.21, 0.31, 0.41, 0.51]}}

    Cdi_data = {0: {"y_span": [0.01, 0.02, 0.03, 0.04, 0.05], "coefficient": [0.011, 0.021, 0.031, 0.041, 0.051]}}

    Cm_data = {
        0: {"y_span": [0.001, 0.002, 0.003, 0.004, 0.005], "coefficient": [0.0011, 0.0021, 0.0031, 0.0041, 0.0051]}
    }

    interpolated_values = interpolate_Cl_Cd_Cm(Cl_data, Cdi_data, Cm_data, 1)
    assert interpolated_values == (
        {0: {"coefficient": 0.51, "y_span": 1}},
        {0: {"coefficient": 0.051, "y_span": 1}},
        {0: {"coefficient": 0.0051, "y_span": 1}},
    )
