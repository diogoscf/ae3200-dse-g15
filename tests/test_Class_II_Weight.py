from HumanAir.Class_II_Weight.Class_II_Weight import Class_II_Weight
import numpy as np
from HumanAir.unit_conversions import m_to_ft, N_to_lbs, m_squared_to_ft_squared, m_s_to_kt, W_to_hp
import math
def test_class_II_weight_init():
    # Mock input data
    aircraft_data = {
        "Weights": {"MTOW_N": 1, "W_L": 2, "WF_N": 3},
        "Aero": {
            "S_Wing": 15,
            "AR": 6,
            "AR_HS": 4,
            "AR_v": 6.9,
            "QuarterChordSweep_Wing_deg": 4.2,
            "HalfChordSweep_Wing_deg": 50,
            "QuarterChordSweep_v_deg": 50,
            "Taper_Wing": 0.5,
            "tc_m_Wing": 0.1,
            "S_h": 13,
            "S_v": 11,
            "b_Wing": 89,
            "b_h": 90,
            "b_v": 91,
            't_root_max_Wing': 121,
            't_root_max_h': 122,
            't_root_max_v': 123,



        },
        "Performance": {"n_ult": 3, "V_H": 5, "Vc_m/s": 6, "N_pax": 4},
        "Stability": {"QCW_to_QCh": 7},
        "Geometry": {
            "l_f_nonosecone": 8,
            "fuselage_max_perimeter": 9,
            "fus_length_m": 10,
            "fus_width_m": 11,
            "fus_height_m": 12,
        },
        "Power_prop": {
            "P_TO": 4200,
            "K_n": 0.24,
            "int": 0.5,
            "N_e": 5.5,
            "N_t": 6.6,
        },
        "Landing_gear": {"l_s_m": 14, "l_s_n": 15},
    }

    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    # Check if the attributes are correctly initialized
    assert weight_class.W_TO == N_to_lbs(1)
    assert weight_class.W_L == N_to_lbs(2)
    assert weight_class.W_F == N_to_lbs(3)

    assert weight_class.S_Wing == m_squared_to_ft_squared(15)
    assert weight_class.S_h == m_squared_to_ft_squared(13)
    assert weight_class.S_v == m_squared_to_ft_squared(11)

    assert weight_class.n_ult == 3
    assert weight_class.n_ult_l == 5.7

    assert weight_class.AR_Wing == 6
    assert weight_class.AR_h == 4
    assert weight_class.AR_v == 6.9

    assert weight_class.QuarterChordSweep_Wing == np.deg2rad(4.2)
    assert weight_class.HalfChordSweep_Wing == np.deg2rad(50)
    assert weight_class.QuarterChordSweep_v == np.deg2rad(50)

    assert weight_class.Taper_Wing == 0.5
    assert weight_class.tc_m_Wing == m_to_ft(0.1)

    assert weight_class.b_Wing == m_to_ft(89)
    assert weight_class.b_h == m_to_ft(90)
    assert weight_class.b_v == m_to_ft(91)

    assert weight_class.t_root_max_Wing == m_to_ft(121)
    assert weight_class.t_root_max_h == m_to_ft(122)
    assert weight_class.t_root_max_v == m_to_ft(123)

    assert weight_class.V_H == m_s_to_kt(5)
    assert weight_class.V_c == m_s_to_kt(6)
    assert weight_class.QCW_to_QCh == m_to_ft(7)
    assert weight_class.l_f_nonosecone == m_to_ft(8)
    assert weight_class.p_max == m_to_ft(9)
    assert weight_class.N_pax == 4
    assert weight_class.l_f == m_to_ft(10)
    assert weight_class.w_f == m_to_ft(11)
    assert weight_class.h_f == m_to_ft(12)

    assert weight_class.P_TO == W_to_hp(4200)
    assert weight_class.K_n == 0.24
    assert weight_class.K_p == 1.1
    assert weight_class.K_pg == 1.16
    assert weight_class.K_fsp == 6.55

    assert weight_class.int == 0.5

    assert weight_class.l_s_m == m_to_ft(14)
    assert weight_class.l_s_n == m_to_ft(15)

    assert weight_class.N_e == 5.5
    assert weight_class.N_t == 6.6


def test_wing_weight():
    aircraft_data = {
        "Weights": {"MTOW_N": 1, "W_L": 2, "WF_N": 3},
        "Aero": {
            "S_Wing": 15,
            "AR": 6,
            "AR_HS": 4,
            "AR_v": 6.9,
            "QuarterChordSweep_Wing_deg": 4.2,
            "HalfChordSweep_Wing_deg": 50,
            "QuarterChordSweep_v_deg": 50,
            "Taper_Wing": 0.5,
            "tc_m_Wing": 0.1,
            "S_h": 13,
            "S_v": 11,
            "b_Wing": 89,
            "b_h": 90,
            "b_v": 91,
            't_root_max_Wing': 121,
            't_root_max_h': 122,
            't_root_max_v': 123,

        },
        "Performance": {"n_ult": 3, "V_H": 5, "Vc_m/s": 6, "N_pax": 4},
        "Stability": {"QCW_to_QCh": 7},
        "Geometry": {
            "l_f_nonosecone": 8,
            "fuselage_max_perimeter": 9,
            "fus_length_m": 10,
            "fus_width_m": 11,
            "fus_height_m": 12,
        },
        "Power_prop": {
            "P_TO": 4200,
            "K_n": 0.24,
            "int": 0.5,
            "N_e": 5.5,
            "N_t": 6.6,
        },
        "Landing_gear": {"l_s_m": 14, "l_s_n": 15},
    }

    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    answer1 = 0.002933 * weight_class.S_Wing ** 1.018 * weight_class.AR_Wing ** 2.473 * weight_class.n_ult ** 0.611

    answer2 = 96.948 * ((weight_class.W_TO * weight_class.n_ult / (10 ** 5)) ** 0.65
               * (weight_class.AR_Wing / np.cos(weight_class.QuarterChordSweep_Wing)) ** 0.57
                  * (weight_class.S_Wing / 100) ** 0.61
                        * ((1 + weight_class.Taper_Wing) / 2 * weight_class.tc_m_Wing) ** 0.36
                        * (1 + weight_class.V_H / 500) ** 0.5) ** 0.993

    answer3 = (0.00125 * weight_class.W_TO * (weight_class.b_Wing / np.cos(weight_class.HalfChordSweep_Wing)) ** 0.75
               * (1 + (6.3 * np.cos(weight_class.HalfChordSweep_Wing) / weight_class.b_Wing) ** 0.5)
               * (weight_class.n_ult) ** 0.55
               * (weight_class.b_Wing * weight_class.S_Wing * weight_class.W_TO / weight_class.t_root_max_Wing * np.cos(weight_class.HalfChordSweep_Wing)) ** 0.30)

    assert math.isclose(weight_class.WingWeight()['Average'], (answer1 + answer2 + answer3)/3, rel_tol=1e-3)


def test_empennage_weight():
    assert False


def test_fuselage_weight():
    assert False


def test_nacelle_weight():
    assert False


def test_landing_gear_weight():
    assert False


def test_structure_weight_total():
    assert False


def test_powerplant_weight_total():
    assert False


def test_fixed_equipment_weight_total():
    assert False


