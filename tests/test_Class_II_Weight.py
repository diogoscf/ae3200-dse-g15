from
def test_class_II_weight_init():
    # Mock input data
    aircraft_data = {
        "Weights": {"MTOW_N": 1, "W_L": 2, "WF_N": 3},
        "Aero": {
            "S_Wing": 1,
            "AR": 6,
            "AR_HS": 4,
            "QuarterChordSweep_Wing_deg": 10,
            "Taper_Wing": 0.5,
            "tc_m_Wing": 0.1,
            "b_Wing": 15,
            "t_root_max_Wing": 0.2,
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
            "P_TO": 13,
            "K_n": 0.24,
            "int": 0.5,
            "N_e": 2,
            "N_t": 1,
        },
        "Landing_gear": {"l_s_m": 14, "l_s_n": 15},
    }

    # Initialize the Class_II_Weight object
    weight_class = Class_II_Weight(aircraft_data)

    # Check if the attributes are correctly initialized
    assert weight_class.W_TO == 2000
    assert weight_class.W_L == 2000
    assert weight_class.W_F == 2000
    assert weight_class.S_Wing == 150
    assert weight_class.S_h == 150
    assert weight_class.S_v == 150
    assert weight_class.n_ult == 3
    assert weight_class.AR_Wing == 6
    assert weight_class.AR_h == 4
    assert weight_class.QuarterChordSweep_Wing == np.deg2rad(10)
    assert weight_class.Taper_Wing == 0.5
    assert weight_class.tc_m_Wing == 50
    assert weight_class.b_Wing == 50
    assert weight_class.t_root_max_Wing == 50
    assert weight_class.V_H == 100
    assert weight_class.V_c == 100
    assert weight_class.QCW_to_QCh == 50
    assert weight_class.l_f_nonosecone == 50
    assert weight_class.p_max == 50
    assert weight_class.N_pax == 4
    assert weight_class.l_f == 50
    assert weight_class.w_f == 50
    assert weight_class.h_f == 50
    assert weight_class.P_TO == 500
    assert weight_class.K_n == 0.24
    assert weight_class.int == 0.5
    assert weight_class.l_s_m == 14
    assert weight_class.l_s_n == 15
    assert weight_class.N_e == 2
    assert weight_class.N_t == 1


def test_wing_weight():
    assert False


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


