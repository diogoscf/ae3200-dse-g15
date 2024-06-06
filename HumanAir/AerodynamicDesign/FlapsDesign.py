import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data

def flaps_design(ac_data=aircraft_data):
    CLmax_clean = ac_data["Aero"]["CLmax_clean"]
    CLmax_Land = ac_data["Aero"]["CLmax_Land"]
    CLmax_TO=ac_data["Aero"]["CLmax_TO"]
    
    deltaCLmax = CLmax_Land-CLmax_clean
    CLalpha_clean=0.1096 # update this if available after CFD simulations

    Swf = ac_data["Aero"]["S_Wing"] * deltaCLmax / 0.9 / 0.9 # plain flap has a deltaClmax=0.9

    if Swf> 0.8 * ac_data["Aero"]["S_Wing"]:
        raise Exception("Not enough space for flaps")
    
    cf_c = 0.25 # value valid only for plain flap

    ac_data["Flaps"]["Swf"]=Swf # update dictionary
    ac_data["Flaps"]["cf_c"]=cf_c

    delta_alphaL0_landing_airfoil = -15 # assumed value taken from ADSEE-II
    delta_alphaL0_takeoff_airfoil = -10 # assumed value taken from ADSEE-II
    alphaL0_airfoil = -7.85

    # calculate with how much the lift curve shifts to the left
    delta_alphaL0_landing = delta_alphaL0_landing_airfoil * Swf / ac_data["Aero"]["S_Wing"]
    delta_alphaL0_takeoff = delta_alphaL0_takeoff_airfoil * Swf / ac_data["Aero"]["S_Wing"]

    # get deflection from dictionary
    deflection_landing = ac_data["Flaps"]["deflection_landing"]
    deflection_takeoff = ac_data["Flaps"]["deflection_takeoff"]

    # define deltac/cf vs flap deflection curve using ADSEE-II
    plain_flap_slope = (0.17 - 0.0)/(45-0) #based on the graph from ADSEE-II summaries 
    starting_point_y = 0.0 # at 10 degree flap deployment
    starting_point_x = 0 # at 10 degree flap deployment

    # get cprime/c ratio
    cprime_c_landing = 1 + cf_c * (starting_point_y + (deflection_landing - starting_point_x) * plain_flap_slope)
    cprime_c_takeoff = 1 + cf_c * (starting_point_y + (deflection_takeoff - starting_point_x) * plain_flap_slope)
    
    Sprime_S_landing = 1 + Swf/ac_data["Aero"]["S_Wing"] * (cprime_c_landing - 1)
    Sprime_S_takeoff = 1 + Swf/ac_data["Aero"]["S_Wing"] * (cprime_c_takeoff - 1)

    CLalpha_flapped_landing = Sprime_S_landing * CLalpha_clean
    CLalpha_flapped_takeoff = Sprime_S_takeoff * CLalpha_clean

    # calculate required AoA in order to land and take-off
    AoA_landing = (CLmax_clean + deltaCLmax)/CLalpha_flapped_landing - abs(alphaL0_airfoil + delta_alphaL0_landing)
    AoA_takeoff = (CLmax_clean - CLmax_Land + CLmax_TO + deltaCLmax)/CLalpha_flapped_takeoff - abs(alphaL0_airfoil + delta_alphaL0_takeoff)

    if AoA_landing > 15:
        raise Exception("AoA for landing is higher than the scrap angle. Try changing the deflection in design.json")
    if AoA_takeoff > 15:
        raise Exception("AoA for take-off is higher than the scrap angle. Try changing the deflection in design.json")

    # updating the dictionary with the remainign calculated values
    ac_data["Flaps"]["cprime_c_landing"]=cprime_c_landing
    ac_data["Flaps"]["cprime_c_takeoff"]=cprime_c_takeoff
    ac_data["Flaps"]["AoA_landing"] = AoA_landing
    ac_data["Flaps"]["AoA_takeoff"] = AoA_takeoff
    return ac_data

if __name__ == "__main__":
    dict_test=flaps_design()




