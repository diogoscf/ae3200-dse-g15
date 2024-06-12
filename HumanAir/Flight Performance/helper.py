import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from isa import isa

def density(h, dT):
    return isa(h, delta_T=dT)[2]

# unused
#def temperature(h, dT):
#    return isa(h, delta_T=dT)[0]