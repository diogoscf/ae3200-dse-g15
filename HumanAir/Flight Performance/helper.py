import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from isa import isa

def density(h, dT):
    return isa(h, dT)[2]