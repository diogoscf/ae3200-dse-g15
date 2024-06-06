import numpy as np
from Functions import chord, import_data2

x_pos = np.array([0.15, 0.5])


def spars(chord, x_pos):
    # x_pos = position of spars in list
    chord = np.array(chord)
    array = []
    for i in range(len(x_pos)):
        array.append(x_pos[i] * chord)
        print(chord)
    return array


def MOI(y, chord, spars):
    return


########## input ###########
Sw = 39  # [m2]
taper_ratio = 0.4
Cl_DATA = import_data2("HumanAir\WingBox\Cl_DATA.txt")
AoA = -6
n = 5


# Assuming chord function takes the arguments as described and returns a chord length
chord_length = chord(Sw, taper_ratio, Cl_DATA, AoA, n)

# Calculate the spars positions based on the chord length
spars_positions = spars(chord_length, x_pos)

# print(spars_positions)
