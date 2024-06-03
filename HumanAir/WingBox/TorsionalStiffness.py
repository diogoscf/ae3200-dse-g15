import numpy as np
from Functions import chord, import_data, import_data2, interp

x_pos = np.array([0.15, 0.5])

def spars(chord, x_pos):
    chord = np.array(chord).reshape((len(chord),1))
    return chord * x_pos 

def MOI(y, chord, spars):
    return

########## input ###########
Sw = 39  # [m2]
taper_ratio = 0.4
Cl_DATA = import_data2('HumanAir\WingBox\Cl_DATA.txt')
AoA = -6
n = 5

######## Execution ########
# Assuming chord function takes the arguments as described and returns a chord length
chord_length, y = chord(Sw, taper_ratio, Cl_DATA, AoA, n)

# Calculate the spars positions based on the chord length
spars_positions = spars(chord_length, x_pos)

# airfoil data
df = import_data('HumanAir/WingBox/airfoil.txt')
print(df)
x_pos_flat = np.array(spars_positions.flatten())
print(x_pos_flat)
y_values = np.array([interp(x, df['x'], df['y']) for x in x_pos_flat])

'''
for i in range(len(spars_positions)):
    for j in range(len(spars_positions[i])):
        spars_positions[i][j] = [spars_positions[i][j], interp(spars_positions[i][j], df['x'], df['y'])]
'''