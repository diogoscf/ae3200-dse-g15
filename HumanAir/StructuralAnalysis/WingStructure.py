import numpy as np
import os
import sys
import copy

# import pandas as pd
# import matplotlib.pyplot as plt
# import re

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from HumanAir.aircraft_data import aircraft_data, airfoil_shape

# from Functions import chord, import_data2


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def interp(x, x_given, y_given):
    i = find_nearest(x_given, x)
    if i == 0:
        return y_given[0]
    elif i == len(x_given):
        return y_given[-1]
    else:
        x1, x2 = x_given[i - 1], x_given[i]
        y1, y2 = y_given[i - 1], y_given[i]
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


class WingStructure:
    # TODO: Make it so we can choose between wing and HT
    def __init__(self, ac_data, airfoil_data, nodes=501):
        self.airfoil_data = airfoil_data
        # self.file_path_y = file_path_y  # file path to span
        self.Sw = ac_data["Aero"]["S_Wing"]
        self.taper_ratio = ac_data["Aero"]["Taper_Wing"]
        # self.AoA = AoA  # [deg]
        self.nodes = nodes  # number of nodes
        self.t1_spar = ac_data["Geometry"]["t_spar_tip"]  # [m] thickness at the tip
        self.t2_spar = ac_data["Geometry"]["t_spar_root"]  # [m] thickness at the root
        self.t_skin = ac_data["Geometry"]["t_skin_wing"]  # [m] thickness of skin
        self.spar_pos = np.array(ac_data["Geometry"]["spar_pos"])  # position of spars
        self.stringer_area = ac_data["Geometry"]["wing_stringer_area_m"]  # [m2] area of stringers
        self.cr = ac_data["Aero"]["c_root_wing"]
        self.ct = self.cr * self.taper_ratio
        self.b = ac_data["Aero"]["b_Wing"]
        self.ypts = np.linspace(-self.b / 2, self.b / 2, self.nodes)

        self.chord_distribution = self.calculate_chord_distribution()
        self.spars = self.chord_distribution.reshape(len(self.chord_distribution), 1) * self.spar_pos

        # Assume spar thickness is uniformly decreasing from root to tip
        self.t_spar_dist = self.calc_spar_dist()

        self.airfoil_division = self.calc_airfoil_division()
        self.y_up, self.y_down = self.y_updown()

    # def import_data(self):
    #     df = pd.read_csv(self.file_path, sep="\s+", header=None, names=["x", "y"], skiprows=1)
    #     return df

    # def import_data2(self):
    #     data = {}
    #     with open(self.file_path_y, "r") as file:
    #         lines = file.readlines()
    #         angle_file = lines[0]
    #         angles = re.findall(r"VLM1 -\s*-?\d+\.?\d*", angle_file)
    #         angles = [float(re.search(r"VLM1 -\s*(-?\d+\.?\d*)", angle).group(1)) for angle in angles]

    #         for line in lines[1:]:
    #             values = line.split()
    #             for i in range(0, len(values), 2):
    #                 angle_index = i // 2
    #                 angle = angles[angle_index]
    #                 y_positions = float(values[i])
    #                 lift_distribution = float(values[i + 1])

    #                 if angle not in data:
    #                     data[angle] = {"y_span": [], "coefficient": []}

    #                 data[angle]["y_span"].append(y_positions)
    #                 data[angle]["coefficient"].append(lift_distribution)
    #     return data

    def calc_spar_dist(self):
        return self.t1_spar + (self.t2_spar - self.t1_spar) / (self.cr - self.ct) * (self.chord_distribution - self.ct)

    def calc_airfoil_division(self):
        df = copy.deepcopy(self.airfoil_data)
        for i in range(len(df)):
            if df["x"][i] == 0:
                idx = i
        df_1 = np.array(df.iloc[:idx, :].values.tolist())
        df_2 = np.array(df.iloc[idx:, :].values.tolist())
        return df_1, df_2

    def calculate_chord_distribution(self):
        chord_length = (
            2
            * self.Sw
            / (1 + self.taper_ratio)
            / self.b
            * (1 - ((1 - self.taper_ratio) / self.b * np.abs(2 * self.ypts)))
        )
        return chord_length

    # def chord(self):
    #     Cl_DATA = import_data2(self.file_path_y)
    #     b = Cl_DATA[self.AoA]['y_span'][-1] * 2
    #     y = np.linspace(Cl_DATA[AoA]['y_span'][0], Cl_DATA[AoA]['y_span'][-1], self.n)

    # def spars(self):
    #     return self.chord_distribution * self.spar_pos

    def y_spar(self, df):
        x_spar = self.spars
        y_spar = []
        for i in range(len(x_spar)):
            for j in range(len(x_spar[i])):
                y_spar.append(
                    interp(x_spar[i][j] / self.chord_distribution[i], df[:, 0], df[:, 1]) * self.chord_distribution[i]
                )
        y_spar = np.array(y_spar)
        return y_spar.reshape((len(x_spar), 2))

    def y_updown(self):
        df_up, df_down = self.airfoil_division
        y_up = self.y_spar(df_up)
        y_down = self.y_spar(df_down)
        return y_up, y_down

    def diff(self, data):
        dx = []
        for i in range(len(data)):
            for j in range(len(data[i])):
                if j == len(data[i]) - 1:
                    continue
                else:
                    dxx = data[i][j + 1] - data[i][j]
                    dx.append(dxx)
        return dx

    def y(self, data):
        dx = []
        for i in range(len(data)):
            for j in range(len(data[i])):
                if j == len(data[i]) - 1:
                    continue
                else:
                    dxx = (data[i][j + 1] + data[i][j]) / 2
                    dx.append(dxx)
        return np.array(dx)

    def d_s1s2(self):
        dx = self.diff(self.spars)
        dy_up = self.diff(self.y_up)
        dy_down = self.diff(self.y_down)
        l_box_up = []
        l_box_down = []
        for i in range(len(dx)):
            l_box_up.append(np.sqrt(dx[i] ** 2 + dy_up[i] ** 2))
            l_box_down.append(np.sqrt(dx[i] ** 2 + dy_down[i] ** 2))
        return np.array(l_box_up), np.array(l_box_down)

    def h_s1s2(self):
        h_mid = np.empty((len(self.y_up), len(self.y_up.T)))
        h_emp = np.empty((len(self.y_up), len(self.y_up.T)))
        for i in range(len(self.y_up)):
            for j in range(len(self.y_up[i])):
                h_mid[i, j] = (self.y_up[i][j] + self.y_down[i][j]) / 2
                h_emp[i, j] = self.y_up[i][j] - self.y_down[i][j]
        return h_mid, h_emp

    def centroid(self):
        y_up, y_down = self.y_updown()
        l_box_up, l_box_down = self.d_s1s2()
        h_mid, h_s1s2 = self.h_s1s2()

        # upper skin
        A_uskin = l_box_up.reshape((len(l_box_up), 1)) * self.t_skin
        y_uskin = self.y(y_up).reshape((len(y_up), 1))

        # lower skin
        A_lskin = l_box_down.reshape((len(l_box_down), 1)) * self.t_skin
        y_lskin = self.y(y_down).reshape((len(y_down), 1))

        # left spar
        A_lspar = self.t_spar_dist * h_s1s2[:, 0].reshape((len(h_s1s2), 1))
        y_lspar = h_mid[:, 0].reshape((len(h_mid), 1))

        # right spar
        A_rspar = self.t_spar_dist * h_s1s2[:, 1].reshape((len(h_s1s2), 1))
        y_rspar = h_mid[:, 1].reshape((len(h_mid), 1))

        return (A_uskin * y_uskin + A_lskin * y_lskin + A_lspar * y_lspar + A_rspar * y_rspar) / (
            A_uskin + A_lskin + A_lspar + A_rspar
        )

    def Ixx(self, no_str):
        h_mid, h_s1s2 = self.h_s1s2()
        l_box_up, l_box_down = self.d_s1s2()
        h_15c = h_s1s2[:, 0].flatten()
        h_50c = h_s1s2[:, 1].flatten()
        htot = h_15c + h_50c
        h_avemax = htot / 4
        I_stringer = self.stringer_area * no_str * h_avemax**2
        I_spar = 1 / 12 * self.t_spar_dist * (h_15c * 3 + h_50c * 3)
        I_skin = self.t_skin * (l_box_down + l_box_up) * h_avemax**2
        return I_stringer + I_spar + I_skin


if __name__ == "__main__":
    # INPUT
    # Sw = 34.56  # [m2]
    # taper_ratio = 0.4
    # AoA = -6  # [deg]
    n = 500
    # t1_spar = 0.010  # [m] thickness at the tip
    # t2_spar = 0.025  # [m] thickness at the root
    # t_skin = 0.007  # [m] thickness of skin
    # file_path = "HumanAir/WingBox/airfoil.txt"
    # file_path_y = "HumanAir\WingBox\Cl_DATA.txt"
    # A_str = 0.02
    # Cr = 2.5  # [m] root chord length
    # b = 19.93
    # x_pos = np.array([0.15, 0.5])

    wing_structure_data = WingStructure(aircraft_data, airfoil_shape, nodes=501)

    df = wing_structure_data.airfoil_data
    chord1 = wing_structure_data.chord_distribution
    y = wing_structure_data.ypts
    # plt.plot(df['x']*chord1[1], df['y']*chord1[1])
    # plt.axis('equal')
    # plt.show()
    idx = find_nearest(y, 3.986)
    print(y[idx])

    # print(torisonal_stiffness.spars())
    h_mid, h_s1s2 = wing_structure_data.h_s1s2()
    print(h_s1s2[idx])

    # plt.plot(df_down[:,0], df_down[:,1])
    # plt.plot(df_up[:,0], df_up[:,1])
    # plt.plot(x_spar[0]/chord_length[0], y_spar_down[0]/chord_length[0])
    # plt.plot(x_spar[0]/chord_length[0], y_spar_up[0]/chord_length[0])
    # plt.axis('equal')
    # plt.show()
