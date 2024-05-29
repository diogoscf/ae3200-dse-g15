import numpy as np
import matplotlib.pyplot as plt
import json


class WeightEstm:
    def __init__(self, dict):
        self.dict = dict

    def OEW_prime(self):
        return self.dict["A"] * self.dict["MTOW"] + self.dict["B"]

    def PowertrainWeight(self, bat):
        return 9.81 * self.dict["MTOW"] / 1000 / self.dict["W/P"] / self.dict["eta_gb"] / self.dict["eta_ptr"] * (
                    1 - bat) / self.dict["ptr_specific_power"]

    def BatteryWeight(self, bat):
        return 9.81 * self.dict["MTOW"] / self.dict["W/P"] * self.dict["endurance"] * bat / self.dict[
            "bat_specific_energy"] / self.dict["eta_bat"] / self.dict["DoD"] / self.dict["eta_EM"]

    def FuelWeight(self, bat):
        return 9.81 * self.dict["MTOW"] / self.dict["W/P"] * (1 - bat) * self.dict["endurance"] / self.dict[
            "fuel_specific_energy"] / self.dict["eta_thermal_engine"]

    def WingWeight(self):
        return self.dict["Aw"] * self.dict["MTOW"] + self.dict["Bw"]

    def Iterations(self, bat):
        MTOW_new = 0
        MTOW_old = self.dict["MTOW"]
        ok = False

        while np.abs((MTOW_new - self.dict["MTOW"]) / self.dict["MTOW"]) > 0.02:
            if ok:
                self.dict["MTOW"] = MTOW_new

            OEW_prime = self.OEW_prime()
            PowertrainWeight = self.PowertrainWeight(bat)
            BatteryWeight = self.BatteryWeight(bat)
            FuelWeight = self.FuelWeight(bat)
            WingWeight = self.WingWeight()

            MTOW_new = OEW_prime + PowertrainWeight + BatteryWeight + FuelWeight + WingWeight + self.dict["Payload"]

            if MTOW_new > 8000:
                break

            ok = True

        if MTOW_new < 8000:
            self.dict["MTOW"] = MTOW_old
            return (
                MTOW_new,
                self.dict["cont"] * MTOW_new,
                self.dict["cont"] * OEW_prime,
                self.dict["cont"] * PowertrainWeight,
                self.dict["cont"] * BatteryWeight,
                self.dict["cont"] * FuelWeight,
                self.dict["cont"] * WingWeight,
                self.dict["cont"] * self.dict["Payload"]
            )
        else:
            self.dict["MTOW"] = MTOW_old
            return (
                0,
                self.dict["cont"] * MTOW_new,
                self.dict["cont"] * OEW_prime,
                self.dict["cont"] * PowertrainWeight,
                self.dict["cont"] * BatteryWeight,
                self.dict["cont"] * FuelWeight,
                self.dict["cont"] * WingWeight,
                self.dict["cont"] * self.dict["Payload"]
            )

    def PolynomialRegression(self, bat):
        lst_P = []
        lst_bat = []

        for pbat in bat:
            row = self.Iterations(pbat)
            if row[0] != 0:
                lst_P.append(9.81 * row[1] / self.dict["W/P"])
                lst_bat.append(pbat)

        lst_bat = np.array(lst_bat)
        lst_P = np.array(lst_P)

        # Filter out non-positive values in lst_P
        valid_indices = lst_P > 0
        lst_bat = lst_bat[valid_indices]
        lst_P = lst_P[valid_indices]

        if len(lst_P) == 0:
            #print("No valid power values to fit.")
            return np.array([10,10]), np.array([10,10])

        coeff_exp = np.polyfit(lst_bat, np.log(lst_P), 1)
        coeff_pol = np.polyfit(lst_bat, lst_P, 2)

        y_pol = coeff_pol[0] * lst_bat ** 2 + coeff_pol[1] * lst_bat + coeff_pol[2]
        y_exp = np.exp(coeff_exp[1]) * np.exp(coeff_exp[0] * lst_bat)

        return coeff_exp, coeff_pol


if __name__ == "__main__":
    # Replace 'file_path' with the correct path to your JSON configuration file
    file_path = 'HumanAir/Configurations/conventional - Nicholas.json'

    data=WeightEstm(file_path)

    bat=0.11
    row=data.Iterations(bat)
    print(row)
    data.PolynomialRegression()


            