import numpy as np
import json
import os

FILE = "config-1.json"


def nmax(MTOW_lbs):
    return min(2.1 + 24000 / (MTOW_lbs + 10000), 3.8)


def nmin(MTOW_lbs):
    return -0.4 * nmax(MTOW_lbs)


def nult(MTOW_lbs):
    return 1.5 * nmax(MTOW_lbs)


if __name__ == "__main__":
    aircraft_data = json.load(open(os.path.join(os.path.dirname(__file__), FILE), "r"))
    MTOW_lbs = aircraft_data["MTOW_lbs"]
    MTOW_kg = MTOW_lbs * 0.45359237

    if MTOW_kg > 5670:
        raise ValueError(
            "MTOW exceeds 5670 kg, so this is not a CS-23 aircraft. This module is only for CS-23 aircraft."
        )
