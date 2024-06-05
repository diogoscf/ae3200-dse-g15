import json
import os

# from types import SimpleNamespace as Namespace

with open(os.path.join(os.path.dirname(__file__), "Configurations", "design.json"), "r", encoding="utf-8") as f:
    aircraft_data = json.load(f)

with open(os.path.join(os.path.dirname(__file__), "Configurations", "c206.json"), "r", encoding="utf-8") as f:
    c206_data = json.load(f)


def save_ac_data_to_json(ac_data=aircraft_data, filename="design.json"):
    with open(os.path.join(os.path.dirname(__file__), "Configurations", filename), "w", encoding="utf-8") as f:
        json.dump(ac_data, f, indent=4)
