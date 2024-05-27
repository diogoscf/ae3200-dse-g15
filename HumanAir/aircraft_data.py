import json
import os
from types import SimpleNamespace as Namespace

with open(os.path.join(os.path.dirname(__file__), "Configurations", "design.json"), "r", encoding="utf-8") as f:
    aircraft_data = json.load(f)

with open(os.path.join(os.path.dirname(__file__), "Configurations", "c206.json"), "r", encoding="utf-8") as f:
    c206_data = json.load(f)
