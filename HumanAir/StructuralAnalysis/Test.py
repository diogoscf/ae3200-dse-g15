import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from HumanAir.StructuralAnalysis.Functions import import_data, import_data2

file_path = "HumanAir/FuselageSizing/DATA_EXAMPLE.txt"
df = import_data(file_path)
print(df["CL"])

file_path2 = "HumanAir\\FuselageSizing\\WING_SPAN.txt"
dictionary = import_data2(file_path2)
print(dictionary[-10]["y_positions"])
