from Functions import import_data
from Functions import import_data2

file_path = 'HumanAir/FuselageSizing/DATA_EXAMPLE.txt'
df = import_data(file_path)

file_path2 = 'HumanAir\FuselageSizing\WING_SPAN.txt'
dict = import_data2(file_path2)
print(dict[-10])
