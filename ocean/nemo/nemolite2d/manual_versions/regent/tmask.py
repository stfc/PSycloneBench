import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

base_data = pd.read_csv('/home/aidan/PSycloneBench/ocean/nemo/nemolite2d/manual_versions/psykal_omp/tmask_00001.dat', delim_whitespace=True,header=None )
base_data = base_data.set_index([0,1])
base_data = base_data.sort_values([0,1])
print(base_data)
new_data = pd.read_csv('tmask_0.dat', delim_whitespace=True,header=None)
new_data = new_data.set_index([0,1])
new_data = new_data.sort_values([0,1])
print(new_data)
print(base_data.equals( new_data))
#base_data['Match'] = np.where(base_data[2] == new_data[2], True, False)
res =new_data[ ~new_data.isin(base_data)]
res.columns = ['A']
res2 = res[ res.A == res.A  ]
print(res2)
