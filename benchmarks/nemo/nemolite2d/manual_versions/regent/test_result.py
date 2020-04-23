import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
import sys

print(sys.argv[1])
val = int(sys.argv[1])

base_data = pd.read_csv(f'go2d_{val:05d}_00000.dat', delim_whitespace=True,header=None )
base_data = base_data.set_index([0,1])
base_data = base_data.sort_values([0,1])
print(base_data)
new_data = pd.read_csv(f'../psykal_omp/go2d_{val:05d}_00001.dat', delim_whitespace=True,header=None)
new_data = new_data.set_index([0,1])
new_data = new_data.sort_values([0,1])
print(new_data)
new_data = new_data.apply(pd.to_numeric)
print(base_data.equals( new_data))

res =new_data[ ~new_data.isin(base_data)]
res.columns = ['A', 'B', 'C', 'D']
res2 = res[ res.A == res.A  ]
print(res2)
res3 = res[ res.B == res.B ]
print(res3)
res4 = res[ res.C == res.C ]
print(res4)
res = res[ res.D == res.D ]
print(res)
