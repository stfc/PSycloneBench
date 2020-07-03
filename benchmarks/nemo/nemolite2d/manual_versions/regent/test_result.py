#Takes a single argument which is the timestep to compare, e.g. to compare results after 3000 steps:
# python3 test_result.py 3000
#Also requires the psykal_omp version to have been run with the same parameters, and is used as comparison.
#For accurate comparison it is recommended to compile the psykal_omp version with only -O2
#The output of the script will be True if the results are identical, and False otherwise

import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
import sys

#print(sys.argv[1])
val = int(sys.argv[1])

base_data = pd.read_csv(f'go2d_{val:05d}_00000.dat', delim_whitespace=True,header=None )
base_data = base_data.set_index([0,1])
base_data = base_data.sort_values([0,1])
#print(base_data)
new_data = pd.read_csv(f'../psykal_omp/go2d_{val:05d}_00001.dat', delim_whitespace=True,header=None)
new_data = new_data.set_index([0,1])
new_data = new_data.sort_values([0,1])
#print(new_data)
new_data = new_data.apply(pd.to_numeric)
print(base_data.equals( new_data))

#res =new_data[ ~new_data.isin(base_data)]
#res.columns = ['A', 'B', 'C', 'D']
#print(res)
#res2 = res[ res.A == res.A  ]
#print(res2)
#res3 = res[ res.B == res.B ]

#res2 = base_data[ ~base_data.isin(new_data)]
#res2.columns =  ['A', 'B', 'C', 'D']
#res4 = res[ res.C == res.C ]
#print(res4)
#res = res[ res.D == res.D ]
#print(res)
