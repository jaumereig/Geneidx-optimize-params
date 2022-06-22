import sys
from matplotlib.pyplot import axis
import pandas as pd
import numpy as np
from sklearn.utils import as_float_array

df = pd.read_csv("infoseq.Adalia_bipunctata-GCA_910591895.1.tbl", sep='\s+')
# print(df.iloc[:,6])
df.columns = df.columns.str.replace('[%]', '')

# df = df.assign(product=df.iloc[:,5:6].product(axis=1))
#df['product'] = df.iloc[:,5] * df.iloc[:,6]
df['product'] = df['Length'] * as_float_array(df['GC'])

sum_prod = df.iloc[:,9].sum(axis=0)
genome_length = df.iloc[:,5].sum(axis=0) 
print(sum_prod, genome_length)
Weight_GC_content = sum_prod/genome_length
print("Weight_GC_content = ", Weight_GC_content)