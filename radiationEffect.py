# %%

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.optimize import curve_fit

# # %%
#
# src = 'analys/анализ излучения'
# out = 'analys/tables'
#
# # %%
#
# testTable = pd.read_csv('analys/анализ излучения/1_1 А контроль АЛА 75мВт/table1.csv', index_col=0)
#
# # %%
#
# # testTable.head()
#
# # testTable.drop(['name'], axis=1, inplace=True)
# # testTable.head()
# testTable.to_csv('test.csv')
#
# # %%
#
# name = ''
# oldTable = pd.DataFrame()
# for group in os.listdir(src):
#     if group == 'graphs.opj':
#         continue
#     if group.endswith('(2)'):
#         table = pd.read_csv(os.path.join(src, group, 'table1.csv'), index_col=0)
#         result = pd.concat([oldTable, table])
#         result.drop(['name'], axis=1, inplace=True)
#         result.to_csv(os.path.join(out, name + '.csv'))
#     else:
#         name = group
#         oldTable = pd.read_csv(os.path.join(src, group, 'table1.csv'), index_col=0)

# %%
src = 'analys/tables_2'
out = 'analys/columns_2'

Table = pd.read_csv(os.path.join(src, 'control_ALA_75mWt.csv'), index_col=0)

for col in Table.columns:
    res = pd.DataFrame()
    for file in os.listdir(src):
        table = pd.read_csv(os.path.join(src, file), index_col=0)
        res = pd.concat([res, pd.DataFrame({file[:-4]: table[col].values})], axis=1)
    res.to_csv(os.path.join(out, col.replace('/', '%')+'.csv'))