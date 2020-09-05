import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.optimize import curve_fit


def bell(x, period, centre, maximum):
    return maximum*(np.exp(-((x - centre) / period) ** 2))


# Надо вписывать кортежи (папка с исходными данными, папка для вывода данных)
# pathes = ["analys/Many", "analys/Easy", "analys/Medium"]
pathes = ["analys/apoptose"]

for dir in pathes:
    table = pd.read_csv(os.path.join(dir, 'table1.csv'),  index_col=0)
    for col in table:
        if table[col].dtype != 'O':
            y, x = np.histogram(table[col], bins=500)
            x = x[:-1]
            popt, conv = curve_fit(bell, x, y, p0=[table[col].mean(), table[col].mean() if col != 'Kurtosis' else 0, y.max()])
            table = table[abs(table[col]-popt[1]) < abs(popt[0])*4]
    table.to_csv(os.path.join(dir, "tableFiltered.csv"))

    table.hist(figsize=(20, 20), bins=50)
    plt.savefig(os.path.join(dir, 'histsFiltered200.png'), dpi=300)
    plt.show()
    plt.close()
    # pd.plotting.scatter_matrix(table, figsize=(50, 50))
    # plt.savefig(os.path.join(dir, 'scatter.png'), dpi=300)
    # plt.close()

    # pd.plotting.parallel_coordinates(table, figsize=(10, 10))
    # plt.savefig(os.pathes.join(dir, 'hists.png'), dpi=300)
    # plt.close()


