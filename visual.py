import numpy as np
import matplotlib.pyplot as plt
import os

dir = 'HeLa апоптоз'

files = os.listdir(dir)

for f in files:
    if f[-4:] == '.txt':
        data = np.loadtxt(os.path.join(dir, f))
        print(f, "opened successfully")
        plt.pcolormesh(data)
        plt.colorbar()
        plt.savefig(os.path.join(dir, f[:-4] + '.png'), dpi=100)
        # plt.show()
        plt.close()
        # plt.show()
