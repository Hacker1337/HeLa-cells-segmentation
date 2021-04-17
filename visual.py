import numpy as np
import matplotlib.pyplot as plt
import os


pathfile = "visual.path"
file = open(pathfile, encoding='utf-8')
pathes = file.read().strip().split('\n')

for dir in pathes:
    files = os.listdir(dir)
    print("\nNew path:", dir)

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
