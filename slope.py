import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit


def flat(x, a, b, c):
    return x[1]*a + x[0]*b + c


dirs = [('cellsEasy', 'slope/Easy'), ('cellsMedium', 'slope/Medium')]

for dir, outdir in dirs:
    print("New path", dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(os.path.join(outdir, 'img')):
        os.makedirs(os.path.join(outdir, 'img'))
    files = os.listdir(dir)

    for f in files:
        if f[-4:] == '.txt' :   # and f[:2] == "19"
            print("Working with file", f)
            data = np.loadtxt(os.path.join(dir, f))
            aver = data.mean()
            backAv = data[data > aver].mean()
            # print("Mean equals to", aver, "backAv to", backAv)
            # vals, bins, _ = plt.hist(data.ravel(), bins=500)
            # plt.show()
            # plt.close()

            flatDiff = 1.2
            used = np.zeros_like(data)

            mask = np.abs(data - backAv) < flatDiff
            XY = np.array(np.where(mask))
            Z = data[mask]
            popt, pcov = curve_fit(flat, XY, Z, p0=[0.01, 0.01, 1])     # , p0=[4.77129651e-01, 9.21051620e+08]

            # plt.pcolormesh(data*mask)
            # plt.colorbar()
            # plt.show()
            # plt.close()

            untilted = data - flat(np.array([np.arange(data.shape[0]).reshape(-1, 1), np.arange(data.shape[1]).reshape(1, -1)]), *popt)


            plt.pcolormesh(untilted)
            plt.colorbar()
            plt.savefig(f'{os.path.join(outdir, "img")}/untilted{f[:-4]}.png', dpi=400)
            # plt.show()
            plt.close()
            #
            np.savetxt(f'{outdir}/untilted{f[:-4]}.txt', untilted)

