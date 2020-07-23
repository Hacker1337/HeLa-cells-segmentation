import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd


# pathes = [("result/Easy", "analys/Easy"), (
#     "result/Medium", "analys/Medium")]  # Надо вписывать кортежи (папка с исходными данными, папка для вывода данных)
pathes = [("result/Easy", "analys/Easy"), ('result/Many', "analys/Many"), ( "result/Medium", "analys/Medium")]
# createPictures = True

alpha = 0.19e-3  # in ci so m^3/kg
pixelSize = 0.165e-6  # in meters
backN = 1.336
cellN = 1.35

for dir, outdir in pathes:
    print("New path:", dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    files = os.listdir(dir)
    M, avSurfDens, vol, SA, SAV, SDM, PAR, ψ, σ, Kurtosis, Skewness, ε, name = ([] for i in range(13))
    for f in files:
        if f.count('.') != 0:
            continue
        cellsdir = os.path.join(dir, f)
        print("\n\tFrame", f, "\n\t\tWorking with", end=' ')
        for cell in os.listdir(cellsdir):
            if cell[-4:] != ".txt":
                continue
            print(cell, end=', ')
            data = np.loadtxt(os.path.join(cellsdir, cell))
            name.append(os.path.join(cellsdir, cell))
            stat = data[data != 0]
            S = stat.size * (pixelSize ** 2)

            M.append(abs(data.sum()) * (pixelSize ** 2) / alpha)

            avSurfDens.append(M[-1] / S)
            # np.savetxt(os.path.join(outdir, cell[:-4] + "surfdens.txt"), data / alpha)  # surfdens
            vol.append(abs(data.sum()) * pixelSize ** 2)
            SA.append(0)
            for i in np.arange(data.shape[0] - 1):
                for j in np.arange(data.shape[1] - 1):
                    if data[i, j] != 0:
                        SA[-1] += pixelSize ** 2
                    if data[i + 1, j] != 0 or data[i, j + 1] != 0 or data[i, j] != 0:
                        SA[-1] += (pixelSize ** 2) * np.sqrt(1 + (data[i + 1, j] - data[i, j]) ** 2 +
                                                           (data[i, j + 1] - data[i, j]) ** 2)
            SAV.append(SA[-1] / vol[-1])
            SDM.append(SA[-1] / M[-1])
            PAR.append(S / vol[-1])
            ψ.append(np.power(np.pi, 1 / 3) * np.power(6 * vol[-1], 2 / 3) / SA[-1])

            µ = stat.mean()
            σ.append(((stat - µ) ** 2).sum() / (stat.size - 1))
            Kurtosis.append(((stat - µ) ** 4).sum() / σ[-1] ** 4)
            Skewness.append(((stat - µ) ** 3).sum() / σ[-1] ** 3)
            cent = np.array(np.where(data != 0)).sum(axis=0) / S
            rmin = data.size ** 2
            rmax = 0
            for i in np.arange(data.shape[0]):
                for j in np.arange(data.shape[1]):
                    if data[i, j] == 0:
                        rmin = min(rmin, (i - cent[0]) ** 2 + (j - cent[1]) ** 2)
                    else:
                        rmax = max(rmax, (i - cent[0]) ** 2 + (j - cent[1]) ** 2)
            rmax = np.sqrt(rmax)
            rmin = np.sqrt(rmin)
            ε.append((rmax - rmin) / (rmax + rmin))

    output = pd.DataFrame(
        {"M": M, "rho": avSurfDens, "V": vol, "SA": SA, "SAV": SAV, "SDM": SDM, "PAR": PAR, "ψ": ψ, "σ": σ,
         "Kurtosis": Kurtosis, "Skewness": Skewness, "ε": ε, "name": name})
    output.to_csv(os.path.join(outdir, "table.csv"))
    output.hist(bins=50, figsize=(10, 10))
    plt.savefig(os.path.join(outdir, 'hists.png'), dpi=400)
    plt.close()
    pd.plotting.scatter_matrix(output, figsize=(50, 50))
    plt.savefig(os.path.join(dir, 'scatter.png'), dpi=300)
    plt.close()

    # for name, value in {"M": M, "\rho": avSurfDens, "V": vol, "SA": SA, "SAV": SAV, "SDM": SDM, "PAR": PAR, "ψ": ψ,
    #                     "σ": σ, "Kurtosis": Kurtosis, "Skewness": Skewness, "ε": ε, "name": name}:
    #     plt.hist(value, bins=50)
    #     plt.title(name)
    #     plt.savefig(os.path.join(outdir, name+'.png'), dpi=300)
