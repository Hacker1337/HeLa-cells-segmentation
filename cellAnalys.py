import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd


# pathes = [("result/Easy", "analys/Easy"), (
#     "result/Medium", "analys/Medium")]  # Надо вписывать кортежи (папка с исходными данными, папка для вывода данных)
# pathes = [('result/Many', "analys/Many"), ("result/Easy", "analys/Easy"), ( "result/Medium", "analys/Medium")]
# pathes = [('result/Many', "analys/Many")]
pathes = [('resultHuge/Many', "analysHuge/Many")]

# createPictures = True

alpha = 0.19  # in mkm^3/pg == ml/g
pixelSize = 0.165  # in mkm
backN = 1.336
cellN = 1.35
waveLen = 0.6328    # in micrometers
vaim = 2425         # in mkm^3

for dir, outdir in pathes:
    print("\nNew path:", dir, end=' ')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    files = os.listdir(dir)
    M, avSurfDens, avPhaseProfile, vol, realV, S, SA, SAV, SDM, PAR, ψ, σ, Kurtosis, Skewness, ε, name = ([] for i in range(16))
    for f in files:
        if f.count('.') != 0:
            continue
        cellsdir = os.path.join(dir, f)
        print("\n\tFrame", f, "\n\t\tWorking with", end=' ')
        for cell in os.listdir(cellsdir):
            if cell[-4:] != ".txt":
                continue
            print(cell, end=', ')
            phaseProfile = np.loadtxt(os.path.join(cellsdir, cell))
            OPD = phaseProfile/(2*np.pi)*waveLen
            name.append(os.path.join(cellsdir, cell))
            stat = OPD[OPD != 0]
            S.append(stat.size * (pixelSize ** 2))

            M.append(abs(OPD.sum()) * (pixelSize ** 2) / alpha)

            avSurfDens.append(M[-1] / S[-1])
            # np.savetxt(os.path.join(outdir, cell[:-4] + "surfdens.txt"), data / alpha)  # surfdens
            vol.append(abs(OPD.sum()) * pixelSize ** 2)
            SA.append(0)
            for i in np.arange(OPD.shape[0] - 1):
                for j in np.arange(OPD.shape[1] - 1):
                    if OPD[i, j] != 0:
                        SA[-1] += pixelSize ** 2
                    if OPD[i + 1, j] != 0 or OPD[i, j + 1] != 0 or OPD[i, j] != 0:
                        SA[-1] += (pixelSize ** 2) * np.sqrt(1 + ((OPD[i + 1, j] - OPD[i, j]) / pixelSize) ** 2 +
                                                             ((OPD[i, j + 1] - OPD[i, j]) / pixelSize) ** 2)
            SAV.append(SA[-1] / vol[-1])
            SDM.append(SA[-1] / M[-1])
            PAR.append(S[-1] / vol[-1])
            ψ.append(np.power(np.pi, 1 / 3) * np.power(6 * vol[-1], 2 / 3) / SA[-1])

            µ = stat.mean()
            avPhaseProfile.append(phaseProfile[phaseProfile!=0].mean())
            σ.append(((stat - µ) ** 2).sum() / (stat.size - 1))
            Kurtosis.append(((stat - µ) ** 4).sum() / σ[-1] ** 4)
            Skewness.append(((stat - µ) ** 3).sum() / σ[-1] ** 3)
            cent = np.array(np.where(OPD != 0)).sum(axis=1) / stat.size
            rmin = OPD.size ** 2
            rmax = 0
            for i in np.arange(OPD.shape[0]):
                for j in np.arange(OPD.shape[1]):
                    if OPD[i, j] == 0:
                        rmin = min(rmin, (i - cent[0]) ** 2 + (j - cent[1]) ** 2)
                    else:
                        rmax = max(rmax, (i - cent[0]) ** 2 + (j - cent[1]) ** 2)
            rmax = np.sqrt(rmax)
            rmin = np.sqrt(rmin)
            ε.append((rmax - rmin) / (rmax + rmin))

    deltaN = sum(vol)/len(vol)/vaim
    print("\ncells refractive index =", backN + deltaN)
    realV = [x/deltaN for x in vol]
    output = pd.DataFrame(
        {"M, pg": M, "av surface density, pg/mkm^2": avSurfDens, "av phase profile, rad": avPhaseProfile,
         "V_ϕ, mkm^3": vol,"realVolume, mkm^3": realV, "S_p, mkm^2": S, "SA, mkm^2": SA, "SAV, mkm^-1": SAV, "SDM, mkm^2/pg": SDM, "PAR, mkm^-1": PAR, "ψ": ψ,
         "σ, mkm^2": σ, "Kurtosis": Kurtosis, "Skewness": Skewness, "ε": ε, "name": name})
    output.to_csv(os.path.join(outdir, "table1.csv"))

    for col in output:
        if output[col].dtype != 'O':
            av = abs(output[col].mean())
            output = output[abs(output[col]) < 5*av]

    output.hist(bins=50, figsize=(15, 15))
    plt.savefig(os.path.join(outdir, 'hists.png'), dpi=400)
    plt.close()
    pd.plotting.scatter_matrix(output, figsize=(50, 50))
    plt.savefig(os.path.join(dir, 'scatter.png'), dpi=300)
    plt.close()


