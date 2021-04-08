import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from natsort import natsorted


'''Версия программы для подбора коэффициентов'''


def bellFixed(x, period):
    return np.exp(-(x / period) ** 2)


def bell(x, period, centre, maximum):
    return maximum*(np.exp(-((x - centre) / period) ** 2))


def flat(x, a, b, c):
    return x[1]*a + x[0]*b + c


def untilt(massive: np.ndarray, ground: float, error: float) -> np.ndarray:
    mask = np.abs(massive - ground) < error
    XY = np.array(np.where(mask))
    Z = massive[mask]
    popt, pcov = curve_fit(flat, XY, Z, p0=[0.01, 0.01, 1])
    return massive - flat(np.array(np.meshgrid(np.arange(massive.shape[1]), np.arange(massive.shape[0]))[::-1]), *popt)


pathfile = "segmentation.path"
file = open(pathfile, encoding='utf-8')
pathes = [i.split('\t') for i in file.read().strip().split('\n')]

minSpaceWithBorders = 1000
minSpaceCentre = 1000
createPictures = True


for dir, outdir in pathes:
    print("\nNew path:", dir, "\nWorking with file ")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    files = os.listdir(dir)
    maxTemp = 5         # Количество файлов из одной папки
    temp = 0
    coefsNumber = 5             # Количество тестируемых коэффициентов
    plt.figure(figsize=(9*(coefsNumber+1), 3 * 3 * maxTemp))
    plt.suptitle(dir)

    for f in natsorted(files):
        if f[-4:] == '.txt':
            name = f[:-4]
            print(f, end=' ')
            data = -np.loadtxt(os.path.join(dir, f))
            data[np.isnan(data)] = 0

            "Выравнивание фона"
            y, x = np.histogram(data.ravel(), bins=500)
            x = x[:-1]

            y = y / y.max()
            centerVal = x[y.argmax()]
            x = x - centerVal
            popt, pcov = curve_fit(bellFixed, x, y)
            counter = 0
            wid = abs(popt[0])

            minWidth = wid
            minCent = centerVal
            bestUntilt = data
            while counter == 0 or (counter < 5):
                data = untilt(data, centerVal, abs(popt[0]))

                y, x = np.histogram(data.ravel(), bins=500)
                x = x[:-1]
                y = y / y.max()
                centerVal = x[y.argmax()]
                x = x - centerVal
                popt, pcov = curve_fit(bellFixed, x, y)
                wid = abs(popt[0])
                if wid < minWidth:
                    minWidth = wid
                    minCent = centerVal
                    bestUntilt = data
                counter += 1

            data = bestUntilt
            data -= minCent
            bordDiff = max(-3 * abs(minWidth), -0.4)

            plt.subplot(maxTemp, 6, 1+(temp)*6)
            plt.pcolormesh(data)
            plt.title(name)
            bords = np.histogram_bin_edges([1.1, 3], coefsNumber - 1)                   # выбор диапазона в формате [начало, конец, количество - 1]
            for i in range(coefsNumber):
                plt.subplot(maxTemp, 6, temp * 6 + 2 + i)
                bord = (1.1 + 0.3 * i)
                bord = bords[i]
                plt.pcolormesh(data < -abs(minWidth) * bord)
                plt.title(f"coef = {bord} ; value = {round(-abs(minWidth) * bord, ndigits=2)}")
            temp += 1
            if temp == maxTemp:
                plt.show()
                break
            else:
                continue


