import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit


def bell(x, period, centre, maximum):
    return maximum*(np.exp(-((x - centre) / period) ** 2))


def near(i, j):
    res = []
    for par in [(i+1, j), (i, j+1), (i-1, j), (i, j-1)]:
        if 0 <= par[0] < data.shape[0] and 0 <= par[1] < data.shape[1]:
            res.append(par)
    return res


def otherNum(i, j, diff):
    res = 0
    for ni, nj in near(i, j):
        if data[ni, nj] >= - diff:
            res += 1
    return res


def distance(cellN, i, j):
    return ((cents[cellN-1][0] - j)**2 + (cents[cellN-1][1] - i)**2)/cells[cellN-1]


pathes = [("slope/Easy", "primitivOutput/Easy"), ("slope/Medium", "primitivOutput/Medium")]
minSpace = 500

for dir, outdir in pathes:
    print("New path:", dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    files = os.listdir(dir)

    for f in files:
        if f[-4:] == '.txt' and f[:8] == "untilted":
            name = f[8:-4]
            print("Working with file", f)
            data = np.loadtxt(os.path.join(dir, f))
            y, x, _ = plt.hist(data.ravel(), bins=500)
            plt.close()
            x = x[:-1]
            popt, pcov = curve_fit(bell, x, y, p0=[1, 0.1, 1000])
            bordDiff = 2 * abs(popt[0])
            centDiff = 1.65

            c = 1
            used = np.zeros_like(data, dtype="int64")

            cells = []  # периметр, площадь     # TODO убрать периметр
            cents = []  # j, i координаты центров
            borders = []  # массив с очертаниями центров клетки
            "Нахождение больших кучностей клеток"
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if data[i, j] < -centDiff and used[i, j] == 0:
                        borders.append([])
                        used[i, j] = c
                        queue = [(i, j)]
                        cells.append(0)
                        cents.append([0, 0])
                        x = 0
                        while x < len(queue):
                            cents[-1][0] += queue[x][1]
                            cents[-1][1] += queue[x][0]
                            if otherNum(*queue[x], centDiff) > 0:  # клетка на границе
                                borders[-1].append(queue[x])
                            cells[-1] += 1
                            for nei in near(*queue[x]):
                                if data[nei[0], nei[1]] < -centDiff and used[nei[0], nei[1]] != c:
                                    if used[nei[0], nei[1]] != 0:
                                        print("Я сломался", i, j)
                                        # exit(1)
                                    used[nei[0], nei[1]] = c
                                    queue.append(nei)
                            x += 1

                        if cells[-1] < minSpace:
                            used[used == c] = 0
                            borders.pop()
                            cells.pop()
                            cents.pop()
                        else:
                            cents[-1][0] /= cells[-1]
                            cents[-1][1] /= cells[-1]
                            c += 1




            # plt.figure()
            # plt.suptitle(name)
            # plt.pcolormesh(used > 0)
            # for i in range(len(cents)):
            #     plt.text(*cents[i], str(cells[i]), {'color':'w'})
            # plt.savefig(os.path.join(outdir, f'spaces{name}diff={centDiff}.png'), dpi=300)
            # plt.show()
            # plt.close()

            # plt.subplots()
            # plt.pcolormesh(used)
            # plt.title(name)
            # plt.savefig(os.path.join(outdir, f'{name}coloring.png'), dpi=300)
            # plt.show()
            # plt.close()

            "Расширение границ клеток"
            for c in range(1, len(borders)+1):
                queue = borders[c-1]
                x = 0
                while x < len(queue):
                    for nei in near(*queue[x]):
                        if data[nei[0], nei[1]] < -bordDiff and used[nei[0], nei[1]] != c:
                            if used[nei[0], nei[1]] != 0 and \
                                 distance(used[nei[0], nei[1]], nei[0], nei[1]) < distance(c, nei[0], nei[1]):
                                continue

                            used[nei[0], nei[1]] = c
                            queue.append(nei)
                    x += 1

            np.savetxt(os.path.join(outdir, f'{name}coloring.txt'), used)
            plt.subplots()
            plt.pcolormesh(used)
            plt.title(name)
            # plt.savefig(os.path.join(outdir, f'{name}coloringAnal.png'), dpi=300)
            plt.show()
            plt.close()

            # outPlace = os.path.join(outdir, name)
            # if not os.path.exists(outPlace):
            #     os.makedirs(outPlace)
            # for cell in range(c-1):
            #     if cells[cell] >= minSpace:
            #         output = data*(used == cell+1)
            #         np.savetxt(os.path.join(outPlace, f'cell{cell}.txt'), output)
            #         fig, ax = plt.subplots()
            #         ax.pcolormesh(output)
            #         ax.set_title(f"{cells[cell]}")
            #         plt.savefig(os.path.join(outPlace, f'im{cell}.png'), dpi=100)
            #         plt.close(fig)


