import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit


def bell(x, period, centre, maximum):
    return maximum*(np.exp(-((x - centre) / period) ** 2))


def near(i, j, ilim, jlim):
    res = []
    for par in [(i+1, j), (i, j+1), (i-1, j), (i, j-1)]:
        if 0 <= par[0] < ilim and 0 <= par[1] < jlim:
            res.append(par)
    return res


def otherNum(i, j, diff, ilim, jlim):
    res = 0
    for ni, nj in near(i, j, ilim, jlim):
        if data[ni, nj] >= diff:
            res += 1
    return res


def distance(cellN, i, j):
    return ((cents[cellN-1][0] - j)**2 + (cents[cellN-1][1] - i)**2)/cells[cellN-1]


def flat(x, a, b, c):
    return x[1]*a + x[0]*b + c


def untilt(massive, ground, error):
    mask = np.abs(massive - ground) < error
    XY = np.array(np.where(mask))
    Z = data[mask]
    popt, pcov = curve_fit(flat, XY, Z, p0=[0.01, 0.01, 1])
    return massive - flat(np.array([np.arange(data.shape[0]).reshape(-1, 1), np.arange(data.shape[1]).reshape(1, -1)]), *popt)


pathes = [("cellsEasy", "result/Easy"), ("cellsMedium", "result/Medium")]       # Надо вписывать кортежи (папка с исходными данными, папка для вывода данных)
minSpaceWithBorders = 1000
minSpaceCentre = 1000
createPictures = True


for dir, outdir in pathes:
    print("New path:", dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    files = os.listdir(dir)

    for f in files:
        if f[-4:] == '.txt':
            name = f[:-4]
            print("Working with file", f)
            data = np.loadtxt(os.path.join(dir, f))
            y, x = np.histogram(data.ravel(), bins=500)
            x = x[:-1]
            popt, pcov = curve_fit(bell, x, y, p0=[1, 0.1, 1000])

            counter = 0
            while counter == 0 or (counter < 6):
                # print("Step", counter, 'popt', popt)
                data = untilt(data, popt[0], abs(popt[1])*2 + 0.3)
                y, x = np.histogram(data.ravel(), bins=500)
                x = x[:-1]
                popt, pcov = curve_fit(bell, x, y, p0=[1, 0.1, 1000])
                counter += 1

            # print("Step", counter, 'popt', popt)


            bordDiff = max(-3 * abs(popt[0]) + popt[1] - 0.3, -1)
            c = 1
            used = np.zeros_like(data, dtype="int64")

            cells = []  # площадь
            cents = []  # j, i координаты центров
            # borders = []  # массив с очертаниями центров клетки
            "Нахождение больших кучностей клеток"
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if data[i, j] < bordDiff and used[i, j] == 0:
                        used[i, j] = c
                        queue = [(i, j)]
                        cells.append(0)
                        x = 0
                        while x < len(queue):

                            cells[-1] += 1
                            for nei in near(*queue[x], data.shape[0], data.shape[1]):
                                if data[nei[0], nei[1]] < bordDiff and used[nei[0], nei[1]] != c:
                                    if used[nei[0], nei[1]] != 0:
                                        print("Я сломался", nei[0], nei[1])
                                        # exit(1)
                                    used[nei[0], nei[1]] = c
                                    queue.append(nei)
                            x += 1

                        if cells[-1] < minSpaceCentre*3:
                            used[used == c] = -1
                            cells.pop()
                            continue
                        cells.pop()

                        y, x = np.histogram(data[used == c], bins=500)

                        x = x[:-1]
                        halfSum = y.sum()/3     # не верить названиям переменных
                        k = 0
                        sum = 0
                        while sum < halfSum:
                            sum += y[k]
                            k += 1
                        centDiff = x[k]
                        used[used == c] = 0
                        borders = []

                        "Нахождение центров в куче"
                        for I, J in queue:
                            if data[I, J] < centDiff and used[I, J] == 0:
                                borders.append([])
                                used[I, J] = c
                                subqueue = [(I, J)]
                                cells.append(0)
                                cents.append([0, 0])
                                x = 0
                                while x < len(subqueue):
                                    cents[-1][0] += subqueue[x][1]
                                    cents[-1][1] += subqueue[x][0]
                                    if otherNum(*subqueue[x], centDiff, *data.shape) > 0:  # клетка на границе
                                        borders[-1].append(subqueue[x])
                                    cells[-1] += 1
                                    for nei in near(*subqueue[x], *data.shape):
                                        if data[nei[0], nei[1]] < centDiff and used[nei[0], nei[1]] != c:
                                            if used[nei[0], nei[1]] != 0:
                                                print("Я сломался", nei[0], nei[1])
                                                # exit(1)
                                            used[nei[0], nei[1]] = c
                                            subqueue.append(nei)
                                    x += 1

                                if cells[-1] < minSpaceCentre:
                                    used[used == c] = -1
                                    borders.pop()
                                    cells.pop()
                                    cents.pop()
                                else:
                                    cents[-1][0] /= cells[-1]
                                    cents[-1][1] /= cells[-1]
                                    c += 1
                        used[used == -1] = 0
                        "Расширение границ центров из кучи"
                        for k in range(len(borders)):
                            subc = k + c - len(borders)
                            queue = borders[k]
                            x = 0
                            while x < len(queue):
                                for nei in near(*queue[x], *data.shape):
                                    if data[nei[0], nei[1]] < bordDiff and used[nei[0], nei[1]] != subc:
                                        if (used[nei[0], nei[1]] != 0 and
                                                distance(used[nei[0], nei[1]], nei[0], nei[1]) < distance(subc, nei[0], nei[1])):
                                            continue

                                        used[nei[0], nei[1]] = subc
                                        queue.append(nei)
                                x += 1
                        if len(borders) == 0:
                            for I, J in queue:
                                used[I, J] = -1

            used[used == -1] = 0

            "Возвращение отрезанных кусков"
            pieces = []     # [cell number, queue, set(neibours)]
            merged = np.zeros_like(used, dtype=bool)
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    if used[i, j] and not merged[i, j]:
                        merged[i, j] = True
                        queue = [(i, j)]
                        pieces.append([used[i, j], [], set()])       # [cell number, queue, set(neibours)]
                        x = 0
                        while x < len(queue):

                            for nei in near(*queue[x], *data.shape):
                                if used[nei[0], nei[1]] == pieces[-1][0]:
                                    if not merged[nei[0], nei[1]]:
                                        merged[nei[0], nei[1]] = True
                                        queue.append(nei)
                                elif used[nei[0], nei[1]] != 0:

                                    pieces[-1][-1].add(used[nei[0], nei[1]])
                            x += 1
                        pieces[-1][1] = queue
            pieces.sort(key=lambda x: len(x[1]), reverse=True)
            pieces.sort(key=lambda x: x[0])
            for k in range(1, len(pieces)):
                if pieces[k][0] == pieces[k - 1][0]:
                    new = 0 if len(pieces[k][-1]) == 0 else pieces[k][-1].pop()
                    for i, j in pieces[k][1]:
                        used[i, j] = new

            # np.savetxt(os.path.join(outdir, f'{name}coloring.txt'), used)
            if createPictures:
                plt.subplots()
                plt.pcolormesh(used)
                plt.title(name)
                plt.savefig(os.path.join(outdir, f'{name}coloring3.png'), dpi=300)
                # plt.show()
                plt.close()

            "Выделение рамок для клеток и их печать"
            outPlace = os.path.join(outdir, name)
            if not os.path.exists(outPlace):
                os.makedirs(outPlace)

            for cell in range(c-1):
                maski, maskj = np.where(used == cell+1)
                imin, imax = maski.min(), maski.max()
                jmin, jmax = maskj.min(), maskj.max()
                if imin == 0 or jmin == 0 or imax == used.shape[0] - 1 or jmax == used.shape[1] - 1:
                    # клетка около границы
                    continue
                map = used[imin-1: imax+2, jmin-1:jmax+2] == cell + 1
                queue = [(0, 0)]

                x = 0
                while x < len(queue):
                    map[0, 0] = True
                    for nei in near(*queue[x], *map.shape):
                        if map[nei[0], nei[1]] == False:
                            map[nei[0], nei[1]] = True
                            queue.append(nei)
                    x += 1
                map = np.bitwise_or((np.bitwise_not(map)), used[imin-1: imax+2, jmin-1:jmax+2] == cell + 1)      # включение внутренних пустот
                output = data[imin-1: imax+2, jmin-1:jmax+2]*map
                np.savetxt(os.path.join(outPlace, str(cell) + '.txt'), output)
                if createPictures:
                    plt.pcolormesh(output)
                    plt.savefig(os.path.join(outPlace, f'im{cell}.png'))
                    plt.close()


