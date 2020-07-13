import numpy as np
import matplotlib.pyplot as plt
import os


def near(i, j):
    res = []
    for par in [(i+1, j), (i, j+1), (i-1, j), (i, j-1)]:
        if 0 <= par[0] < field.shape[0] and 0 <= par[1] < field.shape[1]:
            res.append(par)
    return res


def otherNum(i, j):
    res = 0
    for ni, nj in near(i, j):
        if not field[ni, nj]:
            res += 1
    return res


def distance(cellN, i, j):
    return ((cents[cellN-1][0] - j)**2 + (cents[cellN-1][1] - i)**2)/cells[cellN-1][1]


pathes = [("slope", "primitivOutput/Easy"), ("slopeMedium", "primitivOutput/Medium")]
minSpace = 500

for dir, outdir in pathes:
    print("New path:", dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    files = os.listdir(dir)

    for f in files:
        if f[-7:] == '1.2.txt' and f[:8] == "untilted":
            name = f[8:-4]
            print("Working with file", f)
            data = np.loadtxt(os.path.join(dir, f))
            centDiff = 1.65             # TODO вставить часть с частотным анализом
            bordDiff = 0.8
            field = data < -centDiff        # TODO Надо удалить из тел функций

            c = 1
            used = np.zeros_like(data, dtype="int64")

            cells = []  # периметр, площадь     # TODO убрать периметр
            cents = []  # j, i координаты центров
            borders = []  # массив с очертаниями центров клетки
            "Нахождение больших центров"
            for i in range(field.shape[0]):
                for j in range(field.shape[1]):
                    if data[i, j] < -centDiff and used[i, j] == 0:
                        borders.append([])
                        used[i, j] = c
                        queue = [(i, j)]
                        cells.append([0, 0])
                        cents.append([0, 0])
                        x = 0
                        while x < len(queue):
                            cents[-1][0] += queue[x][1]
                            cents[-1][1] += queue[x][0]
                            if otherNum(*queue[x]) > 0:  # клетка на границе
                                borders[-1].append(queue[x])
                            cells[-1][1] += 1
                            for nei in near(*queue[x]):
                                if data[nei[0], nei[1]] < -centDiff and used[nei[0], nei[1]] != c:
                                    if used[nei[0], nei[1]] != 0:
                                        print("Я сломался", i, j)
                                        # exit(1)
                                    used[nei[0], nei[1]] = c
                                    queue.append(nei)
                            x += 1

                        if cells[-1][1] < minSpace:
                            used[used == c] = 0
                            borders.pop()
                            cells.pop()
                            cents.pop()
                        else:
                            cents[-1][0] /= cells[-1][1]
                            cents[-1][1] /= cells[-1][1]
                            c += 1




            # plt.figure()
            # plt.suptitle(name)
            # plt.pcolormesh(used > 0)
            # for i in range(len(cents)):
            #     plt.text(*cents[i], str(cells[i][1]), {'color':'w'})
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
            plt.savefig(os.path.join(outdir, f'{name}coloringExp.png'), dpi=300)
            plt.show()
            plt.close()

            outPlace = os.path.join(outdir, name)
            if not os.path.exists(outPlace):
                os.makedirs(outPlace)
            for cell in range(c-1):
                if cells[cell][1] >= minSpace:
                    output = data*(used == cell+1)
                    np.savetxt(os.path.join(outPlace, f'cell{cell}.txt'), output)
                    fig, ax = plt.subplots()
                    ax.pcolormesh(output)
                    ax.set_title(f"{cells[cell]} {cells[cell][0]**2/cells[cell][1]}")
                    plt.savefig(os.path.join(outPlace, f'im{cell}.png'), dpi=100)
                    plt.close(fig)


