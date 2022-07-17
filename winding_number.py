import numpy as np
from numpy.linalg import *

def qx(kx, ky):
    v = -1.6
    w = -1.8
    t = -2.2
    m = 1j * np.exp(1j*kx) * (np.exp(1j*ky) * (1 + np.exp(1j*kx)) * v**4 - t**2 * w**2)\
    -1j * (np.exp(1j*ky) * v**4 + np.exp(1j*(kx + ky)) * v**4 - np.exp(1j*kx) * t**2 * w**2)
    n = np.exp(1j*ky) * (1 + np.exp(1j*kx))**2 * v**4 - np.exp(1j*kx) * t**2 * w**2
    return m / n

def qy(kx, ky):
    v = -1.6
    w = -1.8
    t = -2.2
    n = 1j + 4 * v**4 * np.cos(kx/2)**2 * (-1j * np.cos(ky) + np.sin(ky)) / (t**2 * w**2)
    return 1 / n

def main():
    delta_abs = 0.0001
    # kpaths = np.array(
    #     [
    #         [-np.pi, -np.pi],
    #         [0, -np.pi],
    #         [0, np.pi],
    #         [-np.pi, np.pi],
    #         [-np.pi, -np.pi]
    #     ]
    # )

    kpaths = np.array(
        [
            [-np.pi, -np.pi],
            [-np.pi, np.pi],
            [np.pi, np.pi],
            [np.pi, -np.pi],
            [-np.pi, -np.pi]
        ]
    )
    # kpaths = np.array(
    #     [
    #         [-3, -3],
    #         [-3, 3],
    #         [3, 3],
    #         [3, -3],
    #         [-3, -3]
    #     ]
    # )

    # kpaths = np.array(
    #     [
    #         [-np.pi, -0.1],
    #         [0, -0.1],
    #         [0, np.pi],
    #         [-np.pi, np.pi],
    #         [-np.pi, -0.1]
    #     ]
    # )

    W = 0
    for idx in range(len(kpaths) - 1):
        if kpaths[idx][0] == kpaths[idx + 1][0]:
            if kpaths[idx][1] < kpaths[idx + 1][1]:
                delta = delta_abs
            if kpaths[idx][1] > kpaths[idx + 1][1]:
                delta = -delta_abs
            b = np.arange(kpaths[idx][1], kpaths[idx + 1][1], delta)
            a = kpaths[idx][0]
            for i in range(len(b)):
                # print(a, b[i])
                W += qy(a, b[i]) * delta

        if kpaths[idx][1] == kpaths[idx + 1][1]:
            if kpaths[idx][0] < kpaths[idx + 1][0]:
                delta = delta_abs
            if kpaths[idx][0] > kpaths[idx + 1][0]:
                delta = -delta_abs
            a = np.arange(kpaths[idx][0], kpaths[idx + 1][0], delta)
            b = kpaths[idx][1]
            for i in range(len(a)):
                # print(a[i], b)
                W += qx(a[i], b) * delta

    print('Winding number: {}'.format(W/(2*np.pi*1j)))


if __name__ == '__main__':
    main()
