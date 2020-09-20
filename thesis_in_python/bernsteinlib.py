import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb


def bernsteinAntiDerivMat(n, t):
    return np.tril(np.ones((n + 2, n + 1)) * t / (n + 1), -1)


def bernsteinAntiDeriv(p, t, p0):
    return p0 + bernsteinAntiDerivMat(p.shape[0] - 1, t) @ p


def bernsteinBasis(k, n, tau):
    return comb(n, k) * (1 - tau) ** (n - k) * tau ** k


def bernsteinDegrElevMat(n, m):
    return np.array([[comb(i, j) * comb(m - i, n - j) if max(0, i - m + n) <= j <= min(n, i) + 1 else 0
                      for j in range(n + 1)] for i in range(m+1)]) / comb(m, n)


def bernsteinDegrElev(p, m):
    if p.shape[0] - 1 > m:
        return p
    return bernsteinDegrElevMat(p.shape[0] - 1, m) @ p


def bernsteinDerivMat(n, t):
    return n / t * (np.hstack((np.zeros((n, 1)), np.eye(n))) - np.hstack((np.eye(n), np.zeros((n, 1)))))


def bernsteinDeriv(p, t):
    return bernsteinDerivMat(p.shape[0] - 1, t) @ p


def bernsteinDerivElevMat(n, t):
    return bernsteinDerivMat(n + 1, t) @ bernsteinDegrElevMat(n, n + 1)


def bernsteinDerivElev(p, t):
    return bernsteinDerivElevMat(p.shape[0] - 1, t) @ p


def bernsteinEvalMat(n, t, times):
    return np.array([[bernsteinBasis(j, n, ti / t) for j in range(n + 1)] for ti in np.array(times).flatten()])


def bernsteinEval(p, t, times):
    return bernsteinEvalMat(p.shape[0] - 1, t, times) @ p


def bernsteinIntegr(p, t):
    return t / p.shape[0] * np.sum(p, 0)


def bernsteinMul(p1, p2):
    if p1.shape[0] < p2.shape[0]:
        p1, p2 = p2, p1
    m = p1.shape[0] - 1
    n = p2.shape[0] - 1
    return np.array([np.sum(
        [comb(i, j) * comb(m + n - i, m - j) * p1[j, :] * p2[i - j, :] for j in range(max(0, i - n), min(m, i) + 1)], 0)
                     for i in range(m + n + 1)]) / comb(m + n, n)


def bernsteinPow(p, y):
    if y == 0:
        return np.ones((1, p.shape[1]))
    temp_p = bernsteinPow(p, y // 2)
    if y % 2 == 0:
        return bernsteinMul(temp_p, temp_p)
    else:
        return bernsteinMul(p, bernsteinMul(temp_p, temp_p))


def bernsteinSum(p1, p2):
    if p1.shape[0] > p2.shape[0]:
        p2 = bernsteinDegrElev(p2, p1.shape[0] - 1)
    elif p2.shape[0] > p1.shape[0]:
        p1 = bernsteinDegrElev(p1, p2.shape[0] - 1)
    return p1 + p2


def bernsteinToMonMat(n, t):
    return np.flipud([[0 if i > k else comb(n, k) * comb(k, i) * (-1) ** (k - i) for i in range(n + 1)] for k in
                      range(n + 1)]) / np.array([t ** i for i in range(n, -1, -1)]).reshape((-1, 1))


def bernsteinToMon(p, t):
    return bernsteinToMonMat(p.shape[0] - 1, t) @ p


def bernsteinPlot(p, t, ax=None):
    n, dim = p.shape
    times = np.linspace(0, t, 100)
    vals = bernsteinEval(p, t, times)
    curveplot = None
    pointsplot = None
    if ax is None:
        if dim == 3:
            ax = plt.figure().add_subplot(111, projection='3d')
        elif dim == 1 or dim == 2:
            _, ax = plt.subplots()
    if dim != 1 and dim != 2 and dim != 3:
        raise ValueError('Unsupported dim')
    if dim == 1:
        curveplot, = ax.plot(times, vals)
        pointsplot = ax.scatter(np.linspace(0, t, n), p)
    elif dim == 2:
        curveplot, = ax.plot(vals[:, 0], vals[:, 1])
        pointsplot = ax.scatter(p[:, 0], p[:, 1])
    elif dim == 3:
        ax.plot(vals[:, 0], vals[:, 1], vals[:, 2])
    return curveplot, pointsplot
