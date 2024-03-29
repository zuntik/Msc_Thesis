import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb


def antiderivmat(n, t):
    """Returns an integral/antiderivative matrix"""
    return np.tril(np.ones((n + 2, n + 1)) * t / (n + 1), -1)


def antideriv(p, t, p0):
    return p0 + antiderivmat(p.shape[0] - 1, t) @ p


def basis(k, n, tau):
    return comb(n, k) * (1 - tau) ** (n - k) * tau ** k


def degrelevmat(n, m):
    return np.array([[comb(i, j) * comb(m - i, n - j) if max(0, i - m + n) <= j <= min(n, i) + 1 else 0
                      for j in range(n + 1)] for i in range(m+1)]) / comb(m, n)


def degrelev(p, m):
    if p.shape[0] - 1 > m:
        return p
    return degrelevmat(p.shape[0] - 1, m) @ p


def derivmat(n, t):
    return n / t * (np.hstack((np.zeros((n, 1)), np.eye(n))) - np.hstack((np.eye(n), np.zeros((n, 1)))))


def deriv(p, t):
    return derivmat(p.shape[0] - 1, t) @ p


def derivelevmat(n, t):
    return derivmat(n + 1, t) @ degrelevmat(n, n + 1)


def derivelev(p, t):
    return derivelevmat(p.shape[0] - 1, t) @ p


def evalmat(n, t, times):
    return np.array([[basis(j, n, ti / t) for j in range(n + 1)] for ti in np.array(times).flatten()])


def eval(p, t, times):
    return evalmat(p.shape[0] - 1, t, times) @ p


def integr(p, t):
    return t / p.shape[0] * np.sum(p, 0)


def mul(p1, p2):
    if p1.shape[0] < p2.shape[0]:
        p1, p2 = p2, p1
    m = p1.shape[0] - 1
    n = p2.shape[0] - 1
    return np.array([np.sum(
        [comb(i, j) * comb(m + n - i, m - j) * p1[j, :] * p2[i - j, :] for j in range(max(0, i - n), min(m, i) + 1)], 0)
                     for i in range(m + n + 1)]) / comb(m + n, n)


def pow(p, y):
    if y == 0:
        return np.ones((1, p.shape[1]))
    temp_p = pow(p, y // 2)
    if y % 2 == 0:
        return mul(temp_p, temp_p)
    else:
        return mul(p, mul(temp_p, temp_p))


def sum(p1, p2):
    if p1.shape[0] > p2.shape[0]:
        p2 = degrelev(p2, p1.shape[0] - 1)
    elif p2.shape[0] > p1.shape[0]:
        p1 = degrelev(p1, p2.shape[0] - 1)
    return p1 + p2


def tomonmat(n, t):
    return np.flipud([[0 if i > k else comb(n, k) * comb(k, i) * (-1) ** (k - i) for i in range(n + 1)] for k in
                      range(n + 1)]) / np.array([t ** i for i in range(n, -1, -1)]).reshape((-1, 1))


def tomon(p, t):
    return tomonmat(p.shape[0] - 1, t) @ p


def plot(p, t, ax=None):
    n, dim = p.shape
    times = np.linspace(0, t, 100)
    vals = eval(p, t, times)
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
