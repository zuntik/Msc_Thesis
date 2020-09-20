from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from bernsteinlib import *
import numpy as np
import matplotlib.pyplot as plt
from time import time


def main():

    constants = {
        'T': 10,
        'xi': np.array([[0, 0, 0, 1, 0]]),
        'xf': np.array([[5, 5, np.pi / 2, 1, 0]]),
        'statebounds': np.array([[-10000, -20000, -30000, -40000, -50000], [10000, 20000, 30000, 40000, 50000]]),
        'inputbounds': np.array([[]]),
        'N': 50,
        'obstacles_circles': [[5, 0, 3]],
        # 'obstacles_circles': [],
        'obstacles_polygons': [],
        'min_dist_int_veh': 3,
        'min_dist_obs': 0,
    }

    #    constants = {
    #        'N': 30,
    #        'T': 15,
    #        'xi': np.array([[0, 5, 0, 1, 0], [5, 0, np.pi / 2, 1, 0]]),
    #        'xf': np.array([[10, 5, 0, 1, 0], [5, 10, np.pi / 2, 1, 0]]),
    #        'statebounds': None,
    #        'inputboudns': None,
    #        'obstacles_circles': [],
    #        'obstacles_polygons': [],
    #        'min_dist_int_veh': 3,
    #        'min_dist_obs': 0,
    #    }

    #    constants = {
    #        'N': 40,
    #        'T': 15,
    #        'xi': np.array([
    #            [-10, 4, 0, 1, 0],
    #            [-10, -4, 0, 1, 0],
    #            [-10, 0, 0, 1, 0],
    #        ]),
    #        'xf': np.array([
    #            [10, -1, 0, 1, 0],
    #            [10, 1, 0, 1, 0],
    #            [10, 0, 0, 1, 0],
    #        ]),
    #        'obstacles_circles': [[0, 0, 3]],
    #        'obstacles_polygons': [],
    #        'min_dist_obs': 0,
    #        'min_dist_int_veh': 0.9,
    #        'statebounds': None,  # np.array([[-20, -20, -10, -10, -10], [20, 20, 10, 10, 10]]),
    #        'inputbounds': None,
    #    }

    constants = {**constants, **{
        # common parameters
        'DiffMat': bernsteinDerivElevMat(constants['N'], constants['T']),
        'ElevMat': bernsteinDegrElevMat(constants['N'], constants['N'] * 10),
        'numvars': constants['xi'].shape[1],
        'numinputs': 0,
        'Nv': constants['xi'].shape[0],
        'uselogbar': False,
        'usesigma': True,
        # functions
        'costfun_single': costfun_single,
        'dynamics': dynamics5vars,
        'init_guess': init_guess,
        'recoverxy': recoverplot,
    }}

    # noinspection PyTypeChecker
    constants = {**constants, **{
        'obstacles':
            [Circle(c[:-1], c[-1], constants['ElevMat'], constants['min_dist_obs'])
             for c in constants['obstacles_circles']] +
            [Polygon(m) for m in constants['obstacles_polygons']]
    }}

    res, elapsedtime, singtimes = run_problem(constants)
    print('The final cost is ' + str(res.fun))
    plot_xy(res, constants)


################################################################################
# functions
################################################################################
def recoverplot(x, constants):
    def odefunc(t, val, v, w):
        _, y, psi = val
        dx = v(t) * np.cos(psi)
        dy = v(t) * np.sin(psi)
        dpsi = w(t)
        return np.array([dx, dy, dpsi])

    def pol_v(t): return bernsteinEval(x[:, 3], constants['T'], t)

    def pol_w(t): return bernsteinEval(x[:, 4], constants['T'], t)

    sol = solve_ivp(odefunc, [0, constants['T']], x[0, :3], args=(pol_v, pol_w), dense_output=True, vectorized=True)

    return np.linspace(0, constants['T'], 1000), sol.sol(np.linspace(0, constants['T'], 1000))


def dynamics5vars(x, constants):
    diffmat = constants['DiffMat']
    xp = x[:, 0]
    yp = x[:, 1]
    psi = x[:, 2]
    v = x[:, 3]
    w = x[:, 4]

    return np.vstack((
        diffmat @ xp - v * np.cos(psi),
        diffmat @ yp - v * np.sin(psi),
        diffmat @ psi - w,
    )).flatten()


def costfun_single(x, constants):
    v = x[:, 3]
    w = x[:, 4]
    a = constants['DiffMat'] @ v
    # return np.sum((constants['ElevMat']@a)**2)+2*np.sum((constants['ElevMat']@w)**2)
    return np.sum(a ** 2) + 2 * np.sum(w ** 2)


def init_guess(constants):
    return np.random.rand((constants['N'] - 1) * constants['numvars'], constants['Nv'])


################################################################################
# run
################################################################################
def logbarrierfunc(delta, z, usesigma):
    if usesigma:
        z = np.where(z >= delta, np.tanh(z), z)
    k = 2
    return np.where(z > delta, -np.log(np.abs(z)), ((k-1)/k)*(((z-k*delta)/((k-1)*delta))**k-1) - np.log(delta))


def matrify(x, constants):
    x = x.reshape((constants['Nv'], -1))
    x_mat = [
        np.concatenate((
            np.concatenate((
                constants['xi'][i, :].reshape((1, -1)),
                x[i, :(constants['N']-1)*constants['numvars']].reshape((-1, constants['numvars'])),
                constants['xf'][i, :].reshape((1, -1))
            ), axis=0),
            x[i, (constants['N']-1)*constants['numvars']+1:].reshape((constants['N']+1, constants['numinputs']))
        ), axis=1)[:, :, np.newaxis]  # .reshape((constants['N']+1, constants['numvars']+constants['numinputs'], 1))
        for i in range(constants['Nv'])]
    # return np.reshape(x, (constants['N'] - 1, constants['numvars'], constants['Nv']))
    return np.concatenate(x_mat, axis=2)


def costfun(x, constants):
    j = 0
    if constants['uselogbar']:
        c = ineqconstr(x, constants)
        ceq = eqconstr(x, constants)
        j += np.sum(logbarrierfunc(0.1, c, constants['usesigma']))
        j += 1e5 * np.sum(logbarrierfunc(0.01, -(ceq ** 2), constants['usesigma']))

    x = matrify(x, constants)
    j += np.sum([constants['costfun_single'](x[:, :, i], constants) for i in range(constants['Nv'])])

    return j


def eqconstr(x, constants):
    x = matrify(x, constants)
    return np.concatenate([constants['dynamics'](x[:, :, i], constants) for i in range(constants['Nv'])])
    #    # initial and final conditions
    #    constraints = [
    #        (constants['xi'] - x[0, :, :].T).flatten(),
    #        (constants['xf'] - x[-1, :, :].T).flatten()
    #    ]
    #    # dynamics
    #    constraints += [constants['dynamics'](x[:, :, i], constants) for i in range(constants['Nv'])]
    #    return np.concatenate(constraints)


def variablebounds(constants):
    return ([
        (constants['statebounds'][0, var], constants['statebounds'][1, var])
        for _ in range(constants['N']-1)
        for var in range(constants['numvars'])
    ]+[
        (constants['inputbounds'][1, inp], constants['inputbounds'][1, inp])
        for _ in range(constants['N']+1)
        for inp in range(constants['numinputs'])
    ])*constants['Nv'] if constants['statebounds'] is not None else None

    # return [(constants['vehiclebounds'][0, var], constants['vehiclebounds'][1, var])
    #         for _ in range(constants['N'] + 1)
    #         for var in range(constants['numvars'])
    #         for __ in range(constants['Nv'])] if constants['vehiclebounds'] is not None else None


# nonlinear inequality constraints
def ineqconstr(x, constants):
    x = matrify(x, constants)
    c = []

    # inter vehicles
    c += [veh_coll_avoid_isaac(x[:, :2, v1], x[:, :2, v2], constants)
          for v1 in range(constants['Nv']) for v2 in range(v1+1, constants['Nv'])]

    # obstacles
    c += [obs.avoid(x[:, :2, veh]) for obs in constants['obstacles'] for veh in range(constants['Nv'])]
    return np.concatenate(c) if c else np.array([])


def run_problem(constants):
    xin = constants['init_guess'](constants)
    opts = {'disp': True, 'maxiter': 1000}
    xuc = []
    singtimes = []
    for i in range(constants['Nv']):
        constants2 = constants.copy()
        constants2['xi'] = constants['xi'][i, :].reshape((1, -1))
        constants2['xf'] = constants['xf'][i, :].reshape((1, -1))
        constants2['Nv'] = 1
        constants2['obstacles'] = []
        cons2 = {'type': 'eq', 'fun': lambda x: eqconstr(x, constants2)} if not constants['uselogbar'] else []
        t = time()
        if constants['uselogbar']:
            print('Doing alg for vehicle: ' + str(i))
            xuc.append(minimize(costfun, xin[:, i], args=constants2, method='Nelder-Mead', options=opts).x)
        else:
            xuc.append(minimize(costfun, xin[:, i], args=constants2, method='SLSQP', options=opts, constraints=cons2).x)
        singtimes.append(time()-t)
        print('Elapsed time for vehicle '+str(i) + ': ' + str(singtimes[-1]) + ' s.')

    xuc = np.concatenate(xuc)

    if not constants['obstacles'] and constants['Nv'] == 1:
        return xuc

    cons = (
        {'type': 'eq', 'fun': lambda x: eqconstr(x, constants)},
        {'type': 'ineq', 'fun': lambda x: ineqconstr(x, constants)}
    ) if not constants['uselogbar'] else []
    bnds = variablebounds(constants) if constants['uselogbar'] else None
    t = time()
    if constants['uselogbar']:
        res = minimize(costfun, xuc, args=constants, options=opts)
    else:
        res = minimize(costfun, xuc, args=constants, method='SLSQP', bounds=bnds, constraints=cons, options=opts)
    elapsedtime = time()-t
    print('Elapsed time for joined up problem: ' + str(elapsedtime) + ' s.')
    return res, elapsedtime, singtimes


def plot_xy(res, constants):
    x = res.x
    _, ax = plt.subplots()
    ax.axis('equal')
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    x = matrify(x, constants)
    for i in range(constants['Nv']):
        curveplot, pointsplot = bernsteinPlot(np.fliplr(x[:, :2, i]), constants['T'], ax=ax)
        curveplot.set_label('Bernstein Polynomial for vehicle ' + str(i))
        t, xy = constants['recoverxy'](x[:, :, i], constants)
        recoveredplot, = ax.plot(xy[1, :], xy[0, :].T)
        recoveredplot.set_label('ODE solution for vehicle ' + str(i))
        ax.legend(loc='upper right', fontsize='x-small')
    for obs in constants['obstacles']:
        obs.plot(plotinverted=True, ax=ax)
    plt.show()


def veh_coll_avoid_isaac(x1, x2, constants):
    return np.min(np.sqrt(np.sum((constants['ElevMat']@(x1-x2))**2, axis=1))).flatten()-constants['min_dist_int_veh']
    # return np.min(np.sum((constants['ElevMat']@(x1-x2))**2, axis=1)).flatten()-constants['min_dist_int_veh']**2
    # return np.sqrt(np.min(np.sum((constants['ElevMat']@(x1-x2))**2, axis=1))).flatten()-constants['min_dist_int_veh']
    # return np.sqrt(np.sum((constants['ElevMat'] @ (x1-x2))**2, axis=1)).flatten() - constants['min_dist_int_veh']


class Circle:
    def __init__(self, centre, rad, elevmat, mindist):
        self.centre = centre
        self.rad = rad
        self.elevmat = elevmat
        self.mindist = mindist

    def avoid(self, poly):
        return np.sqrt(np.min(np.sum((self.elevmat@(poly-self.centre))**2, axis=1))).flatten()-self.rad-self.mindist
        # return np.sqrt(np.min(np.sum(self.elevmat@(poly-self.centre)**2, axis=1))).flatten()-self.rad-self.mindist

    def plot(self, plotinverted=False, ax=None):
        x = self.centre[1*plotinverted] + self.rad * np.cos(np.linspace(0, 2*np.pi, 100))
        y = self.centre[1*(not plotinverted)] + self.rad * np.sin(np.linspace(0, 2*np.pi, 100))
        if ax is None:
            plt.plot(x, y)
        else:
            ax.plot(x, y)


class Polygon:
    def __init__(self, matrix):
        self.matrix = matrix

    def avoid(self, poly):
        pass

    def obs(self):
        pass


if __name__ == "__main__":
    main()
