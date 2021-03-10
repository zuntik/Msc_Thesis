from trajecoptim import run_problem, plot_xy
from scipy.integrate import solve_ivp
import bernsteinlib as bern
import numpy as np


def main():

    constants = {
        'T': 10,
        'xi': np.array([[0, 0, 0, 1, 0]]),
        'xf': np.array([[5, 5, np.pi / 2, 1, 0]]),
        'statebounds': np.array([[-10000, -20000, -30000, -40000, -50000], [10000, 20000, 30000, 40000, 50000]]),
        #  'inputbounds': None,
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

    constants = {** constants, **{
        'uselogbar': False,
        'usesigma': True,
        # functions
        'costfun_single': costfun_single,
        'dynamics': dynamics5vars,
        'recoverxy': recoverplot,
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

    def pol_v(t): return bern.eval(x[:, 3], constants['T'], t)

    def pol_w(t): return bern.eval(x[:, 4], constants['T'], t)

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


################################################################################
# run
################################################################################

if __name__ == "__main__":
    main()
