import numpy as np
from argparse import Namespace
import datetime
import tensorflow as tf

# Suppress TF warnings
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

from bayes_opt import BayesianOptimization, UtilityFunction
from scipy import optimize

from pyemittance.emit_eval_example import eval_emit_surrogate

# importing lcls surrogate model
from lcls_functions import Lcls


class Opt():
    def __init__(self, init_scan=[-6, -4, -2, 0]):
        self.energy = 0.135
        self.varscan = init_scan
        self.num_points_adapt = 7
        self.pbounds = ((0.46, 0.485), (-0.01, 0.01), (-0.01, 0.01))
        self.plot = False
        self.bsfn = None
        self.uncertainty_lim = 0.25
        self.timestamp = None
        self.noise = False

    def init_run(self):
        lcls_params = Namespace(
            config_bounds=[(self.pbounds[0][0], self.pbounds[0][1]),
                           (self.pbounds[1][0], self.pbounds[1][1]),
                           (self.pbounds[2][0], self.pbounds[2][1])],
            quad_bounds=(-8.0, 1.0),
            beamsizes_bounds=[(0.0, 5e-4), (0.0, 5e-4)],
        )
        lcls = Lcls(params=lcls_params)
        self.bsfn = lcls.beamsizes_list_fn

    def get_beamsizes_model(self, config, val):
        beamsizes_list = self.bsfn(config, [val])[0]
        xrms = beamsizes_list[0]
        yrms = beamsizes_list[1]
        xrms_err = xrms * 0.03
        yrms_err = yrms * 0.03

        return xrms, yrms, xrms_err, yrms_err

    def evaluate(self, varx, vary, varz):
        # Set lcls and functions
        self.init_run()

        # fixed varscan
        quad_init = self.varscan
        config = [varx, vary, varz]

        out_dict, total_num_points = eval_emit_surrogate(self.get_beamsizes_model,
                                                         config,
                                                         quad_init=list(quad_init),
                                                         adapt_ranges=True,
                                                         num_points=self.num_points_adapt,
                                                         check_sym=True,
                                                         infl_check=True,
                                                         add_pnts=True,
                                                         show_plots=self.plot,
                                                         add_noise=self.noise)

        emit = out_dict['nemit'] / 1e-6
        emit_err = out_dict['nemit_err'] / 1e-6

        if np.isnan(emit):
            print("NaN emit")
            return np.nan, np.nan

        # if emit_err / emit < self.uncertainty_lim:
        #     # save total number of points added
        #     timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S")
        #     f= open(f"bo_points_added_iter_10_15_final_noiseless.txt", "a+")
        #     f.write(f'{timestamp},{total_num_points},{emit},{emit_err}\n')
        #     f.close()

        return -emit, -emit_err  # in um

    def run_bo_opt_w_reject(self, rnd_state=11, init_pnts=3, n_iter=120):
        # Set domain
        bounds = {'varx': self.pbounds[0], 'vary': self.pbounds[1], 'varz': self.pbounds[2]}

        # Run BO
        optimizer = BayesianOptimization(
            f=None,
            pbounds=bounds,
            random_state=rnd_state,
            verbose=2
        )

        # utility = UtilityFunction(kind="ucb", kappa=0.1, xi=0.0)
        utility = UtilityFunction(kind="ucb", kappa=2.5, xi=0.0)

        target_list = []

        # init random points
        x = []
        emit_list = []
        emit_err_list = []

        emit_res = (np.nan, np.nan)
        while len(emit_list) < init_pnts:
            x_i = [np.random.uniform(self.pbounds[0][0], self.pbounds[0][1]),
                   np.random.uniform(self.pbounds[1][0], self.pbounds[1][1]),
                   np.random.uniform(self.pbounds[2][0], self.pbounds[2][1])]
            emit_res = self.evaluate(x_i[0], x_i[1], x_i[2])

            if not np.isnan(emit_res[0]) and not np.isnan(emit_res[1]) and abs(emit_res[0]) > 1.0:
                x.append(x_i)
                emit_list.append(emit_res[0])
                emit_err_list.append(emit_res[1])

        # get init points
        for i in range(len(x)):
            # target, error = np.nan, np.nan
            #             while np.isnan(target) or np.isnan(error) or error/target > self.uncertainty_lim:
            next_point = {'varx': x[i][0],
                          'vary': x[i][1],
                          'varz': x[i][2]
                          }
            #                 # evaluate next point
            target = emit_list[i]

            optimizer.register(params=next_point, target=target)
            if target_list and target > np.max(target_list):
                color = '\033[95m', '\033[0m'
            else:
                color = '\u001b[30m', '\033[0m'

            print(
                f"{color[0]}iter {i} | target {-1 * target:.3f} | config {next_point['varx']:.6f}"
                f" {next_point['vary']:.6f} {next_point['varz']:.6f}{color[1]}")
            target_list.append(target)

        # BO iters
        for i in range(n_iter):
            target, error = np.nan, np.nan
            while np.isnan(target) or np.isnan(error) or error / target > self.uncertainty_lim:
                next_point = optimizer.suggest(utility)
                target, error = self.evaluate(**next_point)

            optimizer.register(params=next_point, target=target)
            if target_list and target > np.max(target_list):
                color = '\033[95m', '\033[0m'
            else:
                color = '\u001b[30m', '\033[0m'

            print(
                f"{color[0]}iter {i} | target {-1 * target:.3f} | config {next_point['varx']:.6f} "
                f"{next_point['vary']:.6f} {next_point['varz']:.6f}{color[1]}")
            emit_list.append(target)
            emit_err_list.append(error)
            target_list.append(target)

        # timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S")
        # emit_total = emit_list + emit_err_list
        # np.save(f'results_bo_10_15_final_noiseless/bo_opt_res_emit_list_{rnd_state}_{init_pnts}_{n_iter}_{timestamp}.npy',
        #         emit_total, allow_pickle=True)
        # np.save(f'results_bo_10_15_final_noiseless/optimizer_res/bo_opt_res_{rnd_state}_{init_pnts}_{n_iter}_{timestamp}.npy',
        #         optimizer.res, allow_pickle=True)

        return optimizer

    def eval_simplex(self, x):
        emit, err = self.evaluate(x[0], x[1], x[2])
        if err / emit > self.uncertainty_lim:
            return np.nan
        return -1 * self.evaluate(x[0], x[1], x[2])[0]

    def run_simplex_opt(self, max_iter):
        initial_guess = np.array(
            [np.random.uniform(self.pbounds[0][0], self.pbounds[0][1]),
             np.random.uniform(self.pbounds[1][0], self.pbounds[1][1]),
             np.random.uniform(self.pbounds[2][0], self.pbounds[2][1])
             ])

        # initial_guess1 = self.pbounds[0][0]+ np.random.rand(1) * (self.pbounds[0][1] - self.pbounds[0][0])
        # initial_guess2 = self.pbounds[1][0]+ np.random.rand(1) * (self.pbounds[1][1] - self.pbounds[1][0])
        # initial_guess3 = self.pbounds[2][0]+ np.random.rand(1) * (self.pbounds[2][1] - self.pbounds[2][0])

        # initial_guess = np.array([initial_guess1, initial_guess2, initial_guess3])

        min = optimize.minimize(self.eval_simplex, initial_guess,
                                method='Nelder-Mead', options={'maxiter': max_iter,
                                                               'return_all': True,
                                                               'adaptive': True
                                                               },
                                )
        # timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S")
        # np.save(f'simplex_{timestamp}.npy', min["allvecs"], allow_pickle=True)

        return min

    def run_bo_opt(self, rnd_state=11, init_pnts=3, n_iter=200):
        # Set domain
        bounds = {'varx': self.pbounds[0], 'vary': self.pbounds[1], 'varz': self.pbounds[2]}

        # Run BO
        optimizer = BayesianOptimization(
            f=self.evaluate,
            pbounds=bounds,
            random_state=rnd_state,
        )

        #        optimizer.maximize(init_points=init_pnts, n_iter=n_iter)
        optimizer.maximize(init_points=init_pnts,
                           n_iter=n_iter,
                           kappa=0.01
                           # kappa_decay = 0.8,
                           # kappa_decay_delay = 25
                           )

        return optimizer