from pyemittance.observer import Observer
from pyemittance.data_handler import adapt_range, check_symmetry, find_inflection_pnt, add_measurements_btwn_pnts
from pyemittance.emittance_calc import EmitCalc
from pyemittance.load_json_configs import load_configs


class PyEmittance:

    def __init__(self,
                 config_name='LCLS_OTR2',
                 config_dict=None,  # supersedes json configs
                 meas_type='OTRS',
                 use_model=False,
                 online=False
                 ):

        # if config is not provided, use LCLS-OTR2 as default
        self.config_name = config_name
        self.config_dict = config_dict if config_dict else load_configs(self.config_name)

        self.meas_type = meas_type
        # if running on machine, use_model=False
        self.use_model = use_model
        # only True if setting PVs
        self.online = online
        self.verbose = True

        # injector settings (SOL, CQ, SQ) if optimizing
        self.inj_config = None
        # initial rough quad scan
        self.quad_init = [-6, -4, -2, 0]

        # pyemittance method options
        self.adapt_ranges = True
        self.num_points = 7
        self.check_sym = True
        self.infl_check = True
        self.add_pnts = True
        self.show_plots = True
        self.use_prev_meas = True
        self.quad_tol = 0.05
        self.save_runs = False
        self.calc_bmag = False

        # simulation/model options
        # beamsize function from model
        self.get_bs_model = None
        self.add_noise = False

        # to save total number of points queried
        self.return_num_points = False

    def measure_emittance(self):
        # get initial points from the observer
        o = Observer([], {'x': [], 'y': []}, {'x': [], 'y': []})
        o.use_model = self.use_model
        o.inj_config = self.inj_config
        o.online = self.online
        o.meas_type = self.meas_type
        o.use_prev_meas = self.use_prev_meas
        o.tolerance = self.quad_tol

        # print warning
        if self.online and self.verbose:
            print("Running online!")
        else:
            print("Running offline.")

        # if using sim
        # set beamsize fn
        o.get_beamsizes_model = self.get_bs_model
        o.add_noise = self.add_noise

        o.config_name = self.config_name
        o.config_dict = self.config_dict

        energy = o.config_dict['beamline_info']['energy']
        l_quad = o.config_dict['beamline_info']['l']

        # get initial beamsizes (rough scan)
        bs_x_list, bs_y_list, bs_x_list_err, bs_y_list_err = o.measure_beam(self.quad_init)

        quad_range_x = self.quad_init
        quad_range_y = self.quad_init

        if self.adapt_ranges:
            quad_range_x = adapt_range(quad_range_x,
                                       bs_x_list,
                                       'x',
                                       w=bs_x_list_err,
                                       energy=energy,
                                       l_eff=l_quad,
                                       num_points=self.num_points
                                       )
            quad_range_y = adapt_range(quad_range_y,
                                       bs_y_list,
                                       'y',
                                       w=bs_y_list_err,
                                       energy=energy,
                                       l_eff=l_quad,
                                       num_points=self.num_points
                                       )

            new_beamsize_x = o.measure_beam(quad_range_x)
            bs_x_list, bs_x_list_err = new_beamsize_x[0], new_beamsize_x[2]
            new_beamsize_y = o.measure_beam(quad_range_y)
            bs_y_list, bs_y_list_err = new_beamsize_y[1], new_beamsize_y[3]

        if self.check_sym:
            add_points_x = check_symmetry(quad_range_x, bs_x_list, bs_x_list_err, 'x',
                                          bs_fn=o.measure_beam, add_meas=True)
            add_points_y = check_symmetry(quad_range_y, bs_y_list, bs_y_list_err, 'y',
                                          bs_fn=o.measure_beam, add_meas=True)

            if add_points_x is not None:
                quad_range_x = add_points_x[0]
                bs_x_list = add_points_x[1]
                bs_x_list_err = add_points_x[2]

            if add_points_y is not None:
                quad_range_y = add_points_y[0]
                bs_y_list = add_points_y[1]
                bs_y_list_err = add_points_y[2]

        if self.infl_check:
            left_x, right_x = find_inflection_pnt(quad_range_x,
                                                  bs_x_list,
                                                  show_plots=self.show_plots
                                                  )
            left_y, right_y = find_inflection_pnt(quad_range_y,
                                                  bs_y_list,
                                                  show_plots=self.show_plots
                                                  )

            # truncate data
            quad_range_x = quad_range_x[left_x:right_x]
            bs_x_list = bs_x_list[left_x:right_x]
            bs_x_list_err = bs_x_list_err[left_x:right_x]

            quad_range_y = quad_range_y[left_y:right_y]
            bs_y_list = bs_y_list[left_y:right_y]
            bs_y_list_err = bs_y_list_err[left_y:right_y]

        if self.add_pnts:
            quad_range_x, bs_x_list, bs_x_list_err = add_measurements_btwn_pnts(quad_range_x,
                                                                                bs_x_list,
                                                                                bs_x_list_err,
                                                                                self.num_points,
                                                                                'x',
                                                                                bs_fn=o.measure_beam
                                                                                )
            quad_range_y, bs_y_list, bs_y_list_err = add_measurements_btwn_pnts(quad_range_y,
                                                                                bs_y_list,
                                                                                bs_y_list_err,
                                                                                self.num_points,
                                                                                'y',
                                                                                bs_fn=o.measure_beam
                                                                                )

        # finally get emittance
        ef = EmitCalc({'x': quad_range_x, 'y': quad_range_y},
                      {'x': bs_x_list, 'y': bs_y_list},
                      {'x': bs_x_list_err, 'y': bs_y_list_err},
                      config_dict=o.config_dict,
                      config_name=o.config_name
                      )
        ef.plot = self.show_plots
        ef.save_runs = self.save_runs
        ef.calc_bmag = self.calc_bmag

        # get normalized transverse emittance
        ef.get_emit()

        # save total number of points queried
        if self.return_num_points:
            ef.output["total_points_measured"] = len(o.quad_meas)

        return ef.output
