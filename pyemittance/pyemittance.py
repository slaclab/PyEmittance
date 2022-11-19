

from pyemittance.observer import Observer
from pyemittance.data_handler import adapt_range, check_symmetry, find_inflection_pnt, add_measurements_btwn_pnts
from pyemittance.emittance_calc import EmitCalc
from pyemittance.emittance_calc_multiwire import MultiWireCalc
from pyemittance.load_json_configs import load_configs

# For setting quad back after scan is done
from pyemittance.beam_io import MachineIO

# TODO: set up unit testing for PyEmittance, EmitCalc, MultiWireCalc, Observer
class PyEmittance:

    def __init__(self,
                 emit_calc_type='quadscan',  # 'quadscan' or 'multiwire'
                 config_name='LCLS_OTR2',  # make sure config corresponds to calc type
                 config_dict=None,  # supersedes json configs
                 meas_type='OTRS',  # only relevant if emit_calc_type=='quadscan'
                 use_model=False,
                 online=False
                 ):

        self.emit_calc_type = emit_calc_type
        self.meas_type = meas_type

        # if config is not provided, use LCLS-OTR2 as default
        self.config_name = config_name
        self.config_dict = config_dict if config_dict else load_configs(self.config_name)

        # if running on machine, use_model=False
        self.use_model = use_model
        # only True if setting PVs
        self.online = online
        self.verbose = True

        # injector settings (SOL, CQ, SQ) if optimizing
        # TODO: remove injector settings options from all modules (optimizer needs to be separate)
        self.inj_config = None
        if self.emit_calc_type == 'quadscan':
            # initial rough quad scan
            self.quad_init = [-6, -4, -2, 0]

        # pyemittance method options
        self.save_runs = False
        self.calc_bmag = False
        self.show_plots = True
        if self.emit_calc_type == 'quadscan':
            self.adapt_ranges = True
            self.num_points = 7
            self.check_sym = True
            self.infl_check = True
            self.add_pnts = True
            self.use_prev_meas = True
            self.quad_tol = 0.05

        # simulation/model options
        # beamsize function from model
        self.get_bs_model = None
        self.add_noise = False

        # to save total number of points queried
        self.return_num_points = False

    def measure_emittance(self):
        if self.emit_calc_type == 'quadscan':
            self.measure_emittance_quad_scan()
        elif self.emit_calc_type == 'multiwire':
            self.measure_emittance_multiwire()
        else:
            raise Exception("Cannot perform measurement. 'emit_calc_type' needs to be 'quadscan' or 'multiwire'.")

    def measure_emittance_quad_scan(self):
        # get initial points from the observer
        o = Observer([], {'x': [], 'y': []}, {'x': [], 'y': []})
        o.use_model = self.use_model
        o.inj_config = self.inj_config
        o.online = self.online
        o.meas_type = self.meas_type
        o.use_prev_meas = self.use_prev_meas
        o.tolerance = self.quad_tol
        o.config_name = self.config_name
        o.config_dict = self.config_dict

        # print warning
        if self.online and self.verbose:
            print("Running online!")
        else:
            print("Running offline.")

        # if using sim
        # set beamsize fn
        o.get_beamsizes_model = self.get_bs_model
        o.add_noise = self.add_noise

        energy = o.config_dict['beamline_info']['energy']
        l_quad = o.config_dict['beamline_info']['l']

        # Save init quad value
        io = MachineIO(self.config_name, self.config_dict, self.meas_type)
        io.online = self.online
        quad_init = io.getquad()

        # Get initial beamsizes (rough scan)
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

        # Put quad back
        io.online = self.online
        io.setquad(quad_init)

        # Finally get emittance
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
            ef.out_dict["total_points_measured"] = len(o.quad_meas)

        self.results = ef.out_dict
        return ef.out_dict

    def measure_emittance_multiwire(self):
        # get wire measurements from the observer
        o = Observer()
        o.inj_config = self.inj_config
        o.online = self.online
        o.meas_type = 'WIRE'
        o.emit_calc_type = self.emit_calc_type
        o.config_name = self.config_name
        o.config_dict = self.config_dict

        # print warning
        if self.online and self.verbose:
            print("Running online!")
        else:
            print("Running offline.")

        # Get beamsizes
        # This returns lists now: xrms, yrms, xrms_err, yrms_err
        # TODO: implement error dict to be returned to track multiwire success
        bs_x_list, bs_y_list, bs_x_list_err, bs_y_list_err = o.multiwire_measure_beam()

        # Finally get emittance
        ef = MultiWireCalc(beam_vals={'x': bs_x_list, 'y': bs_y_list},
                           beam_vals_err={'x': bs_x_list_err, 'y': bs_y_list_err},
                           config_dict=o.config_dict,
                           config_name=o.config_name
                           )
        ef.plot = self.show_plots
        ef.save_runs = self.save_runs
        ef.calc_bmag = self.calc_bmag

        # get normalized transverse emittance
        ef.get_emit()
        self.results = ef.out_dict

        return ef.out_dict