from pyemittance.observer import Observer
from pyemittance.data_handler import adapt_range, check_symmetry, find_inflection_pnt, add_measurements_btwn_pnts
from pyemittance.emittance_calc import EmitCalc
from pyemittance.load_json_configs import load_configs

# This file is OBSOLETE, please use the pyemittance.PyEmittance class


def eval_emit_machine(inj_config=None,
                      quad_init=[-6, -4, -2, 0],
                      online=False,
                      config_name='LCLS_OTR2',
                      config_dict=None,  # supersedes json configs
                      meas_type='OTRS',
                      adapt_ranges=True,
                      num_points=7,
                      check_sym=True,
                      infl_check=True,
                      add_pnts=True,
                      show_plots=True,
                      use_prev_meas=True,
                      quad_tol=0.05,
                      save_runs=False,
                      calc_bmag=False
                      ):

    # get initial points from the observer
    o = Observer([], {'x': [], 'y': []}, {'x': [], 'y': []})
    o.use_model = False
    o.inj_config = inj_config
    o.online = online
    o.meas_type = meas_type
    o.use_prev_meas = use_prev_meas
    o.tolerance = quad_tol

    # if config is not provided, use LCLS OTR2 as default
    o.config_name = config_name
    o.config_dict = config_dict if config_dict else load_configs(o.config_name)

    energy = o.config_dict['beamline_info']['energy']
    l_quad = o.config_dict['beamline_info']['l']

    # get initial beamsizes (rough scan)
    bs_x_list, bs_y_list, bs_x_list_err, bs_y_list_err = o.measure_beam(quad_init)

    quad_range_x = quad_init
    quad_range_y = quad_init

    if adapt_ranges:
        quad_range_x = adapt_range(quad_range_x, bs_x_list, 'x', w=bs_x_list_err, energy=energy, l_eff=l_quad, num_points=num_points)
        quad_range_y = adapt_range(quad_range_y, bs_y_list, 'y', w=bs_y_list_err, energy=energy, l_eff=l_quad, num_points=num_points)

        new_beamsize_x = o.measure_beam(quad_range_x)
        bs_x_list, bs_x_list_err = new_beamsize_x[0], new_beamsize_x[2]
        new_beamsize_y = o.measure_beam(quad_range_y)
        bs_y_list, bs_y_list_err = new_beamsize_y[1], new_beamsize_y[3]
    else:
        quad_range_x = quad_init
        quad_range_y = quad_init

    if check_sym:
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

    if infl_check:
        left_x, right_x = find_inflection_pnt(quad_range_x, bs_x_list, show_plots=show_plots)
        left_y, right_y = find_inflection_pnt(quad_range_y, bs_y_list, show_plots=show_plots)

        # truncate data
        quad_range_x = quad_range_x[left_x:right_x]
        bs_x_list = bs_x_list[left_x:right_x]
        bs_x_list_err = bs_x_list_err[left_x:right_x]

        quad_range_y = quad_range_y[left_y:right_y]
        bs_y_list = bs_y_list[left_y:right_y]
        bs_y_list_err = bs_y_list_err[left_y:right_y]

    if add_pnts:
        quad_range_x, bs_x_list, bs_x_list_err = add_measurements_btwn_pnts(quad_range_x, bs_x_list, bs_x_list_err,
                                                                        num_points, 'x', bs_fn=o.measure_beam)
        quad_range_y, bs_y_list, bs_y_list_err = add_measurements_btwn_pnts(quad_range_y, bs_y_list, bs_y_list_err,
                                                                        num_points, 'y', bs_fn=o.measure_beam)

    # finally get emittance
    ef = EmitCalc({'x': quad_range_x, 'y': quad_range_y},
                  {'x': bs_x_list, 'y': bs_y_list},
                  {'x': bs_x_list_err, 'y': bs_y_list_err},
                  config_dict=o.config_dict,
                  config_name=o.config_name
                  )
    ef.plot = show_plots
    ef.save_runs = save_runs
    ef.calc_bmag = calc_bmag

    # get normalized transverse emittance
    ef.get_emit()

    # save total number of points queried
    ef.output["total_points_measured"] = len(o.quad_meas)

    return ef.output


def eval_emit_surrogate(get_bs_model,
                        inj_config=None,
                        quad_init=[-6, -4, -2, 0],
                        config_name='LCLS_OTR2',
                        config_dict=None,
                        adapt_ranges=True,
                        num_points=7,
                        check_sym=True,
                        infl_check=True,
                        add_pnts=True,
                        show_plots=False,
                        add_noise=False,
                        use_prev_meas=False,
                        quad_tol=0.05,
                        save_runs=False,
                        calc_bmag=False):

    # get initial points from the observer
    o = Observer([], {'x': [], 'y': []}, {'x': [], 'y': []})
    o.use_model = True

    # set beamsize fn for MODEL
    o.get_beamsizes_model = get_bs_model
    o.inj_config = inj_config
    o.add_noise = add_noise
    o.use_prev_meas = use_prev_meas
    o.tolerance = quad_tol

    # get initial beamsizes (rough scan)
    bs_x_list, bs_y_list, bs_x_list_err, bs_y_list_err = o.measure_beam(quad_init)

    quad_range_x = quad_init
    quad_range_y = quad_init

    # if config is not provided, use LCLS OTR2 as default
    o.config_name = config_name
    o.config_dict = config_dict if config_dict else load_configs(o.config_name)

    energy = o.config_dict['beamline_info']['energy']
    l_quad = o.config_dict['beamline_info']['l']

    if adapt_ranges:
        quad_range_x = adapt_range(quad_range_x, bs_x_list, 'x', w=bs_x_list_err, energy=energy, l_eff=l_quad, num_points=num_points)
        quad_range_y = adapt_range(quad_range_y, bs_y_list, 'y', w=bs_y_list_err, energy=energy, l_eff=l_quad, num_points=num_points)

        new_beamsize_x = o.measure_beam(quad_range_x)
        bs_x_list, bs_x_list_err = new_beamsize_x[0], new_beamsize_x[2]

        new_beamsize_y = o.measure_beam(quad_range_y)
        bs_y_list, bs_y_list_err = new_beamsize_y[1], new_beamsize_y[3]

    if check_sym:
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

    if infl_check:
        left_x, right_x = find_inflection_pnt(quad_range_x, bs_x_list, show_plots=show_plots)
        left_y, right_y = find_inflection_pnt(quad_range_y, bs_y_list, show_plots=show_plots)

        # truncate data
        quad_range_x = quad_range_x[left_x:right_x]
        bs_x_list = bs_x_list[left_x:right_x]
        bs_x_list_err = bs_x_list_err[left_x:right_x]

        quad_range_y = quad_range_y[left_y:right_y]
        bs_y_list = bs_y_list[left_y:right_y]
        bs_y_list_err = bs_y_list_err[left_y:right_y]

    if add_pnts:
        quad_range_x, bs_x_list, bs_x_list_err = add_measurements_btwn_pnts(quad_range_x, bs_x_list, bs_x_list_err,
                                                                        num_points, 'x', bs_fn=o.measure_beam)
        quad_range_y, bs_y_list, bs_y_list_err = add_measurements_btwn_pnts(quad_range_y, bs_y_list, bs_y_list_err,
                                                                        num_points, 'y', bs_fn=o.measure_beam)

    # finally get emittance
    ef = EmitCalc({'x': quad_range_x, 'y': quad_range_y},
                  {'x': bs_x_list, 'y': bs_y_list},
                  {'x': bs_x_list_err, 'y': bs_y_list_err},
                  config_dict=o.config_dict,
                  config_name=o.config_name
                  )
    ef.plot = show_plots
    ef.save_runs = save_runs
    ef.calc_bmag = calc_bmag

    # get normalized transverse emittance
    ef.get_emit()

    # save total number of points queried
    ef.output["total_points_measured"] = len(o.quad_meas)

    return ef.output