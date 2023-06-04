import os, json

import logging
logger = logging.getLogger(__name__)

this_dir, this_filename = os.path.split(__file__)
CONFIG_PATH = os.path.join(this_dir, "configs")

# default jsons that include the configurations for the emittance scan
json_namelist = [
    "beamline_info",
    "img_proc",
    "meas_pv_info",
  #   "opt_pv_info",
  #  "save_scalar_pvs",
    "savepaths",
]


def load_configs(dir_name="LCLS2_OTR0H04"):
    all_data = {}
    for i in range(len(json_namelist)):
        # load all jsons and save into one dict
        # skip jsons that don't exit and print that it is not found
        # TODO: validate that all jsons/needed keywords configs do exist
        # TODO: validate that all configs are consistent across directories/locations
        # TODO: eg all have rmatx and rmaty
       
        fname = os.path.join(CONFIG_PATH + "/" + dir_name, json_namelist[i] + ".json")
        if not os.path.exists(fname):
            raise FileNotFoundError(
                f"*** File '{json_namelist[i]}.json' does not exist,"
                f" please create appropriate json file for configuration. *** \n"
                f"*** Or alternatively, initialize EmitCalc with dict directly. ***"
            )
        # make dict keys json file names for reference
        data = {json_namelist[i]: json.load(open(fname))}
        all_data = {**all_data, **data}

    # dict of all configs
    return all_data


if __name__ == "__main__":
    all_data = load_configs()
    print(all_data.keys())
    print(all_data.values())
