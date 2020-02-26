# import multifil from forked repository

import warnings
from .. import hs

import os
import ujson as json
import multiprocessing as mp
import time as millis
import uuid

warnings.filterwarnings(action='once')


def run_async(profiles, force=False):
    if len(profiles) >= mp.cpu_count():
        if not force:
            max_threads = max(2, int(0.75 * mp.cpu_count()))
            profiles = profiles[0:max_threads]
        else:
            print("More instances than threads, forcing anyway")

    start = millis.time()

    pool = mp.Pool(len(profiles))

    try:
        processes = pool.map(run_model, profiles)
    except KeyboardInterrupt:
        pool.terminate()

    print(processes)
    end = millis.time()
    print(end - start, "seconds")


def run_model(profile):
    actin = profile['actin']
    if 'output_dir' in profile.keys():
        output_dir = profile['output_dir']
    else:
        cur_dir = os.path.abspath(os.curdir)
        output_dir = cur_dir + "\\_data\\"
        os.makedirs(output_dir, exist_ok=True)
        dir_warning = 'output directory not specified. Defaulting to\n\t' + str(output_dir)
        warnings.warn(dir_warning)

    bar = True
    if 'tm_bar' in profile.keys():
        bar = hs.hs.tm_bar

    every = 100
    if 'every' in profile.keys():
        every = profile['every']

    td = {'pCa': actin}
    sarc = hs.hs(time_dependence=td, timestep_len=0.5)
    ts = len(actin) - 1

    start_model = millis.time()
    run_id = str(uuid.uuid1())
    result, exit_code = sarc.run(time_steps=ts, every=every, bar=bar)
    file_path = output_dir + str(run_id) + ".data.json"
    with open(file_path, 'w') as outputFile:
        json.dump(result, outputFile)
    end_model = millis.time()

    print("model took", end_model - start_model, "seconds")
    return result
