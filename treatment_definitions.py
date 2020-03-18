import multifil.aws.metas as metas
import numpy as np
from numpy import radians


def general_form():
    """This is an example form, for three time_steps without titin"""
    time = [0.0, 0.5, 1.0]
    length = [1250, 1250, 1250]
    ap = [1.0, 1.0, 1.0]
    params = {"ti_a": 0}

    return time, length, ap, params


def t_spring(stiffness=40, ms_p_ts=0.5, duration=500):
    time = metas.time_trace(
        ms_p_ts,  # ms per timestep
        duration,  # ms to run for
    )
    
    length = metas.zline_workloop(
        mean=1200,  # resting hs length
        amp=0,  # peak to peak amp
        freq=1,  # cycle freq in Hz
        time=time  # time_trace
    )
    
    ap = metas.actin_permissiveness_workloop(
        freq=1,  # freq in Hz
        phase=0.01,  # phase mean
        stim_duration=20,  # stimulus duration in ms
        influx_time=2,  # millis it takes for ca to go from 10 to 90% of influx level
        half_life=50,  # half life of Ca decay
        time=time  # time_trace
    )

    params = {"mh_c_ks":stiffness, "mh_c_kw":stiffness}

    return time, length, ap, params


"""Utilities"""


def stepsize_to_tht_mag(step_size):
    c_rest_weak = radians(47.16)
    g_rest_weak = 19.93
    # weak_ls = g_rest_weak * np.sin(c_rest_weak)

    c_rest_strong = radians(73.20)
    g_rest_strong = 16.47
    strong_ls = g_rest_strong * np.sin(c_rest_strong)

    delta_x_weak = g_rest_weak * np.sin(np.pi / 2 - c_rest_weak)
    # delta_x_strong = g_rest_strong * np.sin(np.pi / 2 - c_rest_strong)

    # Begin calculations
    needed = delta_x_weak - step_size  # What delta_x do we need 

    # tan (theta) = needed / strong_ls
    tht = np.arctan(needed / strong_ls)

    tht = np.pi / 2 - tht  # change to mh_converter angle
    mag = np.sqrt(needed * needed + strong_ls * strong_ls)
    return tht, mag
