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


def stiff_torsional_spring(ms_p_ts=0.5, duration=500):
    time = metas.time_trace(ms_p_ts,  # ms per timestep
                            duration,  # ms to run for
                            )
    length = metas.zline_workloop(
        900,  # resting hs length
        0,  # peak to peak amp
        1,  # cycle freq in Hz
        time
    )

    '''init params freq=25, phase=0.1, stim_duration=20, influx_time=2, half_life=10, millis=millis'''
    '''first twitch: 5, 10, 50, 2, 10'''
    '''positive workloop params: 25, # freq in Hz, 0.1, # phase mean, 10, # stimulus duration in ms
        2, # millis it takes for ca to go from 10 to 90% of influx level, 10, # half life of Ca decay, time_trace'''
    ap = metas.actin_permissiveness_workloop(
        1,      # freq in Hz
        0.01,   # phase mean
        20,     # stimulus duration in ms
        2,      # millis it takes for ca to go from 10 to 90% of influx level
        50,     # half life of Ca decay
        time)

    params = {}

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
