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
        timestep_length = ms_p_ts,  # ms per timestep
        run_length_in_ms = duration,  # ms to run for
    )
    
    length = 900  # nanometers
#     length = metas.zline_workloop(
#         mean=900,  # resting hs length
#         amp=0,  # peak to peak amp
#         freq=1,  # cycle freq in Hz
#         time=time  # time_trace
#     )
    
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

def davis_t_spring(stiffness=40):
    time = metas.time_trace(
        timestep_length=0.5,  # ms per timestep
        run_length_in_ms=500,  # ms to run for
    )
    
    length = 900  # nanometers
#     length = metas.zline_workloop(
#         mean=900,  # resting hs length
#         amp=0,  # peak to peak amp
#         freq=1,  # cycle freq in Hz
#         time=time  # time_trace
#     )
    
    ap = davis_calcium_transient()[0:1000]  # trim to 0.5 s
#     ap = metas.actin_permissiveness_workloop(
#         freq=1,  # freq in Hz
#         phase=0.01,  # phase mean
#         stim_duration=20,  # stimulus duration in ms
#         influx_time=2,  # millis it takes for ca to go from 10 to 90% of influx level
#         half_life=50,  # half life of Ca decay
#         time=time  # time_trace
#     )

    params = {"mh_c_ks":stiffness, "mh_c_kw":stiffness}

    return time, length, ap, params


"""Utilities"""


def davis_calcium_transient(treatment="WT"):
    import json
    with open("davis_2016_ca_as_ap.json", 'r') as ca_file:
        return json.load(ca_file)[treatment]


def stepsize_to_tht_mag(step_size, mh_c_rw=radians(47.16), mh_g_rw=19.93):
    """ Parameters:
            step_size: how much the spring should change along the axial direction during a power stroke
            mh_c_rw: the weak state angle - optional. If changed, make sure that the model is changing mh_c_rw as well.
            mh_g_rw: the weak state length - optional. If changed, make sure that the model is changing mh_g_rw as well.
        Returns:
            mh_c_rs: what the strong state angle should be set to so that the power stroke has the specified step size
            mh_g_rs: what the strong state length should be set to so that the power stroke has the specified step size
    """
    c_rest_weak = mh_c_rw
    g_rest_weak = mh_g_rw
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

    mh_c_rs = np.pi / 2 - tht  # change to mh_converter angle
    mh_g_rs = np.sqrt(needed * needed + strong_ls * strong_ls)
    return mh_c_rs, mh_g_rs
