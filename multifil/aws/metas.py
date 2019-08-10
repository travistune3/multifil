#!/usr/bin/env python
# encoding: utf-8
"""
metas.py - create the meta file that will configure a run

This was originally included in the run.py file, but has grown compliated
enough to warrent the creation of a fully separate management system

metas.emit produces a meta file that describes what we want a run to do: the
values of the z_line, lattice spacing, and actin permissiveness through the run
and where it will be stored after completion.

Example
--------
>>> freq, phase = 10, .8
>>> timetrace = metas.time(.1, 10)
>>> zline = metas.zline_workloop(1250, 25, freq, timetrace)
>>> activation = metas.actin_permissiveness_workloop(freq, phase, 10, 3, 3,
...  timetrace)
>>> metas.emit('./', None, timetrace, z_line=zline,
...  actin_permissiveness=activation, write=False, phase=phase, freq=freq)
{'actin_permissiveness': None,
...  'actin_permissiveness_func': None,
...  'comment': None,
...  'lattice_spacing': None,
...  'lattice_spacing_func': None,
...  'name': ...,
...  'path_local': './',
...  'path_s3': None,
...  'timestep_length': 0.1,
...  'timestep_number': 100,
...  'z_line': None,
...  'z_line_func': None}


Created by Dave Williams on 2017-03-08
"""

import os
import uuid
import ujson as json
import numpy as np

from .. import hs

## Define traces to be used in runs
def time(timestep_length, run_length_in_ms):
    """Create a time series in ms. This is easily doable through other methods
    but this documents it a bit.

    timestep_length: float
        Length of a timestep in ms
    run_length_in_ms: int
        Number of timesteps run is simulated for
    """
    return np.arange(0, run_length_in_ms, timestep_length)

def zline_workloop(mean, amp, freq, time):
    """A sinusoidal oscillatory length trace.

    Parameters:
        offset: resting z-line value, will start here
        amp: peak-to-peak amplitude
        freq: frequency of oscillation
        time: time trace in ms to provide length trace for
    """
    period = 1000/freq
    zline = mean + 0.5 * amp * np.cos(2*np.pi*time/period)
    return zline

def zline_forcevelocity(L0, hold_time, L0_per_sec, time):
    """Takes initial length, time to hold there in ms, & shortening in L0/sec"""
    # Things we need to know for the shortening
    number_of_timesteps = len(time) #for ease of reading
    timestep_length = np.diff(time)[0]
    hold_steps = int(hold_time/timestep_length)
    shorten_steps = number_of_timesteps - hold_steps
    nm_per_step = timestep_length * 1/1000 * L0_per_sec * L0
    # Construct the length signal
    zline = [L0 for i in range(hold_steps)]
    for i in range(shorten_steps):
        zline.append(zline[-1] - nm_per_step)
    return zline

def actin_permissiveness_time_independent_workloop(freq, phase, influx_factor, duty_factor, half_life, time, max_signal=1.0, min_signal=0.054):
    steps_per_millis = 0
    timepoint = time[0]
    while timepoint < 1.0 and steps_per_millis < len(time) - 1:
        steps_per_millis += 1
        timepoint = time[steps_per_millis]
    millis_per_cycle = 1 / (freq / 1000)
    steps_per_cycle = int(np.round(steps_per_millis * millis_per_cycle))
    timestep_length = 1 / steps_per_millis
    steps_on = int(duty_factor * steps_per_cycle)
    steps_off = (steps_per_cycle - steps_on)
    
    print("SPC", steps_per_cycle)
    
    out = [min_signal]
    max_step =  int(influx_factor * steps_per_cycle)
    if max_step == 0:
        max_step = 1
    
    growth_step = influx_factor * millis_per_cycle / 4.2
    growth_rate = -np.log(0.5) / growth_step
    
    decay_step = half_life * millis_per_cycle / 4.2
    decay_rate = np.log(0.5) / decay_step
    
    step = 0
    while(len(out) < len(time) + steps_per_cycle):
        if step == max_step - 1:
            out.append(max_signal)
        elif step < steps_on - 1:
            new_val = out[-1] * (np.e ** (growth_rate * timestep_length))
            if new_val >= max_signal:
                out.append(max_signal)
            elif new_val < min_signal:
                out.append(min_signal)
            else:
                out.append(new_val)
        elif step == steps_on - 1:
            out.append(max_signal)
        else:
            new_val = out[-1] * (np.e ** (decay_rate * timestep_length))
            if new_val <= min_signal:
                out.append(min_signal)
            else:
                out.append(new_val)
        step += 1
        if step == steps_per_cycle:
            step = 0
    '''#Things we need to know for smoothing
    sd = 1 #standard deviation of smoothing window in ms
    sw = 3 #smoothing window in ms
    base_normal = np.exp(-np.arange(-sw,sw,timestep_length)**2/(2*sd**2))
    normal = base_normal/sum(base_normal)
    out = np.convolve(normal, out)
    #'''
    
    phase_index = steps_per_cycle * phase
    stt = int(np.round(phase_index))
    end = len(time) + stt
    
    print(phase_index, stt, end)
    return out[stt:end]

def actin_permissiveness_workloop(freq, phase, stim_duration,
                                  influx_time, half_life, time, max_signal=1.0):
    """Requires cycle frequency, phase relative to longest length
    point, duration of on time, time from 10 to 90% influx level, and
    the half-life of the Ca2+ out-pumping.
    """
    # Convert frequency to period in ms
    period = 1000/freq
    # Things we need to know for the shape of a single cycle
    decay_rate = np.log(1/2)/half_life
    growth_rate = 0.5*influx_time
    max_signal = max_signal
    # Things we need to know for the cyclical nature of the signal
    number_of_timesteps = len(time) #for ease of reading
    timestep_length = np.diff(time)[0]
    cycle_step_number = int(period/timestep_length)
    cycle_time_trace = np.arange(0, period, timestep_length)
    try:
        steps_before_stim = np.argwhere(
            cycle_time_trace>=(period*(phase%1)))[0][0]
    except IndexError:
        assert 0 == len(np.argwhere(
            cycle_time_trace>=(period*(phase%1))))
        steps_before_stim = 0 #b/c phase was 0.999 or similar
    stim_step_number = int(stim_duration/timestep_length)
    no_stim_step_number = cycle_step_number - stim_step_number
    # Things we need to know for smoothing
    sd = 1 #standard deviation of smoothing window in ms
    sw = 3 #smoothing window in ms
    base_normal = np.exp(-np.arange(-sw,sw,timestep_length)**2/(2*sd**2))
    normal = base_normal/sum(base_normal)
    # Step through, generating signal
    out = [0.1]
    for i in range(steps_before_stim):
        out.append(out[-1])
    while len(out)<(4*cycle_step_number+number_of_timesteps):
        for i in range(stim_step_number):
            growth = timestep_length * out[-1] * (growth_rate) *\
                    (1-out[-1]/max_signal)
            out.append(out[-1]+growth)
        for i in range(no_stim_step_number):
            decay = timestep_length * out[-1] * decay_rate
            out.append(out[-1]+decay)
    # Smooth signal
    out = np.convolve(normal, out)
    return out[2*cycle_step_number:2*cycle_step_number+number_of_timesteps]


## Configure a run via a saved meta file
def emit(path_local, path_s3, time,  poisson=0.0, ls=None, z_line=None,
         actin_permissiveness=None, comment = None, 
         write = True, titin_params=None, **kwargs):
    """Produce a structured JSON file that will be consumed to create a run

    Import emit into an interactive workspace and populate a directory with
    run configurations to be executed by a cluster.

    Parameters
    ----------
    path_local: string
        The local (absoluzte or relative) directory to which we save both
        emitted files and run output.
    path_s3: string
        The s3 bucket (and optional folder) to save run output to and to which
        the emitted files should be uploaded.
    time: iterable
        Time trace for run, in ms
    poisson: float
        poisson ratio of lattice. 0.5 const vol; 0 default const lattice;
        negative for auxetic
    ls: float, optional
        Specifies the initial starting lattice spacing which will act as a
        zero or offset for the spacing. If not given, the default lattice
        spacing from hs.hs will be used.
    z_line: float or iterable, optional
        If not given, default distance specified in hs.hs is used. If given as
        float, the z-line distance for the run. If given as an iterable, used as
        trace for run, timestep by timestep.
    actin_permissiveness: float or iterable, optional
        Same as for z-line.
    comment: string, optional
        Space for comment on the purpose or other characteristics of the run
    write: bool, optional
        True (default) writes file to path_local/name.meta.json. Other values
        don't. In both cases the dictionary describing the run is returned.
    **kwargs:
        Further keyword args will be included in the output dictionary. These
        are used to sort the resulting runs by their properties of interest.
        For example, where we are varying phase of activation across a series
        of runs we would include the argument, e.g. 'phase=0.2', in order to
        sort over phase when looking at results.

    Returns
    -------
    rund: dict
        Copy of run dictionary saved to disk as json.

    Examples
    --------
    >>> emit('./', None, .1, 100, write=False)
    {'actin_permissiveness': None,
    ...  'actin_permissiveness_func': None,
    ...  'comment': None,
    ...  'lattice_spacing': None,
    ...  'lattice_spacing_func': None,
    ...  'name': ...,
    ...  'path_local': './',
    ...  'path_s3': None,
    ...  'timestep_length': 0.1,
    ...  'timestep_number': 100,
    ...  'z_line': None,
    ...  'z_line_func': None}
    """
    rund = {}
    name = str(uuid.uuid1())
    ## Build dictionary
    rund['name'] = name
    rund['comment'] = comment
    rund['path_local'] = path_local
    rund['path_s3'] = path_s3
    rund['poisson_ratio'] = poisson
    rund['lattice_spacing'] = ls
    rund['z_line'] = z_line
    rund['actin_permissiveness'] = actin_permissiveness
    rund['titin_params'] = titin_params
    rund['timestep_length'] = np.diff(time)[0]
    rund['timestep_number'] = len(time)
    ## Include kwargs
    for k in kwargs:
        rund[k] = kwargs[k]
    ## Write out the run description
    if write is True:
        output_filename = os.path.join(path_local, name+'.meta.json')
        with open(output_filename , 'w') as metafile:
            json.dump(rund, metafile, indent=4)
        return rund, output_filename
    return rund
