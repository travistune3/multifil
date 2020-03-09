import multifil

def general_form():
    time = [0.0, 0.5, 1.0]
    length = [1250, 1250, 1250]
    ap = [1.0, 1.0, 1.0]
    
    return time, length, ap

def alison_trace():
    time = multifil.aws.metas.time(
    .5, # ms per timestep
    500, # ms to run for
                                    )
    length = multifil.aws.metas.zline_workloop(
        900, # resting hs length 
        0, # peak to peak amp 
        1, # cycle freq in Hz 
        time
                                                    )

    '''init params freq=25, phase=0.1, stim_duration=20, influx_time=2, half_life=10, time=time'''
    '''first twitch: 5, 10, 50, 2, 10'''
    '''positive workloop params: 25, # freq in Hz, 0.1, # phase offset, 10, # stimulus duration in ms
        2, # time it takes for ca to go from 10 to 90% of influx level, 10, # half life of Ca decay, time_trace'''
    ap = multifil.aws.metas.actin_permissiveness_workloop(
        1, # freq in Hz 
        0.01, # phase offset 
        20, # stimulus duration in ms
        2, # time it takes for ca to go from 10 to 90% of influx level
        50, # half life of Ca decay 
        time)
    
    return time, length, ap