#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:48:30 2019

@author: tfernandes
"""
from __future__ import division
from math import pi

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(1,'/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')
# os.chdir('/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')
from Sequence.sequence import Sequence
from calc_duration import calc_duration
from make_adc import make_adc
from make_delay import make_delay
from make_sinc_pulse import make_sinc_pulse
from make_trap_pulse import make_trapezoid
from make_arbitrary_grad import make_arbitrary_grad
from make_extended_trapezoid import make_extended_trapezoid
from opts import Opts
sys.path.insert(1,'/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq/Sequence')


def generate_SeqFile_Spiral(gx, gy, tr, n_shots, mg, ms, fA, n_slices, reps, st, tPlot, tReport, rf_offset=0):
    
#%% ===========================================================================
    # --- 1 - Instantiation and gradient limits ---

    system = Opts(max_grad=mg, grad_unit='mT/m', max_slew=ms, slew_unit='T/m/s', rf_ringdown_time=20e-6, rf_dead_time=100e-6,
           adc_dead_time=10e-6)
    seq = Sequence(system)

    # I need to check the manually inputed inverse raster time and the actual value are the same.
    i_raster_time = 100000
    assert 1/i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."

    slice_thickness, slice_gap, TE, TR = st, 15e-3, 5e-3, tr
    delta_z = n_slices * slice_thickness
    z = np.linspace((-delta_z / 2), (delta_z / 2), n_slices) + rf_offset
    
    #%% --- 2 - Gradients
#    scl = 1.05
#    gx  = gx/(scl*np.max(gx)/(max_grad*42576000e-3));
#    gy  = gy/(scl*np.max(gy)/(max_grad*42576000e-3));
    
    G   = gx + 1J*gy
    #%% --- 3 - Slice Selection
    # =========
    # Create 90 degree slice selection pulse and gradient
    # =========

    # RF90
    flip90 = fA * pi / 180

    rf, gz, gzr = make_sinc_pulse(flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                              apodization=0.5, time_bw_product=4, system=system)

    # RF180


    #%% --- 4 - ADCs / Readouts
    adc_samples = math.floor(len(G) / 4) * 4 - 4
    adc = make_adc(num_samples=adc_samples, duration=adc_samples/i_raster_time, system=system)
    

    #%% --- 5 - Spoilers
    pre_time   = 8e-4
    n_pre_time = math.ceil(pre_time*i_raster_time)

    gz_rephz = make_trapezoid(channel='z', area=-gz.area / 2, duration=1e-3, system=system)    
    gz_spoil = make_trapezoid(channel='z', area=-gz.area / 2, duration=3 * n_pre_time/i_raster_time, system=system)


    aux_gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(np.real(G[:, 0])), system=system)

    # Make the spiral in Gx and Gy finish in zero - I use pre_time because I know for sure it's long enough.
    # Furthermore, this is after readout and TR is supposed to be long.
    spoil_time = [0, calc_duration(gz_spoil)]
    gxAmplitud_Spoil = [np.real(G[-1][0]),0]
    gyAmplitud_Spoil = [np.imag(G[-1][0]),0]
    
    gx_spoil = make_extended_trapezoid(channel='x',times=spoil_time,amplitudes=gxAmplitud_Spoil,system=system)
    gy_spoil = make_extended_trapezoid(channel='y',times=spoil_time,amplitudes=gyAmplitud_Spoil,system=system)
    

#    # =========
#    # Prephase and Rephase
#    # =========
#    phase_areas = (np.arange(Ny) - (Ny / 2)) * delta_k
#    kwargs_for_gy_pre = {"channel": 'y', "system": system, "area": phase_areas[-1], "duration": 2e-3}
#    gy_pre = maketrapezoid(kwargs_for_gy_pre)
#
#    kwargs_for_gxpre = {"channel": 'x', "system": system, "area": -gx.area / 2, "duration": 2e-3}
#    gx_pre = maketrapezoid(kwargs_for_gxpre)
#
#    kwargs_for_gz_reph = {"channel": 'z', "system": system, "area": -gz.area / 2, "duration": 2e-3}
#    gz_reph = maketrapezoid(kwargs_for_gz_reph)
    
    #%% --- 6 - Calculate timing/Delays

    n_TE           = math.ceil(TE * i_raster_time)
    n_dur_gz_rephz = math.ceil(calc_duration(gz_rephz) * i_raster_time)
    n_dur_rf       = math.ceil(calc_duration(rf) * i_raster_time)

    n_del_preGRS   = n_TE - (n_dur_gz_rephz + math.ceil(n_dur_rf / 2))                # Time before input in points
    del_preGRS     = n_del_preGRS / i_raster_time
    delay_preGRS   = make_delay(del_preGRS)

    n_TR           = math.ceil(TR * i_raster_time)
    n_dur_aux_gx   = math.ceil(calc_duration(aux_gx) * i_raster_time)
    n_dur_gz_spoil = math.ceil(calc_duration(gz_spoil) * i_raster_time)

    n_delTR        = n_TR - (n_dur_gz_rephz + n_dur_rf - n_dur_aux_gx - n_dur_gz_spoil)  # Time after input to the system in points
    delTR          = n_delTR / i_raster_time
    delayTR        = make_delay(delTR)
        
    #%% --- 7 - Define sequence blocks/Readouts + Create '.seq' file
    for r in range(reps):
        for s in range(n_slices):
            for ns in range(n_shots):
                freq_offset = gz.amplitude * z[s]
                rf.freq_offset = freq_offset

                # RF90
                seq.add_block(rf, gz)
                seq.add_block(gz_rephz)

                # Delay for spiral
                seq.add_block(delay_preGRS)

                # if (r % 2 == 0):
                # Read k-space - Imaging Gradient waveforms
                gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(np.real(G[:, ns])),system=system)
                gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(np.imag(G[:, ns])),system=system)
                seq.add_block(gx, gy, adc)

                # Gradients Spoils
                seq.add_block(gz_spoil,gx_spoil,gy_spoil)
                # seq.add_block(gz_spoil)

                # else: # to allow measure phase with no gradient
                #     # Read k-space - Imaging Gradient waveforms
                #     seq.add_block(adc)
                #
                #     # Gradients Spoils
                #     seq.add_block(gz_spoil)

                # Delay to TR
                seq.add_block(delayTR)

    # %% --- 8 - Plot ADCs, GX, GY, GZ, RFpulse, RFphase
    # #    plt.figure()
    if tPlot:
        seq.plot()

    if tReport:
        print(seq.test_report())
        seq.check_timing()

    return seq
