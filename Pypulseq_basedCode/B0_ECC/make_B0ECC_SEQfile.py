#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wedn Sept  2 14:21:30 2020

@author: tfernandes
"""
from __future__ import division
from math import pi

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

PC = 0 # PC = 1 - Seia OR PC = 0 - My PC

if PC == 1:
    sys.path.insert(1, '/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')
    os.chdir('/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')
    from Sequence.sequence import Sequence
    from calc_duration import calc_duration
    from make_adc import make_adc
    from make_delay import make_delay
    from make_sinc_pulse import make_sinc_pulse
    from make_sinc_pulse_channel import make_sinc_pulse_channel
    from make_trap_pulse import make_trapezoid
    from make_arbitrary_grad import make_arbitrary_grad
    from make_extended_trapezoid import make_extended_trapezoid
    from opts import Opts

    sys.path.insert(1, '/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq/Sequence')

elif PC == 0:
    from pypulseq.Sequence.sequence import Sequence
    from pypulseq.calc_duration import calc_duration
    from pypulseq.make_adc import make_adc
    from pypulseq.make_delay import make_delay
    from pypulseq.make_sinc_pulse import make_sinc_pulse
    from pypulseq.make_trap_pulse import make_trapezoid
    from pypulseq.make_arbitrary_grad import make_arbitrary_grad
    from pypulseq.make_extended_trapezoid import make_extended_trapezoid
    from pypulseq.opts import Opts

    sys.path.insert(1,'D:/Tiago/Trabalho/2021_2025_PhD/Projects/Toolboxes/Python/Pypulseq_addon')
    from make_sinc_pulse_channel import make_sinc_pulse_channel



def generate_B0ECC_SeqFile(iT, tr, te, Npuls, mg, ms, grt, fA, reps, sliceT, tPre, tPos, tX, tY, tZ, tPlot, tReport, testPRE_POST, rf_offset, dummy, Ndummies):
    # %% ===========================================================================


    # %% --- 1 - Instantiation and gradient limits ---
    system = Opts(max_grad=mg, grad_unit='mT/m', max_slew=ms, slew_unit='T/m/s', rf_ringdown_time=20e-6,
                  rf_dead_time=100e-6,
                  adc_dead_time=10e-6)
    seq = Sequence(system)

    # I need to check the manually inputed inverse raster time and the actual value are the same.
    i_raster_time = 100000
    assert 1 / i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."

    slice_thickness, slice_gap, TR, TE, rf_offset_z = sliceT, 15e-3, tr, te, 0
    posAxxis = [slice_thickness + rf_offset, - slice_thickness - rf_offset]

    # %% --- 2 - Gradients
    G = iT
    G_zero = np.zeros(len(iT[:,0]))
    G_pre  = np.zeros(math.floor( tPre * i_raster_time ))
    G_pos  = np.zeros(math.floor( tPos * i_raster_time ))


    # %% --- 3 - Slice Selection
    # =========
    # Create 90 degree slice selection pulse and gradient
    # =========

    # RF90
    flip90 = fA * pi / 180

    rfx, gxRF, gxr = make_sinc_pulse_channel(channel='x', flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4, system=system)
    rfy, gyRF, gyr = make_sinc_pulse_channel(channel='y', flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4, system=system)
    rfz, gzRF, gzr = make_sinc_pulse(flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4, system=system)


    # %% --- 4 - ADCs / Readouts
    adc_samples     = math.floor( len(G) / 4) * 4 - 1
    adc_samples_pre = math.floor( len(G_pre) / 4) * 4 - 2
    adc_samples_pos = math.floor( len(G_pos) / 4) * 4 - 100

    adc = make_adc(num_samples=adc_samples, duration=adc_samples / i_raster_time, system=system)
    adc_pre = make_adc(num_samples=adc_samples_pre, duration=adc_samples_pre / i_raster_time, system=system)
    adc_pos = make_adc(num_samples=adc_samples_pos, duration=adc_samples_pos / i_raster_time, system=system)


    # %% --- 5 - Spoilers
    pre_time = 8e-4
    n_pre_time = math.ceil(pre_time * i_raster_time)

    gx_rephz = make_trapezoid(channel='x', area=-gxRF.area / 2, duration=2e-3, system=system)
    gx_spoil = make_trapezoid(channel='x', area=-gxRF.area / 2, duration=3 * n_pre_time / i_raster_time, system=system)

    gy_rephz = make_trapezoid(channel='y', area=-gyRF.area / 2, duration=2e-3, system=system)
    gy_spoil = make_trapezoid(channel='y', area=-gyRF.area / 2, duration=3 * n_pre_time / i_raster_time, system=system)

    gz_rephz = make_trapezoid(channel='z', area=-gzRF.area / 2, duration=2e-3, system=system)
    gz_spoil = make_trapezoid(channel='z', area=-gzRF.area / 2, duration=3 * n_pre_time / i_raster_time, system=system)

    # double area for all second pulses
    gx_spoil_double = make_trapezoid(channel='x', area=-gxRF.area, duration=3 * n_pre_time / i_raster_time, system=system)
    gy_spoil_double = make_trapezoid(channel='y', area=-gyRF.area, duration=3 * n_pre_time / i_raster_time, system=system)
    gz_spoil_double = make_trapezoid(channel='z', area=-gzRF.area, duration=3 * n_pre_time / i_raster_time, system=system)

    aux_gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:, 0]), system=system)
    aux_gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:, 0]), system=system)
    aux_gz = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:, 0]), system=system)


    # %% --- 6 - Calculate timing/Delays

    n_TE = math.ceil(TE * i_raster_time)
    n_TR = math.ceil(TR * i_raster_time)
    n_dur_adc_pre  = math.ceil(calc_duration(adc_pre) * i_raster_time)
    n_dur_adc_pos  = math.ceil(calc_duration(adc_pos) * i_raster_time)

        # for x
    # delay pre ADC
    n_dur_gx_rephz = math.ceil(calc_duration(gx_rephz) * i_raster_time)
    n_dur_rfx      = math.ceil(calc_duration(rfx) * i_raster_time)

    n_del_preGRS_x = n_TE - (n_dur_gx_rephz + math.ceil(n_dur_rfx / 2))  # Time before input in points
    del_preGRS_x   = n_del_preGRS_x / i_raster_time
    delay_preGRS_x = make_delay(del_preGRS_x)

    # delay pos ADC
    n_dur_aux_gx   = math.ceil(calc_duration(aux_gx) * i_raster_time)
    n_dur_gx_spoil = math.ceil(calc_duration(gx_spoil) * i_raster_time)

    n_delTR_x = n_TR - (n_dur_gx_rephz + math.ceil(n_dur_rfx / 2) + n_dur_adc_pre + n_dur_aux_gx + n_dur_adc_pos + n_dur_gx_spoil)  # Time after input to the system in points
    delTR_x   = n_delTR_x / i_raster_time
    delayTR_x = make_delay(delTR_x)

        # for y
    # delay pre ADC
    n_dur_gy_rephz = math.ceil(calc_duration(gy_rephz) * i_raster_time)
    n_dur_rfy      = math.ceil(calc_duration(rfy) * i_raster_time)

    n_del_preGRS_y = n_TE - (n_dur_gy_rephz + math.ceil(n_dur_rfy / 2))  # Time before input in points
    del_preGRS_y   = n_del_preGRS_y / i_raster_time
    delay_preGRS_y = make_delay(del_preGRS_y)

    # delay pos ADC
    n_dur_aux_gy   = math.ceil(calc_duration(aux_gy) * i_raster_time)
    n_dur_gy_spoil = math.ceil(calc_duration(gy_spoil) * i_raster_time)

    n_delTR_y = n_TR - (n_dur_gy_rephz + math.ceil(n_dur_rfy / 2) + n_dur_adc_pre + n_dur_aux_gy + n_dur_adc_pos + n_dur_gy_spoil)
    delTR_y   = n_delTR_y / i_raster_time
    delayTR_y = make_delay(delTR_y)

        # for z
    # delay pre ADC
    n_dur_gz_rephz = math.ceil(calc_duration(gz_rephz) * i_raster_time)
    n_dur_rfz      = math.ceil(calc_duration(rfz) * i_raster_time)

    n_del_preGRS_z = n_TE - (n_dur_gz_rephz + math.ceil(n_dur_rfz / 2))  # Time before input in points
    del_preGRS_z   = n_del_preGRS_z / i_raster_time
    delay_preGRS_z = make_delay(del_preGRS_z)

    # delay pos ADC
    n_dur_aux_gz   = math.ceil(calc_duration(aux_gz) * i_raster_time)
    n_dur_gz_spoil = math.ceil(calc_duration(gx_spoil) * i_raster_time)

    n_delTR_z = n_TR - (n_dur_gz_rephz + math.ceil(n_dur_rfz / 2) + n_dur_adc_pre + n_dur_aux_gz + n_dur_adc_pos + n_dur_gz_spoil)
    delTR_z   = n_delTR_z / i_raster_time
    delayTR_z = make_delay(delTR_z)


    # %% --- 7 - Define sequence for measuring k_trajectory - blocks/Readouts
    if dummy:
        for nd in range(Ndummies):
            # for x and y rf pulse
            freq_offset_gx = gxRF.amplitude * posAxxis[0]
            rfx.freq_offset = freq_offset_gx
            freq_offset_gy = gyRF.amplitude * posAxxis[0]
            rfy.freq_offset = freq_offset_gy
            # RF90
            seq.add_block(rfx, gxRF)
            seq.add_block(gx_rephz)
            # Delay for spiral
            seq.add_block(delay_preGRS_x)
            # Read k-space - Imaging Gradient waveforms
            gx_zero = make_arbitrary_grad(channel='x', waveform=np.squeeze(G_zero), system=system)
            seq.add_block(adc, gx_zero)
            # Gradients Spoils
            seq.add_block(gx_spoil)
            # Delay to TR
            seq.add_block(delayTR_x)

    n_x = 0
    if tX:
        for ns in range(Npuls):
            for r in range(reps):
                for s in range(2): # for position of axis
                    for a in range(2): # for positive or negative gradient

                        n_x=n_x+1
                        # for x rf pulse
                        freq_offset_gx  = gxRF.amplitude * posAxxis[s]
                        rfx.freq_offset = freq_offset_gx

                        # RF90
                        seq.add_block(rfx, gxRF)
                        seq.add_block(gx_rephz)

                        # Delay for spiral
                        seq.add_block(delay_preGRS_x)


                        # Read - pre GIRF
                        if testPRE_POST:
                            gx_pre = make_arbitrary_grad(channel='x', waveform=np.squeeze(G_pre), system=system)
                            seq.add_block(adc_pre, gx_pre)
                        # Read k-space - Imaging Gradient waveforms
                        if a==1:
                            # %% --- for x-axis negative gradient --- %% #
                            gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:,ns]*-1), system=system)
                        else:
                            # %% --- for x-axis positive gradient --- %% #
                            gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:,ns]), system=system)
                        seq.add_block(adc, gx)
                        # Read - pos GIRF
                        if testPRE_POST:
                            gx_pos = make_arbitrary_grad(channel='x', waveform=np.squeeze(G_pos), system=system)
                            seq.add_block(adc_pos, gx_pos)
                        # Gradients Spoils
                        if (n_x % 2) == 0:
                            # seq.add_block(gx_spoil,gy_spoil,gz_spoil)
                            seq.add_block(gx_spoil_double)
                            # seq.add_block(gy_spoil_double)
                        else:
                            # seq.add_block(gx_spoil,gy_spoil,gz_spoil)
                            seq.add_block(gx_spoil)
                            # seq.add_block(gy_spoil)
                        # Delay to TR
                        seq.add_block(delayTR_x)

    n_y=0
    if tY:
        for ns in range(Npuls):
            for r in range(reps):
                for s in range(2):
                    for a in range(2):
                        n_y=n_y+1
                        # for x rf pulse
                        freq_offset_gy  = gyRF.amplitude * posAxxis[s]
                        rfy.freq_offset = freq_offset_gy

                        # RF90
                        seq.add_block(rfy, gyRF)
                        seq.add_block(gy_rephz)

                        # Delay for spiral
                        seq.add_block(delay_preGRS_y)

                        # Read - pre GIRF
                        if testPRE_POST:
                            gy_pre = make_arbitrary_grad(channel='y', waveform=np.squeeze(G_pre), system=system)
                            seq.add_block(adc_pre, gy_pre)
                        # Read k-space - Imaging Gradient waveforms
                        if a==1:
                            # %% --- for y-axis negative gradient --- %% #
                            gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:,ns]*-1), system=system)
                        else:
                            # %% --- for y-axis positive gradient --- %% #
                            gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:,ns]), system=system)
                        seq.add_block(adc, gy)
                        # Read - pos GIRF
                        if testPRE_POST:
                            gy_pos = make_arbitrary_grad(channel='y', waveform=np.squeeze(G_pos), system=system)
                            seq.add_block(adc_pos, gy_pos)
                        # Gradients Spoils
                        if (n_y % 2) == 0:
                            # seq.add_block(gx_spoil,gy_spoil,gz_spoil)
                            seq.add_block(gy_spoil_double)
                            # seq.add_block(gx_spoil_double)
                        else:
                            # seq.add_block(gx_spoil,gy_spoil,gz_spoil)
                            seq.add_block(gy_spoil)
                            # seq.add_block(gx_spoil)
                        # Delay to TR
                        seq.add_block(delayTR_y)

    if tZ:
        for ns in range(Npuls):
            for r in range(reps):
                for s in range(2):
                    for a in range(2):
                        # for x rf pulse
                        freq_offset_gz  = gzRF.amplitude * posAxxis[s]
                        rfz.freq_offset = freq_offset_gz
                        # RF90
                        seq.add_block(rfz, gz_rephz)
                        # Delay for spiral
                        seq.add_block(delay_preGRS_z)
                        # Read - pre GIRF
                        if testPRE_POST:
                            gz_pre = make_arbitrary_grad(channel='z', waveform=np.squeeze(G_pre), system=system)
                            seq.add_block(adc_pre, gz_pre)
                        # Read k-space - Imaging Gradient waveforms
                        if a==1:
                            # %% --- for y-axis negative gradient --- %% #
                            gz = make_arbitrary_grad(channel='z', waveform=np.squeeze(G[:,ns]*-1), system=system)
                        else:
                            # %% --- for y-axis positive gradient --- %% #
                            gz = make_arbitrary_grad(channel='z', waveform=np.squeeze(G[:,ns]), system=system)
                        seq.add_block(adc, gz)
                        # Read - pos GIRF
                        if testPRE_POST:
                            gz_pos = make_arbitrary_grad(channel='z', waveform=np.squeeze(G_pos), system=system)
                            seq.add_block(adc_pos, gz_pos)
                        # Gradients Spoils
                        seq.add_block(gx_spoil,gy_spoil,gz_spoil)
                        # seq.add_block(gz_spoil)
                        # Delay to TR
                        seq.add_block(delayTR_z)

    # %% --- 8 - Plot ADCs, GX, GY, GZ, RFpulse, RFphase
    # #    plt.figure()
    if tPlot:
        seq.plot()

    if tReport:
        print(seq.test_report())
        seq.check_timing()

    return seq

# seq.write()
