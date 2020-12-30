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

sys.path.insert(1, '/Toolboxes/pypulseq')
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

sys.path.insert(1, '/Toolboxes/pypulseq/Sequence')


def generate_GIRF_SeqFile(iT, tr, te, Npuls, mg, ms, grt, fA, n_slices, reps, sliceT, tPre, tPos, tX, tY, tPlot, tReport, rf_offset, testPRE_POST, dummy):
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
    delta_z = n_slices * slice_thickness
    z = np.linspace((-delta_z / 2), (delta_z / 2), n_slices) + rf_offset_z
    x = np.linspace((-delta_z / 2), (delta_z / 2), n_slices) + rf_offset
    y = np.linspace((-delta_z / 2), (delta_z / 2), n_slices) + rf_offset


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

    rf, gz, gzr = make_sinc_pulse(flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                                                apodization=0.5, time_bw_product=4, system=system)
    rfx, gxRF, gxr = make_sinc_pulse_channel(channel='x', flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                                                apodization=0.5, time_bw_product=4, system=system)
    rfy, gyRF, gyr = make_sinc_pulse_channel(channel='y', flip_angle=flip90, duration=3e-3, slice_thickness=slice_thickness,
                                                apodization=0.5, time_bw_product=4, system=system)


    # %% --- 4 - ADCs / Readouts
    # timeStart = -5e-3
    # timeEnd   = 63e-3
    # n_timeStart = math.ceil(timeStart * i_raster_time)
    # n_timeEnd   = math.ceil(timeEnd * i_raster_time)
    # adc_samples = n_timeEnd - n_timeStart
    adc_samples = math.floor( len(G) / 4) * 4 - 1
    adc_samples_pre = math.floor( len(G_pre) / 4) * 4 - 2
    adc_samples_pos = math.floor( len(G_pos) / 4) * 4 - 100

    adc = make_adc(num_samples=adc_samples, duration=adc_samples / i_raster_time, system=system)
    adc_pre = make_adc(num_samples=adc_samples_pre, duration=adc_samples_pre / i_raster_time, system=system)
    adc_pos = make_adc(num_samples=adc_samples_pos, duration=adc_samples_pos / i_raster_time, system=system)


    # %% --- 5 - Spoilers
    pre_time = 8e-4
    n_pre_time = math.ceil(pre_time * i_raster_time)

    gz_rephz = make_trapezoid(channel='z', area=-gz.area / 2, duration=2e-3, system=system)
    gz_spoil = make_trapezoid(channel='z', area=-gz.area / 2, duration=3 * n_pre_time / i_raster_time, system=system)

    gx_rephz = make_trapezoid(channel='x', area=-gxRF.area / 2, duration=2e-3, system=system)
    gx_spoil = make_trapezoid(channel='x', area=-gxRF.area / 2, duration=3 * n_pre_time / i_raster_time, system=system)

    gy_rephz = make_trapezoid(channel='y', area=-gyRF.area / 2, duration=2e-3, system=system)
    gy_spoil = make_trapezoid(channel='y', area=-gyRF.area / 2, duration=3 * n_pre_time / i_raster_time, system=system)

    aux_gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:, 0]), system=system)
    aux_gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:, 0]), system=system)

    # %% --- 6 - Calculate timing/Delays
        # for z
    # delay pre ADC
    n_TE = math.ceil(TE * i_raster_time)
    n_dur_gz_rephz = math.ceil(calc_duration(gz_rephz) * i_raster_time)
    n_dur_rf = math.ceil(calc_duration(rf) * i_raster_time)

    n_del_preGRS = n_TE - (n_dur_gz_rephz + math.ceil(n_dur_rf / 2))  # Time before input in points
    del_preGRS   = n_del_preGRS / i_raster_time
    delay_preGRS = make_delay(del_preGRS)

    # delay pos ADC
    n_TR = math.ceil(TR * i_raster_time)
    n_dur_aux_gx = math.ceil(calc_duration(aux_gx) * i_raster_time)
    n_dur_adc = math.ceil(calc_duration(adc) * i_raster_time)
    n_dur_gz_spoil = math.ceil(calc_duration(gx_spoil) * i_raster_time)

    n_delTR = n_TR - (n_dur_gz_rephz + n_dur_rf + n_dur_aux_gx + n_dur_gz_spoil)   # Time after input to the system in points
    # n_delTR = n_TR - (n_TE + n_dur_adc + n_dur_gz_spoil)  # Time after input to the system in points
    delTR = n_delTR / i_raster_time
    delayTR = make_delay(delTR)

        # for x
    # delay pre ADC
    n_TE = math.ceil(TE * i_raster_time)
    n_dur_gx_rephz = math.ceil(calc_duration(gx_rephz) * i_raster_time)
    n_dur_rfx = math.ceil(calc_duration(rfx) * i_raster_time)

    n_del_preGRS_x = n_TE - (n_dur_gx_rephz + math.ceil(n_dur_rfx / 2))  # Time before input in points
    del_preGRS_x = n_del_preGRS_x / i_raster_time
    delay_preGRS_x = make_delay(del_preGRS_x)

    # delay pos ADC
    n_TR = math.ceil(TR * i_raster_time)
    n_dur_gx_spoil = math.ceil(calc_duration(gx_spoil) * i_raster_time)
    n_dur_adc_pre  = math.ceil(calc_duration(adc_pre) * i_raster_time)
    n_dur_adc_pos  = math.ceil(calc_duration(adc_pos) * i_raster_time)

    n_delTR_x = n_TR - (n_dur_gx_rephz + math.ceil(n_dur_rfx / 2) + n_dur_adc_pre + n_dur_aux_gx + n_dur_adc_pos + n_dur_gx_spoil)  # Time after input to the system in points
    # n_delTR_x = n_TR - (n_TE + n_dur_adc + n_dur_gx_spoil)  # Time after input to the system in points
    delTR_x = n_delTR_x / i_raster_time
    delayTR_x = make_delay(delTR_x)

        # for y
    # delay pre ADC
    n_TE = math.ceil(TE * i_raster_time)
    n_dur_gy_rephz = math.ceil(calc_duration(gy_rephz) * i_raster_time)
    n_dur_rfy = math.ceil(calc_duration(rfy) * i_raster_time)

    n_del_preGRS_y = n_TE - (n_dur_gy_rephz + math.ceil(n_dur_rfy / 2))  # Time before input in points
    del_preGRS_y = n_del_preGRS_y / i_raster_time
    delay_preGRS_y = make_delay(del_preGRS_y)

    # delay pos ADC
    n_TR = math.ceil(TR * i_raster_time)
    n_dur_aux_gy = math.ceil(calc_duration(aux_gy) * i_raster_time)
    n_dur_gy_spoil = math.ceil(calc_duration(gy_spoil) * i_raster_time)

    n_delTR_y = n_TR - (n_dur_gy_rephz + math.ceil(n_dur_rfy / 2) + n_dur_adc_pre + n_dur_aux_gy + n_dur_adc_pos + n_dur_gy_spoil)
    # n_delTR_y = n_TR - (n_TE + n_dur_adc + n_dur_gy_spoil)  # Time after input to the system in points
    delTR_y = n_delTR_y / i_raster_time
    delayTR_y = make_delay(delTR_y)


    # %% --- 7 - Define sequence for measuring k_trajectory - blocks/Readouts
    if dummy:
        # for x and y rf pulse
        freq_offset_gx = gxRF.amplitude * x[0]
        rfx.freq_offset = freq_offset_gx
        freq_offset_gy = gyRF.amplitude * y[0]
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

    if tX:
        for ns in range(Npuls):
            for r in range(reps):
                for s in range(n_slices):

                    # for x and y rf pulse
                    freq_offset_gx  = gxRF.amplitude * x[s]
                    rfx.freq_offset = freq_offset_gx
                    freq_offset_gy  = gyRF.amplitude * y[s]
                    rfy.freq_offset = freq_offset_gy

                # %% --- for x-axis no gradient --- %% #
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
                    gx_zero = make_arbitrary_grad(channel='x', waveform=np.squeeze(G_zero), system=system)
                    seq.add_block(adc, gx_zero)
                    # Read - pos GIRF
                    if testPRE_POST:
                        gx_pos = make_arbitrary_grad(channel='x', waveform=np.squeeze(G_pos), system=system)
                        seq.add_block(adc_pos, gx_pos)
                    # Gradients Spoils
                    seq.add_block(gx_spoil)
                    # Delay to TR
                    seq.add_block(delayTR_x)

                # %% --- for x-axis w/ gradient --- %% #
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
                    gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:,ns]), system=system)
                    seq.add_block(gx, adc)
                    # Read - pos GIRF
                    if testPRE_POST:
                        gx_pos = make_arbitrary_grad(channel='x', waveform=np.squeeze(G_pos), system=system)
                        seq.add_block(adc_pos, gx_pos)
                    # Gradients Spoils
                    seq.add_block(gx_spoil)
                    # Delay to TR
                    seq.add_block(delayTR_x)

    if tY:
        for ns in range(Npuls):
            for r in range(reps):
                for s in range(n_slices):
                    # %% --- for y-axis no gradient --- %% #
                    seq.add_block(rfy, gyRF)
                    seq.add_block(gy_rephz)
                    # Delay for spiral
                    seq.add_block(delay_preGRS_y)
                    # Read - pre GIRF
                    if testPRE_POST:
                        gy_pre = make_arbitrary_grad(channel='y', waveform=np.squeeze(G_pre), system=system)
                        seq.add_block(adc_pre, gy_pre)
                    # Read k-space - Imaging Gradient waveforms
                    gy_zero = make_arbitrary_grad(channel='y', waveform=np.squeeze(G_zero), system=system)
                    seq.add_block(gy_zero, adc)
                    # Read - pos GIRF
                    if testPRE_POST:
                        gy_pos = make_arbitrary_grad(channel='y', waveform=np.squeeze(G_pos), system=system)
                        seq.add_block(adc_pos, gy_pos)
                    # Gradients Spoils
                    seq.add_block(gy_spoil)
                    # Delay to TR
                    seq.add_block(delayTR_y)

                # %% --- for y-axis w/ gradient --- %% #
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
                    gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:,ns]), system=system)
                    seq.add_block(gy, adc)
                    # Read - pos GIRF
                    if testPRE_POST:
                        gy_pos = make_arbitrary_grad(channel='y', waveform=np.squeeze(G_pos), system=system)
                        seq.add_block(adc_pos, gy_pos)
                    # Gradients Spoils
                    seq.add_block(gy_spoil)
                    # Delay to TR
                    seq.add_block(delayTR_y)



    # %% --- 8 - Plot ADCs, GX, GY, GZ, RFpulse, RFphase
    # #    plt.figure()
    if tPlot:
        seq.plot()

    if tReport:
        print(seq.test_report())
        seq.check_timing()

    return seq

# seq.write()
