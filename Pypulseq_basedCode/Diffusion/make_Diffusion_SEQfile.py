"""
Tiago T. Fernandes
October 2020
LASEEB, Instituto Superior Técnico, Lisbon, Portugal

This function builds a spin echo diffusion weighting sequence with spiral trajectories
optimizing the TE for a desired b-value.
"""
from __future__ import division
import math
import sys
import numpy as np
import scipy.io
import os
import matplotlib.pyplot as plt
from math import pi
import matplotlib.pyplot as plt
import os
import sys

import diff_funcs as difunc

# I do this because apparently the pip pypulseq was crashing on some functions (make_trap_pulse)
sys.path.append('/home/tfernandes/Documents/PYTHON/Toolboxes')
sys.path.insert(1, '/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq')

import make_adc
from make_adc import make_adc
import make_sinc_pulse
from make_sinc_pulse import make_sinc_pulse
import make_gauss_pulse
from make_gauss_pulse import make_gauss_pulse
import make_trap_pulse
from make_trap_pulse import make_trapezoid
import make_extended_trapezoid
from make_extended_trapezoid import make_extended_trapezoid
import make_arbitrary_grad
from make_arbitrary_grad import make_arbitrary_grad
import opts
from opts import Opts
import calc_duration
from calc_duration import calc_duration
import calc_rf_center
from calc_rf_center import calc_rf_center
import make_delay
from make_delay import make_delay
from Sequence.sequence import Sequence
from vds_2d import vds_2d

sys.path.insert(1, '/home/tfernandes/Documents/PYTHON/Toolboxes/pypulseq/Sequence')
import sequence
from sequence import Sequence

import split_gradient_at
from split_gradient_at import split_gradient_at
import align
from align import align
import add_gradients
from add_gradients import add_gradients


def generate_SeqFile_SpiralDiffusion(gx, gy, tr, n_shots, mg, ms, fA, n_slices, reps, st, tPlot, tReport,b_values, n_dirs, fov, Nx):


    #%% --- 1 - Create new Sequence Object + Parameters
    seq = Sequence()

    # =========
    # Parameters
    # =========
    i_raster_time = 100000
    assert 1/i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."

    # =========
    # Code parameters
    # =========
    fatsat_enable   = 0                             # Fat saturation
    kplot           = 0

    # =========
    # Acquisition Parameters
    # =========
    TR              = tr                            # Spin-Echo parameters - TR in [s]
    n_TR            = math.ceil(TR * i_raster_time) # Spin-Echo parameters - number of points TR
    bvalue          = b_values                      # b-value [s/mm2]
    nbvals          = np.shape(bvalue)[0]           # b-value parameters
    ndirs           = n_dirs                        # b-value parameters
    Ny              = Nx
    slice_thickness = st                            # Acquisition Parameters in [m]
    Nshots          = n_shots

    # =========
    # Gradient Scaling
    # =========
    gscl        = np.zeros(nbvals+1)
    gscl[1:]    = np.sqrt(bvalue/np.max(bvalue))
    gdir, nb0s  = difunc.get_dirs(ndirs)

    # =========
    # Create system
    # =========
    system = Opts(max_grad=mg, grad_unit='mT/m', max_slew=ms, slew_unit='T/m/s', rf_ringdown_time=20e-6,
                  rf_dead_time=100e-6, adc_dead_time=10e-6)

    #%% --- 2 - Fat saturation
    if fatsat_enable:
        fatsat_str = "_fatsat"
        b0 = 1.494
        sat_ppm = -3.45
        sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
        rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3, bandwidth=abs(sat_freq),
                                   freq_offset=sat_freq)
        gz_fs = make_trapezoid(channel='z', system=system, delay=calc_duration(rf_fs), area=1 / 1e-4)
    else:
        fatsat_str = ""

    #%% --- 3 - Slice Selection
    # =========
    # Create 90 degree slice selection pulse and gradient
    # =========
    flip90 = fA * pi / 180
    rf, gz, _ = make_sinc_pulse(flip_angle=flip90, system=system, duration=3e-3, slice_thickness=slice_thickness,
                                      apodization=0.5, time_bw_product=4)

    # =========
    # Refocusing pulse with spoiling gradients
    # =========
    rf180, gz180, _ = make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                                apodization=0.5, time_bw_product=4)
    rf180.phase_offset = math.pi/2
    gz_spoil = make_trapezoid(channel='z', system=system, area=6/slice_thickness, duration=3e-3)

    #%% --- 4 - Gradients
    # =========
    # Spiral trajectory
    # =========
    G = gx + 1J * gy

    #%% --- 5 - ADCs / Readouts
    delta_k = 1 / fov
    adc_samples = math.floor(len(G) / 4) * 4 -2 # Apparently, on Siemens the number of samples needs to be divisible by 4...
    adc = make_adc(num_samples=adc_samples, system=system, duration=adc_samples/i_raster_time)

    # =========
    # Pre-phasing gradients
    # =========
    pre_time   = 1e-3
    n_pre_time = math.ceil(pre_time * i_raster_time)
    gz_reph    = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)

    #%% --- 6 - Obtain TE and diffusion-weighting gradient waveform
    # For S&T monopolar waveforms
    # From an initial TE, check we satisfy all constraints -> otherwise increase TE.
    # Once all constraints are okay -> check b-value, if it is lower than the target one -> increase TE
    # Looks time-inefficient but it is fast enough to make it user-friendly.
    # TODO: Re-scale the waveform to the exact b-value because increasing the TE might produce slightly higher ones.

    # Calculate some times constant throughout the process
    # We need to compute the exact time sequence. For the normal SE-MONO-EPI sequence micro second differences
    # are not important, however, if we wanna import external gradients the allocated time for them needs to
    # be the same, and thus exact timing is mandatory. With this in mind, we establish the following rounding rules:
    # Duration of RFs + spoiling, and EPI time to the center of the k-space is always math.ceil().

    # The time(gy) refers to the number of blips, thus we substract 0.5 since the number of lines is always even.
    # The time(gx) refers to the time needed to read each line of the k-space. Thus, if Ny is even, it would take half of the lines plus another half.
    n_duration_center = 0  # The spiral starts right in 0 -- or ADC_dead_time??
    rf_center_with_delay = rf.delay + calc_rf_center(rf)[0]

    n_rf90r = math.ceil((calc_duration(gz) - rf_center_with_delay + pre_time) / seq.grad_raster_time)
    n_rf180r = math.ceil((calc_duration(rf180) + 2 * calc_duration(gz_spoil)) / 2 / seq.grad_raster_time)
    n_rf180l = math.floor((calc_duration(rf180) + 2 * calc_duration(gz_spoil)) / 2 / seq.grad_raster_time)

    # =========
    # Find minimum TE considering the readout times.
    # =========
    n_TE = math.ceil(20e-3 / seq.grad_raster_time)
    n_delay_te1 = -1
    while n_delay_te1 <= 0:
        n_TE = n_TE + 2

        n_tINV = math.floor(n_TE / 2)
        n_delay_te1 = n_tINV - n_rf90r - n_rf180l

    # =========
    # Find minimum TE for the target b-value
    # =========
    bvalue_tmp = 0
    while bvalue_tmp < np.max(bvalue):
        n_TE = n_TE + 2

        n_tINV = math.floor(n_TE / 2)
        n_delay_te1 = n_tINV - n_rf90r - n_rf180l
        delay_te1 = n_delay_te1 / i_raster_time
        n_delay_te2 = n_tINV - n_rf180r - n_duration_center
        delay_te2 = n_delay_te2 / i_raster_time

        # Waveform Ramp time
        n_gdiff_rt = math.ceil(system.max_grad / system.max_slew / seq.grad_raster_time)

        # Select the shortest available time
        n_gdiff_delta = min(n_delay_te1, n_delay_te2)
        n_gdiff_Delta = n_delay_te1 + 2 * math.ceil(calc_duration(gz_spoil) / seq.grad_raster_time) + math.ceil(calc_duration(gz180) / seq.grad_raster_time)

        gdiff = make_trapezoid(channel='x', system=system, amplitude=system.max_grad, duration=n_gdiff_delta / i_raster_time)

        # delta only corresponds to the rectangle.
        n_gdiff_delta = n_gdiff_delta - 2 * n_gdiff_rt

        bv = difunc.calc_bval(system.max_grad, n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time)
        bvalue_tmp = bv * 1e-6

    # =========
    # Show final TE and b-values:
    # =========
    print("TE:", round(n_TE / i_raster_time * 1e3, 2), "ms")
    for bv in range(1, nbvals+1):
        print(round(difunc.calc_bval(system.max_grad * gscl[bv], n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time) * 1e-6, 2), "s/mm2")

    TE = n_TE / i_raster_time
    TR = n_TR / i_raster_time

    #%% --- 7 - Crusher gradients
    gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
    gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

    # TR delay - Takes everything into account
    # Distance between the center of the RF90s must be TR
    # The n_pre_time here is the time used to drive the Gx, and Gy spiral gradients to zero.
    n_spiral_time = adc_samples
    n_tr_per_slice = math.ceil(TR / n_slices * i_raster_time)
    if fatsat_enable:
        n_tr_delay = n_tr_per_slice - (n_TE - n_duration_center + n_spiral_time) \
                            - math.ceil(rf_center_with_delay * i_raster_time) \
                            - n_pre_time \
                            - math.ceil(calc_duration(gx_crush, gz_crush) * i_raster_time) \
                            - math.ceil(calc_duration(rf_fs, gz_fs) * i_raster_time)
    else:
        n_tr_delay = n_tr_per_slice - (n_TE - n_duration_center + n_spiral_time) \
                        - math.ceil(rf_center_with_delay * i_raster_time) \
                        - n_pre_time \
                        - math.ceil(calc_duration(gx_crush, gz_crush) * i_raster_time)
    tr_delay = n_tr_delay / i_raster_time

    #%% --- 8 - Checks
    # =========
    # Check TR delay time
    # =========
    assert n_tr_delay > 0, "Such parameter configuration needs longer TR."

    # =========
    # Delay time,
    # =========
    # Time between the gradient and the RF180. This time might be zero some times, although it is not normal.
    if n_delay_te1 > n_delay_te2:
        n_gap_te1 = n_delay_te1 - n_delay_te2
        gap_te1 = n_gap_te1 / i_raster_time
        gap_te2 = 0
    else:
        n_gap_te2 = n_delay_te2 - n_delay_te1
        gap_te2 = n_gap_te2 / i_raster_time
        gap_te1 = 0

    #%% --- 9 - b-zero acquisition
    for r in range(reps):
        for d in range(nb0s):
            for nshot in range(Nshots):
                for s in range(n_slices):
                    # Fat saturation
                    if fatsat_enable:
                        seq.add_block(rf_fs, gz_fs)

                    # RF90
                    rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
                    seq.add_block(rf, gz)
                    seq.add_block(gz_reph)

                    # Delay for RF180
                    seq.add_block(make_delay(delay_te1))

                    # RF180
                    seq.add_block(gz_spoil)
                    rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
                    seq.add_block(rf180, gz180)
                    seq.add_block(gz_spoil)

                    # Delay for spiral
                    seq.add_block(make_delay(delay_te2))

                    # Read k-space
                    # Imaging Gradient waveforms
                    gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:, nshot].real), system=system)
                    gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:, nshot].imag), system=system)
                    seq.add_block(gx, gy, adc)

                    # Make the spiral finish in zero - I use pre_time because I know for sure it's long enough.
                    # Furthermore, this is after readout and TR is supposed to be long.
                    amp_x = [G[:, nshot].real[-1], 0]
                    amp_y = [G[:, nshot].imag[-1], 0]
                    gx_to_zero = make_extended_trapezoid(channel='x', amplitudes=amp_x, times=[0, pre_time], system=system)
                    gy_to_zero = make_extended_trapezoid(channel='y', amplitudes=amp_y, times=[0, pre_time], system=system)
                    seq.add_block(gx_to_zero, gy_to_zero)

                    seq.add_block(gx_crush, gz_crush)

                    # Wait TR
                    if tr_delay > 0:
                        seq.add_block(make_delay(tr_delay))

    #%% --- 9 - DWI acquisition
    for r in range(reps):
        for bv in range(1, nbvals+1):
            for d in range(ndirs):
                for nshot in range(Nshots):
                    for s in range(n_slices):
                        # Fat saturation
                        if fatsat_enable:
                            seq.add_block(rf_fs, gz_fs)

                        # RF90
                        rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
                        seq.add_block(rf, gz)
                        seq.add_block(gz_reph)

                        # Diffusion-weighting gradient
                        gdiffx = make_trapezoid(channel='x', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 0],
                                                duration=calc_duration(gdiff))
                        gdiffy = make_trapezoid(channel='y', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 1],
                                                duration=calc_duration(gdiff))
                        gdiffz = make_trapezoid(channel='z', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 2],
                                                duration=calc_duration(gdiff))

                        seq.add_block(gdiffx, gdiffy, gdiffz)

                        # Delay for RF180
                        seq.add_block(make_delay(gap_te1))

                        # RF180
                        seq.add_block(gz_spoil)
                        rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
                        seq.add_block(rf180, gz180)
                        seq.add_block(gz_spoil)

                        # Diffusion-weighting gradient
                        seq.add_block(gdiffx, gdiffy, gdiffz)

                        # Delay for spiral
                        seq.add_block(make_delay(gap_te2))

                        # Read k-space
                        # Imaging Gradient waveforms
                        gx = make_arbitrary_grad(channel='x', waveform=np.squeeze(G[:, nshot].real), system=system)
                        gy = make_arbitrary_grad(channel='y', waveform=np.squeeze(G[:, nshot].imag), system=system)
                        seq.add_block(gx, gy, adc)

                        # Make the spiral finish in zero - I use pre_time because I know for sure it's long enough.
                        # Furthermore, this is after readout and TR is supposed to be long.
                        amp_x = [G[:, nshot].real[-1], 0]
                        amp_y = [G[:, nshot].imag[-1], 0]
                        gx_to_zero = make_extended_trapezoid(channel='x', amplitudes=amp_x, times=[0, pre_time], system=system)
                        gy_to_zero = make_extended_trapezoid(channel='y', amplitudes=amp_y, times=[0, pre_time], system=system)
                        seq.add_block(gx_to_zero, gy_to_zero)

                        seq.add_block(gx_crush, gz_crush)

                        # Wait TR
                        if tr_delay > 0:
                            seq.add_block(make_delay(tr_delay))


    if tPlot:
        seq.plot()

    if tReport:
        print(seq.test_report())
        seq.check_timing()



    return seq, TE, TR, fatsat_str



def generate_SeqFile_EPIDiffusion(FOV, nx, ny, ns, mg, ms, reps, st, tr, fA, b_values, n_dirs, partialFourier, tPlot, tReport):


    #%% --- 1 - Create new Sequence Object + Parameters
    seq = Sequence()

    # =========
    # Parameters
    # =========
    i_raster_time = 100000
    assert 1/i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."

    # =========
    # Code parameters
    # =========
    fatsat_enable   = 0                             # Fat saturation
    kplot           = 0

    # =========
    # Acquisition Parameters
    # =========
    fov             = FOV
    Nx              = nx
    Ny              = ny
    n_slices        = ns
    TR              = tr                                # Spin-Echo parameters - TR in [s]
    n_TR            = math.ceil(TR * i_raster_time)     # Spin-Echo parameters - number of points TR
    bvalue          = b_values                          # b-value [s/mm2]
    nbvals          = np.shape(bvalue)[0]               # b-value parameters
    ndirs           = n_dirs                            # b-value parameters
    slice_thickness = st                                # Acquisition Parameters in [m]

    # =========
    # Partial Fourier
    # =========
    pF = partialFourier
    Nyeff = int(pF * Ny)  # Number of Ny samples acquired
    if pF is not 1:
        pF_str = "_" + str(pF) + "pF"
    else:
        pF_str = ""

    # =========
    # Gradient Scaling
    # =========
    gscl        = np.zeros(nbvals+1)
    gscl[1:]    = np.sqrt(bvalue/np.max(bvalue))
    gdir, nb0s  = difunc.get_dirs(ndirs)

    # =========
    # Create system
    # =========
    system = Opts(max_grad=mg, grad_unit='mT/m', max_slew=ms, slew_unit='T/m/s', rf_ringdown_time=20e-6,
                  rf_dead_time=100e-6, adc_dead_time=10e-6)

    #%% --- 2 - Fat saturation
    if fatsat_enable:
        fatsat_str = "_fatsat"
        b0 = 1.494
        sat_ppm = -3.45
        sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
        rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3, bandwidth=abs(sat_freq),
                                   freq_offset=sat_freq)
        gz_fs = make_trapezoid(channel='z', system=system, delay=calc_duration(rf_fs), area=1 / 1e-4)
    else:
        fatsat_str = ""

    #%% --- 3 - Slice Selection
    # =========
    # Create 90 degree slice selection pulse and gradient
    # =========
    flip90 = fA * pi / 180
    rf, gz, _ = make_sinc_pulse(flip_angle=flip90, system=system, duration=3e-3, slice_thickness=slice_thickness,
                                      apodization=0.5, time_bw_product=4)

    # =========
    # Refocusing pulse with spoiling gradients
    # =========
    rf180, gz180, _ = make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                                apodization=0.5, time_bw_product=4)
    rf180.phase_offset = math.pi/2
    gz_spoil = make_trapezoid(channel='z', system=system, area=6/slice_thickness, duration=3e-3)

    #%% --- 4 - Define other gradients and ADC events
    delta_k = 1 / fov
    k_width = Nx * delta_k
    dwell_time = seq.grad_raster_time  # Full receiver bandwith
    readout_time = Nx * dwell_time  # T_acq (acquisition time)
    flat_time = math.ceil(readout_time / seq.grad_raster_time) * seq.grad_raster_time
    gx = make_trapezoid(channel='x', system=system, amplitude=k_width / readout_time, flat_time=flat_time)
    adc = make_adc(num_samples=Nx, duration=readout_time,
                   delay=gx.rise_time + flat_time / 2 - (readout_time - dwell_time) / 2)

    # =========
    # Pre-phasing gradients
    # =========
    pre_time = 1e-3
    gx_pre = make_trapezoid(channel='x', system=system, area=-gx.area / 2, duration=pre_time)
    gz_reph = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)
    gy_pre = make_trapezoid(channel='y', system=system, area=-(Ny / 2 - 0.5 - (Ny - Nyeff)) * delta_k,
                            duration=pre_time)  # Es -0.5 y no +0.5 porque hay que pensar en areas, no en rayas!

    # =========
    # Phase blip in shortest possible time
    # =========
    gy = make_trapezoid(channel='y', system=system, area=delta_k)
    dur = math.ceil(calc_duration(gy) / seq.grad_raster_time) * seq.grad_raster_time

    #%% --- 5 - Obtain TE and diffusion-weighting gradient waveform
    # =========
    # Calculate some times constant throughout the process
    # =========
    duration_center = math.ceil((calc_duration(gx) * (Ny / 2 + 0.5 - (Ny - Nyeff)) + calc_duration(gy) * (
                Ny / 2 - 0.5 - (Ny - Nyeff))) / seq.grad_raster_time) * seq.grad_raster_time
    rf_center_with_delay = rf.delay + calc_rf_center(rf)[0]
    rf180_center_with_delay = rf180.delay + calc_rf_center(rf180)[0]

    # =========
    # Find minimum TE considering the readout times.
    # =========
    TE = 40e-3  # [s]
    delay_te2 = -1
    while delay_te2 <= 0:
        TE = TE + 0.02e-3  # [ms]
        delay_te2 = math.ceil((TE / 2 - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil) - \
                               calc_duration(gx_pre,
                                             gy_pre) - duration_center) / seq.grad_raster_time) * seq.grad_raster_time

    # =========
    # Find minimum TE for the target b-value
    # =========
    bvalue_tmp = 0
    while bvalue_tmp < np.max(bvalue):
        TE = TE + 2 * seq.grad_raster_time  # [ms]
        delay_te1 = math.ceil((TE / 2 - calc_duration(gz) + rf_center_with_delay - pre_time - calc_duration(gz_spoil) - \
                               rf180_center_with_delay) / seq.grad_raster_time) * seq.grad_raster_time
        delay_te2 = math.ceil((TE / 2 - calc_duration(rf180) + rf180_center_with_delay - calc_duration(gz_spoil) - \
                               calc_duration(gx_pre,
                                             gy_pre) - duration_center) / seq.grad_raster_time) * seq.grad_raster_time

        # Waveform Ramp time
        gdiff_rt = math.ceil(system.max_grad / system.max_slew / seq.grad_raster_time) * seq.grad_raster_time

        # Select the shortest available time
        gdiff_delta = min(delay_te1, delay_te2)
        gdiff_Delta = math.ceil(
            (delay_te1 + 2 * calc_duration(gz_spoil) + calc_duration(
                gz180)) / seq.grad_raster_time) * seq.grad_raster_time

        gdiff = make_trapezoid(channel='x', system=system, amplitude=system.max_grad, duration=gdiff_delta)

        # delta only corresponds to the rectangle.
        gdiff_delta = math.ceil((gdiff_delta - 2 * gdiff_rt) / seq.grad_raster_time) * seq.grad_raster_time

        bv = difunc.calc_bval(system.max_grad, gdiff_delta, gdiff_Delta, gdiff_rt)
        bvalue_tmp = bv * 1e-6

    # =========
    # Show final TE and b-values:
    # =========
    print("TE:", round(TE * 1e3, 2), "ms")
    for bv in range(1, nbvals + 1):
        print(round(difunc.calc_bval(system.max_grad * gscl[bv], gdiff_delta, gdiff_Delta, gdiff_rt) * 1e-6, 2),
              "s/mm2")

    # =========
    # Crusher gradients
    # =========
    gx_crush = make_trapezoid(channel='x', area=2 * Nx * delta_k, system=system)
    gz_crush = make_trapezoid(channel='z', area=4 / slice_thickness, system=system)

    #%% --- 6 - Delays
    # =========
    # TR delay - Takes everything into account
    # EPI reading time:
    # Distance between the center of the RF90s must be TR
    # =========
    EPI_time = calc_duration(gx) * Nyeff + calc_duration(gy) * (Nyeff - 1)
    if fatsat_enable:
        tr_delay = math.floor(
            (TR - (TE - duration_center + EPI_time) - rf_center_with_delay - calc_duration(gx_crush, gz_crush) \
             - calc_duration(rf_fs, gz_fs)) \
            / seq.grad_raster_time) * seq.grad_raster_time
    else:
        tr_delay = math.floor(
            (TR - (TE - duration_center + EPI_time) - rf_center_with_delay - calc_duration(gx_crush, gz_crush)) \
            / seq.grad_raster_time) * seq.grad_raster_time

    # =========
    # Check TR delay time
    # =========
    assert tr_delay > 0, "Such parameter configuration needs longer TR."

    # =========
    # Delay time
    # =========

    # =========
    # Time between the gradient and the RF180. This time might be zero some times, although it is not normal.
    # =========
    gap_te1 = math.ceil((delay_te1 - calc_duration(gdiff)) / seq.grad_raster_time) * seq.grad_raster_time

    # =========
    # Time between the gradient and the locate k-space gradients.
    # =========
    gap_te2 = math.ceil((delay_te2 - calc_duration(gdiff)) / seq.grad_raster_time) * seq.grad_raster_time

    #%% --- 9 - b-zero acquisition
    for d in range(nb0s):
        for s in range(n_slices):
            # Fat saturation
            if fatsat_enable:
                seq.add_block(rf_fs, gz_fs)

            # RF90
            rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf, gz)
            seq.add_block(gz_reph)

            # Delay for RF180
            seq.add_block(make_delay(delay_te1))

            # RF180
            seq.add_block(gz_spoil)
            rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
            seq.add_block(rf180, gz180)
            seq.add_block(gz_spoil)

            # Delay for EPI
            seq.add_block(make_delay(delay_te2))

            # Locate k-space
            seq.add_block(gx_pre, gy_pre)

            for i in range(Nyeff):
                seq.add_block(gx, adc)  # Read one line of k-space
                if i is not Nyeff - 1:
                    seq.add_block(gy)  # Phase blip
                gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

            seq.add_block(gx_crush, gz_crush)

            # Wait TR
            if tr_delay > 0:
                seq.add_block(make_delay(tr_delay))

    #%% --- 10 - DWI acquisition
    for bv in range(1, nbvals + 1):
        for d in range(ndirs):
            for s in range(n_slices):
                # Fat saturation
                if fatsat_enable:
                    seq.add_block(rf_fs, gz_fs)

                # RF90
                rf.freq_offset = gz.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
                seq.add_block(rf, gz)
                seq.add_block(gz_reph)

                # Diffusion-weighting gradient
                gdiffx = make_trapezoid(channel='x', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 0],
                                        duration=calc_duration(gdiff))
                gdiffy = make_trapezoid(channel='y', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 1],
                                        duration=calc_duration(gdiff))
                gdiffz = make_trapezoid(channel='z', system=system, amplitude=system.max_grad * gscl[bv] * gdir[d, 2],
                                        duration=calc_duration(gdiff))

                seq.add_block(gdiffx, gdiffy, gdiffz)

                # Delay for RF180
                seq.add_block(make_delay(gap_te1))

                # RF180
                seq.add_block(gz_spoil)
                rf180.freq_offset = gz180.amplitude * slice_thickness * (s - (n_slices - 1) / 2)
                seq.add_block(rf180, gz180)
                seq.add_block(gz_spoil)

                # Diffusion-weighting gradient
                seq.add_block(gdiffx, gdiffy, gdiffz)

                # Delay for EPI
                seq.add_block(make_delay(gap_te2))

                # Locate k-space
                seq.add_block(gx_pre, gy_pre)

                for i in range(Nyeff):
                    seq.add_block(gx, adc)  # Read one line of k-space
                    if i is not Nyeff - 1:
                        seq.add_block(gy)  # Phase blip
                    gx.amplitude = -gx.amplitude  # Reverse polarity of read gradient

                seq.add_block(gx_crush, gz_crush)

                # Wait TR
                if tr_delay > 0:
                    seq.add_block(make_delay(tr_delay))


    if tPlot:
        seq.plot()

    if tReport:
        print(seq.test_report())
        seq.check_timing()



    return seq, TE, TR, fatsat_str