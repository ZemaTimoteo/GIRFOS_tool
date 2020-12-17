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

# from pypulseq.Sequence.sequence import Sequence
# from pypulseq.make_adc import make_adc
# from pypulseq.make_sinc_pulse import make_sinc_pulse
# from pypulseq.make_gauss_pulse import make_gauss_pulse
# from pypulseq.make_trap_pulse import make_trapezoid
# from pypulseq.opts import Opts
# from pypulseq.calc_duration import calc_duration
# from pypulseq.make_delay import make_delay

# Create new Sequence Object
seq = Sequence()
# Due to the floating point uncertainty (https://docs.python.org/3/tutorial/floatingpoint.html),
# sometimes I'm having trouble when '* seq.grad_raster_time',  it returns number like:
# 0.005560000000000001 which give me problems afterwards. However, such problem is solved if I used a manually inputed
# raster time (1/100000). To see this problem do: 556*1e-5 vs 556/1e5. Thus, whenever we have to multiply we divide
# by the equivalent (100000)

# I'm still getting some problems, so the best option to avoid them is to work with integers and only use the actual
# timings by the end. Thus making all intermediate operations (+ and -) with n.

i_raster_time = 100000
# I need to check the manually inputed inverse raster time and the actual value are the same.
assert 1/i_raster_time == seq.grad_raster_time, "Manualy inputed inverse raster time does not match the actual value."

# Code parameters
save_flag          = False                         # Save seq
external_gradients = False                         # External gradients
fatsat_enable      = 0                             # Fat saturation
seqplot            = 0                             # Plot sequence and k-space trajectory
kplot              = 0

# Acquisition Parameters
TR              = 2500e-3                       # Spin-Echo parameters - TR in [s]
n_TR            = math.ceil(TR * i_raster_time) # Spin-Echo parameters - number of points TR
bvalue          = np.array([500])               # b-value [s/mm2]
nbvals          = np.shape(bvalue)[0]           # b-value parameters
ndirs           = 1                             # b-value parameters
Nshots          = 4                             # Spiral parameters - Number of spiral arms
alpha           = 3                             # Spiral parameters - oversampling factor
fov             = 240e-3                        # FOV and resolution - [m]
Nx              = 64                            # FOV and resolution - Number of points.
Ny              = Nx
slice_thickness = 2.5e-3                        # Acquisition Parameters in [m]
n_slices        = 2

# pe_enable = 1
if fatsat_enable:
    fatsat_str = "_fatsat"
else:
    fatsat_str = ""

# Gradient Scaling
gscl = np.zeros(nbvals+1)
gscl[1:] = np.sqrt(bvalue/np.max(bvalue))
gdir, nb0s = difunc.get_dirs(ndirs)

# Set system limits
system = Opts(max_grad=28, grad_unit='mT/m', max_slew=125, slew_unit='T/m/s', rf_ringdown_time=20e-6,
              rf_dead_time=100e-6, adc_dead_time=10e-6)

# Fat saturation
if fatsat_enable:
    b0 = 1.494
    sat_ppm = -3.45
    sat_freq = sat_ppm * 1e-6 * b0 * system.gamma
    rf_fs, _, _ = make_gauss_pulse(flip_angle=110 * math.pi / 180, system=system, duration=8e-3, bandwidth=abs(sat_freq),
                               freq_offset=sat_freq)
    gz_fs = make_trapezoid(channel='z', system=system, delay=calc_duration(rf_fs), area=1 / 1e-4)

# Create 90 degree slice selection pulse and gradient
rf, gz, _ = make_sinc_pulse(flip_angle=math.pi / 2, system=system, duration=3e-3, slice_thickness=slice_thickness,
                                  apodization=0.5, time_bw_product=4)

# Refocusing pulse with spoiling gradients
rf180, gz180, _ = make_sinc_pulse(flip_angle=math.pi, system=system, duration=5e-3, slice_thickness=slice_thickness,
                            apodization=0.5, time_bw_product=4)
rf180.phase_offset = math.pi/2
gz_spoil = make_trapezoid(channel='z', system=system, area=6/slice_thickness, duration=3e-3)

""" K-space reaing with a spiral """
# Spiral trajectory
ktraj, G, _ = vds_2d(fov, Nx, Nshots, alpha, system)

# Define ADC events
delta_k = 1 / fov
adc_samples = math.floor(len(G) / 4) * 4  # Apparently, on Siemens the number of samples needs to be divisible by 4...
# adc = make_adc(num_samples=adc_samples, dwell=system.grad_raster_time, system=system)
adc = make_adc(num_samples=adc_samples, system=system, duration=adc_samples/i_raster_time)

# Remove extra points of the spiral trajectory
# ktraj = ktraj[:adc_samples+1, :]
# G = G[:adc_samples+1, :]

# Pre-phasing gradients
pre_time = 1e-3
n_pre_time = math.ceil(pre_time * i_raster_time)
gz_reph = make_trapezoid(channel='z', system=system, area=-gz.area / 2, duration=pre_time)

"""" Obtain TE and diffusion-weighting gradient waveform """
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

# Find minimum TE considering the readout times.
n_TE = math.ceil(20e-3 / seq.grad_raster_time)
n_delay_te1 = -1
while n_delay_te1 <= 0:
    n_TE = n_TE + 2

    n_tINV = math.floor(n_TE / 2)
    n_delay_te1 = n_tINV - n_rf90r - n_rf180l

# Find minimum TE for the target b-value
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

# Show final TE and b-values:
print("TE:", round(n_TE / i_raster_time * 1e3, 2), "ms")
for bv in range(1, nbvals+1):
    print(round(difunc.calc_bval(system.max_grad * gscl[bv], n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time) * 1e-6, 2), "s/mm2")

# Crusher gradients
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

# Check TR delay time
assert n_tr_delay > 0, "Such parameter configuration needs longer TR."

# In case we would like to use external gradients, we need to obtain them from MATLAB for what we need to know the
# exact parameters. (I print it as the MATLAB code of ODGD - https://doi.org/10.1002/mrm.27462).
# There are two time resolutions printed because usually in MATLAB we work with a grad_raster_time = 0.5e-3 [s]
print("----------------------------------------------------------------")
print("-----------------------MATLAB DATA------------------------------")
print("----------------------------------------------------------------")
print("G_Max:", str(system.max_grad/system.gamma*1e3), "mT/m.")
print("S_Max:", str(system.max_slew/system.gamma), "T/m/s.")
print("GAMMA:", str(system.gamma/1e6), "MHz/T.")
print("T_90:", str(round(n_rf90r / i_raster_time * 1e3, 2)), "ms.")
print("T_RF:", str(round((n_rf180r + n_rf180l) / i_raster_time * 1e3, 2)), "ms.")
print("dt:", str(system.grad_raster_time), "s - Waveform length - n:", str(int(n_TE - n_duration_center)), "points.")
print("dt:", str(50 * system.grad_raster_time), "s - Waveform length - n:", str(int(n_TE - n_duration_center)/50), "points.")
print("TEn", str(int(n_TE)), " points.")
print("b-value:", str(round(bvalue_tmp, 2)), "s/mm2.")
print("----------------------------------------------------------------")


# Delay time
# Time between the gradient and the RF180. This time might be zero some times, although it is not normal.
if n_delay_te1 > n_delay_te2:
    n_gap_te1 = n_delay_te1 - n_delay_te2
    gap_te1 = n_gap_te1 / i_raster_time
    gap_te2 = 0
else:
    n_gap_te2 = n_delay_te2 - n_delay_te1
    gap_te2 = n_gap_te2 / i_raster_time
    gap_te1 = 0

""" b-zero acquisition """
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

""" DWI acquisition """
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


print(seq.test_report())
if seqplot:
    seq.plot()

seq.plot()


if kplot:
    ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc = seq.calculate_kspace()

    time_axis = np.arange(1, ktraj.shape[1] + 1) * system.grad_raster_time
    plt.plot(time_axis, ktraj.T)
    plt.plot(t_adc, ktraj_adc[0, :], '.')
    plt.figure()
    plt.plot(ktraj[0, :], ktraj[1, :], 'b')
    plt.axis('equal')
    plt.plot(ktraj_adc[0, :], ktraj_adc[1, :], 'r.')
    plt.show()

# Save the sequence
if save_flag:
    os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/Diffusion')
    TE = n_TE / i_raster_time
    TR = n_TR / i_raster_time
    seqfname = "se_SPIRAL_dwi_" + str(n_slices) + "slices_" + str(bvalue) + "bvalues_" + str(ndirs) + "dirs_" + str(
        round(TE * 1e3, 2)) + "TE_" + str(round(TR * 1e3)) + "TR_" + str(system.max_grad/system.gamma*1e3) + "G_Max_" +\
               str(system.max_slew/system.gamma) + "SR_Max_" + str(Nshots) + "NShots_" + str(alpha) + "alpha" + fatsat_str
    # os.mkdir("tests/" + seqfname)
    # seq.write("tests/" + seqfname + "/" + seqfname + ".seq")
    seq.write(seqfname)
    # scipy.io.savemat("tests/" +  seqfname + "/sequence_info.mat", mdict={'gradients':G, 'Nshots':Nshots, 'alpha':alpha, \
    #                                                                   'adc_samples':adc_samples, 'fov':fov, 'Nx':Nx, \
    #                                                                   'n_slices':n_slices, 'bvalue':bvalue, \
    #                                                                   'ndirs':ndirs, 'nb0s':nb0s, 'dt':system.grad_raster_time})
    # scipy.io.savemat("/sequence_info.mat", mdict={'gradients':G, 'Nshots':Nshots, 'alpha':alpha, \
    #                                                     'adc_samples':adc_samples, 'fov':fov, 'Nx':Nx, \
    #                                                     'n_slices':n_slices, 'bvalue':bvalue, \
    #                                                     'ndirs':ndirs, 'nb0s':nb0s, 'dt':system.grad_raster_time})
    print("Seq file saved --", seqfname)