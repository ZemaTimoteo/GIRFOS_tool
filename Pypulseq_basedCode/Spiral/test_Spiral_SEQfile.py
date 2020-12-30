#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Name:
    - Test_Spiral_SEQfile

Description:
    - This code creates spiral (with external gradients) '.seq' file

@author: tfernandes, IST 2020
"""

# =============================================================================

#%% --- 0 - Import functions
#import make_t1_mprage as mtp
import os
import scipy.io
import matplotlib
import numpy as np
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import sys

sys.path.append('/Toolboxes')
sys.path.append('/Pypulseq_basedCode/Diffusion')
# os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/Spiral/pulses')

from make_Spiral_SEQfile import generate_SeqFile_Spiral
from make_Diffusion_SEQfile import generate_SeqFile_SpiralDiffusion
from scipy.io import savemat
from pypulseq.calc_duration import calc_duration



#%% --- 0.2 - Settings
test       = 85   #test for scanner
testSPI    = 4    #spiral used in this test

save_flag  = False  # save file
plotTest   = True   # Plot sequence
reportTest = False  # Report the sequence
diffusion  = True   # Sequence w/ diffusion - 'True' OR Sequence w/out diffusion - 'False'
getGRADs   = True   # Get gradients for Matlab analyis - 'True' OR - 'False'
seq2jem    = False   # convert '.seq' to '.xml'

#%% --- 0.5 - Import data

    # --- load data
TR      = 2000e-3 # [s]

if diffusion:
    ndirs   = 2                # b-value parameters
    bvalue  = np.array([500])  # b-value [s/mm2]

mat = scipy.io.loadmat('vectSpiral_test{}.mat'.format(testSPI))
# mat = scipy.io.loadmat('vectSpiral.mat')
aux_vectSpiral = mat['vectSpiral']

vectSpiral = aux_vectSpiral[0, 0][0]
nshots     = 1         #aux_vectSpiral[0, 0][1][0,0]
nslices    = 1
nreps      = 4
alpha      = 2
N_x        = 133       # FOV and resolution - Number of points.
FOV        = 200e-3    # FOV and resolution - [m]


    # --- time & gradients
t   = vectSpiral[:,2] # in s
gx  = vectSpiral[:,0] # in T/m
gy  = vectSpiral[:,1] # in T/m
ADC = len(t)

    # --- reshape data
t  = np.reshape(t,(int(len(t)/nshots),nshots),'F')
gx = np.reshape(gx,(int(len(gx)/nshots),nshots),'F')
gy = np.reshape(gy,(int(len(gy)/nshots),nshots),'F')
#TR = round(t[-1,0])
# TR = 1
DT = t[1,0] - t[0,0]

    # --- select index for ADCs
ind = np.nonzero(gx[:,0])[0]
ind = np.concatenate((ind[0]-1,ind), axis=None)

gx = gx[ind,:]
gy = gy[ind,:]
t  = t[ind,:]

del aux_vectSpiral, mat, vectSpiral, ind

#%% --- 0.75 - Plot loaded data

# plt.figure('F')
# plt.subplot(211)
# plt.plot(np.reshape(gx,(len(gx)*nshots,1),'F'))
# plt.title('Gx')
# plt.grid(True)
#
# plt.subplot(212)
# plt.plot(np.reshape(gy,(len(gy)*nshots,1),'F'))
# plt.title('Gy')
# plt.grid(True)
# plt.show()


#%% --- 1 - Run function to create struct 'seq' for Spiral
# parameters for spiral sequence
max_grad, max_slew, grdrastTime, flipAngle, sliceT, RFoffset = 28, 125, DT, 90, 3e-3, 0

if diffusion:
    seq, TE, TR, fatsat_str = generate_SeqFile_SpiralDiffusion(gx, gy, tr=TR, n_shots=nshots, mg=max_grad, ms=max_slew, fA=flipAngle, n_slices=nslices, reps=nreps, st=sliceT, tPlot=plotTest,
                                           tReport=reportTest, b_values=bvalue, n_dirs=ndirs, fov=FOV, Nx=N_x)
    diff = "_Diff"
else:
    seq = generate_SeqFile_Spiral(gx, gy, tr=TR, n_shots=nshots, mg=max_grad, ms=max_slew, fA=flipAngle, n_slices=nslices, reps=nreps, st=sliceT, tPlot=plotTest, tReport=reportTest, rf_offset=RFoffset)
    fatsat_str, diff, TE, ndirs = "", "", 0, 0
    bvalue = np.array([0])

#%% --- 2 - Get Gradients in vector
if getGRADs:
    # xx-axis
    fig_xx = plt.figure("xx")
    fig_xx_sp_list = [fig_xx.add_subplot(111)]  # , fig1.add_subplot(212)]

    t_factor = [1]
    t0 = 0
    for iB in range(1, len(seq.block_events) + 1):
        block = seq.get_block(iB)
        if hasattr(block, 'gx'):
            grad = getattr(block, 'gx')
            if grad.type == 'grad':
                # In place unpacking of grad.t with the starred expression
                t = grad.delay + [0, *(grad.t + (grad.t[1] - grad.t[0]) / 2),
                                  grad.t[-1] + grad.t[1] - grad.t[0]]
                waveform_x = np.array([grad.first, grad.last])
                waveform_x = 1e-3 * np.insert(waveform_x, 1, grad.waveform)
            else:
                t = np.cumsum([0, grad.delay, grad.rise_time, grad.flat_time, grad.fall_time])
                waveform_x = 1e-3 * grad.amplitude * np.array([0, 0, 1, 1, 0])
            fig_xx_sp_list[0].plot(t_factor * (t0 + t), waveform_x)
        t0 += calc_duration(block)

    ax = plt.gca()
    j = 0
    xd_xx = []
    yd_xx = []
    while j < len(ax.lines):
        aux_xd = ax.lines[j].get_xdata()
        aux_yd = ax.lines[j].get_ydata()
        xd_xx.append(aux_xd)
        yd_xx.append(aux_yd)
        j += 1
        del aux_xd, aux_yd

    # plt.figure('F')
    # plt.plot(xd_xx,yd_xx)

    # yy-axis
    fig_yy = plt.figure("yy")
    fig_yy_sp_list = [fig_yy.add_subplot(111)]
    t_factor = [1]
    t0 = 0
    for iB in range(1, len(seq.block_events) + 1):
        block = seq.get_block(iB)
        if hasattr(block, 'gy'):
            grad = getattr(block, 'gy')
            if grad.type == 'grad':
                # In place unpacking of grad.t with the starred expression
                t = grad.delay + [0, *(grad.t + (grad.t[1] - grad.t[0]) / 2),
                                  grad.t[-1] + grad.t[1] - grad.t[0]]
                waveform_y = np.array([grad.first, grad.last])
                waveform_y = 1e-3 * np.insert(waveform_y, 1, grad.waveform)
            else:
                t = np.cumsum([0, grad.delay, grad.rise_time, grad.flat_time, grad.fall_time])
                waveform_y = 1e-3 * grad.amplitude * np.array([0, 0, 1, 1, 0])
            # gradY = np.concatenate((gradY, waveform_y), axis=0)
            fig_yy_sp_list[0].plot(t_factor * (t0 + t), waveform_y)
        t0 += calc_duration(block)

    ay = plt.gca()
    i = 0
    xd_yy = []
    yd_yy = []
    while i < len(ay.lines):
        aux_xd = ay.lines[i].get_xdata()
        aux_yd = ay.lines[i].get_ydata()
        xd_yy.append(aux_xd)
        yd_yy.append(aux_yd)
        i += 1
        del aux_xd, aux_yd

    del j, i
# plt.show()
# fig1_sp_list = [fig1.add_subplot(211), fig1.add_subplot(212)]
#
# line = plt.gca(fig1_sp_list).get_lines(fig1_sp_list)[n]
# line = plt.gca().get_lines()[n]
#
# xd_yy = line.get_xdata(fig1_sp_list[0])
# yd_yy = line.get_ydata(fig1_sp_list[0])


# plt.figure('F')
# plt.subplot(121)
# plt.plot(np.reshape(gradX,(len(gradX),1),'F'))
# plt.title('Gx')
# plt.grid(True)
#
# plt.subplot(122)
# plt.plot(t)
# plt.title('Gy')
# plt.grid(True)
# plt.show()


#%% --- 3 - save
if save_flag:
    if reportTest:
        print("Test Report: " + seq.test_report())
    # os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/1_Spiral_test')
    standName = 'test%s_Spiraltest%s_Nreps-%s_Nslices-%s_TR-%ss_TE-%sms_sm-%s_RFbw-1_RFoff-%smm%s_Ndirs-%s_bvalue-%s%s' % (test, testSPI, nreps, nslices, round(TR), round(TE*1e3), max_slew, round(RFoffset * 1e3), diff, ndirs, bvalue[0], fatsat_str)
    if os.path.isdir("Test/" + standName):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs("Test/" + standName)

    if os.path.isdir("Test/" + standName + "/cfl_hdr_files/"):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs("Test/" + standName + "/cfl_hdr_files/")

    # save gradients
    if getGRADs:
        os.chdir("Test/" + standName + "/cfl_hdr_files/")
        gradX = {'time': xd_xx, 'Grad': yd_xx}
        savemat("GradX.mat", gradX)
        gradY = {'time': xd_yy, 'Grad': yd_yy}
        savemat("GradY.mat", gradY)

    # save '.seq' file
    # os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/1_Spiral_test')
    os.chdir("Test/" + standName)
    seq.write(standName + '.seq')

    # save parameters
    parametersMAT = {'test': test, 'testSPI':testSPI, 'TR': TR, 'DT': DT, 'nslices': nslices, 'nreps': nreps, 'Nshots':nshots,
                     'st': sliceT, 'adc_samples': ADC, 'Nx': N_x, 'alpha': alpha, 'x0': RFoffset, 'max_grad': max_grad,
                     'max_slew': max_slew, 'flipAngle': flipAngle, 'TE': TE, 'ndirs': ndirs, 'diff': diff, 'bvalue': bvalue}
    savemat("sequence_info.mat", parametersMAT)

    print("\n\n ----- Seq file saved -----  ", standName)


if seq2jem:
    seqName = seq
    xmlName = standName
    out_folder ="/Test/" + standName + "/Reults"
    if os.path.isdir(out_folder):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs(out_folder)
    seq2xml(seqName, xmlName, out_folder)
    print("\n\n ----- Converted to .xml -----  ", standName)
