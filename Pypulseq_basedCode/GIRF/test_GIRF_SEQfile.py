#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Name:
    - test_GIRF_SEQfile
Description:
    - This code creates GIRF (with external gradients) '.seq' file to run in scanner

by: TTFernandes, IST - 2020
"""

# =============================================================================

#%% --- 0 - Import functions
#import make_t1_mprage as mtp
import os
import matplotlib
import scipy.io
import numpy as np
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import sys

sys.path.append('/Toolboxes')

os.chdir('GIRFOS_tool/Pypulseq_basedCode/GIRF')
from make_GIRF_SEQfile import generate_GIRF_SeqFile
from scipy.io import savemat

#%% --- 0.2 - Settings

test       = 64      # define test for trace in google drive

save_flag  = False    # save file
testX      = True    # test X-axis
testY      = True    # test Y-axis
plotTest   = True   # Plot sequence
reportTest = False   # Report the sequence
tPrePost   = False   # Test ADC pre gradient and post gradient
Dum        = True    # Dummy for first pulse


    # --- load data
dir = 'Pypulseq_basedCode/GIRF/pulses'

Npulses   = 12      # Number of triangular pulses
nreps     = 3       # Number of repetitions
nslices   = 1       # Number of slices
st        = 3e-3    # slice thickness (m)
RFoffset  = 50e-3   # RF offset - Distance for measuring (m)
type      = 'SCAN'  # Type: 'SCAN' or 'PULS' ('0_4')
tPulTotal = '10'    # Time total of scan - '0_4'; '10'; '30'; '40'
p         = 115

TR        = 2000e-3           # [s]
TE        = 10e-3             # [s]
DT        = np.float(10e-6)   # [s]
ndirs     = 1                 # number of directions

os.chdir(dir)
mat = scipy.io.loadmat('girf_{}_pulse_{}_tTotal_{}ms_DT_1_1e-5_p{}.mat'.format(type,Npulses,tPulTotal,p))
# mat = scipy.io.loadmat('girf_pulse_{}_DT_1_1e-5_p{}.mat'.format(Npulses,p))
# mat = scipy.io.loadmat('girf_pulse_SCAN_{}.mat'.format(Npulses))
i_t = mat['i_t']
# i_t = i_t[0:-1:8,:]

shape_i_t = i_t.shape
ADC       = shape_i_t[0]

del mat

#%% --- 0.75 - Plot loaded data
#
# if plotTest:
#     # for i in range(shape_i_t[1]):
#     #     plt.figure('F')
#     #     plt.plot(np.reshape(i_t[:,i],(len(i_t[:,i]),1),'F'))
#     #     plt.title('GIRF')
#     #     plt.grid(True)
#     # plt.show()

#%% --- 1 - Run function to create struct 'seq' for Spiral
max_grad, max_slew, grdrastTime, flipAngle, Tpre, Tpos = 30, 125, DT, 90, 0, 0
if tPrePost:
    Tpre, Tpos = 5e-3, 61e-3
seq = generate_GIRF_SeqFile(iT=i_t, tr=TR, te=TE, Npuls=Npulses, mg=max_grad, ms=max_slew, grt=grdrastTime, fA=flipAngle, n_slices=nslices, reps=nreps, sliceT=st, tPre=Tpre, tPos=Tpos, tX=testX, tY=testY, tPlot=plotTest, tReport=reportTest,rf_offset=RFoffset, testPRE_POST=tPrePost, dummy=Dum)

#%% --- 2 - save file '.seq'
if save_flag:
    if testX:
        ttX = 1
    else:
        ttX = 0
    if testY:
        ttY = 1
    else:
        ttY = 0
    # print("Test Report: " + seq.test_report())
    # os.chdir('/GIRF_test')
    standName = 'test%s_GIRF_tX-%s_tY-%s_Npulses-%s_Nreps-%s_Nslices-%s_TR-%ss_sm-%s_RFbw-1_RFoff-%smm_tPulTotal-%sms' % (test, ttX, ttY, Npulses, nreps, nslices, round(TR), p, round(RFoffset*1e3),tPulTotal)
    if os.path.isdir("Test/" + standName):
        print("\n\n - Directory already exist - ")
    else:
        os.makedirs("Test/" + standName)
    os.chdir("Test/" + standName)
    seq.write(standName + '.seq')
    parametersMAT = {'test':test,'TR':TR,'TE':TE,'DT':DT,'Npulses':Npulses,'nslices':nslices,'nreps':nreps,'st':st,
                     'x0':RFoffset,'max_grad':max_grad,'max_slew':p,'flipAngle':flipAngle,'Tpre':Tpre,
                     'Tpos':Tpos,'ndirs':ndirs,'adc_samples':ADC,'type':type,'tPulTotal':tPulTotal}
    savemat("sequence_info.mat", parametersMAT)

    print("\n\n ----- Seq file saved -----  ", standName)
