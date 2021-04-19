#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sept 1 18:41:25 2020
Run B0ECC test
@author: tfernandes
"""

# =============================================================================
PC = 0 # PC = 1 - Seia OR PC = 0 - My PC

#####
#%% --- 0 - Import functions
#####
import os
import matplotlib
import sys
import scipy.io
import numpy as np
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from scipy.io import savemat

if PC == 1:
    sys.path.append('/home/tfernandes/Documents/PYTHON/Toolboxes')
    os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/B0_ECC')
elif PC == 0:
    os.chdir('D:/Tiago/Trabalho/2021_2025_PhD/Projects/GIRFOS_tool/Code/1_pythonCodes/B0_ECC')

from make_B0ECC_SEQfile import generate_B0ECC_SeqFile

#####
#%% --- 0.2 - Settings
#####
test       = 131    # define test for trace in google drive

save_flag  = True  # save file
testX      = True   # test X-axis
testY      = True   # test Y-axis
testZ      = False  # test Z-axis
tPrePost   = False  # Test ADC pre gradient and post gradient
plotTest   = True
# Plot sequence
reportTest = False  # Report the sequence
dum        = True   #

#####
#%% --- 0.5 - Import data
#####

    # --- load data
Npulses   = 12       # Number of triangular pulses
nreps     = 1
nslices   = 1
st        = 2e-3    # slice thickness (m)
RFoffset  = 10e-3   # RF offset - distance to isocenter (mm)
type      = 'SCAN'  # Type: 'SCAN' or 'PULS' ('0_4')
tPulTotal = '10'    # Time total of scan - '0_4'; '10'; '30'; '40'
p         = 115

flipAngle = 74                # [degrees] - see flip angle vs TR - http://www.mritoolbox.com/ErnstAngle.html
TR        = 400e-3            # [s]
TE        = 10e-3             # [s]
DT        = float(10e-6)   # [s]
ndirs     = 1
dummies   = 5                 # number of dummies


if PC == 1:
    dir = '/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/B0_ECC/pulses'
elif PC == 0:
    dir = 'D:/Tiago/Trabalho/2021_2025_PhD/Projects/GIRFOS_tool/Code/1_pythonCodes/B0_ECC/pulses'

os.chdir(dir)
mat = scipy.io.loadmat('girf_{}_pulse_{}_tTotal_{}ms_DT_1_1e-5_p{}.mat'.format(type,Npulses,tPulTotal,p))
# mat = scipy.io.loadmat('girf_pulse_{}.mat'.format(Npulses))
# mat = scipy.io.loadmat('girf_pulse_SCAN_{}.mat'.format(Npulses))
i_t = mat['i_t']
# i_t = i_t[0:-1:8,:]

shape_i_t = i_t.shape
ADC       = shape_i_t[0]

del mat

#####
#%% --- 0.75 - Plot loaded data
#####
#
# if plotTest:
#     # for i in range(shape_i_t[1]):
#     #     plt.figure('F')
#     #     plt.plot(np.reshape(i_t[:,i],(len(i_t[:,i]),1),'F'))
#     #     plt.title('GIRF')
#     #     plt.grid(True)
#     # plt.show()

#####
#%% --- 1 - Run function to create struct 'seq' for Sequence
#####
max_grad, max_slew, grdrastTime = 32, 130, DT,
if tPrePost:  # Tpre - time preADC, Tpos - time posADC
    Tpre, Tpos = 5e-3, 61e-3
    testPrePost = 1
else:
    Tpre, Tpos = 0, 0
    testPrePost = 0

seq = generate_B0ECC_SeqFile(iT=i_t, tr=TR, te=TE, Npuls=Npulses, mg=max_grad, ms=max_slew, grt=grdrastTime, fA=flipAngle, reps=nreps, sliceT=st, tPre=Tpre, tPos=Tpos, tX=testX, tY=testY, tZ=testZ, tPlot=plotTest, tReport=reportTest, testPRE_POST=tPrePost,rf_offset=RFoffset, dummy=dum, Ndummies=dummies)

# Duration of sequence in ms
seqDur = seq.duration()
print("\n\n .. seq duration .. ", seqDur[0]/60, "in min")

#####
#%% --- 2 - save file '.seq'
#####
if save_flag:
    Naxes = 0
    if testX:
        ttX = 1
        Naxes = Naxes + 1
    else:
        ttX = 0
    if testY:
        ttY = 1
        Naxes = Naxes + 1
    else:
        ttY = 0
    if testZ:
        ttZ = 1
        Naxes = Naxes + 1
    else:
        ttZ = 0


    standName  = 'test%s_B0-ECC_tX-%s_tY-%s_tZ-%s_Npulses-%s_Nreps-%s_TR-%ss_adc_5ms-PULSE-61ms_sm-%s_RFbw-1_RFoff-%smm_PrePostTest_%s' % (test, ttX, ttY, ttZ, Npulses, nreps, round(TR), max_slew, round(RFoffset*1e3), testPrePost)
    folderName = 'test%s_B0-ECC' % (test)

    if PC == 1:
        os.chdir('/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/3_B0_ECC_test')
        if os.path.isdir("Test/" + standName):
            print("\n\n - Directory already exist - ")
        else:
            os.makedirs("Test/" + standName)
        os.chdir("Test/" + standName)
    elif PC == 0:
        os.chdir('C:/Users/tiago/OneDrive/Documentos/PhD/Projetos/qMRI/Sequences/B0ECC')
        if os.path.isdir(folderName):
            print("\n\n - Directory already exist - ")
        else:
            os.makedirs(folderName)
        os.chdir(folderName)

    seq.write(standName + '.seq')
    parametersMAT = {'test':test,'TR':TR,'TE':TE,'DT':DT,'Npulses':Npulses,'nslices':nslices,'nreps':nreps,'st':st,
                     'x0':RFoffset,'max_grad':max_grad,'max_slew':p,'flipAngle':flipAngle,'Tpre':Tpre,
                     'Tpos':Tpos,'ndirs':ndirs,'adc_samples':ADC,'Naxes':Naxes, 'dummy':dum, 'Ndummies':dummies, 'Dr':RFoffset,
                     'tPulTotal':tPulTotal, 'duration':seqDur[0]}
    savemat("sequence_info.mat", parametersMAT)

    print("\n\n ----- Seq file saved -----  ", standName)