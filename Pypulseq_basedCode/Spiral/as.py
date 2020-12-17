

import os
import scipy.io
import matplotlib
import numpy as np
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import sys


sys.path.append('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/py2jemris')
sys.path.append('/home/tfernandes/Documents/PYTHON/Toolboxes')

from seq2xml import seq2xml

os.chdir("/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/py2jemris/try_seq2xml")
seqName = "spgr_gspoil_N15_Ns1_TE5ms_TR10ms_FA30deg.seq"
xmlName = "testSeq_xml"
out_folder = "/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/py2jemris/try_seq2xml"
seq2xml(seqName, xmlName, out_folder)




#
#
#
#
#
#
#
#
#
#
# if getGRADs:
#     fig2 = plt.figure(2)
#     fig2_sp_list = [fig2.add_subplot(311, sharex=sp11), fig2.add_subplot(312, sharex=sp11),
#                     fig2.add_subplot(313, sharex=sp11)]
#
#     t_factor = [1]
#     time_range = (0, np.inf)
#     grad_channels = ['gx', 'gy', 'gz']
#     t0 = 0
#     for iB in range(1, len(seq.block_events) + 1):
#         block = seq.get_block(iB)
#         for x in range(0, len(grad_channels)):
#             if hasattr(block, grad_channels[x]):
#                 grad = getattr(block, grad_channels[x])
#                 if grad.type == 'grad':
#                     # In place unpacking of grad.t with the starred expression
#                     t = grad.delay + [0, *(grad.t + (grad.t[1] - grad.t[0]) / 2),
#                                       grad.t[-1] + grad.t[1] - grad.t[0]]
#                     waveform = np.array([grad.first, grad.last])
#                     waveform = 1e-3 * np.insert(waveform, 1, grad.waveform)
#                 else:
#                     t = np.cumsum([0, grad.delay, grad.rise_time, grad.flat_time, grad.fall_time])
#                     waveform = 1e-3 * grad.amplitude * np.array([0, 0, 1, 1, 0])
#                 fig2_sp_list[x].plot(t_factor * (t0 + t), waveform)
#         t0 += calc_duration(block)
#
#     plt.show()