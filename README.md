# GIRFOS_tool  (Gradient Impulse Response Function - Open-Source Toolbox)
# 

A major challenge in MRI, particularly when aiming to reduce equipment costs, is to be able to cope with imperfect hardware performance. Even when using current clinical systems, calibration is often required to reduce image artifacts; an example is the impact of eddy-currents in readouts (i.e. spiral). 

- This tool helps control for such effect by determining the gradient impulse response function (GIRF) using a calibration pulse sequence.
- This tool allows to design MRI sequences to measure the Impulse Response Function of a specific scan.
- This tool work with a couple of packages: in python and Matlab.

Prior knowledge to pypulseq is helpful to understand the design of the pulse sequence.

This implementation will contribute towards the development of open source tools for spiral imaging.


Disclaim:
In this tool we used the following toolbox:

- https://mrirecon.github.io/bart/

- https://github.com/ndwork/griddingCodes

- https://github.com/mriphysics/verse-mb

- https://github.com/imr-framework/pypulseq

- https://github.com/openjournals/joss-reviews/issues/2478



To run the matlab codes add to the Tools folder the following toolboxes:
- bart - https://mrirecon.github.io/bart/
- Gridding - https://github.com/ndwork/griddingCodes


** To run the python codes create a Toolbox folder with the following toolboxes:
- code from pypulseq - https://github.com/imr-framework/pypulseq       
- code from pulseqDiffusion - https://github.com/ritagnunes/PulseqDiffusion
