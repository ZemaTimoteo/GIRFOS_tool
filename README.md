# GIRFOS_tool  (Gradient Impulse Response Function - Open-Source Toolbox)
# 

A major challenge in MRI, particularly when aiming to reduce equipment costs, is to be able to cope with imperfect hardware performance. This is particularly important when aiming to reduce equipment costs as this may imply relaxing on hardware performance. Even when using current clinical systems, calibration is often required to reduce image artifacts; an example is the impact of eddy-currents in fast imaging readouts (i.e. spiral).

  - To help control for eddy current effects, it is important to determine the gradient impulse response function (GIRF) using a calibration pulse sequence.
  - GIRFOS allows to design MRI pulse sequences to measure the Impulse Response Function of a specific MRI scanner.
  - This tool enables pulse sequence design in both Python and Matlab. Prior knowledge of Pypulseq/Pulseq, respectively is thus helpful to understand the pulse sequence design

The present implementation aims to contribute towards the development of open-source tools for MRI, particularly of sequences based on spiral imaging 

Disclaimer: we use the following toolboxes:
- https://mrirecon.github.io/bart/

- https://github.com/ndwork/griddingCodes

- https://github.com/mriphysics/verse-mb

- https://github.com/imr-framework/pypulseq

- https://github.com/openjournals/joss-reviews/issues/2478

To run the matlab codes, download and add the following toolboxes to the Tools folder:
- bart - https://mrirecon.github.io/bart/
- Gridding - https://github.com/ndwork/griddingCodes

** To run the python codes create a Toolbox folder with the following toolboxes:
- code from pypulseq - https://github.com/imr-framework/pypulseq       
- code from pulseqDiffusion - https://github.com/ritagnunes/PulseqDiffusion
