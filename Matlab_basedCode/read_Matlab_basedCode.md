This file helps describing Matlab_based Code.

 
### Folder: Example
  - dat_file - folder with '.dat' file for example (splitted in three parts - use: https://pinetools.com/join-files to combine it in one file)
  - Spiral_dat_file: folder with spiral '.dat' file for example (splitted in three parts - use: https://pinetools.com/join-files to combine it in one file)
  - test64_GIRF_tX-1_tY-1_Npulses-12_Nreps-3_Nslices-1_TR-2s_sm-115_RFbw-1_RFoff-50mm_tPulTotal-10ms - folder for example
  - girf_SCAN_Input_freq_12_tTotal_10ms.mat - file used with 'get_GIRF.m'
  - girf_SCAN_Input_freq_sincFunct_12.mat - file used with 'get_GIRF.m'
  - girf_pulse_12_DT_1_1e-5_p110.mat - file used with 'get_GIRF.m'
  - vectSpiral_test4.mat - file for example of spiral implementation with GIRF
  - test64_IRF.mat - file for example of spiral implementation with GIRF
  
### Folder: GIRF
  - get_GIRF.m (Saves the GIRF from scanner of the respective GIRF sequence from pypulseq.seq )
  
  - To run:
      - Set your path
      - Add the following tools to the Tool folder:
          - bart - https://mrirecon.github.io/bart/
          - Gridding - https://github.com/ndwork/griddingCodes
      - As an exemple choose file '.dat' in /Example folder
            - 'dat_file' folder contains the '.dat' file in three parts
            - use https://pinetools.com/join-files to join the three files
            - copy that file to the path: 'GIRFOS_tool/Matlab_basedCode/Example/test64_GIRF_tX-1_tY-1_Npulses-12_Nreps-3_Nslices-1_TR-2s_sm-115_RFbw-1_RFoff-         
                                              50mm_tPulTotal-10ms/'
            - when running the script, select the '.dat' file you copy to the path in the last step.                                              
  
### Folder: Spiral
  - pypulseqSPIRAL_seqRecon.m - Obtain external gradients in matlab matrix from pypulseq sequence
  - recon_spiral_IRF.m - Reconstruct Spiral with GIRF
  - recon_spiral_IRF_both_IRF_THEOR.m - Run for each test of SCAN the reconstruction with Theoric ktraj and GIRF ktraj
  
  - To run 'recon_spiral_IRF.m':
      - Set your path
      - Add the following tools to the Tool folder:
          - bart - https://mrirecon.github.io/bart/
          - Gridding - https://github.com/ndwork/griddingCodes
      - As an exemple choose file Spiral '.dat' in /Example folder
            - 'Spiral_dat_file' folder contains the '.dat' file in three parts
            - use https://pinetools.com/join-files to join the three files
            - copy that file to the path: 'GIRFOS_tool/Matlab_basedCode/Example/test84_Spiraltest4_Nreps-4_Nslices-1_TR-2s_TE-0ms_sm-125_RFbw-1_RFoff-0mm_Ndirs-
                                              _bvalue-0/'
            - when running the script, select the '.dat' file you copy to the path in the last step.
  
  - To run 'recon_spiral_IRF_both_IRF_THEOR.m':
      - Set your path
      - Add the following tools to the Tool folder:
          - bart - https://mrirecon.github.io/bart/
          - Gridding - https://github.com/ndwork/griddingCodes
      - As an exemple choose file Spiral '.dat' in /Example folder
            - 'Spiral_dat_file' folder contains the '.dat' file in three parts
            - use https://pinetools.com/join-files to join the three files
            - copy that file to the path: 'GIRFOS_tool/Matlab_basedCode/Example/test84_Spiraltest4_Nreps-4_Nslices-1_TR-2s_TE-0ms_sm-125_RFbw-1_RFoff-0mm_Ndirs-
                                              _bvalue-0/'
            - when running the script, select the '.dat' file you copy to the path in the last step.
            
### Folder: Tools
  - Input_sequences (Folder)
  - mapVBVD (toolbox from https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD)
  - pulses (Folder)
  - GIRF_shape_script.m (Generate TrianglePulse in '.mat' file, for generate with external gradients a '.seq' file w/ 'test_GIRF_SEQfile.py')
  - sharp_imProfile.m (his function returns a struct 'bound' which contains the bourders of distinct regions for a imProfile matlab function)
  - gradient_create_GIRF.m (function adapted from: https://github.com/mriphysics/verse-mb/tree/master/bin)
  - gradient_distort_GIRF.m (function adapted from: https://github.com/mriphysics/verse-mb/tree/master/bin)
  - AuxSpiralTF_grid_kb.m (function to 'recon_spiral_IRF.m')


* To run this add to your 'Tool' folder the following toolboxes:
  - bart - https://mrirecon.github.io/bart/
  - Gridding - https://github.com/ndwork/griddingCodes
 
  
  Also you need to install JEMRIS (http://www.jemris.org/) to run 
