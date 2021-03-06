This file helps describing pypulseq_based Code.

### Folder: GIRF
  - Pulses (folder with files for GIRF pulses)
  - test_GIRF_SEQfile.py (This code creates GIRF (with external gradients) '.seq' file to run in scanner)
  - make_GIRF_SEQfile.py (aux function)
  
  - Steps:
      - Create a Toolboxes folder with:
          - code from pypulseq - https://github.com/imr-framework/pypulseq
          - code from pulseqDiffusion - https://github.com/ritagnunes/PulseqDiffusion
      - Set paths
      - To define different GIRF pulses run code from folder ('GIRF_shape_script.m') - https://github.com/ZemaTimoteo/GIRFOS_tool/tree/main/Matlab_basedCode/Tools
      - Run the test_GIRF_SEQfile.py
      
### Folder: GRE
  - test_GRE_SEQfile.py (This code creates GRE '.seq' file to run in scanner)
  
  - Steps:
      - Set paths
      - Run the test_GRE_SEQfile.py
  
### Folder: Spiral
  - Pulses (folder with files for SPIRAL pulses)
  - test_Spiral_SEQfile.py (This code creates SPIRAL (with external gradients) '.seq' file to run in scanner)
  - make_Spiral_SEQfile.py (aux function)
  
    
  - Steps:
      - Create a Toolboxes folder with:
          - code from pypulseq - https://github.com/imr-framework/pypulseq
          - code from pulseqDiffusion - https://github.com/ritagnunes/PulseqDiffusion      
      - Set paths
      - Set parameters
      - Run the test_Spiral_SEQfile.py
      
### Folder: Diffusion (To be runned by test_Spiral_SEQfile w/ Diffusion)
  - diff_funcs.py
  - k2g.py
  - make_Diffusion_SEQfile.py (This function builds a spin echo diffusion weighting sequence with spiral trajectories and EPI trajectory optimizing the TE for a desired b-value.)
  - vds_2d.py 
  
