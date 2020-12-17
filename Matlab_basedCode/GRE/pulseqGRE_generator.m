%run & save '.seq' Spiral file from pulseq 
clear all
clc
close all
cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/5_Spiral/')
addpath(genpath('/home/tfernandes/Documents/MATLAB/Toolboxes/pulseq-master/'))
%% 0 - Parameters
test = 69;
saveFile = 'Fals';
%% 1 - Run pulseq function
gre

%% 2 - Save file
if saveFile == 'True'
    seqName = ['test',num2str(test),'_pulseqGRE']
    cd_dir    = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/1_Spiral_test/Test/',seqName];
    mkdir(cd_dir)
    cd(cd_dir)
    seq.write([seqName,'.seq']);   % Output sequence for scanner
end