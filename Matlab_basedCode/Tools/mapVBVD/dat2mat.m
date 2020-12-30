function [kspace_data,dimsData]= dat2mat(varargin)


%% Author - Sairam Geethanath
% Date 20130824
% Input - *.dat file
% Output - Kspace, Image space
% We may have to have a relook at this code after ISMRM 2014 due to its
% rigidity for N x N acquisitions


%% Modification history
% to make it generic - remove hardcoding
% 06/09/2014
%% Obtain the file
if nargin<1
    [Filename,Pathname] = uigetfile('*.dat','Pick the raw data file');
else
    [Pathname,name,ext]=fileparts(varargin{1});
    Filename=strcat(name,ext);
end

%% Read data using mapVBVD
% image_obj = mapVBVD(fullfile(Pathname,Filename));
twix_obj = mapVBVDVE(fullfile(Pathname,Filename));


if(iscell(twix_obj))
           image_obj = twix_obj{2}.image;
else
    image_obj = twix_obj.image;
end
image_objD = squeeze(image_obj()); %very rare usage - interesting
% sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
dimsData = image_obj.sqzDims;
% dimsallData = image_obj.dataDims;%all data dimensions
kspace_data = squeeze(image_objD);












