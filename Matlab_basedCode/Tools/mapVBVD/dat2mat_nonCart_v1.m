function [kspace_data,dimsData]= dat2mat_nonCart_v1()


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
[Filename,Pathname] = uigetfile('*.dat','Pick the raw data file');


%% Read data using mapVBVD
% image_obj = mapVBVD(fullfile(Pathname,Filename));
twix_obj = mapVBVDVE(fullfile(Pathname,Filename));

% image_obj = twix_obj.image; %Body coil
if(~iscell(twix_obj))
% if ( size(twix_obj.image) ==1)
    image_obj = twix_obj.image;
else
    image_obj = twix_obj{2}.image; %change it back to 2
end

sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
dimsData = image_obj.sqzDims;
dimsallData = image_obj.dataDims;


%%
 kspace_data = squeeze(image_obj(:,:,:,:));
 kspace_data = squeeze(kspace_data);
