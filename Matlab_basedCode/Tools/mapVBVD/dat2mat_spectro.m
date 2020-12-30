function [kspace_data,image_data]= dat2mat_spectro(recon)


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


switch recon
    case 'remodel'
        image_obj = twix_obj.image;
    case 'spectro'
           image_obj = twix_obj.image;
    otherwise
        image_obj = twix_obj.image;
end
image_objD = squeeze(image_obj()); %very rare usage - interesting
sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
dimsData = image_obj.sqzDims;
dimsallData = image_obj.dataDims;



%% Concatenate slices and partitions together and has a N x N acquisition set
kspace_data = zeros([sizeData(1)*sizeData(3) sizeData(2)]);
% kspace_data = permute(kspace_data,[1 3  4  5  2]); %kx, ky, slices, partitions, nch
image_data = zeros(size(kspace_data));




tic;
for nch=1:size(kspace_data,2)
        temp = squeeze(image_objD(:,nch,:));
        temp2 = cat(1, squeeze(temp(:,1)), squeeze(temp(:,2))); %This will work only if you have 2 sets!
%        temp2 = reshape(temp, [512*2 ,1 ]);
        kspace_data(:,nch) = temp2;
        image_data(:,nch) = fftshift(fft2(temp2));
 
end
toc;

image_data = squeeze(image_data);
kspace_data = squeeze(kspace_data);











