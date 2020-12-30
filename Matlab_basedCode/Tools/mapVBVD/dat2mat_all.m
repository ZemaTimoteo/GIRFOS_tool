function [kspace_data,dimsData,image_data]= dat2mat_all(process_kspace)


%% Author - Sairam Geethanath
% Date 20130824
% Input - *.dat file
%       - Process_kspace: for averaged/sos kspace
% Output - Kspace, Image space
% We may have to have a relook at this code after ISMRM 2014 due to its
% rigidity for N x N acquisitions


%% Modification history
% to make it generic - remove hardcoding
% 06/09/2014 - requires specific handling by user for recon
% Will get only kspace_data here if no processing is specified
% Else will get complete image data of kx,ky,nch 
%% Obtain the file
[Filename,Pathname] = uigetfile('*.dat','Pick the raw data file');


%% Read data using mapVBVD
image_obj = mapVBVD(fullfile(Pathname,Filename));
sizeData = image_obj.sqzSize; %
dimsData = image_obj.sqzDims; % which attributes


kspace_data = image_obj(:);
kspace_data = reshape(kspace_data, sizeData);
image_data =0;%Dummy output for placeholder

%% Process kspace - keep adding switch cases as and when, easily extendable now

if (process_kspace)

for k =1:length(dimsData)
    
    switch dimsData{k}
        
%         case 'Col' - DO NOT implement
            
            
%         case  'Cha' - DO NOT implement
            
            
%         case 'Lin' - DO NOT implement
            
            
        case 'Ave'
            
            kspace_data = sum(kspace_data,k).* (1/sizeData(k)); 
            
            
            
            
            
            
            
    end
end


        %% The image space will have only x, y, z, nch at max
        
        imszData = size(kspace_data);
        image_data = zeros(size(kspace_data));
        
        
        switch length(imszData)
            
            case 3 %x, y, nch
                
               for nch =1:size(kspace_data,2)
                  
                   image_data(:,nch,:) = fftshift(fft2(squeeze(kspace_data(:,nch,:))));
                   
% imagesc(abs(fftshift(fft2(squeeze(kspace_data(:,nch,:)))))); %DEBUG code
                   
               end
                image_data = permute(image_data, [1 3 2]);
                figure;imagesc(abs(image_data(:,:)));% DEBUG code
                image_data = squeeze(sqrt(sum(image_data.^2,3))); %SOS recon
                kspace_data = ifft2(image_data); %Generate kspace data on the SOS image
            case 4 %x,y,z,nch
        
        
        end
        

end
