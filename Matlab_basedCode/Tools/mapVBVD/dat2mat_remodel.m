function [kspace_data,image_data]= dat2mat(recon)


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


  image_obj = twix_obj{2}.image;


sizeData = image_obj.sqzSize; %kx, nch, ky, slices, partitions
dimsData = image_obj.sqzDims;
dimsallData = image_obj.dataDims;


%% Create kspace data matrix based on selected Data
if(numel(sizeData)==3)
    sizeData(4)=1;
    sizeData(5)=1;
end


%% Put it into matrices - this has the blocks seperate
% kspace_data = zeros(sizeData);
% kspace_data = permute(kspace_data,[1 3  4  5 2]); %kx, ky, slices, partitions, nch
% image_data = zeros(size(kspace_data));
% tic;
% for nch=1:size(kspace_data,5)
%    for partitions=1:size(kspace_data,4)
%        for nslices=1:size(kspace_data,3)
%        
%         temp = squeeze(image_obj(:,nch,:,nslices,partitions));
%         kspace_data(:,:,nslices,partitions,nch) =temp;
%         image_data(:,:,nslices,partitions,nch) =fft2(fftshift(temp));
%         
%         
%         imagesc(abs(fft2(fftshift(temp))));drawnow;
%        end
%    end
% end
% 
% toc;

%% Concatenate slices and partitions together
% kspace_data = zeros([sizeData(1) sizeData(3) sizeData(4)*sizeData(5) sizeData(2)]);
% % kspace_data = permute(kspace_data,[1 3  4  5  2]); %kx, ky, slices, partitions, nch
% image_data = zeros(size(kspace_ata));
% sl = 0;
% partitions=1;
% 
% 
% tic;
% for nch=1:size(kspace_data,4)
% 
%        for nslices=1:size(kspace_data,3)
%            
%        sl = sl +1;
%            
%         temp = squeeze(image_obj(:,nch,:,sl,partitions));
%         kspace_data(:,:,nslices,nch) = temp;
%         image_data(:,:,nslices,nch) = fftshift(fft2(temp));
%         
%         if(sl == sizeData(4))
%             sl=0;
%             if(partitions == sizeData(5))
%                 partitions =0;
%             end
%             partitions = partitions+1;
%         end
%         
% %         imagesc(abs(fftshift(fft2(temp))));drawnow;
%      imagesc(abs(((temp))));drawnow;
%         pause;
%         
%         
%         
%        end
% 
% end
% toc;

%% Concatenate slices and partitions together and has a N x N acquisition set
kspace_data = zeros([sizeData(1) sizeData(3) sizeData(4)*sizeData(5) sizeData(2)]);
% kspace_data = permute(kspace_data,[1 3  4  5  2]); %kx, ky, slices, partitions, nch
image_data = zeros(size(kspace_data));
sl = 0;
partitions=1;





%% Determine k-space shift to the center

temp = squeeze(image_obj(:,1,:,round(sizeData(4)/2),round(sizeData(5)/2)));
[~,idx] = max(abs(temp(:)));
[idx,idy] = ind2sub(size(temp),idx);

%%


tic;
for nch=1:size(kspace_data,4)

       for nslices=1:size(kspace_data,3)
           
       sl = sl +1;
           
        temp = squeeze(image_obj(:,nch,:,sl,partitions));
        temp = circshift(temp, [round(size(temp,1)/2) - idx,round(size(temp,2)/2) - idy ]);
        trunc_part =round(0.5 .*(size(temp,1) - size(temp,2)));
        temp = squeeze(temp(trunc_part+1:end-trunc_part,:));
        
        
        kspace_data(:,:,nslices,nch) = temp;
        image_data(:,:,nslices,nch) = fftshift(fft2(temp));
        
        if(sl == sizeData(4))
            sl=0;
            if(partitions == sizeData(5))
                partitions =0;
            end
            partitions = partitions+1;
        end
        
%         imagesc(abs(fftshift(fft2(temp))));drawnow;
%      imagesc(abs(fftshift(fft2(temp))));drawnow;
%         pause;
%         
        
        
       end

end
toc;

image_data = squeeze(image_data);
kspace_data = squeeze(kspace_data);











