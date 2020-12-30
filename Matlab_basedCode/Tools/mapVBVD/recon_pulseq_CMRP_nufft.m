%% This file reconstructs data from MBBI "Eve" Prisma 
addpath(genpath('.'));
ktraj = 'nonCartesian';
spiral =1;

switch ktraj
    case 'Cartesian'
        [kspace_data,image_data]= dat2mat();
        %% Sum of squares combination
        sos=abs(sum(image_data.^2,3).^(1/2));
        sos=sos./max(sos(:));
    case 'nonCartesian'
         [kspace_data]= dat2mat_nonCart();
         if(spiral)
             kspace_data = permute(kspace_data, [2 1 3]);
             kspace_data = squeeze(kspace_data(2:end,:,:)); %First data is tricky
             
             %% Create par
  par.recon.psd = 'sim';
             
              [out,par] = get_sens(data,par);
             [Img] = recon_nufft_irt(kspace_data, ktrajs);
             
             for shot = 1:size(kspace_data,2)
                 plot(squeeze(real(kspace_data(:,shot,1))));hold on;
             end
         end
         
    case 'Remodel'
         [kspace_data]= dat2mat_remodel();
end





% figure;
% imshow(sos,[]);
% imwrite(sos, ['img_combined.png'])