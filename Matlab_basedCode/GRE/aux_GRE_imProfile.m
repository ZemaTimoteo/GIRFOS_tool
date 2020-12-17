% read signal
% % clear all
% % clc

cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/5_GRE/test_Andreia_pGRE');
load('pulseq_molli_bSSFP_5RF_Linear200.mat')

imTHEO_GRE = imrotate(sos(:,:,5),-90);

clear sos
%%
% 1.1 - slice lines
x_im(1,:)     = [size(imTHEO_GRE,1)/10*8  size(imTHEO_GRE,1)/10*8];
x_im(2,:)     = [size(imTHEO_GRE,1)/10*3  size(imTHEO_GRE,1)/10*7];
x_im(3,:)     = [size(imTHEO_GRE,1)/10*0  size(imTHEO_GRE,1)/10*10];

y_im(1,:)     = [size(imTHEO_GRE,1)/10*0 size(imTHEO_GRE,1)/10*10];
y_im(2,:)     = [size(imTHEO_GRE,1)/10*0 size(imTHEO_GRE,1)/10*10];
y_im(3,:)     = [size(imTHEO_GRE,1)/10*4 size(imTHEO_GRE,1)/10*4];

n = 133;  %size(Norm_im_crop,1);

% 1.2 - plots
figure()
for ji=1:size(y_im,1)
    imProfGRE_THEO(:,ji) = improfile(imTHEO_GRE,x_im(ji,:),y_im(ji,:),n);
    subplot(2,3,ji)
    plot(imProfGRE_THEO(:,ji),'r')
    hold on
    xlabel('Nx')
    ylabel('Pixel Intensity')
    title(['ImProfile - GRE'])
    ji+3
    subplot(2,3,ji+3)
    imshow(imTHEO_GRE,[]);
    hold on
    line(x_im(ji,:),y_im(ji,:))    
    title('Imag Recon - slice')
end

