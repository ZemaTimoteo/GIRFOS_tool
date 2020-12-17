%% Comparison from JEMRIS Simulation IRF with Theor
clear all
clc

recon_spiral_IRF_both_IRF_THEOR

% image initialization

imGIRF = imGIRF/max(max(imGIRF));
imTHEO = imTHEO/max(max(imTHEO));


%% Initialization

% image initialization
aux_GRE_imProfile


imUPsize = 35;
save_imTHEO_GRE = imTHEO_GRE;
% imTHEO_GRE = save_imTHEO_GRE;

   
new_imTHEO_GRE = zeros(size(imTHEO_GRE,2)+imUPsize*2,size(imTHEO_GRE,2)+imUPsize*2);
new_imTHEO_GRE(imUPsize:size(imTHEO_GRE,2)+imUPsize-1,imUPsize:size(imTHEO_GRE,2)+imUPsize-1) = imTHEO_GRE;
figure()
subplot(121)
imshow(new_imTHEO_GRE,[])
clear imTHEO_GRE
nGRE=133;
imTHEO_GRE = imresize(new_imTHEO_GRE,nGRE/size(new_imTHEO_GRE,1));
imTHEO_GRE = imTHEO_GRE/max(max(imTHEO_GRE));
subplot(122)
imshow(imTHEO_GRE,[])



%% 2nd variaty - plot


x_im(1,:)     = [size(imTHEO,1)/10*7  size(imTHEO,1)/10*7];
x_im(2,:)     = [size(imTHEO,1)/10*3  size(imTHEO,1)/10*7];
x_im(3,:)     = [size(imTHEO,1)/10*0  size(imTHEO,1)/10*10];

y_im(1,:)     = [size(imTHEO,1)/10*0 size(imTHEO,1)/10*10];
y_im(2,:)     = [size(imTHEO,1)/10*0 size(imTHEO,1)/10*10];
y_im(3,:)     = [size(imTHEO,1)/10*4 size(imTHEO,1)/10*4];


figure()
for ii = 1:size(imGIRF,3)
    for ji=1:size(y_im,1)
        imProf_GIRF(:,ji,ii)     = improfile(imGIRF(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_THEO(:,ji,ii)     = improfile(imTHEO(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_THEO_GRE(:,ji,ii) = improfile(imTHEO_GRE(:,:,ii),x_im(ji,:),y_im(ji,:),n);

        
        subplot(3,4,4*ji-3)
        plot(imProf_THEO_GRE(:,ji,ii),'r')
        hold on
        plot(imProf_THEO(:,ji,ii),'b')
        plot(imProf_GIRF(:,ji,ii),'g')        
        xlabel('Nx')
        ylabel('Pixel Intensity')
        if ji==1
            title(['ImProfile: bSSFP(r) THEOR(b) GIRF(g)'])
        end
        hold on
        ylim([0 1])                        
        xlim([0 133])
        subplot(3,4,4*ji-2)
        imshow(imTHEO_GRE(:,:,ii),[]);
        hold on
        if ji==1
            title('Image Recon bSSFP')        
        end
        hold on
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 1])
        
        subplot(3,4,4*ji-1)
        imshow(imTHEO(:,:,ii),[]);
        hold on
        if ji==1
            title('Image Recon THEOR')        
        end
        hold on
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 1])
        
        subplot(3,4,4*ji-0)
        imshow(imGIRF(:,:,ii),[]);
        hold on
        if ji==1
            title('Image Recon GIRF')
        end
        hold on
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 1])
    end
end


%% edge-shaprness GRE

for ji=1:size(y_im,1)
    bounds = sharp_imProfile(imProf_THEO_GRE(:,ji),'True');
    resSHARP_GRE(1,ji,ii) = (bounds.firstUP   - bounds.firstLOW ) * res; % in mm
    resSHARP_GRE(2,ji,ii) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
    clear bounds    
end

mean(mean(resSHARP_GRE))


