%% Compare for each test of SCAN the reconstruction with Theoric ktraj and B0ECC ktraj anc B0ECC correction
%   - uses recon_spiral_B0ECC_multipleTests.m

%  TTFernandes, Dez 2020
% IST - ISR

%% part 0 - Addpath & Get Folders
tic
clear all;
addpath(genpath('GIRFOS_tool/Matlab_basedCode/Tools/'));  % To reconstruct        
% close all
clc


%% part 1 - Get Files


% get file
str ={'select file'}
s = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str);
if s==1
    [file,folder]=uigetfile;
else
    folder=uigetdir
end
strTest1 = '.dat';
auxA     = strfind(file,strTest1);
filename = file(1:auxA-1);

 
file_folder_cfl_hdr = [folder '/cfl_hdr_files'];
file_folder_results = [folder, '/Results'];

mkdir([file_folder_cfl_hdr]);
mkdir([file_folder_results]);
cd(folder)

[files]   = dir(folder);
fn        = files(5).name;       % get '.dat' file
fn_info   = files(6).name;       % get '.mat' file with sequence parameters from Python
tmp       = load(fn_info);       % get values from '.mat' file with sequence parameters

cd('GIRFOS_tool/Matlab_basedCode/Example/')
name_seq_template = ['vectSpiral_test',num2str(tmp.testSPI),'.mat'];

copyfile(name_seq_template,folder)
cd(folder)

fn_MATLAB = files(8).name;       % get '.mat' file with sequence parameters from MATLAB

tmpMATL   = load(fn_MATLAB);
clear strTest1 filename auxA fn_info fn_MATLAB

fprintf('\n\n 1 - Sucessfully finished \n\n')


%% part 2 - select recon with Theoric ktraj with recon_spiral_B0ECC
test      = 120;                % test to get B0ECC

testKtraj = 'THEO';             % Correct with BO_ECC H1 theoretical ktrajct - B0EC or THEO (no test)
H0_corr   = 'Fals';             % Correct with B0_ECC H0 - 'True' or 'Fals'
testInv   = 'AboS';             % 'Freq' or 'Time' or 'AboS'

Diffusion = 'Fals';

plotTest2        = 'Fals';
plotTest4        = 'Fals';
saveData         = 'Fals';       % Save raw Data
saveResults      = 'True';       % Save Results
test_B0ECC_multi = 'True';

recon_spiral_B0ECC


%% part 3 - select recon with H1 correct ktraj recon_spiral_B0ECC
clearvars -except tmpMATL  name_seq_template tmp fn file_folder_results ...
    file_folder_cfl_hdr myCD test testInv plotTest2 plotTest4 folder ...
    saveData saveResults test_B0ECC_multi test PC testKtraj

testKtraj = 'B0EC';              % Correct with BO_ECC H1 theoretical ktrajct - B0EC or THEO (no test)
H0_corr   = 'Fals';              % Correct with B0_ECC H0 - 'True' or 'Fals'

recon_spiral_B0ECC


%% part 4 - select recon with H0 + H1 correct ktraj recon_spiral_B0ECC
clearvars -except tmpMATL  name_seq_template tmp fn file_folder_results ...
    file_folder_cfl_hdr myCD test testInv plotTest2 plotTest4 folder ...
    saveData saveResults test_B0ECC_multi test PC testKtraj

testKtraj = 'B0EC';              % Correct with BO_ECC H1 theoretical ktrajct - B0EC or THEO (no test)
H0_corr   = 'True';              % Correct with B0_ECC H0 - 'True' or 'Fals'

recon_spiral_B0ECC


%% part 5 - obtain differences

close all

cd(file_folder_results)
clear im_sos imTHEO
load(['test_',num2str(tmp.test),'_imRecon_Theor.mat']);
imTHEO = im_sos;
clear imGIRF im_sos
if Diffusion == 'True'
    load(['test_',num2str(tmp.test),'_imRecon_GIRF_Diff.mat']);
    im_BOEC_H1 = im_sos;    
else
    load(['test_',num2str(tmp.test),'_imRecon_BOEC_H1.mat']);
    im_BOEC_H1 = im_sos;  
    clear im_sos
    load(['test_',num2str(tmp.test),'_imRecon_BOEC_H0_H1.mat']);
    im_BOEC_H0_H1 = im_sos;  
end

maxValueCompar = max([max(max(max(im_BOEC_H1))) max(max(max(im_BOEC_H0_H1))) max(max(max(imTHEO)))]);
figure()
for ii=1:size(imTHEO,3)
    subplot(3,size(imTHEO,3),ii)
    imshow(imTHEO(:,:,ii),[]);
    hold on
    title(['Rep ',num2str(ii)])
    ylabel(['Theor'])
    caxis([0 maxValueCompar])
    subplot(3,size(im_BOEC_H1,3),ii+size(im_BOEC_H1,3))
    imshow(im_BOEC_H1(:,:,ii),[]);
    hold on
    title(['Rep ',num2str(ii)])
    ylabel(['H1'])
    caxis([0 maxValueCompar])
    subplot(3,3,3)
    subplot(3,size(im_BOEC_H0_H1,3),ii+3+size(im_BOEC_H0_H1,3))
    imshow(im_BOEC_H0_H1(:,:,ii),[]);
    hold on
    title(['Rep ',num2str(ii)])
    ylabel(['H0 & H1'])
    caxis([0 maxValueCompar])      
end


if tmp.nreps >1
    imTHEO_orig  = imTHEO;
    imH1_orig    = im_BOEC_H1;
    imH0_H1_orig = im_BOEC_H0_H1;
    
% %     imTHEO = imGIRF_orig;
% %     imGIRF = imGIRF_orig;
    clear imTHEO imGIRF
    if Diffusion == 'True'        
        % ---------------- b=0 -------------------
        imTHEO_noDiff = imTHEO_orig(:,:,2:tmp.nreps);
        imH1_noDiff   = imH1_orig(:,:,2:tmp.nreps);
        
        % all plots b=0
        figure()
        for l=1:size(imTHEO_noDiff,3)
            subplot(2,size(imTHEO_noDiff,3)+1,l+1)
            imshow(imTHEO_noDiff(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l)])
            ylabel(['Theor'])
            caxis([0 maxValueCompar])
            subplot(2,size(imH1_noDiff,3)+1,l+2+size(imH1_noDiff,3))
            imshow(imH1_noDiff(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l)])
            ylabel(['GIRF'])
            caxis([0 maxValueCompar])
        end
        imTHEO_noDiff_mean = mean(imTHEO_noDiff,3);
        imGIRF_noDiff_mean = mean(imH1_noDiff,3);
        
        % plot for mean
        subplot(2,size(imH1_noDiff,3)+1,1)
        imshow(imTHEO_noDiff_mean,[]);
        hold on
        title('Avg - b=0')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,size(imH1_noDiff,3)+1,5)
        imshow(imGIRF_noDiff_mean,[]);
        hold on
        title('Avg - b=0')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        
        
        % ---------------- b=500 -------------------
        imTHEO_xx = imTHEO_orig(:,:,tmp.nreps+1:2:tmp.nreps*(tmp.ndirs+1));
        imTHEO_yy = imTHEO_orig(:,:,tmp.nreps+2:2:tmp.nreps*(tmp.ndirs+1));   
        imGIRF_xx = imH1_orig(:,:,tmp.nreps+1:2:tmp.nreps*(tmp.ndirs+1));
        imGIRF_yy = imH1_orig(:,:,tmp.nreps+2:2:tmp.nreps*(tmp.ndirs+1));        
        
        % all plots b=500 - x ---
        figure()
        for l=1:size(imTHEO_xx,3)
            subplot(2,size(imTHEO_xx,3)+1,l+1)
            imshow(imTHEO_xx(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l), ' - xx - b=500 s/mm^2'])
            ylabel(['Theor'])
            caxis([0 maxValueCompar])
            subplot(2,size(imGIRF_xx,3)+1,l+2+size(imGIRF_xx,3))
            imshow(imGIRF_xx(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l), ' - xx - b=500 s/mm^2'])
            ylabel(['GIRF'])
            caxis([0 maxValueCompar])
        end
        imTHEO_xx_mean = mean(imTHEO_xx,3);        
        imGIRF_xx_mean = mean(imGIRF_xx,3);
               
        % plot for mean
        subplot(2,size(imH1_noDiff,3)+1,1)
        imshow(imTHEO_xx_mean,[]);
        hold on
        title('Avg - b=500 s/mm^2')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,size(imH1_noDiff,3)+1,5)
        imshow(imGIRF_yy_mean,[]);
        hold on
        title('Avg - b=500 s/mm^2')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        
        % all plots b=500 - y ---
        figure()
        for l=1:size(imTHEO_yy,3)
            subplot(2,size(imTHEO_yy,3)+1,l+1)
            imshow(imTHEO_yy(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l), ' - yy - b=500 s/mm^2'])
            ylabel(['Theor'])
            caxis([0 maxValueCompar])
            subplot(2,size(imGIRF_yy,3)+1,l+2+size(imGIRF_yy,3))
            imshow(imGIRF_yy(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l), ' - yy - b=500 s/mm^2'])
            ylabel(['GIRF'])
            caxis([0 maxValueCompar])
        end
        imTHEO_yy_mean = mean(imTHEO_yy,3);        
        imGIRF_yy_mean = mean(imGIRF_yy,3);
               
        % plot for mean
        subplot(2,size(imH1_noDiff,3)+1,1)
        imshow(imTHEO_yy_mean,[]);
        hold on
        title('Avg - yy - b=500 s/mm^2')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,size(imH1_noDiff,3)+1,5)
        imshow(imGIRF_yy_mean,[]);
        hold on
        title('Avg - yy - b=500 s/mm^2')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        
        % ---- save file for ADC test ---
        clear im_sos imTest
        nreps = tmp.nreps;        
        im_sos(:,:,1) = imTHEO_noDiff_mean;na
        im_sos(:,:,2) = imTHEO_xx_mean;
        im_sos(:,:,3) = imTHEO_yy_mean;
        imTest = 'THEO'
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')
        
        clear im_sos imTest
        im_sos(:,:,1) = imGIRF_noDiff_mean;
        im_sos(:,:,2) = imGIRF_xx_mean;
        im_sos(:,:,3) = imGIRF_yy_mean;
        imTest = 'GIRF'
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')
        
        
    else
        imTHEO        = mean(imTHEO_orig,3);
        im_BOEC_H1    = mean(imH1_orig,3);
        im_BOEC_H0_H1 = mean(imH0_H1_orig,3);
                
        figure()
        subplot(3,1,1)
        imshow(imTHEO(),[]);
        hold on
        title('Avg')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(3,1,2)
        imshow(im_BOEC_H1(),[]);
        hold on
        title('Avg')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        subplot(3,1,3)
        imshow(im_BOEC_H0_H1(),[]);
        hold on
        title('Avg')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        
        % --- save for ADC test ---
        clear im_sos imTest
        im_sos = imTHEO;
        nreps = tmp.nreps;
        imTest = 'THEO';
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')
        clear im_sos imTest
        im_sos = im_BOEC_H1;
        imTest = 'B0EC_H1';
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')   
        clear im_sos imTest
        im_sos = im_BOEC_H0_H1;
        imTest = 'B0EC_H0_H1';
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')          
    end
end

if Diffusion == 'Fals'
    % test for H0&H1 from BO_ECC    
    imDiff_H1     = zeros(size(im_BOEC_H1));
    imDiff_H1_ABS = zeros(size(im_BOEC_H1));
    figure()
    for ii=1:size(im_BOEC_H1,3)
        imDiff_H1(:,:,ii) = imTHEO(:,:,ii) - im_BOEC_H1(:,:,ii);
        imDiff_H1_ABS(:,:,ii) = imDiff_H1(:,:,ii) / max(max(imDiff_H1(:,:,ii)));
        
        % Simple plot
        subplot(1,size(im_BOEC_H1,3),ii)
        imshow(imDiff_H1_ABS(:,:,ii),[]);
        hold on
        title(['Diff between Recon Theor - Recon B0_E_C_C H1 ',num2str(ii)])
        caxis([-1 1])
    end
    
    % test for H0&H1 from BO_ECC
    imDiff_H0_H1     = zeros(size(im_BOEC_H1));
    imDiff_H0_H1_ABS = zeros(size(im_BOEC_H1));
    figure()
    for ii=1:size(im_BOEC_H1,3)
        imDiff_H0_H1(:,:,ii) = imTHEO(:,:,ii) - im_BOEC_H0_H1(:,:,ii);
        imDiff_H0_H1_ABS(:,:,ii) = imDiff_H0_H1(:,:,ii) / max(max(imDiff_H0_H1(:,:,ii)));
        
        % Simple plot
        subplot(1,size(im_BOEC_H0_H1,3),ii)
        imshow(imDiff_H0_H1_ABS(:,:,ii),[]);
        hold on
        title(['Diff between Recon Theor - Recon B0_E_C_C H0&H1 ',num2str(ii)])
        caxis([-1 1])
    end
    
        % test for H0&H1 from H1
    imDiff_H0_H1     = zeros(size(im_BOEC_H1));
    imDiff_H0_H1_ABS = zeros(size(im_BOEC_H0_H1));
    figure()
    for ii=1:size(im_BOEC_H1,3)
        imDiff_H0_H1_and_H1(:,:,ii)     = im_BOEC_H0_H1(:,:,ii) - im_BOEC_H1(:,:,ii);
        imDiff_H0_H1_and_H1_ABS(:,:,ii) = imDiff_H0_H1_and_H1(:,:,ii) / max(max(imDiff_H0_H1_and_H1(:,:,ii)));
        maxDif_H0_H1_and_H1 = max(max(imDiff_H0_H1_and_H1))
        % Simple plot
        subplot(1,size(im_BOEC_H1,3),ii)
        imshow(imDiff_H0_H1_and_H1(:,:,ii),[]);
        hold on
        title(['Diff between Recon Theor - Recon B0_E_C_C H0&H1 ',num2str(ii)])
        caxis([-0.1 0.1])
    end
    
end
fprintf('\n\n 5 - Sucessfully finished \n\n')

%% part6 - Analysis of Image Profile

% 5.1 - slice lines
x_im(1,:)     = [size(imTHEO,1)/10*7  size(imTHEO,1)/10*7];
x_im(2,:)     = [size(imTHEO,1)/10*3  size(imTHEO,1)/10*7];
x_im(3,:)     = [size(imTHEO,1)/10*0  size(imTHEO,1)/10*10];

y_im(1,:)     = [size(imTHEO,1)/10*0 size(imTHEO,1)/10*10];
y_im(2,:)     = [size(imTHEO,1)/10*0 size(imTHEO,1)/10*10];
y_im(3,:)     = [size(imTHEO,1)/10*4 size(imTHEO,1)/10*4];

n = size(imTHEO,1);  %size(Norm_im_crop,1);

% 5.2 - plots
for ii = 1:size(im_BOEC_H1,3)
    for ji=1:size(y_im,1)
        imProf_H1(:,ji,ii)    = improfile(im_BOEC_H1(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_H0_H1(:,ji,ii) = improfile(im_BOEC_H0_H1(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_THEO(:,ji,ii)  = improfile(imTHEO(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        if tmp.nreps==1
            figure(100*ii+ji)
        else
            figure(100+ji)
        end
        subplot(231)
        plot(imProf_H1(:,ji,ii),'g')
        hold on
        xlabel('Nx')
        ylabel('Pixel Intensity')
        title(['ImProfile_x_x: ', num2str(x_im(ji,1)),', - B0_E_C_C H1'])
        ylim([0 maxValueCompar])
        subplot(232)
        plot(imProf_H0_H1(:,ji,ii),'g')
        hold on
        xlabel('Nx')
        ylabel('Pixel Intensity')
        title(['ImProfile_x_x: ', num2str(x_im(ji,1)),', - B0_E_C_C H0&H1'])
        ylim([0 maxValueCompar])
        subplot(233)
        plot(imProf_THEO(:,ji,ii),'r')
        hold on
        xlabel('Nx')
        ylabel('Pixel Intensity')
        title(['ImProfile_x_x: ', num2str(x_im(ji,1)),', - THEOR'])
        ylim([0 maxValueCompar])
        subplot(234)
        title('Image Recon w/H')
        imshow(im_BOEC_H1(:,:,ii),[]);
        hold on
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
        subplot(235)
        title('Image Recon w/H')
        imshow(im_BOEC_H0_H1(:,:,ii),[]);
        hold on
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
        subplot(236)
        title('Image Recon w/H')
        imshow(imTHEO(:,:,ii),[]);
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
    end
end


% 2nd variaty - plot 
for ii = 1:size(im_BOEC_H1,3)
    for ji=1:size(y_im,1)
        imProf_H1(:,ji,ii) = improfile(im_BOEC_H1(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_THEO(:,ji,ii) = improfile(imTHEO(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        if tmp.nreps==1
            figure(100*ii+ji)
        else
            figure(100+ji)
        end
        subplot(141)
        plot(imProf_H1(:,ji,ii),'g')
        hold on
        plot(imProf_H0_H1(:,ji,ii),'b')
        hold on
        xlabel('Vector Points')
        ylabel('Pixel Intensity')
        plot(imProf_THEO(:,ji,ii),'r')
        hold on
        title(['ImProf - recon H1 (g), H0&H1 (b), recon THEOR (r)'])
        ylim([0 maxValueCompar])
        subplot(142)
        imshow(im_BOEC_H1(:,:,ii),[]);
        hold on
%         title('Image Recon B0_E_C_C H1')        
        title('Image Recon Nominal')        
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
        subplot(143)
        imshow(im_BOEC_H0_H1(:,:,ii),[]);
        hold on
        title('Image Recon B0_E_C_C H0&H1')        
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
        subplot(144)
        imshow(imTHEO(:,:,ii),[]);
        hold on
        title('Image Recon THEOR')
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
    end
end


%% 7 - Numerical analysis

whichcols = 2;                              % how many bourders do I want
numTraces = size(x_im,1)*whichcols;         % Number of boundaries values   
c_H1      = zeros(round(n/2),size(whichcols,2));
c_THEO    = zeros(round(n/2),size(whichcols,2));

FOV = 200;               % mm
Nx  = double(tmp.Nx);
res = FOV/Nx;            % in mm

for ii = 1:size(im_BOEC_H1,3)    
    for ji=1:size(y_im,1)
        bounds = sharp_imProfile(imProf_H1(:,ji,ii),'Fals');
        resSHARP_H1(1,ji,ii) = (bounds.firstUP   - bounds.firstLOW ) * res; % in mm
        resSHARP_H1(2,ji,ii) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds

        bounds = sharp_imProfile(imProf_H0_H1(:,ji,ii),'Fals');
        resSHARP_H0_H1(1,ji,ii) = (bounds.firstUP   - bounds.firstLOW ) * res; % in mm
        resSHARP_H0_H1(2,ji,ii) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds
        
        bounds = sharp_imProfile(imProf_THEO(:,ji,ii),'Fals');
        resSHARP_THEO(1,ji,ii) = (bounds.firstUP  - bounds.firstLOW ) * res; % in mm
        resSHARP_THEO(2,ji,ii) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds
    end

    
    % 6.1 - plots
    
    figure()
    maxCaxis = max([max(max(resSHARP_H1(:,:,ii))) max(max(resSHARP_H0_H1(:,:,ii))) max(max(resSHARP_THEO(:,:,ii)))]);
% %     maxCaxis = max([max(max(resSHARP_H1(:,:,ii))) max(max(resSHARP_THEO(:,:,ii)))]);
% %     vectY = [1 ; 2];
% %     vectX = [1 ; 2 ; 3];
    subplot(311)
    imshow(resSHARP_H1(:,:,ii),[])
    hold on
    title(['Res. Sharpness recon H1'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_H1,1)*2)
    xticks(1:1:size(resSHARP_H1,2)*2)
% %     yticklabels(vectY)
% %     xticklabels(vectX)
    reconH1_mean   = mean(mean(resSHARP_H1))
    maxresSHARP_H1 = max(max(resSHARP_H1))
    minresSHARP_H1 = min(min(resSHARP_H1))    
    
    subplot(312)
    imshow(resSHARP_H0_H1(:,:,ii),[])
    hold on
    title(['Res. Sharpness recon H0&H1'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_H0_H1,1)*2)
    xticks(1:1:size(resSHARP_H0_H1,2)*2)
% %     yticklabels(vectY)
% %     xticklabels(vectX)
    reconH0_H1_mean =mean(mean(resSHARP_H0_H1))
    
    subplot(313)
    imshow(resSHARP_THEO(:,:,ii),[])
    hold on
    title(['Res. Sharpness recon THEO'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_THEO,1)*2)
    xticks(1:1:size(resSHARP_THEO,2)*2)
    yticklabels(['1st Boarder' , '2nd Boarder'])
    xticklabels(['1st cut','2nd cut','3rd cut'])
    
    reconTheo_mean = mean(mean(resSHARP_THEO))
    maxresSHARP_THEO = max(max(resSHARP_THEO))
    minresSHARP_THEO = min(min(resSHARP_THEO))    
end

toc