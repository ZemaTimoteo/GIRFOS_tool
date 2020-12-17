%% Run for each test of SCAN the reconstruction with Theoric ktraj and GIRF ktraj
%   - uses recon_spiral_IRF.m

%  TTFernandes, Dez 2020

%% part 0 - Addpath & Get Folders
clear all;
computer = 2; % 1 - MyPC OR 2 - Seia PC
if computer == 1
    cd('D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\4_Reconstruction')
    addpath('mapVBVD\') % To read
    addpath(genpath('Aux_Functions'));  % To reconstruct
else
    cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/3_Reconstruction_Code')
    addpath('mapVBVD/') % To read
    addpath(genpath('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/bart-master')) % To read
    addpath(genpath('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode'));
    addpath(genpath('Aux_Functions'));  % To reconstruct
end
% close all
clc


%% part 1 - Get Files
if exist('test_GIRF_&_THEO','var') == 1
    fprintf('\n\n 1 - Sucessfully finished \n\n')
else
    if computer == 1
        folder = 'D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\3_Data_scan\SPIRAL_measurement\Test\';
        folder = [folder, 'SPIRAL_test_3\'];
    else
        %     folder = '/home/tfernandes/Documents/Projetos/Project_lfStroke/preProcResults/Data/SPIRAL_measurement/';
        %     folder = [folder, 'SPIRAL_test_3/'];
        myCD = '/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/1_Spiral_test/Test';
        cd(myCD)
        str ={'select file'}
        s = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str)
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
    end
    
    [files]   = dir(folder);
    fn        = files(5).name;       % get '.dat' file
    fn_info   = files(6).name;       % get '.mat' file with sequence parameters from Python
    tmp       = load(fn_info);       % get values from '.mat' file with sequence parameters
    
    cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/Spiral/pulses')
    name_seq_template = ['vectSpiral_test',num2str(tmp.testSPI),'.mat'];
    % name_seq_template = ['vectSpiral.mat'];
    
    copyfile(name_seq_template,folder)
    cd(folder)
    
    fn_MATLAB = files(8).name;       % get '.mat' file with sequence parameters from MATLAB
    
    tmpMATL   = load(fn_MATLAB);
    clear strTest1 filename auxA fn_info fn_MATLAB
    
    fprintf('\n\n 1 - Sucessfully finished \n\n')
end

%% part 2 - select recon with Theoric ktraj with recon_spiral_IRF

plotTest2 = 'Fals';
testInv   = 'AboS';              % 'Freq' or 'Time' or 'AboS'
testKtraj = 'ELSE';              % Correct theoretical ktrajct - GIRF or B0EC or ELSE (no test)
test      = 64;                  % test to get IRF

test_GIRF_THEO = 'True';

recon_spiral_IRF


%% part 3 - select recon with GIRF correct ktraj recon_spiral_IRF
clearvars -except tmpMATL  name_seq_template tmp fn file_folder_results ...
    file_folder_cfl_hdr myCD test testInv plotTest2 test_GIRF_THEO folder

testKtraj      = 'GIRF';              % Correct theoretical ktrajct - GIRF or B0EC or ELSE (no test)

recon_spiral_IRF
    
%% part 4 - obtain differences

close all

cd(file_folder_results)
clear im_sos imTHEO
load(['test_',num2str(tmp.test),'_imRecon_Theor.mat']);
imTHEO = im_sos;
clear imGIRF im_sos
if Diffusion == 'True'
    load(['test_',num2str(tmp.test),'_imRecon_GIRF_Diff.mat']);
    imGIRF = im_sos;    
else
    load(['test_',num2str(tmp.test),'_imRecon_GIRF.mat']);
    imGIRF = im_sos;  
end

maxValueCompar = max([max(max(max(imGIRF))) max(max(max(imTHEO)))]);
figure()
for ii=1:size(imTHEO,3)
    subplot(2,size(imTHEO,3),ii)
    imshow(imTHEO(:,:,ii),[]);
    hold on
    title(['Rep ',num2str(ii)])
    ylabel(['Theor'])
    caxis([0 maxValueCompar])
    subplot(2,size(imGIRF,3),ii+size(imGIRF,3))
    imshow(imGIRF(:,:,ii),[]);
    hold on
    title(['Rep ',num2str(ii)])
    ylabel(['GIRF'])
    caxis([0 maxValueCompar])        
end


if tmp.nreps >1
    imTHEO_orig = imTHEO;
    imGIRF_orig = imGIRF;
% %     imTHEO = imGIRF_orig;
% %     imGIRF = imGIRF_orig;
    clear imTHEO imGIRF
    if Diffusion == 'True'        
        % ---------------- b=0 -------------------
        imTHEO_noDiff = imTHEO_orig(:,:,2:tmp.nreps);
        imGIRF_noDiff = imGIRF_orig(:,:,2:tmp.nreps);
        
        % all plots b=0
        figure()
        for l=1:size(imTHEO_noDiff,3)
            subplot(2,size(imTHEO_noDiff,3)+1,l+1)
            imshow(imTHEO_noDiff(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l)])
            ylabel(['Theor'])
            caxis([0 maxValueCompar])
            subplot(2,size(imGIRF_noDiff,3)+1,l+2+size(imGIRF_noDiff,3))
            imshow(imGIRF_noDiff(:,:,l),[]);
            hold on
            title(['Rep: ',num2str(l)])
            ylabel(['GIRF'])
            caxis([0 maxValueCompar])
        end
        imTHEO_noDiff_mean = mean(imTHEO_noDiff,3);
        imGIRF_noDiff_mean = mean(imGIRF_noDiff,3);
        
        % plot for mean
        subplot(2,size(imGIRF_noDiff,3)+1,1)
        imshow(imTHEO_noDiff_mean,[]);
        hold on
        title('Avg - b=0')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,size(imGIRF_noDiff,3)+1,5)
        imshow(imGIRF_noDiff_mean,[]);
        hold on
        title('Avg - b=0')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        
        
        % ---------------- b=500 -------------------
        imTHEO_xx = imTHEO_orig(:,:,tmp.nreps+1:2:tmp.nreps*(tmp.ndirs+1));
        imTHEO_yy = imTHEO_orig(:,:,tmp.nreps+2:2:tmp.nreps*(tmp.ndirs+1));   
        imGIRF_xx = imGIRF_orig(:,:,tmp.nreps+1:2:tmp.nreps*(tmp.ndirs+1));
        imGIRF_yy = imGIRF_orig(:,:,tmp.nreps+2:2:tmp.nreps*(tmp.ndirs+1));        
        
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
        subplot(2,size(imGIRF_noDiff,3)+1,1)
        imshow(imTHEO_xx_mean,[]);
        hold on
        title('Avg - b=500 s/mm^2')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,size(imGIRF_noDiff,3)+1,5)
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
        subplot(2,size(imGIRF_noDiff,3)+1,1)
        imshow(imTHEO_yy_mean,[]);
        hold on
        title('Avg - yy - b=500 s/mm^2')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,size(imGIRF_noDiff,3)+1,5)
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
        imTHEO = mean(imTHEO_orig,3);
        imGIRF = mean(imGIRF_orig,3);
                
        figure()
        subplot(2,1,1)
        imshow(imTHEO(),[]);
        hold on
        title('Avg')
        ylabel(['Theor'])
        caxis([0 maxValueCompar])
        subplot(2,1,2)
        imshow(imGIRF(),[]);
        hold on
        title('Avg')
        ylabel(['GIRF'])
        caxis([0 maxValueCompar])
        
        % --- save for ADC test ---
        clear im_sos imTest
        im_sos = imTHEO;
        nreps = tmp.nreps;
        imTest = 'THEO'
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')
        clear im_sos imTest
        im_sos = imGIRF;
        imTest = 'GIRF'
        cd(file_folder_results)
        nameSave = ['test_',num2str(tmp.test),'_ADC_mean_imRecon_',imTest,'.mat'];
        delete(nameSave)
        save(nameSave,'im_sos','nreps')        
    end
end

if Diffusion == 'Fals'
    imDiff     = zeros(size(imGIRF));
    imDiff_ABS = zeros(size(imGIRF));
    figure()
    for ii=1:size(imGIRF,3)
        imDiff(:,:,ii) = imTHEO(:,:,ii) - imGIRF(:,:,ii);
        imDiff_ABS(:,:,ii) = imDiff(:,:,ii) / max(max(imDiff(:,:,ii)));
        
        % Simple plot
        subplot(1,size(imGIRF,3),ii)
        imshow(imDiff_ABS(:,:,ii),[]);
        hold on
        title(['Diff between Recon Theor - Recon GIRF ',num2str(ii)])
        caxis([-1 1])
    end
end
fprintf('\n\n 7 - Sucessfully finished \n\n')

%% part 5 - Analysis of Image Profile

% 5.1 - slice lines
x_im(1,:)     = [size(imTHEO,1)/10*7  size(imTHEO,1)/10*7];
x_im(2,:)     = [size(imTHEO,1)/10*3  size(imTHEO,1)/10*7];
x_im(3,:)     = [size(imTHEO,1)/10*0  size(imTHEO,1)/10*10];

y_im(1,:)     = [size(imTHEO,1)/10*0 size(imTHEO,1)/10*10];
y_im(2,:)     = [size(imTHEO,1)/10*0 size(imTHEO,1)/10*10];
y_im(3,:)     = [size(imTHEO,1)/10*4 size(imTHEO,1)/10*4];

n = size(imTHEO,1);  %size(Norm_im_crop,1);

% 5.2 - plots
for ii = 1:size(imGIRF,3)
    for ji=1:size(y_im,1)
        imProf_GIRF(:,ji,ii) = improfile(imGIRF(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_THEO(:,ji,ii) = improfile(imTHEO(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        if tmp.nreps==1
            figure(100*ii+ji)
        else
            figure(100+ji)
        end
        subplot(221)
        plot(imProf_GIRF(:,ji,ii),'g')
        hold on
        xlabel('Nx')
        ylabel('Pixel Intensity')
        title(['ImProfile_x_x: ', num2str(x_im(ji,1)),', - GIRF'])
        ylim([0 maxValueCompar])
        subplot(222)
        plot(imProf_THEO(:,ji,ii),'r')
        hold on
        xlabel('Nx')
        ylabel('Pixel Intensity')
        title(['ImProfile_x_x: ', num2str(x_im(ji,1)),', - THEOR'])
        ylim([0 maxValueCompar])
        subplot(223)
        title('Image Recon w/H')
        imshow(imGIRF(:,:,ii),[]);
        hold on
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
        subplot(224)
        title('Image Recon w/H')
        imshow(imTHEO(:,:,ii),[]);
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])

    end
end


% 2nd variaty - plot 
for ii = 1:size(imGIRF,3)
    for ji=1:size(y_im,1)
        imProf_GIRF(:,ji,ii) = improfile(imGIRF(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        imProf_THEO(:,ji,ii) = improfile(imTHEO(:,:,ii),x_im(ji,:),y_im(ji,:),n);
        if tmp.nreps==1
            figure(100*ii+ji)
        else
            figure(100+ji)
        end
        subplot(131)
        plot(imProf_GIRF(:,ji,ii),'g')
        hold on
        xlabel('Vector Points')
        ylabel('Pixel Intensity')
        plot(imProf_THEO(:,ji,ii),'r')
        hold on
        title(['ImProf - recon GIRF (g) recon THEOR (r)'])
        ylim([0 maxValueCompar])
        subplot(132)
        imshow(imGIRF(:,:,ii),[]);
        hold on
        title('Image Recon GIRF')        
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
        subplot(133)
        imshow(imTHEO(:,:,ii),[]);
        hold on
        title('Image Recon THEOR')
        line(x_im(ji,:),y_im(ji,:))
        caxis([0 maxValueCompar])
    end
end


%% 6 - Numerical analysis

whichcols = 2;                              % how many bourders do I want
numTraces = size(x_im,1)*whichcols;         % Number of boundaries values   
c_GIRF    = zeros(round(n/2),size(whichcols,2));
c_THEO    = zeros(round(n/2),size(whichcols,2));

FOV = 200;               % mm
Nx  = double(tmp.Nx);
res = FOV/Nx;            % in mm

for ii = 1:size(imGIRF,3)    
    for ji=1:size(y_im,1)
        bounds = sharp_imProfile(imProf_GIRF(:,ji,ii),'Fals');
        resSHARP_GIRF(1,ji,ii) = (bounds.firstUP   - bounds.firstLOW ) * res; % in mm
        resSHARP_GIRF(2,ji,ii) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds
        
        bounds = sharp_imProfile(imProf_THEO(:,ji,ii),'Fals');
        resSHARP_THEO(1,ji,ii) = (bounds.firstUP  - bounds.firstLOW ) * res; % in mm
        resSHARP_THEO(2,ji,ii) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds
    end

    
    % 6.1 - plots
    
    figure()
    maxCaxis = max([max(max(resSHARP_GIRF(:,:,ii))) max(max(resSHARP_THEO(:,:,ii)))]);
% %     vectY = [1 ; 2];
% %     vectX = [1 ; 2 ; 3];
    subplot(211)
    imshow(resSHARP_GIRF(:,:,ii),[])
    hold on
    title(['Res. Sharpness recon GIRF'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_GIRF,1)*2)
    xticks(1:1:size(resSHARP_GIRF,2)*2)
% %     yticklabels(vectY)
% %     xticklabels(vectX)
    mean(mean(resSHARP_GIRF))

    subplot(212)
    imshow(resSHARP_THEO(:,:,ii),[])
    hold on
    title(['Res. Sharpness recon THEO'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_THEO,1)*2)
    xticks(1:1:size(resSHARP_THEO,2)*2)
    yticklabels(['1st Boarder' , '2nd Boarder'])
    xticklabels(['1st cut','2nd cut','3rd cut'])
    
    mean(mean(resSHARP_THEO))
end
