%% Reconstruct Spiral with GIRF or B0_ECC
%  by TTFernandes, Out 2020



%% part 0 - Addpath & Get Folders
if exist('test_GIRF_THEO','var') == 1
    fprintf('\n\n 0 - Sucessfully finished \n\n')
else
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
end


%% part 0.5 - Set parameters
% --- 0.5.1 - tests  ---
if exist('test_GIRF_THEO','var') == 1
    fprintf('\n\n 0.5 - Sucessfully finished \n\n')
else
    plotTest2 = 'True';
    testInv   = 'AboS';              % 'Freq' or 'Time' or 'AboS'
    testKtraj = 'GIRF';              % Correct theoretical ktrajct - GIRF or B0EC or ELSE (no test)
    
    test      = 64;                  % test to get IRF
end


%% part 1 - Get Files
if exist('test_GIRF_THEO','var') == 1
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


%% part 2 - Define Parametres

% --- 2.1 - parameters ---
Nshots       = double(tmp.Nshots);          % number of interleaves - number of trajectories (SPIRAL)
alpha        = double(tmp.alpha);           % Oversampling - greater than 1 greater the density in center than outside the spiral. <1 the other way around (SPIRAL)
N            = double(tmp.Nx);              % matrix size (SPIRAL)
nsl          = tmp.nslices;                 % number of slices
Nadc         = tmp.adc_samples;             % number of ADCs
dt           = 1e-5;                        % dt (in seconds)
gamma        = 42.576e6;                    % Gyro-magnetic ratio Hz/T


if testKtraj == 'GIRF'
    if tmp.ndirs >=1
        Diffusion = 'True';
    else
        Diffusion = 'Fals';
    end
else
    Diffusion = 'Fals';
end

% Get Gradients
Gx      = tmpMATL.vectSpiral.data(1:Nadc,1); % [Hz/m]
Gy      = tmpMATL.vectSpiral.data(1:Nadc,2); % [Hz/m]
if Diffusion == 'True'
    sizeSpiral = size(Gx,1); 
    nP_before  = 2400;
    pypulseqSPIRAL_seqRecon 
    clear Gx Gy
    
    Gx = gradForH_x';
    Gy = gradForH_y';
end
% % figure()
% % subplot(211)
% % plot(Gx/gamma)
% % subplot(212)
% % plot(Gy/gamma,'r')

fprintf('\n\n 2 - Sucessfully finished \n\n')


%% part 3 - Get Theoric & Experimental k-Trajectories

% --- 3.1 - load IRF function ---
cdResults = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/2_GIRF_test/Test/Results_IRF'];
cd(cdResults)
load(['test',num2str(test),'_IRF.mat'])
Hx      = H.Hx_w(:,1);
Hy      = H.Hy_w(:,1);
Npulses = H.pulses;
clear H

% --- 3.2 - Obtain Experimental k-space trajectory ---
if testKtraj == 'GIRF' %test according Vannesjo et al, 2013   
    
    
    % ================== 3.2.1 - Test w/ conv. in freq =====================
    if testInv == 'Freq'
        % 3.2.3.1 obtain time domain of Hx & Hy
        Ht_x = ifft(Hx);
        Ht_y = ifft(Hy);
    
        % 3.2.4.2 - Compute the length of the convolution
        n_H    = size(Hx,1);
        n_Gx   = size(Gx,1);
        n_conv = n_Gx + n_H - 1;        

        % 3.2.4.3 - obtain FFT of Hx & Hy (size of conv.) - The FFT will be padded with zeros to match the length of the convolution       
        Hx_w = fft(Ht_x,n_conv);
        Hy_w = fft(Ht_y,n_conv);
        
        % 3.2.4.4 - obtain FFT of kx & ky - get Ix & Iy (size of conv.)
        Ix_w = fft(Gx,n_conv);
        Iy_w = fft(Gy,n_conv);
        
        % 3.2.4.5 - obtained measured k_traj in freq_space - get Ox & Oy
        Ox_w = Ix_w .* (abs(Hx_w))/max(abs(Hx_w));
        Oy_w = Iy_w .* (abs(Hy_w))/max(abs(Hy_w));
        if plotTest2 == 'True'
            figure()
            subplot(131)
            plot(abs(Hx_w)/max(abs(Hx_w)))
            hold on
            plot(abs(Hy_w)/max(abs(Hy_w)))
            title('Hw')
            subplot(132)
            plot(abs(Ix_w))
            hold on
            plot(abs(Iy_w))
            title('Iw')
            subplot(133)
            plot(abs(Ox_w))
            hold on
            plot(abs(Oy_w))
            title('Ow')
        end
        
        % 3.2.4.6 - Perform time-domain convolution by frequency domain multiplication
        Ox_t = ifft(Ox_w, n_conv)./n_conv;
        Oy_t = ifft(Oy_w, n_conv)./n_conv;
        
        O_t = Ox_t + Oy_t;
        Ox_test = real((O_t));
        Oy_test = imag((O_t));
        % plots
        % % %     figure()
        % % %     plot(abs((Ox_t)))
        % % %     hold on
        % % %     plot(abs((Oy_% --- 3.1 - load IRF function ---
        cdResults = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/2_GIRF_test/Test/Results_IRF'];
        cd(cdResults)
        load(['test',num2str(test),'_IRF.mat'])
        Hx      = H.Hx_w(:,1);
        Hy      = H.Hy_w(:,1);
        Npulses = H.pulses;
        clear Ht)))
        
        if plotTest2 == 'True'
            figure()
            subplot(221)
            plot(Gx)
            hold on
            title('Iy_t')
            subplot(222)
            plot(Gy,'r')
            hold on
            title('Iy_t')
            subplot(223)
            plot(Ox_test)
            hold on
            title('Ox_t')
            subplot(224)
            plot(Oy_test,'r')
            title('Oy_t')
        end
        
        % Truncate the result of convolution to obtain a time-series that is the same length as the original signal
        Ox = Ox_test(floor(n_H/2) : end - floor(n_H/2) - 1,:);
        Oy = Oy_test(floor(n_H/2) : end - floor(n_H/2) - 1,:);
        %     clear Ox_test Oy_test
        
        % Obtain Output% --- 3.1 - load IRF function ---
        cdResults = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/2_GIRF_test/Test/Results_IRF'];
        cd(cdResults)
        load(['test',num2str(test),'_IRF.mat'])
        Hx      = H.Hx_w(:,1);
        Hy      = H.Hy_w(:,1);
        Npulses = H.pulses;
        clear H k-space trajectory from the gradients (O_t)
        int_mat = triu(ones(Nadc-1));
        O_kx = zeros(Nshots,Nadc*Nshots-1);
        O_ky = zeros(Nshots,Nadc*Nshots-1);
        for iNs=1:Nshots
            O_kx = Ox(:,iNs)' * int_mat * dt; % vem em unidades gamma
            O_ky = Oy(:,iNs)' * int_mat * dt;
        end
        O_kx    = O_kx(Nshots,1:end-1)';
        O_ky    = O_ky(Nshots,1:end-1)';
        O_ktraj = O_kx + 1i*O_ky;
        
        % plot comparison between Ktraj Theoric & Ktraj Output
        if plotTest2 == 'True'
            figure()
            subplot(121)
            plot(O_ktraj)
            hold on
            subplot(122)
            plot(ktraj,'r')
        end
    
    
    % ================== 3.2.2 - Test w/ conv. in time ====================
    elseif testInv == 'Time'        
            % Normalized Hx & Hy
            norm_Ht_x = Ht_x;
            norm_Ht_y = Ht_y;
            
            % Obtain O(t) according to Vannesjmo, 2013  - in time domain
            Ot_kx = conv(Gx,norm_Ht_x);
            Ot_ky = conv(Gy,norm_Ht_y);
            %         O_ktraj = O_kx + 1i*O_ky;
            
            n_Ht = size(Ht_x,1);
            Ot_kx = Ot_kx(floor(n_Ht/2) : end - floor(n_Ht/2) - 1,:);
            Ot_ky = Ot_ky(floor(n_Ht/2) : end - floor(n_Ht/2) - 1,:);
    
    % ================== 3.2.3 - Test w/ Abo Seaba Implementation =========
    elseif testInv == 'AboS'        
        Npad    = 0;           
        T       = linspace(0.05,0.16,Npulses);             % time range - (ms)            
        maxFreq = 1./T(1)*1e3;
        ff      = linspace(-maxFreq,maxFreq,size(Hx,1));   % Frequency in Hz
        
        for ii=1:size(Gx,2)
            Ginput = Gx(:,ii);       % input in time domain
            [Ot_x(:,ii),~,~] = gradient_distort_GIRF(Ginput,ff,Hx,dt,Npad);
            clear Ginput
        end
        for ii=1:size(Gy,2)
            Ginput = Gy(:,ii);       % input in time domain
            [Ot_y(:,ii),~,~] = gradient_distort_GIRF(Ginput,ff,Hy,dt,Npad);
            clear Ginput
        end    
        
        % figures of G theor & G measured
        if plotTest2 == 'True' 
            timeV = [0:dt:(size(Gx,1)-1)*dt]*1e3;
            figure()
            subplot(211)
            plot(timeV,Gx(:,1)/gamma*1e3,'b')
            hold on
            plot(timeV,Ot_x(:,1)/gamma*1e3,'g')
            hold on
            ylabel(['Gradient Amplitude (mT/m)'])
            xlabel(['time (ms)'])            
            title('Gx THEOR (b) & Gx GIRF (g)')
            subplot(212)
            plot(timeV,Gy(:,1)/gamma*1e3,'b')
            hold on
            plot(timeV,Ot_y(:,1)/gamma*1e3,'g')
            hold on
            ylabel(['Gradient Amplitude (mT/m)'])
            xlabel(['time (ms)'])            
            title('Gy THEOR (b) & Gy GIRF (g)')
        end
    end
    
    
    % --- 3.3 - Obtain Experimental k-space trajectory ---
    if Diffusion == 'True'
        % Get original Gradients
        Gx      = tmpMATL.vectSpiral.data(1:Nadc,1); % [Hz/m]
        Gy      = tmpMATL.vectSpiral.data(1:Nadc,2); % [Hz/m]
        
        Ot_x_aux = Ot_x(nP_before:end-1,:)*1e3;
        Ot_y_aux = Ot_y(nP_before:end-1,:)*1e3;        
                
        Ot_x = Ot_x_aux;
        Ot_y = Ot_y_aux;   
        
        
% %         figure()
% %         plot(Gx,'b')
% %         hold on
% %         plot(Ot_x(:,2),'r')
    end
    
    
    int_mat = triu(ones(Nadc));
    % Have to split the integral because each of the gradients is independent.
    % Otherwise we would have a continuous integral
    for ii=1:size(Ot_x,2)        
        O_kx_aux = Ot_x(:,ii)'  * int_mat * dt; % vem em unidades gamma
        O_ky_aux = Ot_y(:,ii)'  * int_mat * dt;
        O_kx(ii,:) = O_kx_aux(:,1:end-2);
        O_ky(ii,:) = O_ky_aux(:,1:end-2);
        clear O_kx_aux O_ky_aux
        
        O_ktraj(ii,:) = O_kx(ii,:) + 1i*O_ky(ii,:);   % final K-space measured
    end
    
    ntrials = size(O_ktraj,1); % set number of trials for recon
    
    % Theorical K_traj
    for ii=1:size(Gx,2)
        kx_aux(ii,:) = Gx(:,ii)'  * int_mat * dt;    % vem em unidades gamma
        ky_aux(ii,:) = Gy(:,ii)'  * int_mat * dt;
        kx(ii,:) = kx_aux(ii,1:end-2)';
        ky(ii,:) = ky_aux(ii,1:end-2)';
        clear kx_aux ky_aux
        
        ktraj(ii,:) = kx(ii,:) + 1i*ky(ii,:);
    end
    
    
    % plot comparison between Ktraj Theoric & Ktraj Output
    if plotTest2 == 'True'
        figure()
        plot(ktraj,'b')
        hold on
        plot(O_ktraj(2,:),'r')
        title('Ktraj theor (b) & Ktraj Measured (r)')
        ylabel(['1/m'])
        xlabel(['1/m'])
    end
    
    % obtain final Ktraj, kx and ky
    clear ktraj kx ky
    ktraj = O_ktraj;
    kx = O_kx';
    ky = O_ky';
    
elseif testKtraj == 'B0EC' %test according to Robison et al 2019, 

    
    
elseif testKtraj == 'ELSE'     % --- 3.1 - Obtain theoretic k-space trajectory from the gradients ---   
    int_mat = triu(ones(Nadc));

    % Have to split the integral because each of the gradients is independent.
    % Otherwise we would have a continuous integral
    kx = Gx' * int_mat * dt; % vem em unidades gamma
    ky = Gy' * int_mat * dt;
    ntrials = size(kx,1); % set number of trials for recon
    
    kx    = kx(1,1:end-2)';
    ky    = ky(1,1:end-2)';
    
    ktraj = kx + 1i*ky;
end

clc
fprintf('\n\n 3 - Sucessfully finished \n\n')


%% part 4 - Read data
cd(folder)
kdata=dat2mat(fn);                      % [Nadc, nCoils, nsl&nbv]
nimgs=floor(size(kdata,3)/nsl/Nshots);  % Number of images acquired, counting the calibration.
kdata=double(reshape(kdata(:,:,1:nsl*nimgs*Nshots),[size(kdata,1)*Nshots size(kdata,2) nsl nimgs])); % Reshape info data

nc=size(kdata,2); % Number of coils
im=zeros([N N nc nsl nimgs]);
size(kdata)

fprintf('\n\n 4 - Sucessfully finished \n\n')


% -- adjust for pypulseq ADC
% % % kx = kx(1:size(kx,1)-4+1,:);
% % % ky = ky(1:size(ky,1)-4+1,:);

kx = kx(1:size(kdata,1),:);
ky = ky(1:size(kdata,1),:);


%% part 5 - Reconstruction throw the k-space with different parameters

% part 5.1 - preCompensation + Gridding + PosCompensation + iFFT + plot

% ... preparing gridding ...
gridsize = double(N*2);

% controlar para que o meu valor do espaço-k não seja maior do que 0.5
% Requirement of voronoi.... ??
for ji=1:ntrials
    if max(kx(:,ji))>0.5 | max(ky(:,ji))>0.5 % Normalize
        nkx(:,ji)=(kx(:,ji)*0.5)/max(kx(:,ji));
        nky(:,ji)=(ky(:,ji)*0.5)/max(ky(:,ji));
    else
        nkx(:,ji) = kx(:,ji);
        nky(:,ji) = ky(:,ji);
    end
    
    nkTraj(:,:,ji) = [nkx(:,ji) nky(:,ji)];
    ntraj(:,ji)  = (nkx(:,ji)+nky(:,ji)*1i); % Normalized k-space trajectory
    
    % ... selecao de parametros ...
    overgridfactor = 2;
    kwidth         = 3; % kwidth = 50;  % kernel width
    kbeta          = pi*sqrt(kwidth^2/overgridfactor^2*(overgridfactor-0.5)^2-0.8 );  % kernel beta
    
    % ... pre-compensation weights (independent of simulated data) ...
    [acf(:,ji),flag(:,ji),res(:,ji)] = makePrecompWeights_2D_VORONOI( nkTraj(:,:,ji), [gridsize/2 gridsize/2] );
end

for ii=1:nimgs
    for c=1:nc
        for sl=1:nsl
            if Diffusion == 'True'
                ji = ii;
            else
                ji = 1;
            end
            % ... keiser-bessel gridding + fft ...
            %     trimming       -- 'y' or 'n'
            %     apodize        -- 'y' or 'n'
            %     postcompensate -- 'y' or 'n'
            [m_postcomp,kwx,kwy] = AuxSpiralTF_grid_kb(kdata(:,c,sl,ii),ntraj(:,ji)',acf(:,ji),gridsize,overgridfactor,kwidth,kbeta,'n','n','y');
            im_large             = fftshift(fft2(ifftshift(m_postcomp)));

            % ... crop figures ...
            m_postcomp_crop = m_postcomp(gridsize-round(N/2)+1:gridsize+round(N/2),gridsize-round(N/2)+1:gridsize+round(N/2));
            im(:, :, c, sl, ii) = im_large(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2)-1,round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2)-1);

        end
    end
end
im_sos=sqrt(squeeze(sum(im.*conj(im),3)));

fprintf('       .... fim da reconstrucao da imagem - Gridding Pre & Pos Compensation !   ....    \n');



fprintf('\n\n 5 - Sucessfully finished \n\n')

%% part 6 - Plot images

% Simple plot
maxValue = max(max(max(im_sos)));
figure()
for ii=1:size(im_sos,3)
    subplot(1,size(im_sos,3),ii)
    imshow(im_sos(:,:,ii),[]);
    hold on
    title(['Rep ',num2str(ii)])
    caxis([0 maxValue])    
end
fprintf('\n\n 6 - Sucessfully finished \n\n')

% close all

%% part 7 - save data

if testKtraj == 'GIRF'
    if Diffusion == 'True'
        imTest = 'GIRF_Diff';
    else
        imTest = 'GIRF';
    end
elseif testKtraj == 'ELSE'
    imTest = 'Theor';    
end

cd(file_folder_results)
nameSave = ['test_',num2str(tmp.test),'_imRecon_',imTest,'.mat'];
delete(nameSave)
save(nameSave,'im_sos')


