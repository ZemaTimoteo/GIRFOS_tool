%% Reconstruct Spiral with B0_ECC or with Theo
%       Instructions:
%            - run this code after running 'getBOECC.m'
%  by: TTFernandes, Jan 2021
%        @IST_ISR_Lisboa


%% part 0 - Addpath & Get Folders
if exist('test_B0ECC_multi','var') == 1
    fprintf('\n\n 0 - Sucessfully finished \n\n')
else
    clear all;    
    addpath(genpath('GIRFOS_tool/Matlab_basedCode/Tools/'));  % To reconstruct        
    % close all
    clc
end


%% part 0.5 - Set parameters
% --- 0.5.1 - tests  ---
if exist('test_B0ECC_multi','var') == 1
    fprintf('\n\n 0.5 - Sucessfully finished \n\n')
else
    test       = 120;                 % test to get B0_ECC
    testKtraj  = 'B0EC';              % Correct with BO_ECC H1 theoretical ktrajct - B0EC or THEO (no test)
    H0_corr    = 'True';              % Correct with B0_ECC H0 - 'True' or 'Fals'
    testInv    = 'AboS';              % 'Time' or 'AboS'

    saveData    = 'Fals';             % Save raw Data
    saveResults = 'True';             % Save results
    plotTest2   = 'True';
    plotTest4   = 'True';
end


%% part 1 - Get Files
if exist('test_B0ECC_multi','var') == 1
    fprintf('\n\n 0 - Sucessfully finished \n\n')
else
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
end


%% part 2 - Define Parametres

% --- 2.1 - parameters ---
Nshots       = double(tmp.Nshots);          % number of interleaves - number of trajectories (SPIRAL)
alpha        = double(tmp.alpha);           % Oversampling - greater than 1 greater the density in center than outside the spiral. <1 the other way around (SPIRAL)
% N            = 91;%double(tmp.Nx);              % matrix size (SPIRAL)
N            = double(tmp.Nx);              % matrix size (SPIRAL)
nsl          = tmp.nslices;                 % number of slices
Nadc         = tmp.adc_samples;             % number of ADCs
dt           = 1e-5;                        % dt (in seconds)
gamma        = 42.576e6;                    % Gyro-magnetic ratio Hz/T

if tmp.ndirs >=1
    Diffusion = 'True';
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
if PC == 1
cd('GIRFOS_tool/Matlab_basedCode/Example/')
load(['test',num2str(test),'_B0ECC_H.mat'])




% --- 3.2 - Obtain Experimental k-space trajectory throw B0ECC ---
if testKtraj == 'B0EC' %test according Robison et al, 2019
    H1x      = H.H1_w.x(:,1);
    H1y      = H.H1_w.y(:,1);
    Npulses  = H.H0_w.pulses;
    % ================== 3.2.2 - Test w/ conv. in time ====================
    if testInv == 'Time'
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
        % --- parameters ---
        Npad    = 0;
        T       = linspace(0.05,0.16,Npulses);              % time range - (ms)
        maxFreq = 1./T(1)*1e3;
        ff1     = linspace(-maxFreq,maxFreq,size(H1x,1));   % Frequency in Hz
        
        % --- apply H1 & get Gout ---
        for ii=1:size(Gx,2)
            Ginput = Gx(:,ii);       % input in time domain
            [Ot_x(:,ii),~,~] = gradient_distort_GIRF(Ginput,ff1,H1x,dt,Npad);
            clear Ginput
            Ginput = Gy(:,ii);       % input in time domain
            [Ot_y(:,ii),~,~] = gradient_distort_GIRF(Ginput,ff1,H1y,dt,Npad);
            clear Ginput
        end
        
        % --- apply H0 & get Pha_B0 ---        
        if H0_corr == 'True'
            clear phas_B0
            H0x      = H.H0_w.x(:,1);
            H0y      = H.H0_w.y(:,1);            
            
            ff0      = linspace(-maxFreq,maxFreq,size(H0x,1));   % Frequency in Hz            
            for ii=1:size(Gx,2)
                Ginput = Gx(:,ii)/gamma*1e3;       % input in time domain
                [phas_B0x(:,ii),~,~] = gradient_distort_GIRF(Ginput,ff0,H0x,dt,Npad);
                clear Ginput
                Ginput = Gy(:,ii)/gamma*1e3;              % input in time domain

                
                [phas_B0y(:,ii),~,~] = gradient_distort_GIRF(Ginput,ff0,H0y,dt,Npad);
                clear Ginput

                phas_B0(:,ii) = phas_B0x(:,ii) + phas_B0y(:,ii) ; 
                
                if plotTest2 == 'True'
                    timeV = [0:dt:(size(Gx,1)-1)*dt]*1e3;
                    figure()
                    plot(timeV,phas_B0(:,1))
                    hold on
                    ylabel(['phase (rad)'])
                    xlabel(['time (ms)'])
                    title('Phase of B_0')
                end           
            end
        end
        clear H

        % --- Figures of G theor & G measured ---
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
    % ... get diffusion theor gradients ...
    if Diffusion == 'True'
        % Get original Gradients
        Gx      = tmpMATL.vectSpiral.data(1:Nadc,1); % [Hz/m]
        Gy      = tmpMATL.vectSpiral.data(1:Nadc,2); % [Hz/m]
        
        Ot_x_aux = Ot_x(nP_before:end-1,:)*1e3;
        Ot_y_aux = Ot_y(nP_before:end-1,:)*1e3;
        
        Ot_x = Ot_x_aux;
        Ot_y = Ot_y_aux;

    end
    
    % ... integrate Gout ...
    int_mat = triu(ones(Nadc));
    
    % Have to split the integral because each of the gradients is independent.
    % Otherwise we would have a continuous integral
    for ii=1:size(Ot_x,2)
        O_kx_aux   = Ot_x(:,ii)'  * int_mat * dt; % vem em unidades gamma
        O_ky_aux   = Ot_y(:,ii)'  * int_mat * dt;
        O_kx(ii,:) = O_kx_aux(:,1:end-2);
        O_ky(ii,:) = O_ky_aux(:,1:end-2);
        clear O_kx_aux O_ky_aux
        
        O_ktraj(ii,:) = O_kx(ii,:) + 1i*O_ky(ii,:);   % final K-space measured
    end
    
    ntrials = size(O_ktraj,1); % set number of trials for recon
     
    % --- 3.4 - Theorical K_traj ---
    clear kx ky ktraj
    for ii=1:size(Gx,2)
        kx_aux(ii,:) = Gx(:,ii)'  * int_mat * dt;    % vem em unidades gamma
        ky_aux(ii,:) = Gy(:,ii)'  * int_mat * dt;
        kx(ii,:) = kx_aux(ii,1:end-2)';
        ky(ii,:) = ky_aux(ii,1:end-2)';
        clear kx_aux ky_aux
        
        ktraj(ii,:) = kx(ii,:) + 1i*ky(ii,:);
    end
    
    
    % --- 3.5 - Plot comparison: Ktraj Theoric & Ktraj Output ---
    if plotTest2 == 'True'
        figure()
        plot(ktraj,'b')
        hold on
        plot(O_ktraj,'r')
        title('Ktraj theor (b) & Ktraj Measured (r)')
        ylabel(['1/m'])
        xlabel(['1/m'])
    end
    
    % --- 3.6 - Obtain final Ktraj, kx and ky ---
    clear ktraj kx ky
    ktraj = O_ktraj;
    kx = O_kx';
    ky = O_ky';
    
elseif testKtraj == 'THEO'     % --- 3.1 - Obtain theoretic k-space trajectory from the gradients ---
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
% ... 4.1 - Load data in Matlab ...
clear kdata
cd(folder)
kdata=dat2mat(fn);                      % [Nadc, nCoils, nsl&nbv]
nimgs=floor(size(kdata,3)/nsl/Nshots);  % Number of images acquired, counting the calibration.
kdata=double(reshape(kdata(:,:,1:nsl*nimgs*Nshots),[size(kdata,1)*Nshots size(kdata,2) nsl nimgs])); % Reshape info data

nc=size(kdata,2); % Number of coils
im=zeros([N N nc nsl nimgs]);

% ... 4.2 - Adjust for pypulseq ADC ...
kx = kx(1:size(kdata,1),:);
ky = ky(1:size(kdata,1),:);

% ... 4.3 - Correct with Phas_B0 ...
if H0_corr == 'True'
    clear kdata_H0 kdata_HOv2
    phas_B0  = phas_B0(1:size(kdata,1),:);
    kdata_H0 = zeros(size(kdata));
    kdata_HOv2 = zeros(size(kdata));
    
    for nc=1:size(kdata,2) % n coils
        for sli=1:size(kdata,3) % n slices
            for nrep=1:size(kdata,4) % n repetitions
                auxKdat = kdata(:,nc,sli,nrep);
                kdata_H0(:,nc,sli,nrep) = auxKdat .* exp(-1i*phas_B0);
                clear auxKdat
            end
        end
    end
    if plotTest4 == 'True'
        figure()
        subplot(131)
        plot(abs(kdata_H0(:,1,1,1)),'b')
        hold on
        title('method 1')
        subplot(132)
        plot(abs(kdata(:,1,1,1)),'g')
        hold on
        title('original k_data')
        
        figure()
        coil = 4;
        slice = 1;
        rep = 1;
        plot(angle(kdata_H0(:,coil,slice,rep)),'r')
        hold on
        title('method exp(phas) (r) & original k_data (g)')
        plot(angle(kdata(:,coil,slice,rep)),'g')
        
        maxAngl = max((angle(kdata(:,coil,slice,rep))));
        percAngle_Kdata = abs(( angle(kdata_H0(:,coil,slice,rep))-angle(kdata(:,coil,slice,rep)))./ (maxAngl*2) *100);
        figure()
        plot(percAngle_Kdata,'r')
        title('% of Kdata_H0')        
        
    end
    clear kdata
    kdata = kdata_H0;
    clear kdata_H0 kdata_HOv2
end


fprintf('\n\n 4 - Sucessfully finished \n\n')


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
            [m_postcomp,kwx,kwy] = AuxSpiral_grid_kb(kdata(:,c,sl,ii),ntraj(:,ji)',acf(:,ji),gridsize,overgridfactor,kwidth,kbeta,'n','n','y');
            im_large             = fftshift(fft2(ifftshift(m_postcomp)));

            % ... crop figures ...
            m_postcomp_crop = m_postcomp(gridsize-round(N/2)+1:gridsize+round(N/2),gridsize-round(N/2)+1:gridsize+round(N/2));
%             im(:, :, c, sl, ii) = im_large(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2),round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2));
            im(:, :, c, sl, ii) = im_large(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2)-1,round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2)-1);

        end
    end
end
im_sos=sqrt(squeeze(sum(im.*conj(im),3)));

fprintf('       .... fim da reconstrucao da imagem - Gridding Pre & Pos Compensation !   ....    \n');



fprintf('\n\n 5 - Sucessfully finished \n\n')


%% part 5.5 - Save data information

if saveData == 'True'
    if PC == 1
        mkdir([file_folder_results '\Data'])
    else
        addpath(genpath('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/3_Reconstruction_Code'))
        mkdir([file_folder_results '/Data'])
    end
    % --- 5.5.1 - save kdata -
    cd(folder)
    kdata_save = dat2mat(fn);                 % [Nadc, nCoils, nsl&nbv]
    nimgs      = floor(size(kdata_save,3)/Nshots);  % Number of images acquired, counting the calibration.
    kdata_save = double(reshape(kdata_save(:,:,1:nsl*nimgs),[size(kdata_save,1) size(kdata_save,2) Nshots nimgs])); % Reshape info data
    dat        = permute(kdata_save(:,:,:,1),[1 3 2 4]);
    
    
    aux_dat = dat;
    dat = aux_dat(:,:,3)
    
    
    %size(kdata_save)
    nameSave_kdata = ['rawdata_spiral.mat'];
    delete(nameSave_kdata);
    cd([file_folder_results '/Data']);
    save(nameSave_kdata,'dat');
    
    % --- 5.5.2 - save k-trajectory -
    %mat2np(ntraj,'ktraj.npy','float64');
    nameSave_ntraj = ['ktraj.mat'];
    delete(nameSave_ntraj);
    cd([file_folder_results '/Data']);
    save(nameSave_ntraj,'ntraj');
    
    % --- 5.5.3 - save Density correction factor - 
    nameSave_dcf = ['dcf.mat']; 
    dcf_out = acf;
    save(nameSave_dcf,'dcf_out');

    % --- 5.5.4 - save field map -
    getFieldMap    
    
    cd(folder)
end


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

% Plot all coils for 1 rep
maxValue = max(max(max(im_sos)));
figure()
for ii=1:size(im,3)
    subplot(4,size(im,3)/4,ii)
    imshow(im(:,:,ii,:,1),[]);
    hold on
    title(['Coil ',num2str(ii)])
    caxis([0 maxValue])    
end
fprintf('\n\n 6 - Sucessfully finished \n\n')


%% part 7 - save Results

if saveResults == 'True'
    if testKtraj == 'B0EC'
        if Diffusion == 'True'
            imTest = 'B0EC_Diffus';
        else
            if H0_corr == 'True'
                imTest = 'BOEC_H0_H1';                
            else
                imTest = 'BOEC_H1';
            end
        end
    else
        imTest = 'Theor';
    end
    
    cd(file_folder_results)
    nameSave = ['test_',num2str(tmp.test),'_imRecon_',imTest,'.mat'];
    delete(nameSave)
    save(nameSave,'im_sos')
    fprintf('\n\n 7- Save Sucessfully finished \n\n')
    
end