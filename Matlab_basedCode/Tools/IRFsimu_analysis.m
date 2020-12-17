%% Simulate the impact of IRF in reconstructing an image using JEMRIS
% Code by: TTFernandes - Nov, 2020

%% --- 1 - Addpaths
clear all 
clc
myFolder = '/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF';
addpath(genpath(myFolder))
addpath(genpath('/home/tfernandes/Documents/MATLAB/Toolboxes/griddingCodes-master'));
addpath(genpath('/home/tfernandes/Documents/MATLAB/Toolboxes/nifti_matlab'));
addpath(genpath('/home/tfernandes/Documents/MATLAB/Toolboxes/pulseq-master'));
addpath(genpath('/home/tfernandes/Documents/MATLAB/Toolboxes/NUFFT'));
addpath(genpath('/home/tfernandes/Documents/MATLAB/Toolboxes/Reconstruction_Gridding'));
addpath(genpath('/home/tfernandes/Documents/MATLAB/Tests/18_10_Tests/Tiago_code/Codigo_principal'))
addpath(genpath('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF'))
addpath('/usr/local/share/jemris/matlab/')

fprintf('\n\n 1 - Sucessfully finished \n\n')
cd(myFolder)
%% --- 2 - Parameters
% 2.1 - Type of tests
test      = 26;
testSPI   = 24;

tissue    = '2Sp';          % Type of tissue for sample: 'LeS' - only Lesion & 'All' - Full Brain & '2Sp' - Shpere
noiseTest = 'Fals';         % 
commonH   = 'Fals';         % Create a common H for xx and yy - 'True' OR 'Fals' - different H for xx and yy

InputTest = 'SPI';          % Type of input: 'SPI' - Spiral traj OR 'TRI' Triangular traj OR 'EPI' - EPI traj
testConv  = 'AboS';         % 'True' - test with convolution; 'AboS' - With implementation of AboSeada; 'NoTs' - No test with IRF

testHScan = 58;             % test to get IRF
Htype     = 'HScan';        % IRF from Scan - 'HScan' or IRF from Gaussian function - 'HGaus'

Diffusion = 'Fals';         % 'True' Test Diffusion Gradient OR 'Fals' do not test

plotTest3  = 'True';         % Plot Test - 'True' or 'Fals'
plotTest4  = 'True';         % Plot Test - 'True' or 'Fals'
plotTest5  = 'True';         % Plot Test - 'True' or 'Fals'
plotTest6  = 'True';
plotTest11 = 'True';
plotTest12 = 'True';

% 2.2 - Test Parameters
meanGAUS  = 1;              % Mean of Gaussian for H_t
FOV = 200     % mm
Kx  = 100;   
res = FOV/Kx; % mm

fprintf('\n\n 2 - Sucessfully finished \n\n')

%% --- 3 - Build Input
% ======================================================================= %
if InputTest == 'SPI' 
    % load data
    cd([myFolder '/Input_sequences'])
% %     tmpMATL = load('sequence_SPIRAL.mat'); % test107
    load(['sequence_SPIRAL_t',num2str(testSPI),'.mat']); % test207
    tmp = load('sequence_SPIRALinfo.mat'); % test107
    cd(myFolder)
    % get I_t - parametros
    Nshots       = double(tmpMATL.Nshots);                          % number of interleaves - number of trajectories (SPIRAL)
    alpha        = double(tmpMATL.alpha);                           % Oversampling - greater than 1 greater the density in center than outside the spiral. <1 the other way around (SPIRAL)
    N            = double(tmpMATL.Nx);                              % matrix size (SPIRAL)
    nsl          = tmpMATL.n_slices;                                % number of slices
    Nadc         = tmpMATL.adc_samples;                             % number of ADCs
    dt           = tmpMATL.dt*1e3;                                  % dt (s)
    gamma        = 42.576e6;                                        % Gyro-magnetic ratio (Hz/T)
    tConvert     = [0:dt:998*dt-dt];                                % time range (ms) - according to H from scan.
    tt           = [0:dt:size(tmpMATL.gradients,1)*dt-dt];          % time range (ms) - according to H from scan.
    kl           = tmpMATL.kl;
    lenT         = 12;                                              % Length of T
    T            = linspace(0.05,0.16,lenT);                        % time range - (ms)
    
    if Diffusion == 'True'
        % 1 - input parameters
        b     = 400;       % b-value (s*mm^-2)
        smMax = 200;       % unidades (T/m/s)
        gmMax = 0.030;     % unidades (T/m)
        
        % 2 - Get TE
        convertB = b*1e6; % s/m^2       
        TE_ideal = ((12 * convertB) / ((gamma*2*pi)^2 * gmMax^2))^(1/3);       % Echo Time (s) for trapezoidical gradients
        TEr = TE_ideal*1e3;   % in ms
        TEh = TEr/2;          % Half TE in ms
        
        % 3 - Time of slope
        DT_slope      = gmMax / smMax * 1e3;         % in ms
        nPoints_slope = DT_slope/dt;   
        
        % 4 - Get Gdiff Pulse
        vectorTimeDiff = [0:1:round(TEh/dt)]*dt;     % in ms
        for iDiff=1:size(vectorTimeDiff,2)
            if iDiff == 1
                gDiff(iDiff) = 0;
            elseif iDiff <= nPoints_slope
                gDiff(iDiff) = 200*iDiff*dt*1e-3;
            elseif  nPoints_slope < iDiff && iDiff <= size(vectorTimeDiff,2)-nPoints_slope
                gDiff(iDiff) = gmMax;
            elseif iDiff > size(vectorTimeDiff,2)-nPoints_slope
                gDiff(iDiff) = gmMax - 200*(iDiff-size(vectorTimeDiff,2)+nPoints_slope)*dt*1e-3;
            end
        end
        
        gDiff = gDiff * gamma * 1e-6;   % units mT/mm * gamma
        
        % 5 - Plots    
        if plotTest3 == 'True'            
            figure()
            plot(gDiff*1e6/gamma)
            hold on
            grid
        end
    end
    
    % obtain Gx & Gy
    Ix_t = zeros(Nadc*Nshots,1);
    Iy_t = zeros(Nadc*Nshots,1);
    for iNs=1:Nshots
        start    = (iNs-1)*Nadc + 1;
        Ix_t(start:start + Nadc - 1)       = real(tmpMATL.gradients(1:Nadc,iNs)); % [Hz/m]
        Iy_t(start:start + Nadc - 1)       = imag(tmpMATL.gradients(1:Nadc,iNs)); % [Hz/m]
        origSize = size(Ix_t,1);
    end          

    if plotTest3 == 'True'        
        figure()
        subplot(211)
        plot(tt,Ix_t*1e6/gamma)
        hold on
        ylabel(['Gradiente Amplitude (T/m)'])
        xlabel(['time (ms)'])
        title('Ix')
        subplot(212)
        plot(tt,Iy_t*1e6/gamma,'r')
        hold on
        ylabel(['Gradiente Amplitude (T/m)'])
        xlabel(['time (ms)'])
        title('Iy')
    end

    
    % obtain theoretic k-space trajectory from the gradients ---
    int_mat = triu(ones(Nadc));
    kx = zeros(Nadc*Nshots,1);
    ky = zeros(Nadc*Nshots,1);
    for iNs=1:Nshots     % Have to split the integral because each of the gradients is independent. Otherwise we would have a continuous integral
        start    = (iNs-1)*Nadc + 1;
        kx(start:start + Nadc - 1) = real(tmpMATL.gradients(1:Nadc,iNs))' * int_mat * dt; % vem em unidades gamma
        ky(start:start + Nadc - 1) = imag(tmpMATL.gradients(1:Nadc,iNs))' * int_mat * dt;
    end
    kx = kx(1:end-2,1);
    ky = ky(1:end-2,1);
    traj = kx + 1i*ky;
    if plotTest3 == 'True'
        figure()
        subplot(121)
        plot(traj*2*pi)
        subplot(122)
        plot(kl/2/pi)
    end    
    
   
% ======================================================================= %
elseif InputTest == 'TRI'
    
    test  = 'SCAN';
    saveT = 'Fals';
    
    lenT   = 12;                        % Length of T
    T      = linspace(0.05,0.16,lenT);  % time range - (ms)
    p      = 115;                       % slope - mT/(m ms)
    DT     = 0.01;                      % raster time (ms)
    tTotal = 10;                        % Total time (ms)
    tpre   = 5;                         % time Pre Pulse (ms)
    
    if test == 'PULS'
        t    = [-0.2:DT:0.2]; % time range - ms
        tTotal = 2*0.2;
    elseif test == 'SCAN'
        t    = [-tpre:DT:tTotal-tpre]; % time range - ms
    end
    
    
    % convertion
    tConvert = t *1e-3;     % in s
    Tconvert = T *1e-3;     % in s
    pConvert = p;           % in T/(m s)
    gammabar = 42.576e6;    % Gyro-magnetic ratio Hz/T
    
    % Build vector & Plot
    i_t = zeros(length(tConvert),1);
    maxAmp    = zeros(1,lenT);
    maxMoment = zeros(1,lenT);
    
    if plotTest3 == 'True'
        figure()
    end
    for ii=1:lenT
        for i=1:length(tConvert)
            if abs(tConvert(i)) <= Tconvert(ii)
                i_t(i,ii) = (-1) * pConvert*(Tconvert(ii)-abs(tConvert(i))) * gammabar ;
            elseif abs(tConvert(i)) > Tconvert(ii)
                i_t(i,ii) = (-1) * 0 * gammabar;
            end
        end
        maxAmp(ii)    = (-1) * pConvert * Tconvert(ii);  % Maximum Amplitude of the Pulse
        maxMoment(ii) = pConvert * Tconvert(ii)^2;       % Moment of the Pulse
        if plotTest3 == 'True'            
            plot(tConvert,i_t(:,ii))
            hold on
        end
    end
    if plotTest3 == 'True'        
        ylabel(['Gradient (T/m) * Gamma'],'fontsize',15)
        xlabel(['time (s)'],'fontsize',15)
        title('GIRF pulses')
        xlim([-0.2 0.2]*1e-3)
    end
    Ix_t = i_t;
    Iy_t = i_t;
% ======================================================================= %    
elseif InputTest == 'EPI'
    
    
end

if Diffusion == 'True'
    % add Gdiff with K_trajc
    Ix_t = [gDiff Ix_t']';
    Iy_t = [gDiff Iy_t']';
    
    tConvert_All = [0:dt:size(Ix_t,1)*dt-dt];
    
    if plotTest3 == 'True'
        figure()
        subplot(211)
        plot(tConvert_All,Ix_t*1e6/gamma)
        hold on
        ylabel(['Gradient (T/m)'],'fontsize',15)
        title('Gradient Diff + Traject - Ix')
        subplot(212)
        plot(tConvert_All,Iy_t*1e6/gamma,'r')
        hold on
        ylabel(['Gradient (T/m)'],'fontsize',15)
        xlabel(['time (ms)'],'fontsize',15)
        title('Gradient Diff + Traject - Iy')
    end
    
end

fprintf('\n\n 3 - Sucessfully finished \n\n')

%% --- 4 - Build Template H
% 4.1 - H(t)
if Htype == 'HGaus'
    if commonH == 'True'
        Hx_t = normpdf(tConvert,0,0.01);
        Hy_t = Hx_t;
    elseif commonH == 'Fals'
        Hx_t = normpdf(tConvert,0,0.001);
        Hy_t = normpdf(tConvert,0,0.01);
    end
    
    if plotTest4 == 'True'
        figure()
        plot(tConvert,Hx_t/max(abs(Hx_t)))
        hold on
        title(['IRF from Gaussian Function in time w/ mean of Gauss = ',num2str(meanGAUS)])
    end
    
    % 4.2 - H(w)
    Hx_w = fftshift(abs(fft(Hx_t)/numel(tConvert)));
    Hy_w = fftshift(abs(fft(Hy_t)/numel(tConvert)));
    
    if plotTest4 == 'True'
        figure()
        plot(Hx_w)
        hold on
        title('IRF from Gaussian Function in frequence')
    end
    
    % 4.3 - Add noise to H(t) & H(w)
    if noiseTest == 'True'
        nPoints   = size(tConvert,2);
        if commonH == 'True'
            noiseLev   = 10;
            noise      = noiseLev * randn(nPoints,1);
            Hx_t_noisy = Hx_t .* noise';
            Hy_t_noisy = Hy_t .* noise';
        elseif commonH == 'Fals'
            noiseLev_x = 5;
            noise_x    = noiseLev * randn(nPoints,1);
            Hx_t_noisy = Hx_t .* noise_x';
            
            noiseLev_y = 10;
            noise_y    = noiseLev * randn(nPoints,1);
            Hy_t_noisy = Hy_t .* noise_y';
        end
        
        % ==== plot figures ====
        if plotTest4 == 'True'
            figure()
            %     subplot(121)
            plot(tConvert,H_t_noisy')
            hold on
            ylabel(['Gradient (T/m) * Gamma'],'fontsize',15)
            xlabel(['time (s)'],'fontsize',15)
            title('GIRF Distorced')
            %     subplot(122)
            plot(tConvert,H_t')
            
            %     ylabel(['Gradient (T/m) * Gamma'],'fontsize',15)
            %     xlabel(['time (s)'],'fontsize',15)
            %     title('GIRF Original')
        end
        clear noise
        
        % %     input            = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_pulse_', num2str(Npulses),'_DT_1_1e-5_p', num2str(p),'.mat']);
        % %     input_freq       = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_',type,'_Input_freq_', num2str(Npulses),'.mat']);
        % %     input_sincF_freq = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_',type,'_Input_freq_sincFunct_', num2str(Npulses),'.mat']);
    end
elseif Htype == 'HScan'
    cdResults = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/2_GIRF_test/Test/Results_IRF'];
    cd(cdResults)
    load(['test',num2str(testHScan),'_IRF.mat'])   
    Hx_w    = H.Hx_w(:,1)';
    Hy_w    = H.Hy_w(:,1)';
    Npulses = H.pulses;

    maxFreq   = 1./T(1)*1e3;
    ff        = linspace(-maxFreq,maxFreq,size(Hx_w,2))*1e-3; % frequencies vector (Hz)
    
    clear H
    cd(myFolder)
    
    if plotTest4 == 'True'
        figure()
        plot(ff,Hx_w,'b')
        hold on
        plot(ff,Hy_w,'r')        
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        title('GIRF x-blue & y-red')        
    end
end
fprintf('\n\n 4 - Sucessfully finished \n\n')


%% --- 5 - Obtain O according to Vannesjmo, 2013

if testConv == 'True'
    for ii=1:size(Ix_t,2)
        Ox(:,ii) = conv(Ix_t(:,ii),Hx_t)/gamma*1e9;
        Oy(:,ii) = conv(Iy_t(:,ii),Hy_t)/gamma*1e9;
    end
    if plotTest5 == 'True'  
        % Input - I_t
        figure()
        subplot(121)
        plot(Ix_t/gamma)
        hold on
        title('Ix')        
        subplot(122)
        plot(Iy_t/gamma)
        hold on
        title('Iy')            
                                
        % Input - O_t        
        figure()
        subplot(121)
        plot(Ox(1:Nadc,:))
        hold on
        title(['Ox w/ mean of Gauss = ',num2str(meanGAUS)])        
        subplot(122)
        plot(Oy(1:Nadc,:))
        hold on
        title(['Oy w/ mean of Gauss = ',num2str(meanGAUS)])        
    end
        
    % generate '.h5' files
    gan        = Ox(1:Nadc,:) + 1*i * Oy(1:Nadc,:);       
    NPTOTAL    = size(tConvert,2);
    interleave = 1;
    type       = 'SIMU';    
    
    cd_dir    = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/',InputTest,'_test_',num2str(test)];
    dir_inter = [cd_dir '/interleaves'];
    mkdir(cd_dir)
    mkdir(dir_inter)
    
    generate_h5_grad( gan', tt, NPTOTAL, interleave, dir_inter, InputTest , type)
    
    clc
    clear gan NPTOTAL interleave type
    
    
elseif testConv == 'AboS'
    if Diffusion == 'Fals'
        maxT    = max(abs(Hx_w));
        Htest   = Hx_w' / maxT; % normalized freq domain
        % % Htest   = H_w;
        DT      = dt*1e-3;     % in seconds
        Npad    = 0;           % testar
        maxFreq = 1./T(1)*1e3;
        ff      = linspace(-maxFreq,maxFreq,size(Htest,1));    % Frequency in Hz
        
        if commonH == 'True'
            for ii=1:size(Ix_t,2)
                Gx_input = Ix_t(:,ii);       % input in time domain
                Gy_input = Iy_t(:,ii);       % input in time domain
                
                [Ox(:,ii),~,~] = gradient_distort_GIRF(Gx_input,ff,Htest,DT,Npad);
                [Oy(:,ii),~,~] = gradient_distort_GIRF(Gy_input,ff,Htest,DT,Npad);
                
                clear Gx_input Gy_input
            end
        elseif commonH == 'Fals'
            maxT_x  = max(abs(Hx_w));
            Htest_x = Hx_w' / maxT_x;  % normalized freq domain
            maxT_y  = max(abs(Hy_w));
            Htest_y = Hy_w' / maxT_y;  % normalized freq domain
            
            for ii=1:size(Ix_t,2)
                Gx_input = Ix_t(:,ii);       % input in time domain
                Gy_input = Iy_t(:,ii);       % input in time domain
                
                [Ox(:,ii),~,~] = gradient_distort_GIRF(Gx_input,ff,Htest_x,DT,Npad);
                [Oy(:,ii),~,~] = gradient_distort_GIRF(Gy_input,ff,Htest_y,DT,Npad);
                
                clear Gx_input Gy_input
            end
            clear maxT_x Htest_x maxT_y Htest_y
        end
        
        
    % ===== Test with diffusion ======
    elseif Diffusion == 'True'
        maxT    = max(abs(Hx_w));
        Htest   = Hx_w' / maxT; % normalized freq domain
        DT      = dt*1e-3;     % in seconds
        Npad    = 0;           % testar
        maxFreq = 1./T(1)*1e3;
        ff      = linspace(-maxFreq,maxFreq,size(Htest,1));    % Frequency in Hz
        
        if commonH == 'True'
            for ii=1:size(Ix_t,2)
                Gx_input = Ix_t(:,ii);       % input in time domain
                Gy_input = Iy_t(:,ii);       % input in time domain
                
                [Ox(:,ii),~,~] = gradient_distort_GIRF(Gx_input,ff,Htest,DT,Npad);
                [Oy(:,ii),~,~] = gradient_distort_GIRF(Gy_input,ff,Htest,DT,Npad);
                
                clear Gx_input Gy_input
            end
        elseif commonH == 'Fals'
            maxT_x  = max(abs(Hx_w));
            Htest_x = Hx_w' / maxT_x;  % normalized freq domain
            maxT_y  = max(abs(Hy_w));
            Htest_y = Hy_w' / maxT_y;  % normalized freq domain
            
            for ii=1:size(Ix_t,2)
                Gx_input = Ix_t(:,ii);       % input in time domain
                Gy_input = Iy_t(:,ii);       % input in time domain
                
                [Ox(:,ii),~,~] = gradient_distort_GIRF(Gx_input,ff,Htest_x,DT,Npad);
                [Oy(:,ii),~,~] = gradient_distort_GIRF(Gy_input,ff,Htest_y,DT,Npad);
                
                clear Gx_input Gy_input
            end
            clear maxT_x Htest_x maxT_y Htest_y
        end
    end
    
    
    if plotTest5 == 'True'
        % Input - I_t & O_t
        tt = [0:dt:size(Ix_t)*dt-dt];                                % time range (ms) - according to H from scan.        
        figure()
        subplot(211)
        plot(tt,Ix_t*1e6/gamma,'b')
        hold on
        plot(tt,Ox*1e6/gamma,'g')
        hold on
        title(['Ix - blue & Ox - green'])        
        ylabel(['gm (T/m)'])
        xlabel(['time (ms)'])
        subplot(212)
        plot(tt,Iy_t*1e6/gamma,'b')
        hold on
        plot(tt,Oy*1e6/gamma,'g')
        hold on
        title(['Iy - blue & Oy - green'])        
        ylabel(['gm (T/m)'])
        xlabel(['time (ms)'])
    end
    
    % generate '.h5' files
    if Diffusion == 'True'
        Ox_test = Ox;
        Oy_test = Oy;               
        Ox      = Ox_test(size(Ox,1)-origSize+1:end);
        Oy      = Oy_test(size(Oy,1)-origSize+1:end);
        tt      = [0:dt:size(tmpMATL.gradients,1)*dt-dt];          % time range (ms) - according to H from scan.
        
        if plotTest5 == 'True'            
            figure()
            subplot(211)
            plot(tt,Ox)
            subplot(212)
            plot(tt,Oy,'r')
        end
    end
    
    gan        = Ox(1:Nadc,:) + 1*i * Oy(1:Nadc,:);
    NPTOTAL    = size(tt,2);
    interleave = 1;
    type       = 'SIMU';
    
    cd_dir    = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/',InputTest,'_test_',num2str(test)];
    dir_inter = [cd_dir '/interleaves'];
    mkdir(cd_dir)
    mkdir(dir_inter)
    
    generate_h5_grad( gan', tt, NPTOTAL, interleave, dir_inter, InputTest , type)
    
    clc
    clear gan NPTOTAL interleave type
    
    
    
elseif testConv == 'NoTs'
    for ii=1:size(Ix_t,2)
        Ox(:,ii) = Ix_t(:,ii);
        Oy(:,ii) = Iy_t(:,ii);
    end
    if plotTest5 == 'True'
        % Input - I_t
        figure()
        subplot(121)
        plot(Ix_t/gamma*1e6)
        hold on
        grid on
        ylabel(['gm (mT/m)'])
        title('Ix')        
        subplot(122)
        plot(Iy_t/gamma*1e6,'r')
        hold on
        grid on        
        title('Iy')            
                                
        % Input - O_t        
        figure()
        subplot(121)
        plot(Ox/gamma*1e6)
        hold on
        title(['tissueOx w/ mean of Gauss = ',num2str(meanGAUS)])        
        subplot(122)
        plot(Oy/gamma*1e6)
        hold on
        title(['Oy w/ mean of Gauss = ',num2str(meanGAUS)])        
    end
    
    % === generate '.h5' files ===
    gan        = (Ox + 1*i * Oy)';
    NPTOTAL    = size(tt,2);
    interleave = 1;
    type       = 'SIMU';    
    
    cd_dir    = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/',InputTest,'_test_',num2str(test)];
    dir_inter = [cd_dir '/interleaves'];
    mkdir(cd_dir)
    mkdir(dir_inter)
    
    generate_h5_grad( gan, tt, NPTOTAL, interleave, dir_inter, InputTest , type)
    
    clc
    clear NPTOTAL type    
    
    
end

% % % % % --- read sample ---
if plotTest5 == 'True'    
    cd(dir_inter)
    sample_dataY = h5read('mGY1.h5','/extpulse');
    sample_dataX = h5read('mGX1.h5','/extpulse');
    
    figure();
    subplot(211)
    plot(sample_dataX(:,2)/2/pi/gamma*1e6,'g')
    hold on
    ylabel(['Gradiente Amplitude (T/m)'])
    title('Gradient X')
    subplot(212)
    plot(sample_dataY(:,2)/2/pi/gamma*1e6,'r')
    title('Gradient Y')    
    ylabel(['Gradiente Amplitude (T/m)'])
    cd(cd_dir)
end

cd(cd_dir)
clc
fprintf('\n\n 5 - Sucessfully finished \n\n')


%% --- 6 - Select Sample
dir_Orig          = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences'];

if tissue == 'LeS'
    cd(dir_Orig)
    name_seq_template = 'sample.h5';
    copyfile(name_seq_template,cd_dir)
    cd(cd_dir)
    clear name_seq_template    
elseif tissue == 'All'
    % ======= Create new sample ==============
    dir_seq   = cd_dir;
    widenSz   = 10;
    origFOV   = 180;
    low_field = 0.5;     % magnitude of scanner field (T)
    
    Sample_Phantom_generator_BRAIN_all_4regions
    
elseif tissue == '2Sp'
    cd(dir_Orig)
    widenSz = 10;
    Sample_Phantom_generator_2SPHERES
    cd(cd_dir)
    clear name_seq_template   
end

% % % % --- read sample ---
sample_data = h5read('sample.h5','/sample/data');
resolutioN = h5read('sample.h5','/sample/resolution');
offsetS = h5read('sample.h5','/sample/offset');

if plotTest6 == 'True'    
    figure();
    imagesc(reshape(sample_data(2,:,:),size(sample_data,2),size(sample_data,3)))
    colormap gray
end

fprintf('\n\n 6 - Sucessfully finished \n\n')

%% --- 7 - Create 'seq.xml'
% 7.1 - create file folder and '.xml'
% JEMRIS_seq
name    = ['seq_',InputTest,'_test_',num2str(test),'.xml'];
dir_aux = struct2cell(dir);
any(ismember(dir_aux(1,:),name))
if ans == 0
    cd(dir_Orig)
    name_seq_template = 'seq.xml';
    copyfile(name_seq_template,cd_dir)
    cd(cd_dir)
    movefile(name_seq_template,name);        % Rename the file
end
clear name_seq_template ans


fprintf('\n\n 7 - Sucessfully finished \n\n')

%% --- 8 - Run on Jemris
% 8.1 - Run on Jemris
! jemris seq_SPI_test_26.xml
% JEMRIS_seq
% clc

name_seq     = 'seq.h5'; % Nome do ficheiro seqclc

ninterleave  = 1;
Treal        = h5read(name_seq,'/seqdiag/T');
RXP          = h5read(name_seq,'/seqdiag/RXP');       % RF Receiver phase; unit: radiants; if negative, the TPOI was not an ADC
RXP          = logical(RXP);RXP = imcomplement(RXP);  % separa os ADC (=1) de todos os TPOIs (=0)
KX           = h5read(name_seq,'/seqdiag/KX');        % K_trajectory for X-Gradient
KY           = h5read(name_seq,'/seqdiag/KY');        % K_trajectory for Y-Gradient
KX2          = KX(RXP);
KY2          = KY(RXP);

kreal_aux = KX2 + 1i*KY2;                   % Gradient on ADCs
nT        = (length(Treal)-1)/ninterleave; % Total number of points per interleave
kreal     = kreal_aux';

clear kreal_aux KX KY KX2 KY2 RXP Treal ninterleaves

if plotTest6 == 'True'
    figure();
    %     plot(kreal/2/pi*1e3,'r')
    plot(kreal,'r')
    hold on
    ylabel(['rad/mm'])
    xlabel(['rad/mm'])
    grid()
    %     plot(traj)
end

fprintf('\n\n 8 - Sucessfully finished \n\n')

%% --- 9 - Create 'simu.xml'
% 9.1 - copy simu.xml file to test folder...
name_simu_template = 'simu.xml';
cd(cd_dir)
any(ismember(dir_aux(1,:),name_simu_template));
if ans == 0
    cd(dir_Orig)
    copyfile(name_simu_template,cd_dir)
end
clear name_simu_template ans name dir_aux ans
cd(cd_dir)

fprintf('\n\n 9 - Sucessfully finished \n\n')

%% --- 10 - Simulate in JEMRIS
% 1 - Command line change directory to match cd_dir (in ssh -x asia.local)
% 2 - Run the follwing code:
% cd /home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/SPI_test_8
% mpirun -np 15 pjemris simu.xml
close all

fprintf('\n\n 10 - Sucessfully finished \n\n')

%% --- 11 - Read Data from simulation in JEMRIS
% 11.1 - load signals ...
cd(cd_dir)
name_simu = ['Data_',InputTest,'_test',num2str(test),'.mat'];
read_h5signals(name_simu,'signals.h5')
load(name_simu)

% 11.2 - Organize data data obtained with simulation
clear kdata
kdata    = complex(M(:,1),M(:,2));   

fprintf('\n\n 11 - Sucessfully finished \n\n')

%% --- 12 - Reconsctruct Data and Obtaining Results
clear m_postcomp

% ... preparing gridding ...
gridsize = 270;
% gridsize = round((2*kmax)/(1/FOVX));
p2r      = real(kreal); % parte da imagem real
p2i      = imag(kreal); % parte da imagem imaginária
kloc     = p2r+p2i*1i;  % k_trajectory created

% controlar para que o meu valor do espaço-k não seja maior do que 0.5
if max(p2r)>0.5 | max(p2i)>0.5
    p2r=(p2r*0.5)/max(p2r);
    p2i=(p2i*0.5)/max(p2i);
end

kTraj = [imag(kreal)' real(kreal)'];
kloc  = (p2r+p2i*1i);   % k_trajectory created

% [dcf,k] = vecdcf(G2,kloc,1,res,FOVX);
% % figure,
% % plot(kloc, 'b')
% % hold on
% % plot(kdata, 'r')
% % hold on
% % grid()

% ... pre-compensation weights (independent of simulated data) ...
[acf,flag,res] = makePrecompWeights_2D_VORONOI( kTraj, [gridsize/2 gridsize/2] );

% ... selecao de parametros ...
overgridfactor = 2;
kwidth         = 3; % kwidth = 50;  % kernel width
kbeta          = pi*sqrt(kwidth^2/overgridfactor^2*(overgridfactor-0.5)^2-0.8 );  % kernel beta
w              = acf; % w = weigth(kloc);

% ... keiser-bessel gridding + fft ...
%     trimming       -- 'y' or 'n'
%     apodize        -- 'y' or 'n'
%     postcompensate -- 'y' or 'n'
[m_postcomp,kwx,kwy] = grid_kb_teste141118(kdata,kloc,w,gridsize,overgridfactor,kwidth,kbeta,'n','n','y');
im                   = fftshift(fft2(ifftshift(m_postcomp)));

% ... crop figures ...
m_postcomp_crop = m_postcomp(gridsize-round(N/2)+1:gridsize+round(N/2),gridsize-round(N/2)+1:gridsize+round(N/2));
im_crop = im(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2),round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2));

% ... plot das figuras ...
if plotTest12 == 'True'
    figure;
    % subplot(121)
    imshow(abs(m_postcomp_crop),[]);    
    figure;
    % subplot(121)
    imshow(abs(im_crop),[]);
end


fprintf('\n\n 12 - Sucessfully finished \n\n')

%% 13 - RUN without H

cd(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/SPI_test_',num2str(testSPI)])
% 13.1 - Obtain Ktrajct
! jemris seq_SPI_test_24.xml
% JEMRIS_seq
% clc

name_seq     = 'seq.h5'; % Nome do ficheiro seqclc

RXP          = h5read(name_seq,'/seqdiag/RXP');       % RF Receiver phase; unit: radiants; if negative, the TPOI was not an ADC
RXP          = logical(RXP);RXP = imcomplement(RXP);  % separa os ADC (=1) de todos os TPOIs (=0)
KX           = h5read(name_seq,'/seqdiag/KX');        % K_trajectory for X-Gradient
KY           = h5read(name_seq,'/seqdiag/KY');        % K_trajectory for Y-Gradient
KX2          = KX(RXP);
KY2          = KY(RXP);

kreal_aux     = KX2 + 1i*KY2;                   % Gradient on ADCs
kreal_noH     = kreal_aux';

clear kreal_aux KY2 KX2 KY KX RXP Treal
cd(cd_dir)

% 13.2 - Recon image
clear m_postcomp

% ... preparing gridding ...
p2r_noH      = real(kreal_noH); % parte da imagem real
p2i_noH      = imag(kreal_noH); % parte da imagem imaginária
kloc_noH     = p2r_noH+p2i_noH*1i;  % k_trajectory created

% controlar para que o meu valor do espaço-k não seja maior do que 0.5
if max(p2r_noH)>0.5 | max(p2i_noH)>0.5
    p2r_noH=(p2r_noH*0.5)/max(p2r_noH);
    p2i_noH=(p2i_noH*0.5)/max(p2i_noH);
end

kTraj_noH = [imag(kreal_noH)' real(kreal_noH)'];
kloc_noH  = (p2r_noH+p2i_noH*1i);   % k_trajectory created

% [dcf,k] = vecdcf(G2,kloc,1,res,FOVX);
% % figure,
% % plot(kloc, 'b')
% % hold on
% % plot(kdata, 'r')
% % hold on
% % grid()

% ... pre-compensation weights (independent of simulated data) ...
[acf_noH,flag_noH,res_noH] = makePrecompWeights_2D_VORONOI( kTraj_noH, [gridsize/2 gridsize/2] );

% ... selecao de parametros ...
w_noH              = acf_noH; % w = weigth(kloc);

% ... keiser-bessel gridding + fft ...
%     trimming       -- 'y' or 'n'
%     apodize        -- 'y' or 'n'
%     postcompensate -- 'y' or 'n'
[m_postcomp_noH,kwx,kwy] = grid_kb_teste141118(kdata,kloc_noH,w_noH,gridsize,overgridfactor,kwidth,kbeta,'n','n','y');
im_noH                   = fftshift(fft2(ifftshift(m_postcomp_noH)));

% ... crop figures ...
m_postcomp_noH_crop = m_postcomp_noH(gridsize-round(N/2)+1:gridsize+round(N/2),gridsize-round(N/2)+1:gridsize+round(N/2));
im_noH_crop = im_noH(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2),round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2));

% ... plot das figuras comparaçao entre c/ e s/ H ...
if plotTest12 == 'True'
    figure;
    subplot(131)    
    imshow(abs(im_crop),[]);    
    title('Recon with Optimal Traj - H')
    subplot(132)
    imshow(abs(im_noH_crop),[]);
    title('Recon with Theor Traj') 
    subplot(133)
    Diff_im_crop = abs(im_crop) - abs(im_noH_crop);
    imshow(Diff_im_crop,[]);
    colormap jet
    title('Diff between images')     
end

% ... plot das figuras comparaçao entre c/ e s/ H com sample ...
NormSample       = sample_data(2,:,:)./max(max(sample_data(2,:,:)));
NormSample       = reshape(NormSample,size(NormSample,2),size(NormSample,2));
Norm_im_crop     = abs(im_crop)/max(max(abs(im_crop)));
Norm_im_noH_crop = abs(im_noH_crop)/max(max(abs(im_noH_crop)));

reSiz            = size(Norm_im_crop,1)/size(NormSample,1);
NormSample_reS   = imresize(NormSample,reSiz);

CompSample_H     = NormSample_reS-Norm_im_crop;
CompSample_noH   = NormSample_reS-Norm_im_noH_crop;

if plotTest12 == 'True'
    figure;
    subplot(121)    
    imshow(CompSample_H,[]);  
    caxis([-1 1])
    title('Sample - H')
    subplot(122)
    imshow(CompSample_noH,[]);  
    hold on
    caxis([-1 1])
    title('Sample - noH') 
end
fprintf('\n\n 13 - Sucessfully finished \n\n')


%% 14 - Image profile intensities

% 1 - 
clear c_im c_im_noH
x_im(1,:)     = [size(Norm_im_crop,1)/10*3  size(Norm_im_crop,1)/10*3];
x_im(2,:)     = [size(Norm_im_crop,1)/10*3  size(Norm_im_crop,1)/10*7];
x_im(3,:)     = [size(Norm_im_crop,1)/10*0  size(Norm_im_crop,1)/10*10];

y_im(1,:)     = [size(Norm_im_crop,1)/10*0 size(Norm_im_crop,1)/10*10];
y_im(2,:)     = [size(Norm_im_crop,1)/10*0 size(Norm_im_crop,1)/10*10];
y_im(3,:)     = [size(Norm_im_crop,1)/10*6 size(Norm_im_crop,1)/10*6];

n = 186;  %size(Norm_im_crop,1);

for ji=1:size(y_im,1)
    c_im(:,ji)     = improfile(Norm_im_crop,x_im(ji,:),y_im(ji,:),n);
    c_im_noH(:,ji) = improfile(Norm_im_noH_crop,x_im(ji,:),y_im(ji,:),n);
    figure()
    subplot(221)
    plot(c_im(:,ji),'g')
    hold on
    title(['ImProfile_x_x w/ H'])
    subplot(222)
    plot(c_im_noH(:,ji),'r')
    hold on
    title(['ImProfile_x_x: w/out H'])
    subplot(223)
    title('Image Recon w/H')
    imshow(Norm_im_crop,[]);
    hold on
    line(x_im(ji,:),y_im(ji,:))
    subplot(224)
    title('Image Recon w/H')
    imshow(Norm_im_noH_crop,[]);
    line(x_im(ji,:),y_im(ji,:))
end

% save
clear nameSave
nameSave=['test_',num2str(test),'_imProfile']
save(nameSave,'c_im','c_im_noH','Norm_im_crop','Norm_im_noH_crop','testSPI')

% clear x_im y_im

fprintf('\n\n 14 - Sucessfully finished \n\n')

%% 15 - Numerical analysis
if tissue == '2Sp'
    whichcols = 2;                              % how many bourders do I want
    numTraces = size(x_im,1)*whichcols;         % Number of boundaries values
    c_GIRF    = zeros(round(n/2),size(whichcols,2));
    c_THEO    = zeros(round(n/2),size(whichcols,2));
    
    FOV = 200;               % mm
    Nx  = N;
    res = FOV/Nx;            % in mm
    
    for ji=1:size(y_im,1)
        bounds = sharp_imProfile(c_im(:,ji),'Fals');
        resSHARP_H(1,ji) = (bounds.firstUP   - bounds.firstLOW ) * res; % in mm
        resSHARP_H(2,ji) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds
        
        bounds = sharp_imProfile(c_im_noH(:,ji),'Fals');
        resSHARP_noH(1,ji) = (bounds.firstUP  - bounds.firstLOW ) * res; % in mm
        resSHARP_noH(2,ji) = (bounds.secondLOW - bounds.secondUP) * res; % in mm
        clear bounds
    end
    
    
    % 6.1 - plots
    
    figure()
    maxCaxis = max([max(max(resSHARP_H)) max(max(resSHARP_noH))]);
    vectY = [1 ; 2];
    vectX = [1 ; 2 ; 3];
    subplot(211)
    imshow(resSHARP_H,[])
    hold on
    title(['Res. Sharpness recon GIRF'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_H,1)*2)
    xticks(1:1:size(resSHARP_H,2)*2)
    yticklabels(vectY)
    xticklabels(vectX)
    
    subplot(212)
    imshow(resSHARP_noH,[])
    hold on
    title(['Res. Sharpness recon THEO'])
    caxis([0 maxCaxis])
    yticks(1:1:size(resSHARP_noH,1)*2)
    xticks(1:1:size(resSHARP_noH,2)*2)
    yticklabels(['1st Boarder' , '2nd Boarder'])
    xticklabels(['1st cut','2nd cut','3rd cut'])
    
    
    fprintf('\n\n 15 - Sucessfully finished \n\n')
end