%% Recon of images from scanner with pulseq.seq for BO_ECC sequence
% Robison, 2019
% TTFernandes, IST 2021

%% part 0 - Addpath & Get Folders
clear all;

yourpath = ['...']; % fill out with yourpath were you saved your GIRFOS_tool
addpath(genpath([yourpath, 'GIRFOS_tool/Matlab_basedCode/Tools/']));  % To reconstruct
addpath(genpath([yourpath, 'GIRFOS_tool/Matlab_basedCode/B0_ECC/']));  % To reconstruct

% close all
clc


%% part 1 - Get Files

str ={'select file'}
s = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str);
if s==1
    [file,folder]=uigetfile;
else
    folder=uigetdir;
end
strTest1 = '.dat';
auxA     = strfind(file,strTest1);
filename = file(1:auxA-1);

file_folder = folder;
file_folder_cfl_hdr = [folder '/cfl_hdr_files'];
file_folder_results = [folder, '/Results'];
mkdir([file_folder_cfl_hdr]);
mkdir([file_folder_results]);
cd(folder)

[files] = dir(file_folder);
fn = files(5).name;       % get '.dat' file
fn_info = files(6).name;  % get '.mat' file with sequence parameters
tmp = load(fn_info);      % get values from '.mat' file with sequence parameters

clear strTest1 filename auxA fn_info folder s computer str fn_info
clc
fprintf('\n\n 1 - Sucessfully finished - Get Files\n\n')

%% part 1.5 - Main parameters
MATLABfft   = 'Fals';                      % Save - 'True' (According to Matlab) or 'Fals' (acording to Abo Seaba)


%% part 2 - Define Parametres

% ... part 2.1 - set test ...
remChannels = 'Fals';                      % remove channels Test - 'True' or 'Fals'
testDummy   = 'True';                      % Test Dummy 1st pulse - 'True' - test42 or 'Fals' - test27
tMean       = 'True';                      % Test Mean for many reps
syncFun     = 'Fals';                      % Use SyncFunction for freq theorical Gradient
filter      = 'Fals';                      % Filter for H0

saveT       = 'True';                      % Save - 'True' or 'Fals' 
type        = 'SCAN';                      %'PULS' - test102 & test45 or 'SCAN' - test109

plotTest2 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest3 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest4 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest5 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest6 = 'True';                      % Plot Test - 'True' or 'Fals'
plotTest7 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest8 = 'True';                      % Plot Test - 'True' or 'Fals'
test6     = 'Fals';                      % Test for debbuging

% ... part 2.2 - parameters ....
test         = tmp.test;                    % test
Npulses      = double(tmp.Npulses);         % Number of triangular pulses pulses
Nreps        = double(tmp.nreps);           % get Repetitions 
nsl          = tmp.nslices;                 % number of slices
Ny           = 1;                           % for reshape of k vector
dirs         = tmp.ndirs;                   % number of directions
Nadc         = tmp.adc_samples;             % number of ADCs
dt           = tmp.DT;                      % dt (s)
Naxes        = tmp.Naxes;                   % Number of axes
Nobs         = 4;                           % Number of observation - 2 per axes
gamma        = 42.576e6;                    % Gyro-magnetic ratio (Hz/T)
sliceT       = tmp.st;                      % Slice thickness (m)
D            = tmp.x0;                      % Distance of the slice to gradient isocenter (m)
p            = 115;%tmp.max_slew;           % Max_slew rate (mT/s)
T            = linspace(0.05,0.16,Npulses); % time range - (ms)
gammabar     = gamma;                       % Gyro-magnetic ratio Hz/T
nDummies     = 1;%tmp.Ndummies;                % Number of dummies to achieve equilibrium
% nDummies     = tmp.Ndummies;                % Number of dummies to achieve equilibrium
x0           = tmp.x0;                      % RF off-set in m

extra     = 2;  % retirar os dois pontos que o .seq tem a menos em ADC


if testDummy == 'True' % if test dummy is present
    tDum = 1;
else
    tDum = 0;
end

% ... part 2.3 - get theoretical Gx and Gy ...
cd([yourpath, 'GIRFOS_tool/Matlab_basedCode/Example/'])

if type == 'PULS'
    tmp.gradients = load(['girf_pulse_' num2str(Npulses),'_DT_1_1e-5_p', num2str(p),'.mat']);
elseif type == 'SCAN'
    tTotal        = tmp.tPulTotal;               % total time of pulse    
    tmp.gradients = load(['girf_',type,'_pulse_' num2str(Npulses),'_tTotal_',num2str(tTotal),'ms_DT_1_1e-5_p', num2str(p)]);
end

tGx           = zeros((Nadc-extra),Npulses);
tGy           = zeros((Nadc-extra),Npulses);

for iNs=1:Npulses
    tGx(:,iNs) = tmp.gradients.i_t(1:(Nadc-extra),iNs); % [Hz/m]
    tGy(:,iNs) = tmp.gradients.i_t(1:(Nadc-extra),iNs); % [Hz/m]
end

if plotTest2 == 'True'
    figure()
    plot(-1*tGx/gamma*1e3,'b')
    ylabel(['Gradient (mT/m)'],'fontsize',15)    
    hold on
    plot(-1*tGy/gamma*1e3,'r')
    ylabel(['Gradient (mT/m)'],'fontsize',15)    
end

cd(file_folder)

% ... part 2.4 - get k-space trajectory ktraj_x e ktraj_y from the gradients ...
int_mat  = triu(ones((Nadc-extra)));
tKtraj_x = zeros((Nadc-extra),Npulses);
tKtraj_y = zeros((Nadc-extra),Npulses);
% Have to split the integral because each of the gradients is independent. 
% Otherwise we would have a continuous integral
for iNs=1:Npulses 
    tKtraj_x(:,iNs) = tmp.gradients.i_t(1:(Nadc-extra),iNs)' * int_mat * dt;
    tKtraj_y(:,iNs) = tmp.gradients.i_t(1:(Nadc-extra),iNs)' * int_mat * dt;
end

if plotTest2 == 'True'
    figure()
    plot(tKtraj_x,'b')
end

fprintf('\n\n 2 - Sucessfully finished - Define Parametres\n\n')


%%  part 3 - Read Input

% ... part 3.1 - read input in frequency domain ...

input            = load(['GIRFOS_tool/Matlab_basedCode/Example/girf_pulse_', num2str(Npulses),'_DT_1_1e-5_p', num2str(p),'.mat']);
input_freq       = load(['GIRFOS_tool/Matlab_basedCode/Example/girf_',type,'_Input_freq_', num2str(Npulses),'_tTotal_',num2str(tTotal),'ms.mat']);        
input_sincF_freq = load(['GIRFOS_tool/Matlab_basedCode/Example/girf_',type,'_Input_freq_sincFunct_', num2str(Npulses),'.mat']);

if plotTest3 == 'True'
    figure()
    plot(abs(input_freq.i_tFFT))
end

fprintf('\n\n 3 - Sucessfully finished - Read Input\n\n')

%%  part 4 - Read Output
cd(file_folder)
kdata = dat2mat(fn);          % [Nadc, nCoils, nsl&nbv&Nreps]

[NadcK,ncoils,Nreps_Naxes_Nobs] = size(kdata);        % Size - [Nadc, nCoils, N_repetitions & N_axes & N_observations + N_dummies]
nimgs = Nreps*Naxes*Nobs;
aux_k = reshape(kdata(:,:,1+tDum*nDummies:end),[NadcK,ncoils,Ny,nsl,nimgs*Npulses]);

if Naxes == 1
    kx = aux_k(:,:,:,:,:);                              % get ktraj_x
    
    k.kx_Pgrad_Px0 = kx(:,:,:,:,2:4:(nimgs)*Npulses);   % get ktraj_x GxPos & X0pos
    k.kx_Ngrad_Px0 = kx(:,:,:,:,1:4:(nimgs)*Npulses);   % get ktraj_x GxNeg & X0pos
    k.kx_Pgrad_Nx0 = kx(:,:,:,:,4:4:(nimgs)*Npulses);   % get ktraj_x GxPos & X0neg
    k.kx_Ngrad_Nx0 = kx(:,:,:,:,3:4:(nimgs)*Npulses);   % get ktraj_x GxNeg & X0neg 
    
elseif Naxes == 2
    kx = aux_k(:,:,:,:,1:nimgs*Npulses/2);                  % get ktraj_x
    ky = aux_k(:,:,:,:,nimgs*Npulses/2+1:end);              % get ktraj_y
    
    k.kx_Pgrad_Px0 = kx(:,:,:,:,2:4:(nimgs/2)*Npulses);     % get ktraj_x GxPos & X0pos
    k.kx_Ngrad_Px0 = kx(:,:,:,:,1:4:(nimgs/2)*Npulses);     % get ktraj_x GxNeg & X0pos
    k.kx_Pgrad_Nx0 = kx(:,:,:,:,4:4:(nimgs/2)*Npulses);     % get ktraj_x GxPos & X0neg 
    k.kx_Ngrad_Nx0 = kx(:,:,:,:,3:4:(nimgs/2)*Npulses);     % get ktraj_x GxNeg & X0neg 
    
    k.ky_Pgrad_Px0 = ky(:,:,:,:,2:4:(nimgs/2)*Npulses);     % get ktraj_y GyPos & X0pos
    k.ky_Ngrad_Px0 = ky(:,:,:,:,1:4:(nimgs/2)*Npulses);     % get ktraj_y GyNeg & X0pos
    k.ky_Pgrad_Nx0 = ky(:,:,:,:,4:4:(nimgs/2)*Npulses);     % get ktraj_y GyPos & X0neg 
    k.ky_Ngrad_Nx0 = ky(:,:,:,:,3:4:(nimgs/2)*Npulses);     % get ktraj_y GyNeg & X0neg     
    
elseif Naxes == 3
    kx = aux_k(:,:,:,:,1:nimgs*Npulses/3);                      % get ktraj_x
    ky = aux_k(:,:,:,:,nimgs*Npulses/3+1:2*nimgs*Npulses/3);    % get ktraj_y
    kz = aux_k(:,:,:,:,2*nimgs*Npulses/3+1:end);                % get ktraj_z
    
    k.kx_Pgrad_Px0 = kx(:,:,:,:,2:4:(nimgs/3)*Npulses);         % get ktraj_x GxPos & X0pos
    k.kx_Ngrad_Px0 = kx(:,:,:,:,1:4:(nimgs/3)*Npulses);         % get ktraj_x GxNeg & X0pos
    k.kx_Pgrad_Nx0 = kx(:,:,:,:,4:4:(nimgs/3)*Npulses);         % get ktraj_x GxPos & X0neg 
    k.kx_Ngrad_Nx0 = kx(:,:,:,:,3:4:(nimgs/3)*Npulses);         % get ktraj_x GxNeg & X0neg     
    
    k.ky_Pgrad_Px0 = ky(:,:,:,:,2:4:(nimgs/3)*Npulses);         % get ktraj_y GyPos & X0pos
    k.ky_Ngrad_Px0 = ky(:,:,:,:,1:4:(nimgs/3)*Npulses);         % get ktraj_y GyNeg & X0pos
    k.ky_Pgrad_Nx0 = ky(:,:,:,:,4:4:(nimgs/3)*Npulses);         % get ktraj_y GyPos & X0neg 
    k.ky_Ngrad_Nx0 = ky(:,:,:,:,3:4:(nimgs/3)*Npulses);         % get ktraj_y GyNeg & X0neg     
    
    k.kz_Pgrad_Px0 = kz(:,:,:,:,2:4:(nimgs/3)*Npulses);         % get ktraj_z GzPos & X0pos
    k.kz_Ngrad_Px0 = kz(:,:,:,:,1:4:(nimgs/3)*Npulses);         % get ktraj_z GzNeg & X0pos
    k.kz_Pgrad_Nx0 = kz(:,:,:,:,4:4:(nimgs/3)*Npulses);         % get ktraj_z GzPos & X0neg 
    k.kz_Ngrad_Nx0 = kz(:,:,:,:,3:4:(nimgs/3)*Npulses);         % get ktraj_z GzNeg & X0neg         
end


if plotTest4 == 'True'
    testPlt1 = double(abs(reshape(k.kx_Pgrad_Px0(:,12,:,:,:),size(k.kx_Pgrad_Px0,1),size(k.kx_Pgrad_Px0,5))));
    testPlt2 = double(abs(reshape(k.ky_Pgrad_Px0(:,12,:,:,:),size(k.ky_Pgrad_Px0,1),size(k.ky_Pgrad_Px0,5))));   
    phaseX   = unwrap(angle(reshape(k.kx_Pgrad_Nx0(:,12,:,:,:),size(k.kx_Pgrad_Px0,1),size(k.kx_Pgrad_Px0,5))));
    phaseY   = unwrap(angle(reshape(k.ky_Pgrad_Nx0(:,12,:,:,:),size(k.ky_Pgrad_Px0,1),size(k.ky_Pgrad_Px0,5))));
    
    figure()
    subplot(2,2,1)
    plot(testPlt1)
    hold on
    title('kx B0ECC-Mag')        
    subplot(2,2,2)
    plot(testPlt2)
    title('ky B0ECC-Mag')   

    subplot(2,2,3)
    plot(phaseX)
    hold on
    title('kx B0ECC-Phase')
    subplot(2,2,4)
    plot(phaseY)    
    title('ky B0ECC-Phase')    
end

fprintf('\n\n 4 - Sucessfully finished -  Read Output\n\n')


%% part 5 - Coil Compression using BART toolbox

% ... 5.0.1 - Compress coils  ...
cd(file_folder_cfl_hdr)
nc = 1;  % number of coils desired
k_cc = compressCoils_BOECC(k,NadcK,Ny,nsl,nc,Npulses,Nreps,Naxes,file_folder);

% ... 5.0.2 - save data compressed ...
cd(file_folder_results)
save(['Data_test',num2str(test),'.mat'],'k_cc','Naxes')
cd(file_folder)

% ... 5.1 - Figures ...
if plotTest5 == 'True'
    % ... 5.1.1 - Magnitude ...
    figure()
    subplot(2,2,1)
    plot(abs(k_cc.kx_Pgrad_Nx0))
    hold on
    title('kx_c_c B0ECC-Mag')
    subplot(2,2,2)
    plot(abs(k_cc.ky_Pgrad_Nx0))
    hold on
    title('ky_c_c B0ECC-Mag')
    
    % ... 5.1.2 - Angle/Phase ...
    subplot(2,2,3)
    plot(unwrap(angle(k_cc.kx_Pgrad_Nx0)))
    hold on
    title('kx_c_c B0ECC-Phase')
    subplot(2,2,4)
    plot(unwrap(angle(k_cc.ky_Pgrad_Nx0)))
    hold on
    title('ky_c_c B0ECC-Phase')
   
end

fprintf('\n\n 5 - Sucessfully finished - Coil Compression using BART toolbox \n\n')


%% part 6 - Measuring H0 & H1 according to Robison, 2019

% ... part 6.1 - Get Phase, Ktraj & Gout of system - figure 1 - Robison, 2019 ...
[~,~,~, phase_B0ECC, ktraj_out_B0ECC, Gout_B0ECC] = auxCalc_BOECC(k_cc,D,Npulses,Nreps,Naxes,dt);

% ... part 6.3 - Figures phase & Trajectory ...
if plotTest6 == 'True'
    figure();
    for jj = 1: Nreps*Npulses
        % ... 6.2.1 - Phase ...
        subplot(2,3,1)
        plot(phase_B0ECC.x(:,jj),'b')
        hold on
        plot(phase_B0ECC.y(:,jj),'r')
        hold on
        ylabel(['Phase_B_0 (rad)'],'fontsize',15)
        title('Phase B0ECC - x (b) y (r)')
        % ... 6.2.2 - Ktraj ...
        subplot(2,3,2)
        plot(-1*ktraj_out_B0ECC.x(:,jj),'b')
        hold on
        plot(-1*ktraj_out_B0ECC.y(:,jj),'r')
        hold on
        ylabel(['Gradient (Hz/m)'],'fontsize',15)
        title('ktraj B0ECC - x (b) y (r)')
        % ... 6.2.3 - Gout ...
        subplot(2,3,3)
        timeG = [0:dt:dt*(size(Gout_B0ECC.x,1)-1)]*1e3;
        plot(timeG,-1*Gout_B0ECC.x(:,jj),'b')
        hold on
        plot(timeG,-1*Gout_B0ECC.y(:,jj),'r')
        hold on
        title('Gout B0ECC - x (b) y (r)')
        ylabel(['Gradient (mT/m)'],'fontsize',15)
        xlabel(['time (ms)'],'fontsize',15)
    end
    for jj = 1: Npulses        
        % ... 6.2.4 - theoretical plots ...
        subplot(2,3,6)
        plot(tGx(:,jj)/gammabar*1e3,'r')
        hold on
        title('theor G')
        ylabel(['Gradient (mT/m)'],'fontsize',15)
        subplot(2,3,5)
        plot(tKtraj_x,'b')
        hold on
        title('theor K')
        ylabel(['Gradient (Hz/m)'],'fontsize',15)
    end
end

fprintf('\n\n 6 - Sucessfully finished - Measuring ktraj & Gout for xx and yy\n\n')


%% part 7 - Obtaining the impulse response function of the system 

Gin    = tGx/gammabar*1e3;  % mT/m
filter = 'True';
% % % % filter = 'Fals';

% ... part 7.1 - Reshape phase_B0ECC & Gout_B0ECC ...
phase_B0ECC_resp.x = zeros(size(phase_B0ECC.x,1),Npulses,Nreps);
phase_B0ECC_resp.y = zeros(size(phase_B0ECC.y,1),Npulses,Nreps);
Gout_B0ECC_resp.x  = zeros(size(Gout_B0ECC.x,1),Npulses,Nreps);
Gout_B0ECC_resp.y  = zeros(size(Gout_B0ECC.y,1),Npulses,Nreps);
for nr=1:Nreps
    phase_B0ECC_resp.x(:,:,nr) = phase_B0ECC.x(:,nr:double(Nreps):end);
    phase_B0ECC_resp.y(:,:,nr) = phase_B0ECC.y(:,nr:double(Nreps):end);
    Gout_B0ECC_resp.x(:,:,nr)  = Gout_B0ECC.x(:,nr:double(Nreps):end);
    Gout_B0ECC_resp.y(:,:,nr)  = Gout_B0ECC.y(:,nr:double(Nreps):end);
end

% ... part 7.2 - Get H0 & H1 ...
for nr=1:Nreps
    if filter == 'True'
        phase_B0ECC_test.order = 3;
        phase_B0ECC_test.wsize = 51;
    end
    
    phase_B0ECC_test.x = phase_B0ECC_resp.x(:,:,nr);
    phase_B0ECC_test.y = phase_B0ECC_resp.y(:,:,nr);
    
    Gout_B0ECC_test.x  = Gout_B0ECC_resp.x(:,:,nr);
    Gout_B0ECC_test.y  = Gout_B0ECC_resp.y(:,:,nr);
    
    % --- 7.1.1 - Get H0 ...
    H0_w = getH0_BOECC(phase_B0ECC_test,input_freq,Gin,Npulses,Naxes,dt,T,plotTest7,MATLABfft,filter);
    % --- 7.1.2 - Get H1 ...
    H1_w = getH1_BOECC(Gout_B0ECC_test,input_freq,Gin,Npulses,Nreps,Naxes,dt,T,plotTest7,MATLABfft);
    
    clear phase_B0ECC_test Gout_B0ECC_test
    
    H0_w_x(:,nr) = H0_w.x;
    H0_w_y(:,nr) = H0_w.y;
    
    H1_w_x(:,nr) = H1_w.x;
    H1_w_y(:,nr) = H1_w.y;
    
    
end

% ... part 7.3 - save Structure
H.H0_w = H0_w;
H.H1_w = H1_w;

H.H0_w.x = mean(H0_w_x,2);
H.H0_w.y = mean(H0_w_y,2);
H.H1_w.x = mean(H1_w_x,2);
H.H1_w.y = mean(H1_w_y,2);

H_plot.H0_w = H0_w;
H_plot.H1_w = H1_w;

H_plot.H0_w.x = H0_w_x;
H_plot.H0_w.y = H0_w_y;
H_plot.H1_w.x = H1_w_x;
H_plot.H1_w.y = H1_w_y;

fprintf('\n\n 7 - Sucessfully finished -  Obtaining the IRF of the system \n\n')


%% part 8 - plots of H0 & H1
if plotTest8 == 'True'
% % %     % Plot H0
% % %     figure()
% % %     for nr=1:Nreps
% % %         subplot(5,3,nr)
% % %         hold on
% % %         plot(H_plot.H0_w.ff*1e-3,abs(H_plot.H0_w.x(:,nr)),'b')
% % %         ylabel(['B_0 EC Response'],'fontsize',15)
% % %         xlabel(['frequency (kHz)'],'fontsize',15)
% % %         hold on
% % %         if filter == 'True'
% % %             title(['H0_x & H0_y (filter - order:',num2str(H_plot.H0_w.FilterOrder),' Wsize:',num2str(H_plot.H0_w.FilterWsize),')']);
% % %         else
% % %             title(['H0_x & H0_y (no filter)']);
% % %         end
% % %         plot(H_plot.H0_w.ff*1e-3,abs(H_plot.H0_w.y(:,nr)),'r')
% % %         legend('H_0_x','H_0_y')
% % %     end
% % %     
% % %     % Plot H1    
% % %     figure()
% % %     for nr=1:Nreps
% % %         subplot(5,3,nr)        
% % %         hold on
% % %         plot(H_plot.H1_w.ff*1e-3,abs(H_plot.H1_w.x(:,nr)),'b')
% % %         ylabel(['Gradient Field Response'],'fontsize',15)
% % %         xlabel(['frequency (kHz)'],'fontsize',15)
% % %         hold on
% % %         plot(H_plot.H1_w.ff*1e-3,abs(H_plot.H1_w.y(:,nr)),'r')
% % %         legend('H_1_x','H_1_y')
% % %     end
    
    % --- mean plots ----
    % Plot H0
    figure()
    hold on
    plot(H_plot.H0_w.ff*1e-3,abs(mean(H_plot.H0_w.x,2)),'b')
    ylabel(['B_0 EC Response'],'fontsize',15)
    xlabel(['frequency (kHz)'],'fontsize',15)
    hold on
    if filter == 'True'
        title(['H0_x & H0_y (filter - order:',num2str(H_plot.H0_w.FilterOrder),' Wsize:',num2str(H_plot.H0_w.FilterWsize),')']);
    else
        title(['H0_x & H0_y (no filter)']);
    end
    plot(H_plot.H0_w.ff*1e-3,abs(mean(H_plot.H0_w.y,2)),'r')
% %     plot(H_plot.H0_w.ff*1e-3,abs(H_plot.H0_w.x(:,1)),'g')
    legend('H_0_x','H_0_y')
% %     legend('H_0_x','H_0_y','H_0_x_-_r_e_p_1')
    
    % Plot H1
    figure()
%     title(['H1_x & H1_y']);
    title(['GIRF M_2']);
    hold on
    plot(H_plot.H1_w.ff*1e-3,abs(mean(H_plot.H1_w.x,2)),'b')
    ylabel(['GIRF magnitude'],'fontsize',15)
    xlabel(['frequency (kHz)'],'fontsize',15)
    hold on
    plot(H_plot.H1_w.ff*1e-3,abs(mean(H_plot.H1_w.y,2)),'r')
    hold on
%     plot(H_plot.H1_w.ff*1e-3,abs(H_plot.H1_w.x(:,1)),'g')
    legend('H_1_x','H_1_y')
% %     legend('H_1_x','H_1_y','H_1_x_-_r_e_p_1')
    
end

% clear H_plot
fprintf('\n\n 8 - Sucessfully finished - Plots BO_ECC H0&H1\n\n')

%% part 9 - Save IRF
if saveT == 'True'
    if PC == 1
        cd(file_folder_results)
        save(['test',num2str(test),'_B0ECC_H.mat'],'H')
        save(['test',num2str(test),'_IRF.mat'],'H')
        cd(['/Results_IRF' ])
        cd(file_folder)        
        fprintf('\n\n 8 - Sucessfully finished - Save IRF\n\n')
    end
else
    fprintf('\n\n 9 - Sucessfully finished - Did NOT Save BO_ECC H0&H1\n\n')
end

% close all
