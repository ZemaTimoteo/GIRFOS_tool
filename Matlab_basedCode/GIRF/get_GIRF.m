%% Recon of images from scanner with pulseq.seq for GIRF sequence
% TTFernandes, IST 2020

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
    addpath(genpath('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/bart-master'))           % To read
    addpath(genpath('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF'))
    addpath(genpath('GIRF')) % To read    
    addpath(genpath('Aux_Functions'));  % To reconstruct
end


% close all
clc


%% part 1 - Get Files

if computer == 1
    folder = 'D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\3_Data_scan\SPIRAL_measurement\';
    folder = [folder, 'SPIRAL_test_3\'];
else
    myCD = '/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/2_GIRF_test/Test';
    cd(myCD);
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
end

[files] = dir(file_folder);
fn = files(5).name;       % get '.dat' file
fn_info = files(6).name;  % get '.mat' file with sequence parameters
testName = files(7).name; % get '.seq' name
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

saveT       = 'True';                      % Save - 'True' or 'Fals' 
type        = 'SCAN';                      %'PULS' - test27 & test45 or 'SCAN' - test42

testPCA      = 1;                           % Test PCA independent xx and yy - '1' OR all concatenated - '2' OR all independent - '0'
meanPCA      = 0;                           % Test PCA with mean - '1' OR without mean - '0'
nc           = 1;                           % number of coils for PCA

plotTest2 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest4 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest5 = 'Fals';                      % Plot Test - 'True' or 'Fals'
plotTest6 = 'True';                      % Plot Test - 'True' or 'Fals'
plotTest7 = 'True';                      % Plot Test - 'True' or 'Fals'
test6     = 'Fals';                      % Test for debbuging

% ... part 2.2 - parameters ....
test         = tmp.test                     % test
Npulses      = double(tmp.Npulses);         % Number of triangular pulses pulses
Nreps        = tmp.nreps;                   % get Repetitions 
nsl          = tmp.nslices;                 % number of slices
Ny           = 1;                           % for reshape of k vector
dirs         = tmp.ndirs;                   % number of directions
Nadc         = tmp.adc_samples;             % number of ADCs
dt           = tmp.DT;                      % dt (s)
Naxes        = 2;                           % Number of axes
Nobs         = 2;                           % Number of observation - 2: one w/ gradient one without gradient
gamma        = 42.576e6;                    % Gyro-magnetic ratio (Hz/T)
sliceT       = tmp.st;                      % Slice thickness (m)
D            = tmp.x0;                      % Distance of the slice to gradient isocenter (m)
p            = tmp.max_slew;                % Max_slew rate (mT/s)
T            = linspace(0.05,0.16,Npulses); % time range - (ms)
gammabar     = gamma;                       % Gyro-magnetic ratio Hz/T

moreplots = 1;
ADCest    = 0;
extra     = 2;  % retirar os dois pontos que o .seq tem a menos em ADC


if testDummy == 'True'
    tDum = 1;
else
    tDum = 0;
end

% ... part 2.3 - get theoretical Gx and Gy ...
cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/GIRF/pulses');
if type == 'PULS'
    tmp.gradients = load(['girf_pulse_' num2str(Npulses),'_DT_1_1e-5_p', num2str(p)]);
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
    plot(tGx(:,1),'b')
    hold on
    plot(tGy(:,2),'r')
end

cd(file_folder)

% ... part 2.4 - get k-space trajectory ktraj_x e ktraj_y from the gradients ...
int_mat = triu(ones((Nadc-extra)));
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
    plot(tKtraj_x(:,1),'b')
    hold on
    plot(tKtraj_y(:,12),'r')
end

fprintf('\n\n 2 - Sucessfully finished - Define Parametres\n\n')


%%  part 3 - Read Input

% ... part 3.1 - read input in frequency domain ...

if computer == 1
    input            = load(['D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\5_GIRF\1_matlab_code\pulses\girf_pulse_', num2str(Npulses),'_DT_1_1e-5_p', num2str(p),'.mat']);
    input_freq       = load(['D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\5_GIRF\1_matlab_code\pulses\girf_Input_freq_', num2str(Npulses),'.mat']);
    input_sincF_freq = load(['D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\5_GIRF\1_matlab_code\pulses\girf_Input_freq_sincFunct_', num2str(Npulses),'.mat']);
else
    input            = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_pulse_', num2str(Npulses),'_DT_1_1e-5_p', num2str(p),'.mat']);
    if type == 'PULS'
        input_freq   = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_',type,'_Input_freq_', num2str(Npulses),'.mat']);
    elseif type == 'SCAN'
        input_freq   = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_',type,'_Input_freq_', num2str(Npulses),'_tTotal_',num2str(tTotal),'ms.mat']);        
    end
    input_sincF_freq = load(['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF/pulses/girf_',type,'_Input_freq_sincFunct_', num2str(Npulses),'.mat']);
end

fprintf('\n\n 3 - Sucessfully finished - Read Input\n\n')

%%  part 4 - Read Output
cd(file_folder)
k     = dat2mat(fn);          % [Nadc, nCoils, nsl&nbv&Nreps]

[NadcK,ncoils,Nreps_Naxes_Nobs] = size(k);        % Size - [Nadc, nCoils, N_repetitions & N_axes & N_observations = 2]
nreps = floor(size(k,3)/nsl/Npulses/Naxes/Nobs);  % Number of images acquired, counting the calibration. Dividir por dois pq para cada imagem é preciso uma com e sem gradiente
nimgs = nreps*Naxes*Nobs;
aux_k = reshape(k(:,:,1+tDum:end),[NadcK,ncoils,Ny,nsl,nimgs*Npulses]);

if Naxes == 1
    kx = aux_k(:,:,:,:,1:(nimgs*nreps)/Naxes);          % get ktraj_x
    kx_noGrad = aux_k(:,:,:,:,1:2:end);                 % get ktraj_x with no gradient
    kx_grad   = aux_k(:,:,:,:,2:2:end);                 % get ktraj_x with gradient
elseif Naxes == 2
    kx = aux_k(:,:,:,:,1:nimgs*Npulses/2);                              % get ktraj_x
    ky = aux_k(:,:,:,:,nimgs*Npulses/2:end);                          % get ktraj_y
    kx_noGrad = aux_k(:,:,:,:,1:2:(nimgs/2)*Npulses);         % get ktraj_x with no gradient
    kx_grad   = aux_k(:,:,:,:,2:2:(nimgs/2)*Npulses);         % get ktraj_x with gradient
    ky_noGrad = aux_k(:,:,:,:,(nimgs/2)*Npulses+1:2:end);     % get ktraj_x with no gradient
    ky_grad   = aux_k(:,:,:,:,(nimgs/2)*Npulses+2:2:end);     % get ktraj_x with gradient
elseif Naxes == 3
    kx = aux_k(:,:,:,:,1:nimgs/Naxes);                           % get ktraj_x
    ky = aux_k(:,:,:,:,nimgs/Naxes+1:(nimgs/Naxes)*2);           % get ktraj_y
    kz = aux_k(:,:,:,:,(nimgs/Naxes)*2+1:end);                   % get ktraj_z
    kx_noGrad = aux_k(:,:,:,:,1:2:nimgs/Naxes);                  % get ktraj_x with no gradient
    kx_grad   = aux_k(:,:,:,:,2:2:nimgs/Naxes);                  % get ktraj_x with gradient
    ky_noGrad = aux_k(:,:,:,:,nimgs/Naxes+1:2:(nimgs/Naxes)*2);  % get ktraj_x with no gradient
    ky_grad   = aux_k(:,:,:,:,nimgs/Naxes+2:2:(nimgs/Naxes)*2);  % get ktraj_x with gradient
    kz_noGrad = aux_k(:,:,:,:,(nimgs/Naxes)*2+1:2:end);          % get ktraj_x with no gradient
    kz_grad   = aux_k(:,:,:,:,(nimgs/Naxes)*2+2:2:end);          % get ktraj_x with gradient
end


if plotTest4 == 'True'
               
% %     subplot(221)
% %     plot(abs(kx_noGrad(:,9,1,1,1)))
% %     subplot(222)
% %     plot(abs(ky_noGrad(:,9,1,1,1)))
% %     subplot(223)
% %     plot(abs(kx_grad(:,9,1,1,1)))
% %     subplot(224)
% %     plot(abs(ky_grad(:,9,1,1,1)))
% %     
% %     figure()
% %     subplot(221)
% %     plot(angle(kx_noGrad(:,1,1,1,1)))
% %     subplot(222)
% %     plot(angle(ky_noGrad(:,1,1,1,1)))
% %     subplot(223)
% %     plot(angle(kx_grad(:,1,1,1,1)))
% %     subplot(224)
% %     plot(angle(ky_grad(:,1,1,1,1)))
% %     
    % set all figures

    
% %     figure('name',['Figure MAG - kx_noGrad']);
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(abs(kx_noGrad(:,nc,1,1,1)),'r')
% %         hold on
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure MAG - ky_noGrad']);
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(abs(ky_noGrad(:,nc,1,1,1)),'g')
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure MAG - kx_grad'])
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(abs(kx_grad(:,nc,1,1,1)),'b')
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure MAG - ky_grad'])
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(abs(ky_grad(:,nc,1,1,1)),'y')
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure Phase - kx_noGrad'])
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(unwrap(angle(kx_noGrad(:,nc,1,1,1))),'r')
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure Phase - ky_noGrad'])
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(unwrap(angle(ky_noGrad(:,nc,1,1,1))),'g')
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure Phase - kx_grad'])
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(unwrap(angle(kx_grad(:,nc,1,1,1))),'b')
% %         title(['channel ',num2str(nc)])
% %     end
% %     figure('name',['Figure Phase - ky_gra    mean_kx_noGrad = mean(kx_noGrad,5);
    mean_ky_noGrad = mean(ky_noGrad,5);
    mean_kx_grad   = mean(kx_grad,5);
    mean_ky_grad   = mean(ky_grad,5);
% %     for  nc = 1:ncoils
% %         subplot(2,ncoils/2,nc)
% %         plot(unwrap(angle(ky_grad(:,nc,1,1,1))),'y')
% %         title(['channel ',num2str(nc)])
% %     end

    % ... test ktraj from coils ...
    figure('name',['Figure ktraj_measure_coils kx']);
    for  nc = 1:ncoils
        ddd = unwrap(angle(mean(kx_noGrad(:,nc,:,:,:),5)));
        eee = unwrap(angle(mean(kx_grad(:,nc,:,:,:),5)));
        fff = (eee - ddd)/D;
        subplot(2,ncoils/2,nc)
        plot(fff,'r')
        title(['channel ',num2str(nc)])
        clear ddd eee fff
    end
    figure('name',['Figure ktraj_measure_coils ky']);
    for  nc = 1:ncoils
        ddd = unwrap(angle(mean(kx_noGrad(:,nc,:,:,:),5)));
        eee = unwrap(angle(mean(kx_grad(:,nc,:,:,:),5)));
        fff = (eee - ddd)/D;
        subplot(2,ncoils/2,nc)
        plot(fff,'b')
        title(['channel ',num2str(nc)])
        clear ddd eee fff
    end

end

fprintf('\n\n 4 - Sucessfully finished -  Read Output\n\n')


%% part 5 - Coil Compression using BART toolbox

% ... part 5.0.0 - select number of coils after PCA ...
cd(file_folder_cfl_hdr)
nv=size(aux_k,5);          % number of Volumes in the scanner

if meanPCA == 0    
    kx_ccswp_noGrad = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_noGrad = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_grad   = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_grad   = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    kx_swp_noGrad = permute(kx_noGrad,[1 3 4 2 5]); %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_noGrad = permute(ky_noGrad,[1 3 4 2 5]); %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_grad   = permute(kx_grad,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_grad   = permute(ky_grad,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    
    % ... part 5.1 - compression + concatenation ...
    kx_l_conc = [kx_swp_noGrad(:,:,:,:,1); kx_swp_grad(:,:,:,:,1)];
    ky_l_conc = [ky_swp_noGrad(:,:,:,:,1); ky_swp_grad(:,:,:,:,1)];
    
    % create ccoef for x-axis
    writecfl('kx_l_conc',kx_l_conc);
    system(sprintf('bart cc -p %d -M kx_l_conc x_cccoef_conc',nc));
    
    % create ccoef for y-axis
    writecfl('ky_l_conc',ky_l_conc);
    system(sprintf('bart cc -p %d -M ky_l_conc y_cccoef_conc',nc));
    
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        kx_tmp_noGrad = kx_swp_noGrad(:,:,:,:,v);
        writecfl('kx_tmp_noGrad',kx_tmp_noGrad);
        system(sprintf('bart ccapply -p %d kx_tmp_noGrad x_cccoef_conc kx_tmpcc_noGrad',nc));
        kx_ccswp_noGrad(:,:,:,:,v)=readcfl('kx_tmpcc_noGrad');
        
        ky_tmp_noGrad = ky_swp_noGrad(:,:,:,:,v);
        writecfl('ky_tmp_noGrad',ky_tmp_noGrad);
        system(sprintf('bart ccapply -p %d ky_tmp_noGrad y_cccoef_conc ky_tmpcc_noGrad',nc));
        ky_ccswp_noGrad(:,:,:,:,v)=readcfl('ky_tmpcc_noGrad');
        
        kx_tmp_grad=kx_swp_grad(:,:,:,:,v);
        writecfl('kx_tmp_grad',kx_tmp_grad);
        system(sprintf('bart ccapply -p %d kx_tmp_grad x_cccoef_conc kx_tmpcc_grad',nc));
        kx_ccswp_grad(:,:,:,:,v)=readcfl('kx_tmpcc_grad');
        
        ky_tmp_grad=ky_swp_grad(:,:,:,:,v);
        writecfl('ky_tmp_grad',ky_tmp_grad);
        system(sprintf('bart ccapply -p %d ky_tmp_grad y_cccoef_conc ky_tmpcc_grad',nc));
        ky_ccswp_grad(:,:,:,:,v)=readcfl('ky_tmpcc_grad');
        
        clear kx_tmpcc_noGrad ky_tmpcc_noGrad kx_tmpcc_grad ky_tmpcc_grad
    end
    
    cd(file_folder)
    
    % ... part 5.3 - invert permutation ...
    kx_cc_noGrad = ipermute(kx_ccswp_noGrad,[1 3 4 2 5]);
    ky_cc_noGrad = ipermute(ky_ccswp_noGrad,[1 3 4 2 5]);
    kx_cc_grad   = ipermute(kx_ccswp_grad,[1 3 4 2 5]);
    ky_cc_grad   = ipermute(ky_ccswp_grad,[1 3 4 2 5]);
    
    kx_cc_noGrad = reshape(kx_cc_noGrad,[NadcK nc Npulses*Nreps]);
    ky_cc_noGrad = reshape(ky_cc_noGrad,[NadcK nc Npulses*Nreps]);
    kx_cc_grad   = reshape(kx_cc_grad,[NadcK nc Npulses*Nreps]);
    ky_cc_grad   = reshape(ky_cc_grad,[NadcK nc Npulses*Nreps]);
    
elseif meanPCA == 1
    % ... part 5.1 - compr elseif meanPCA == 1
    mean_kx_noGrad = mean(kx_noGrad,5);
    mean_ky_noGrad = mean(ky_noGrad,5);
    mean_kx_grad   = mean(kx_grad,5);
    mean_ky_grad   = mean(ky_grad,5);
    
    kx_ccswp_noGrad = zeros(NadcK, Ny, nsl, nc, Npulses);
    ky_ccswp_noGrad = zeros(NadcK, Ny, nsl, nc, Npulses);
    kx_ccswp_grad   = zeros(NadcK, Ny, nsl, nc, Npulses);
    ky_ccswp_grad   = zeros(NadcK, Ny, nsl, nc, Npulses);
    
    % ... part 5.0.1 - permute for bart ...
    kx_swp_noGrad = permute(mean_kx_noGrad,[1 3 4 2]); %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_noGrad = permute(mean_ky_noGrad,[1 3 4 2]); %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_grad   = permute(mean_kx_grad,[1 3 4 2]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_grad   = permute(mean_ky_grad,[1 3 4 2]);   %Nx Ny Nslice Ncoils Nvolumes
    
    % ... part 5.1 - compression + concatenation ...
    kx_l_conc = [kx_swp_noGrad(:,:,:,:); kx_swp_grad(:,:,:,:)];
    ky_l_conc = [ky_swp_noGrad(:,:,:,:); ky_swp_grad(:,:,:,:)];
    
    % create ccoef for x-axis
    writecfl('kx_l_conc',kx_l_conc);
    system(sprintf('bart cc -p %d -M kx_l_conc x_cccoef_conc',nc));
    
    % create ccoef for y-axis
    writecfl('ky_l_conc',ky_l_conc);
    system(sprintf('bart cc -p %d -M ky_l_conc y_cccoef_conc',nc));
    
    
    % ... part 5.2 - apply compression ...
    kx_tmp_noGrad = kx_swp_noGrad;
    writecfl('kx_tmp_noGrad',kx_tmp_noGrad);
    system(sprintf('bart ccapply -p %d kx_tmp_noGrad x_cccoef_conc kx_tmpcc_noGrad',nc));
    kx_ccswp_noGrad=readcfl('kx_tmpcc_noGrad');
    
    ky_tmp_noGrad = ky_swp_noGrad;
    writecfl('ky_tmp_noGrad',ky_tmp_noGrad);
    system(sprintf('bart ccapply -p %d ky_tmp_noGrad y_cccoef_conc ky_tmpcc_noGrad',nc));
    ky_ccswp_noGrad=readcfl('ky_tmpcc_noGrad');
    
    kx_tmp_grad=kx_swp_grad;
    writecfl('kx_tmp_grad',kx_tmp_grad);
    system(sprintf('bart ccapply -p %d kx_tmp_grad x_cccoef_conc kx_tmpcc_grad',nc));
    kx_ccswp_grad=readcfl('kx_tmpcc_grad');
    
    ky_tmp_grad=ky_swp_grad;
    writecfl('ky_tmp_grad',ky_tmp_grad);
    system(sprintf('bart ccapply -p %d ky_tmp_grad y_cccoef_conc ky_tmpcc_grad',nc));
    ky_ccswp_grad=readcfl('ky_tmpcc_grad');
    
    clear kx_tmpcc_noGrad ky_tmpcc_noGrad kx_tmpcc_grad ky_tmpcc_grad
    
    cd(file_folder)
    
    % ... part 5.3 - invert permutation ...
    kx_cc_noGrad = ipermute(kx_ccswp_noGrad,[1 3 4 2]);
    ky_cc_noGrad = ipermute(ky_ccswp_noGrad,[1 3 4 2]);
    kx_cc_grad   = ipermute(kx_ccswp_grad,[1 3 4 2]);
    ky_cc_grad   = ipermute(ky_ccswp_grad,[1 3 4 2]);
    
    kx_cc_noGrad = reshape(kx_cc_noGrad,[NadcK nc Nreps]);
    ky_cc_noGrad = reshape(ky_cc_noGrad,[NadcK nc Nreps]);
    kx_cc_grad   = reshape(kx_cc_grad,[NadcK nc Nreps]);
    ky_cc_grad   = reshape(ky_cc_grad,[NadcK nc Nreps]);
    
end


if plotTest5 == 'True'
    % ... Magnitude ...
    figure('name',['Figure MAG Component - kx_noGrad']);
    for  nr = 1:Nreps
        if Nreps == 1
            plot(abs(kx_cc_noGrad),'r')
        else
            subplot(2,Nreps/2,nr)
            plot(abs(kx_cc_noGrad(:,:,nr)),'r')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    
    figure('name',['Figure MAG Component - ky_noGrad']);
    for  nr = 1:Nreps
        if Nreps == 1v
            plot(abs(ky_cc_noGrad),'g')
        else
            subplot(2,Nreps/2,nr)
            plot(abs(ky_cc_noGrad(:,:,nr)),'g')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    figure('name',['Figure MAG Component - kx_grad']);
    for  nr = 1:Nreps
        if Nreps == 1
            plot(abs(kx_cc_grad),'b')
        else
            subplot(2,Nreps/2,nr)
            plot(abs(kx_cc_grad(:,:,nr)),'b')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    figure('name',['Figure MAG Component - ky_grad']);
    for  nr = 1:Nreps
        if Nreps == 1    % ... part 5.1 - compression + concatenation ...
    kx_l_conc = [kx_swp_noGrad(:,:,:,:); kx_swp_grad(:,:,:,:)];
    ky_l_conc = [ky_swp_noGrad(:,:,:,:); ky_swp_grad(:,:,:,:)];
            plot(abs(ky_cc_grad),'y')
        else
            subplot(2,Nreps/2,nr)
            plot(abs(ky_cc_grad(:,:,nr)),'y')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    
    % ... Angle/Phase ...
    figure('name',['Figure PHA Component - kx_noGrad']);
    for  nr = 1:Nreps
        if Nreps == 1
            plot(unwrap(angle(kx_cc_noGrad)),'r')
        else    % ... part 5.1 - compression + concatenation ...
    kx_l_conc = [kx_swp_noGrad(:,:,:,:); kx_swp_grad(:,:,:,:)];
    ky_l_conc = [ky_swp_noGrad(:,:,:,:); ky_swp_grad(:,:,:,:)];
            subplot(2,Nreps/2,nr)
            plot(unwrap(angle(kx_cc_noGrad(:,:,nr))),'r')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    figure('name',['Figure PHA Component - ky_noGrad']);
    for  nr = 1:Nreps
        if Nreps == 1
            plot(unwrap(angle(ky_cc_noGrad)),'g')
        else
            subplot(2,Nreps/2,nr)
            plot(unwrap(angle(ky_cc_noGrad(:,:,nr))),'g')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    figure('name',['Figure PHA Component - kx_grad']);
    for  nr = 1:Nreps
        if Nreps == 1
            plot(unwrap(angle(kx_cc_grad)),'b')
        else       subplot(2,Nreps/2,nr)
            plot(unwrap(angle(kx_cc_grad(:,:,nr))),'b')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
    figure('name',['Figure PHA Component - ky_grad']);
    for  nr = 1:Nreps
        if Nreps == 1
            plot(unwrap(angle(ky_cc_grad)),'y')
        else
            subplot(2,Nreps/2,nr)
            plot(unwrap(angle(ky_cc_grad(:,:,nr))),'y')
            hold on
            title(['Rep ',num2str(nr)])
        end
    end
end

cd(file_folder_results)
save(['Data_test',num2str(test),'.mat'],'kx_cc_noGrad','ky_cc_noGrad','kx_cc_grad','ky_cc_grad')
cd(file_folder)

fprintf('\n\n 5 - Sucessfully finished - Coil Compression using BART toolbox \n\n')
if computer == 1
    load('')
end


%% part 6 - Measuring ktraj for xx and yy
% ... part 6.1 - Get Phase angle and its difference point-by-point ...
pkx_noGrad = zeros(NadcK,Npulses*Nreps);
pky_noGrad = zeros(NadcK,Npulses*Nreps);
pkx_grad   = zeros(NadcK,Npulses*Nreps);
pky_grad   = zeros(NadcK,Npulses*Nreps);
diff_pkx   = zeros(NadcK,Npulses*Nreps);
diff_pky   = zeros(NadcK,Npulses*Nreps);


% ----------- test for debbugging ------------------
if test6 == 'True'
    % phase
    for a=1:Npulses
        if plotTest6 == 'True'           
            figure(a)
        end
        for v=(a-1)*Nreps + 1 : a*Nreps
            test_kxNoGrad_phase(:,v) = unwrap(angle(kx_cc_noGrad(:,:,v)));
            test_kxGrad_phase(:,v)   = unwrap(angle(kx_cc_grad(:,:,v)));
            test_kyNoGrad_phase(:,v) = unwrap(angle(ky_cc_noGrad(:,:,v)));
            test_kyGrad_phase(:,v)   = unwrap(angle(ky_cc_grad(:,:,v)));
            
            if plotTest6 == 'True'
                subplot(2,2,1)
                hold on
                title('Measured test k_xNoGrad phase')
                plot(test_kxNoGrad_phase(:,v))
                subplot(2,2,2)
                hold on
                title('Measured test k_xGrad phase')
                plot(test_kxGrad_phase(:,v))
            end
            
            if plotTest6 == 'True'
                subplot(2,2,3)
                hold on
                title('Measured test k_yNoGrad phase')
                plot(test_kyNoGrad_phase(:,v))
                subplot(2,2,4)
                hold on
                title('Measured test k_yGrad phase')
                plot(test_kyGrad_phase(:,v))
            end
        end
    end
    
    % magnitude
    for a=1:Npulses
        if plotTest6 == 'True'
            figure(a+12)
        end
        for v=(a-1)*Nreps + 1 : a*Nreps
            test_kxNoGrad_magn(:,v) = abs(kx_cc_noGrad(:,:,v));
            test_kxGrad_magn(:,v)   = abs(kx_cc_grad(:,:,v));
            test_kyNoGrad_magn(:,v) = abs(ky_cc_noGrad(:,:,v));
            test_kyGrad_magn(:,v)   = abs(ky_cc_grad(:,:,v));
            
            if plotTest6 == 'True'
                subplot(2,2,1)
                hold on
                title('Measured test k_xNoGrad Magnitude')
                plot(test_kxNoGrad_magn(:,v))
                subplot(2,2,2)
                hold on
                title('Measured test k_xGrad Magnitude')
                plot(test_kxGrad_magn(:,v))
            end
            
            if plotTest6 == 'True'
                subplot(2,2,3)
                hold on
                title('Measured test k_yNoGrad Magnitude')
                plot(test_kyNoGrad_magn(:,v))
                subplot(2,2,4)
                hold on
                title('Measured test k_yGrad Magnitude')
                plot(test_kyGrad_magn(:,v))
            end
        end
    end
    
    % measuring T2*
    vectorT = [0:dt:size(kx_cc_noGrad,1)*dt]; % in s
    pointA = 1;
    pointB = 35;
    deltaT = vectorT(pointB) - vectorT(pointA);
    for a=1:Npulses
        for v=(a-1)*Nreps + 1 : a*Nreps            
            T2_starNoGrad_X(a) = (log(test_kxNoGrad_magn(pointA,v))-log(test_kxNoGrad_magn(pointB,v)))/deltaT;
            T2_starNoGrad_Y(a) = (log(test_kyNoGrad_magn(pointA,v))-log(test_kyNoGrad_magn(pointB,v)))/deltaT;
            T2_starGrad_X(a)   = (log(test_kxGrad_magn(pointA,v))-log(test_kxGrad_magn(pointB,v)))/deltaT;
            T2_starGrad_Y(a)   = (log(test_kyGrad_magn(pointA,v))-log(test_kyGrad_magn(pointB,v)))/deltaT;
        end
    end
    
    figure()
    plot(T2_starNoGrad_X,'b')
    hold on
    plot(T2_starNoGrad_Y,'g')
    title(['T2* values for noGrad - X & Y - pointA: ', num2str(vectorT(pointA)*1e3),'ms pointB: ',num2str(vectorT(pointB)*1e3),'ms'])
    ylabel('Time in (s)')
    xlabel('Pulses')

end
% ----------- FINNISH test for debbugging ------------------



for v=1:Npulses*Nreps
    pkx_noGrad(:,v) = unwrap(angle(kx_cc_noGrad(:,:,v)));
    pky_noGrad(:,v) = unwrap(angle(ky_cc_noGrad(:,:,v)));
    pkx_grad(:,v)   = unwrap(angle(kx_cc_grad(:,:,v)));
    pky_grad(:,v)   = unwrap(angle(ky_cc_grad(:,:,v)));
    
    diff_pkx(:,v) = pkx_grad(:,v) - pkx_noGrad(:,v);
    diff_pky(:,v) = pky_grad(:,v) - pky_noGrad(:,v);
end

% ... part 6.2 - Get ktra ...
ktraj_xCorr = (diff_pkx./abs(D))';
ktraj_yCorr = (diff_pky./abs(D))';

% ... part 6.3 - Get Gout ...
Gout_x = zeros(Nreps*Npulses,size(ktraj_xCorr,2)-1);
Gout_y = zeros(Nreps*Npulses,size(ktraj_yCorr,2)-1);
for v=1:Nreps*Npulses
    Gout_x(v,:) = (ktraj_xCorr(v,2:end)-ktraj_xCorr(v,1:end-1))/dt/(2*pi)/gammabar*1e3;
    Gout_y(v,:) = (ktraj_yCorr(v,2:end)-ktraj_yCorr(v,1:end-1))/dt/(2*pi)/gammabar*1e3;
end

% ... align with I_t ...
% % Gout_x = [zeros(Nreps*Npulses,1) Gout_x(:,1:end-1)];
% % Gout_y = [zeros(Nreps*Npulses,1) Gout_y(:,1:end-1)];

% ... part 6.4 - Plot Ktraj corrections ...
if plotTest6 == 'True'
% %     figure()
% %     subplot(2,1,1)
% %     plot(ktraj_xCorr','r')
% %     hold on
% %     title('Measured Ktraj_x')
% %     subplot(2,1,2)
% %     plot(ktraj_yCorr','b')
% %     hold on
% %     title('Measured Ktraj_y')    
% %     
% %     figure()    
% %     subplot(211)
% %     plot(Gout_x,'r')
% %     hold on
% %     title('Measured G_x')
% %     ylabel(['Gradient (mT/m)'],'fontsize',15)
% %     subplot(212)
% %     plot(Gout_y,'b')
% %     hold on
% %     title('Measured G_y')
% %     ylabel(['Gradient (mT/m)'],'fontsize',15)

    % --- theoretical plots ---
% %     figure()
% %     subplot(211)
% %     hold on
% %     title('theor Ktraj_x')
% %     plot((-1)*tKtraj_x,'r')
% %     subplot(212)
% %     plot((-1)*tKtraj_y,'b')
% %     hold on
% %     title('theor Ktraj_y')
% %     
% %     figure()
% %     subplot(211)
% %     hold on
% %     title('theor G_x')
% %     ylabel(['Gradient (mT/m)'],'fontsize',15)
% %     plot((-1)*tGx/gammabar*1e3,'r')
% %     subplot(212)
% %     plot((-1)*tGy/gammabar*1e3,'b')
% %     hold on
% %     title('theor G_y')
% %     ylabel(['Gradient (mT/m)'],'fontsize',15)    
    
% %     figure()
% %     for iii=1:Npulses
% %         subplot(1,Npulses,iii)
% %         plot(ktraj_xCorr(iii*Nreps-3+1,:),'r')
% %         hold on
% %         title('Measured Ktraj_x')
% %     end
% %     figure()
% %     for iii=1:Npulses        
% %         subplot(1,Npulses,iii)
% %         plot(ktraj_yCorr(iii*Nreps-3+1,:),'b')
% %         hold on
% %         title('Measured Ktraj_y')
% %     end

end


% Test - remove first pulse and 12 pulse & Average them out
if tMean == 'True'
    Gout_xAUX = Gout_x;
    Gout_yAUX = Gout_y;
    
    % Gout_x = Gout_xAUX(1:(Npulses-1)*Nreps,:);
    if testDummy == 'Fals'
        Gout_x = Gout_xAUX(2:end,:);
        Gout_y = Gout_yAUX(1:(Npulses-1)*Nreps,:);
    end
    
    mGout_x = zeros(size(Gout_x,1)/Nreps,size(Gout_y,2));
    mGout_y = zeros(size(Gout_y,1)/Nreps,size(Gout_y,2));
    
    for jj=1:Npulses
        if testDummy == 'Fals'            
            % mean for Gx
            if jj==1
                mGout_x(jj,:) = mean(Gout_x(1:2,:),1);
            else
                mGout_x(jj,:) = mean(Gout_x(Nreps*jj-Nreps:Nreps*jj-1,:),1);
            end
            % mean for Gy
            if jj<12
                mGout_y(jj,:) = mean(Gout_y(Nreps*jj-Nreps+1:Nreps*jj,:),1);
            end
        elseif testDummy == 'True'   
            mGout_x(jj,:) = mean(Gout_x(Nreps*jj-Nreps+1:Nreps*jj,:),1);   % mean for Gx (mT/m whitout Gamma)         
            mGout_y(jj,:) = mean(Gout_y(Nreps*jj-Nreps+1:Nreps*jj,:),1);   % mean for Gy (mT/m whitout Gamma)
        end            
    end
    
    % mean plots
    if plotTest6 == 'True' 
        timeV = [0:dt:dt*(size(mGout_x,2)-1)]*1e3;
        mGout_x_plot=mGout_x; %*1e-3*gammabar;
        mGout_y_plot=mGout_y; %*1e-3*gammabar;
        figure()
        subplot(211)
        plot(timeV,(-1)*mGout_x_plot','r')
        hold on
        title('Measured G_x')
        ylabel(['Gradient (mT/m)'],'fontsize',15)
        xlabel(['time (ms)'],'fontsize',15)        
        subplot(212)
        plot(timeV,(-1)*mGout_y_plot','b')
        hold on
        title('Measured G_y')
        ylabel(['Gradient (mT/m)'],'fontsize',15)        
        xlabel(['time (ms)'],'fontsize',15)        
    end
    clear Gout_x Gout_y
    Gout_x = mGout_x'*1e-3*gammabar;  % Gx_out (T/m * Gamma)
    Gout_y = mGout_y'*1e-3*gammabar;  % Gy_out (T/m * Gamma)
end

fprintf('\n\n 6 - Sucessfully finished - Measuring ktraj & Gout for xx and yy\n\n')


%% part 7 - Obtaining the impulse response function of the system 
% ... part 7.1 - Output ...
% Ox_w = zeros(size(Gout_x,2),Npulses);
% Oy_w = zeros(size(Gout_y,2),Npulses);
Ox_w = zeros(size(Gout_x));
Oy_w = zeros(size(Gout_y));

% get O(w)
for ii=1:Npulses
    Ox_w(:,ii) = fftshift(fft(Gout_x(:,ii))/gammabar);
end    
for ii=1:Npulses  % test27 put: 'Npulses-1'
    Oy_w(:,ii) = fftshift(fft(Gout_y(:,ii))/gammabar);
end

% plot O(w)
if plotTest7 == 'True'
    figure()    
    subplot(131)
    plot(abs(Ox_w));
    hold on
    title('O_w x-axis')
    subplot(132)
    plot(abs(Oy_w));
    hold on
    title('O_w y-axis')    
end


% ... part 7.2 - Input ...
if syncFun == 'True'
    I_w = input_sincF_freq.I_w;
    I_w = I_w(2:1+size(Ox_w,1),:);
else
    I_w = input_freq.i_tFFT;
    I_w = I_w(2:1+size(Ox_w,1),:);    
end
if plotTest7 == 'True'
    subplot(133)  
    plot(abs(I_w))
    title('I_w all-axis')
end

% ... part 7.3 - Transfer function (Impulse Response Function) ...
maxFreq = 1./T(1)*1e3;
ff      = linspace(-maxFreq,maxFreq,size(Gout_y,1)); % frequencies vector (Hz)

if MATLABfft == 'True' % Test H from FFT from matlab 
    Hx_w_factorA = zeros(size(I_w,1),1);
    Hx_w_factorB = zeros(size(I_w,1),1);
    Hy_w_factorA = zeros(size(I_w,1),1);
    Hy_w_factorB = zeros(size(I_w,1),1);
    
    for ii=1:Npulses
        Hx_w_factorA    = Hx_w_factorA  +  (conj(I_w(:,ii)).*Ox_w(:,ii)) ;
        Hx_w_factorB    = Hx_w_factorB  +  (abs(I_w(:,ii)).^2) ;
    end
    for ii=1:Npulses-1
        Hy_w_factorA    = Hy_w_factorA  +  (conj(I_w(:,ii)).*Oy_w(:,ii)) ;
        Hy_w_factorB    = Hy_w_factorB  +  (abs(I_w(:,ii)).^2) ;
    end
    
    Hx_w = Hx_w_factorA ./ Hx_w_factorB;
    Hy_w = Hy_w_factorA ./ Hy_w_factorB;

    %save H
    H.Hx_w   = abs(Hx_w)/max(abs(Hx_w));
    H.Hy_w   = abs(Hy_w)/max(abs(Hy_w));
    H.pulses = Npulses;
    
    if plotTest7 == 'True'        
        % Magnitude Plots
        figure()
        subplot(121)
        hold on
        title('Magnitude Hx_w (b) & Magnitude Hy_w (r)')
        plot(ff*1e-3,abs(Hx_w)/max(abs(Hx_w)),'b')
        hold on
        plot(ff*1e-3,abs(Hy_w)/max(abs(Hy_w)),'r')
        % xlim([-20 20])
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
                
        % Phase Plots
        subplot(122)
        hold on
        title('Phase Hx_w & Phase Hy_w')
        plot(ff*1e-3,(angle(Hx_w)),'b')
        hold on
        plot(ff*1e-3,(angle(Hy_w)),'r')    
        ylabel(['GIRF Phase (rad)'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
    end  
    
    
    clear Hx_w Hy_w
elseif MATLABfft == 'Fals'
    % Test H with FFT from: Abo Seda
    tGx     = tGx(1:size(Gout_x,1),:);
    tGy     = tGy(1:size(Gout_y,1),:);
    Npad    = 0; % Add points

    Hx_w_test = gradient_create_GIRF(tGx,ff,Gout_x,dt,Npad);
    Hy_w_test = gradient_create_GIRF(tGy,ff,Gout_y,dt,Npad);
    
    %save H
    H.Hx_w   = abs(Hx_w_test)/max(abs(Hx_w_test));
    H.Hy_w   = abs(Hy_w_test)/max(abs(Hy_w_test));
    H.pulses = Npulses;
    
    if plotTest7 == 'True'
        % Magnitude Plots        
        figure()
        subplot(121)
        hold on
        title('Magnitude Hx (b) & Magnitude Hy (r)')
        plot(ff*1e-3,abs(Hx_w_test)/max(abs(Hx_w_test)),'b')
        hold on
        plot(ff*1e-3,abs(Hy_w_test)/max(abs(Hy_w_test)),'r')
        % xlim([-20 20])
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        
        % Phase Plots
        subplot(122)
        hold on
        title('Phase Hx (b) & Phase Hy (r)')
        plot(ff*1e-3,(angle(Hx_w_test)),'b')
        hold on
        plot(ff*1e-3,(angle(Hy_w_test)),'r')  
        ylabel(['GIRF Phase (rad)'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
    end
    clear Hx_w_test Hy_w_test    
end



fprintf('\n\n 7 - Sucessfully finished -  Obtaining the IRF of the system \n\n')


%% part 8 - Save IRF
if saveT == 'True'
    cd(file_folder_results)
    save(['test',num2str(test),'_IRF.mat'],'H')
    cd([myCD, '/Results_IRF' ])
    save(['test',num2str(test),'_IRF.mat'],'H')
    cd(file_folder)    
    fprintf('\n\n 8 - Sucessfully finished - Save IRF\n\n')    
else
    fprintf('\n\n 8 - Sucessfully finished - Did NOT Save IRF\n\n')    
end

% close all
