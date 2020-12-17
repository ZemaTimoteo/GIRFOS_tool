%% TrianglePulse function for GIRF
%  by TTFernandes, July 2020

%% Initialize
clc
% close all
clear all

FFT      = 1;  % '1' - True, '0' - False
computer = 2;  % 'Seia' - 2 OU 'Casa' - 1
cdTest = '/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF';
cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/1_GIRF')
%% Dirac delta Function

test  = 'SCAN';
saveT = 'Fals';

lenT   = 12;                        % Length of T
T      = linspace(0.05,0.16,lenT);  % time range - (ms)
p      = 115;                       % slope - mT/(m ms)
DT     = 0.01;                      % raster time (ms)
tTotal = 10;                         % Total time (ms)
tpre   = tTotal/2;                         % time Pre Pulse (ms)

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

figure()
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
    
    plot(tConvert,i_t(:,ii))
    hold on
end
ylabel(['Gradient (Hz/m)'],'fontsize',15)
xlabel(['time (s)'],'fontsize',15)
title('GIRF pulses')
xlim([-0.2 0.2]*1e-3)



%% Save GIRF pulse
if saveT == 'True'
    if computer == 1
        cd('./pulses')
    else
        cd('./pulses')
    end
    if test == 'PULS'
        save(['girf_',test,'_pulse_', num2str(lenT),'_tTotal_',num2str(tTotal),'ms_DT_',num2str(DT*1e2),'_1e-5_p',num2str(p),'.mat'],'i_t')
        cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/GIRF/pulses')
        save(['girf_',test,'_pulse_', num2str(lenT),'_tTotal_',num2str(tTotal),'ms_DT_',num2str(DT*1e2),'_1e-5_p',num2str(p),'.mat'],'i_t')
        cd(cdTest)
    elseif test == 'SCAN'
        save(['girf_',test,'_pulse_', num2str(lenT),'_tTotal_',num2str(tTotal),'ms_DT_',num2str(DT*1e2),'_1e-5_p',num2str(p),'.mat'],'i_t')
        cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/1_pythonCodes/GIRF/pulses')
        save(['girf_',test,'_pulse_', num2str(lenT),'_tTotal_',num2str(tTotal),'ms_DT_',num2str(DT*1e2),'_1e-5_p',num2str(p),'.mat'],'i_t')
    end
    if computer == 1
        cd('D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\5_GIRF\1_matlab_code')
    else
        cd(cdTest)
    end
end


%% FFT of the signal & Frequency domain
if FFT == 1
    
    % part 1 - FFT
    figure()
    freq = 1./tConvert*1e-3;
    for ii=1:length(T)
        i_tFFT(:,ii) = fftshift(fft(i_t(:,ii))/gammabar);
%         plot(freq,abs(i_tFFT(:,ii)));
        plot(abs(i_tFFT(:,ii)));
        hold on
    end
%     xlim([-20 20])
    title('Frequency Domain')
    xlabel(['frequency (kHz)'],'fontsize',15)
    ylabel(['Spectral Density (ms*mT/m)'],'fontsize',15)
    
    % part 2 - Frequency domain
    clear I_w
    w = (1./t) * (2*pi);
    freq = 1./tConvert*1e-3;   
    testI_w = 0;
    figure()
    for ii=1:length(T)
        I_w(:,ii) = (pConvert * ( ( sin( w * T(ii)./2 ) .^2 ) ./ (w.^2) ));
        testI_w = testI_w + sqrt(abs(I_w(:,ii)).^2);
        plot(freq,I_w(:,ii))
%         plot(I_w(:,ii))        
        hold on        
    end
    
    xlim([-20 20])
    xlabel(['frequency (kHz)'],'fontsize',15)
    ylabel(['Spectral Density (s*T/m)'],'fontsize',15)     
    maxAmpli_W = p*T(end)^2    
    
    figure()
    test1 = pConvert./(w.^2);
    plot(freq,test1)
    hold on
    xlim([-20 20])
    xlabel(['frequency (kHz)'],'fontsize',15)
    ylabel(['Spectral Density (s*T/m)'],'fontsize',15)    
    
    figure()
    plot(freq,testI_w)
    hold on
    xlim([-20 20])
    xlabel(['frequency (kHz)'],'fontsize',15)
    ylabel(['Spectral Density (s*T/m)'],'fontsize',15)    
    
    
    
 % part 3 - Save
    if saveT == 'True'
        if computer == 1
            cd('./pulses')
            save(['girf_Input_freq_', num2str(lenT),'.mat'],'i_tFFT')
            cd('D:\Tiago\Trabalho\2020_IST\lfStroke\Spiral_Optimization\5_GIRF\1_matlab_code')
        else
            if test == 'PULS'
                cd('./pulses')
                save(['girf_',test,'_Input_freq_', num2str(lenT),'.mat'],'i_tFFT')
                save(['girf_',test,'_Input_freq_sincFunct_', num2str(lenT),'.mat'],'I_w')
                cd(cdTest)
            elseif test == 'SCAN'
                cd('./pulses')
                save(['girf_',test,'_Input_freq_', num2str(lenT),'_tTotal_',num2str(tTotal),'ms.mat'],'i_tFFT')
                save(['girf_',test,'_Input_freq_sincFunct_', num2str(lenT),'.mat'],'I_w')
                cd(cdTest)
            end
        end
    end
end

close all

