%% Get H1 according to Robison, 2019
%  -  works for 'getBOECC.m' implementation
%  -  outputs the H1 (Impulse Response Function) for each axxis
%  see Robison, 2019 paper (10.1002/mrm.27583)
%
% TTFernandes, IST, Jan 2021

%%

function [H1_w] = getH1_BOECC(Gout_B0ECC,input_freq,Gin,Npulses,Nreps,Naxes,dt,T,plotTest7,MATLABfft)

gammabar = 42.576e6;  % Gyro-magnetic ratio (Hz/T)

if Naxes == 1
    % ... initialization ...
    maxFreq = 1./T(1)*1e3;
    ff      = linspace(-maxFreq,maxFreq,size(Gout_B0ECC.x,1)); % frequencies vector (Hz)
    
    if MATLABfft == 'True' % Test H from FFT from matlab
        % ... part 7.1 - Output ...
        O_w.x = zeros(size(Gout_B0ECC.x));
        for ii=1:Npulses*Nreps % get O(w)
            O_w.x(:,ii) = fftshift(fft(Gout_B0ECC.x(:,ii))/gammabar);
        end
        if plotTest7 == 'True' % plot O(w)
            figure()
            plot(abs(O_w.x));
            hold on
        end
        
        % ... part 7.2 - Input ...
        I_w = input_freq.i_tFFT;
        I_w = I_w(2:1+size(O_w,1),:);        
        if plotTest7 == 'True'
            plot(abs(I_w))
            title('I_w all-axis')
        end
        
        % ... part 7.3 - H1 - Transfer function (Impulse Response Function) ...
        Hx_w_factorA = zeros(size(I_w,1),1);
        H_w_factorB  = zeros(size(I_w,1),1);
        for ii=1:Npulses
            Hx_w_factorA   = Hx_w_factorA  +  (conj(I_w(:,ii)).*O_w.x(:,ii)) ;
            H_w_factorB    = H_w_factorB   +  (abs(I_w(:,ii)).^2) ;
        end        
        Hx_w = Hx_w_factorA ./ H_w_factorB;
        
        %save H1
        H1_w.x      = abs(Hx_w)/max(abs(Hx_w));
        H1_w.pulses = Npulses;
        clear Hx_w Hx_w_factorA H_w_factorB
        
    elseif MATLABfft == 'Fals' % Test H with FFT from: Abo Sema
        Gin       = Gin(1:size(Gout_B0ECC.x,1),:);
        Npad      = 0; % Add points
        Hx_w_test = gradient_create_GIRF(Gin,ff,Gout_B0ECC.x,dt,Npad);
        
        %save H1
        H1_w.x      = abs(Hx_w_test)/max(abs(Hx_w_test));
        H1_w.pulses = Npulses;
        H1_w.ff     = ff;
        clear Hx_w_test
    end
    
    % ... part 7.4 - Plots ...    
    if plotTest7 == 'True'
        figure() 
        subplot(121)
        hold on
        title('Magnitude H_1_x(b)')
        plot(ff*1e-3,abs(H1_w.x),'b')
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)               
        legend('H_0_x')
        
        % Phase plot
        subplot(122)
        hold on
        plot(ff*1e-3,(angle(H1_w.x)),'b')
        hold on
        ylabel(['GIRF Phase (rad)'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        legend('H_1_x')   
    end

    
elseif Naxes == 2
    % ... initialization ...
    maxFreq = 1./T(1)*1e3;
    ff      = linspace(-maxFreq,maxFreq,size(Gout_B0ECC.x,1)); % frequencies vector (Hz)
    
    if MATLABfft == 'True' % Test H from FFT from matlab
        % ... part 7.1 - Output ...
        O_w.x = zeros(size(Gout_B0ECC.x));
        O_w.y = zeros(size(Gout_B0ECC.y));
        for ii=1:Npulses*Nreps % get O(w)
            O_w.x(:,ii) = fftshift(fft(Gout_B0ECC.x(:,ii))/gammabar);
            O_w.y(:,ii) = fftshift(fft(Gout_B0ECC.y(:,ii))/gammabar);
        end
        if plotTest7 == 'True' % plot O(w)
            figure()
            plot(ff*1e-3, abs(O_w.x),'b');
            plot(ff*1e-3, abs(O_w.y),'r');
            hold on
%                 xlim([-20 20])
            title('Frequency Domain')
            xlabel(['frequency (kHz)'],'fontsize',15)
            ylabel(['Spectral Density (ms*mT/m)'],'fontsize',15)
        end
        
        % ... part 7.2 - Input ...
        I_w = input_freq.i_tFFT;
        I_w = I_w(2:1+size(O_w,1),:);       
        if plotTest7 == 'True'
            plot(abs(I_w))
            title('I_w all-axis')
        end
        
        % ... part 7.3 - H1 - Transfer function (Impulse Response Function) ...
        Hx_w_factorA = zeros(size(I_w,1),1);
        Hy_w_factorA = zeros(size(I_w,1),1);
        H_w_factorB  = zeros(size(I_w,1),1);
        for ii=1:Npulses
            Hx_w_factorA   = Hx_w_factorA  +  (conj(I_w(:,ii)).*O_w.x(:,ii)) ;
            Hy_w_factorA   = Hy_w_factorA  +  (conj(I_w(:,ii)).*O_w.y(:,ii)) ;
            H_w_factorB    = H_w_factorB   +  (abs(I_w(:,ii)).^2) ;
        end        
        Hx_w = Hx_w_factorA ./ H_w_factorB;
        Hy_w = Hy_w_factorA ./ H_w_factorB;
        
        %save H1
        H1_w.x      = abs(Hx_w)/max(abs(Hx_w));
        H1_w.y      = abs(Hy_w)/max(abs(Hy_w));
        H1_w.pulses = Npulses;
        clear Hx_w Hx_w_factorA Hy_w_factorA H_w_factorB
                
    elseif MATLABfft == 'Fals' % Test H with FFT from: Abo Sema
        Gin       = Gin(1:size(Gout_B0ECC.x,1),:);
        Npad      = 0; % Add points
        Hx_w_test = gradient_create_GIRF(Gin,ff,Gout_B0ECC.x,dt,Npad);
        Hy_w_test = gradient_create_GIRF(Gin,ff,Gout_B0ECC.y,dt,Npad);
        
        %save H1
        H1_w.x      = abs(Hx_w_test)/max(abs(Hx_w_test));
        H1_w.y      = abs(Hy_w_test)/max(abs(Hy_w_test));
        H1_w.pulses = Npulses;
        H1_w.ff     = ff;
        
        clear Hx_w_test Hy_w_test
    end
    
    % ... part 7.4 - Plots ...
    if plotTest7 == 'True'
        figure()
%         subplot(121)        
        hold on
        plot(ff*1e-3,abs(H1_w.x),'b')
        ylabel(['Gradient Field Response'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        hold on
        plot(ff*1e-3,abs(H1_w.y),'r')
        legend('H_1_x','H_1_y')
        
        % Phase Plots
% %         subplot(122)
% %         hold on
% %         plot(ff*1e-3,(angle(H1_w.x)),'b')
% %         hold on
% %         plot(ff*1e-3,(angle(H1_w.y)),'r')
% %         ylabel(['GIRF Phase (rad)'],'fontsize',15)
% %         xlabel(['frequency (kHz)'],'fontsize',15)
% %         legend('H_1_x','H_1_y')           
    end
    
    
elseif Naxes == 3
    % ... initialization ...
    maxFreq = 1./T(1)*1e3;
    ff      = linspace(-maxFreq,maxFreq,size(Gout_B0ECC.x,1)); % frequencies vector (Hz)
    
    if MATLABfft == 'True' % Test H from FFT from matlab
        % ... part 7.1 - Output ...
        O_w.x = zeros(size(Gout_B0ECC.x));
        O_w.y = zeros(size(Gout_B0ECC.y));
        O_w.z = zeros(size(Gout_B0ECC.z));
        for ii=1:Npulses*Nreps  % get O(w)
            O_w.x(:,ii) = fftshift(fft(Gout_B0ECC.x(:,ii))/gammabar);
            O_w.y(:,ii) = fftshift(fft(Gout_B0ECC.y(:,ii))/gammabar);
            O_w.z(:,ii) = fftshift(fft(Gout_B0ECC.z(:,ii))/gammabar);
        end
        if plotTest7 == 'True' % plot O(w)
            figure()
            plot(abs(O_w.x),'b');
            plot(abs(O_w.y),'r');
            plot(abs(O_w.z),'g');
            hold on
        end
        
        % ... part 7.2 - Input ...
        I_w = input_freq.i_tFFT;
        I_w = I_w(2:1+size(O_w,1),:);        
        if plotTest7 == 'True'
            plot(abs(I_w))
            title('I_w all-axis')
        end
        
        % ... part 7.3 - H1 - Transfer function (Impulse Response Function) ...
        Hx_w_factorA = zeros(size(I_w,1),1);
        Hy_w_factorA = zeros(size(I_w,1),1);
        Hz_w_factorA = zeros(size(I_w,1),1);
        H_w_factorB  = zeros(size(I_w,1),1);
        for ii=1:Npulses
            Hx_w_factorA   = Hx_w_factorA  +  (conj(I_w(:,ii)).*O_w.x(:,ii)) ;
            Hy_w_factorA   = Hy_w_factorA  +  (conj(I_w(:,ii)).*O_w.y(:,ii)) ;
            Hz_w_factorA   = Hz_w_factorA  +  (conj(I_w(:,ii)).*O_w.z(:,ii)) ;
            H_w_factorB    = H_w_factorB   +  (abs(I_w(:,ii)).^2) ;
        end        
        Hx_w = Hx_w_factorA ./ H_w_factorB;
        Hy_w = Hy_w_factorA ./ H_w_factorB;
        Hz_w = Hz_w_factorA ./ H_w_factorB;
        
        %save H1
        H1_w.x      = abs(Hx_w)/max(abs(Hx_w));
        H1_w.y      = abs(Hy_w)/max(abs(Hy_w));
        H1_w.z      = abs(Hy_w)/max(abs(Hy_w));
        H1_w.pulses = Npulses;
        clear Hx_w Hx_w_factorA Hy_w_factorA Hz_w_factorA H_w_factorB
        
    elseif MATLABfft == 'Fals' % Test H with FFT from: Abo Sema
        Gin       = Gin(1:size(Gout_B0ECC.x,1),:);
        Npad      = 0; % Add points
        Hx_w_test = gradient_create_GIRF(Gin,ff,Gout_B0ECC.x,dt,Npad);
        Hy_w_test = gradient_create_GIRF(Gin,ff,Gout_B0ECC.y,dt,Npad);
        Hz_w_test = gradient_create_GIRF(Gin,ff,Gout_B0ECC.z,dt,Npad);
        
        %save H1
        H1_w.x      = abs(Hx_w_test)/max(abs(Hx_w_test));
        H1_w.y      = abs(Hy_w_test)/max(abs(Hy_w_test));
        H1_w.z      = abs(Hz_w_test)/max(abs(Hz_w_test));
        H1_w.pulses = Npulses;
        H1_w.ff     = ff;        
        clear Hx_w_test Hy_w_test Hz_w_test
    end
    
    % ... part 7.4 - Plots ...
    if plotTest7 == 'True'
        figure()
        subplot(121)        
        hold on
        title('Magnitude H_1_x(b) & H_1_y(r) & H_1_z(g)')
        plot(ff*1e-3,abs(H1_w.x),'b')
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        hold on
        plot(ff*1e-3,abs(H1_w.y),'r')
        hold on
        plot(ff*1e-3,abs(H1_w.z),'g')
        legend('H_0_x','H_0_y','H_0_z')
        
        % Phase Plots
        subplot(122)
        hold on
        plot(ff*1e-3,(angle(H1_w.x)),'b')
        hold on
        plot(ff*1e-3,(angle(H1_w.y)),'r')
        hold on
        plot(ff*1e-3,(angle(H1_w.z)),'g')
        ylabel(['GIRF Phase (rad)'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        legend('H_0_x','H_0_y','H_0_z')
    end       
end

end
