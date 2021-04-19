%% Get H0 according to Robison, 2019
%  -  works for 'getBOECC.m' implementation
%  -  outputs the H0 (impact of B0 phase) for each axxis
%  see Robison, 2019 paper (10.1002/mrm.27583)
% 
% TTFernandes, IST, Jan 2021

%%
function [H0_w] = getH0_BOECC(phase_B0ECC,input_freq,Gin,Npulses,Naxes,dt,T,plotTest7,MATLABfft,filter)

gammabar = 42.576e6;  % Gyro-magnetic ratio (Hz/T)

if filter == 'True'
    % filter parameters
    order = phase_B0ECC.order;   % order
    wsize = phase_B0ECC.wsize;  % window size
end

if Naxes == 1
    % ... initialization ...
    maxFreq = 1./T(1)*1e3;
    ff      = linspace(-maxFreq,maxFreq,size(phase_B0ECC.x,1)); % frequencies vector (Hz)
    
    if MATLABfft == 'True' % Test H from FFT from matlab 
        % ... part 7.1 - Output ...
        O_w.x = zeros(size(phase_B0ECC.x));
        for ii=1:Npulses*Nreps % get O(w)
            O_w.x(:,ii) = fftshift(fft(phase_B0ECC.x(:,ii))/gammabar);
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
        % Test H with FFT from: Abo Sema
        Gin       = Gin(1:size(phase_B0ECC.x,1),:);
        Npad      = 0; % Add points
        Hx_w_test = gradient_create_GIRF(Gin,ff,phase_B0ECC.x,dt,Npad);
        
        % ... filter with Savitzky-Golay filter (3rd order, window size:51) ...
        if filter == 'True'
            Hx_w_test = sgolayfilt(Hx_w_test,order,wsize);
            H0_w.FilterOrder = order;
            H0_w.FilterWsize = wsize;
        end
        
        % ... save H0 ...
        H0_w.x      = abs(Hx_w_test)/max(abs(Hx_w_test));
        H0_w.pulses = Npulses;
        H0_w.ff     = ff;
        clear Hx_w_test
    end
    
    % ... Plots ...
    if plotTest7 == 'True'
        figure()
        hold on
        plot(ff*1e-3,abs(H0_w.x),'b')
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)               
        legend('H_0_x')
    end

    
elseif Naxes == 2
    % ... initialization ...
    maxFreq = 1./T(1)*1e3;
    ff      = linspace(-maxFreq,maxFreq,size(phase_B0ECC.x,1)); % frequencies vector (Hz)
    
    % ... Test H with FFT from: Abo Sema ...
    Gin  = Gin(1:size(phase_B0ECC.x,1),:);
%     Gin  = unwrap(angle(Gin));
    
    Npad      = 0; % Add points
    
    Hx_w_test = gradient_create_GIRF(Gin,ff,phase_B0ECC.x,dt,Npad);
    Hy_w_test = gradient_create_GIRF(Gin,ff,phase_B0ECC.y,dt,Npad);
    
    % ... filter with Savitzky-Golay filter (3rd order, window size:51) ...
    if filter == 'True'
        Hx_w_test = sgolayfilt(Hx_w_test,order,wsize);
        Hy_w_test = sgolayfilt(Hy_w_test,order,wsize);
        H0_w.FilterOrder = order;
        H0_w.FilterWsize = wsize;        
    end
    
    % ... save H0 ...
    H0_w.x      = Hx_w_test;
    H0_w.y      = Hy_w_test;
    H0_w.pulses = Npulses;
    H0_w.ff     = ff;

%     clear Hx_w_test Hy_w_test Hx_test_filter Hy_test_filter
    
    % ... Plots ...
    if plotTest7 == 'True'
        figure()
        hold on
        plot(ff*1e-3,abs(H0_w.x),'b')
        ylabel(['B_0 EC Response'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        hold on
        if filter == 'True'
            title(['H0_x & H0_y (filter - order:',num2str(order),' Wsize:',num2str(wsize),')']);
        else
            title(['H0_x & H0_y (no filter)']);            
        end
        plot(ff*1e-3,abs(H0_w.y),'r')
        legend('H_0_x','H_0_y')
    end
    
    
elseif Naxes == 3
    % ... initialization ...
    maxFreq = 1./T(1)*1e3;
    ff      = linspace(-maxFreq,maxFreq,size(phase_B0ECC.x,1)); % frequencies vector (Hz)
    
    % ... Test H with FFT from: Abo Sema ...
    Gin       = Gin(1:size(phase_B0ECC.x,1),:);
    Npad      = 0; % Add points
    Hx_w_test = gradient_create_GIRF(Gin,ff,phase_B0ECC.x,dt,Npad);
    Hy_w_test = gradient_create_GIRF(Gin,ff,phase_B0ECC.y,dt,Npad);
    Hz_w_test = gradient_create_GIRF(Gin,ff,phase_B0ECC.z,dt,Npad);

    % ... filter with Savitzky-Golay filter (3rd order, window size:51) ...
    if filter == 'True'
        Hx_w_test = sgolayfilt(Hx_w_test,order,wsize);
        Hy_w_test = sgolayfilt(Hy_w_test,order,wsize);
        Hz_w_test = sgolayfilt(Hz_w_test,order,wsize);
        H0_w.FilterOrder = order;
        H0_w.FilterWsize = wsize;        
    end    
    
    % ... save H0 ...
    H0_w.x      = abs(Hx_w_test)/max(abs(Hx_w_test));
    H0_w.y      = abs(Hy_w_test)/max(abs(Hy_w_test));
    H0_w.z      = abs(Hz_w_test)/max(abs(Hz_w_test));
    H0_w.pulses = Npulses;
    H0_w.ff     = ff;
    clear Hx_w_test Hy_w_test Hz_w_test
    
    % ... Plots ...
    if plotTest7 == 'True'
        figure()
        hold on
        plot(ff*1e-3,abs(H0_w.x),'b')
        ylabel(['GIRF magnitude'],'fontsize',15)
        xlabel(['frequency (kHz)'],'fontsize',15)
        hold on
        plot(ff*1e-3,abs(H0_w.y),'r')
        hold on
        plot(ff*1e-3,abs(H0_w.z),'g')
        if filter == 'True'
            title(['H0_x, H0_y & H0_z (filter - order:',num2str(order),' Wsize:',num2str(wsize),')']);
        else
            title(['H0_x, H0_y & H0_z (no filter)']);            
        end
        legend('H_0_x','H_0_y','H_0_z')
    end
        
end

end
