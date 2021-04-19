%% Calculations according to Robison, 2019
%  -  works for 'getBOECC.m' implementation
%  -  outputs the kdata compressed for each of the particular situation
%  see Robison, 2019 paper (10.1002/mrm.27583)
% 
% TTFernandes, IST, Jan 2021

%%

function [angle_kcc,A,B,phase_B0ECC,ktraj_out_B0ECC, Gout_B0ECC] = auxCalc_BOECC(k_cc,D,Npulses,Nreps,Naxes,dt)

gammabar = 42.576e6;  % Gyro-magnetic ratio (Hz/T)

if Naxes == 1    
    for  nr = 1:Nreps*Npulses        
        angle_kcc.I_x(:,nr)   = unwrap(angle(k_cc.kx_Pgrad_Px0(:,nr)));
        angle_kcc.II_x(:,nr)  = unwrap(angle(k_cc.kx_Ngrad_Px0(:,nr)));
        angle_kcc.III_x(:,nr) = unwrap(angle(k_cc.kx_Pgrad_Nx0(:,nr)));
        angle_kcc.IV_x(:,nr)  = unwrap(angle(k_cc.kx_Ngrad_Nx0(:,nr)));  
        
        A.x(:,nr) = angle_kcc.I_x(:,nr)  - angle_kcc.II_x(:,nr);
        B.x(:,nr) = angle_kcc.IV_x(:,nr) - angle_kcc.III_x(:,nr);
        
        % Phase from B0 eddy currents
        phase_B0ECC.x(:,nr)     = ( A.x(:,nr) - B.x(:,nr) ) / 4;        
        % Actual k-space trajectory in the presence of eddy currents
        ktraj_out_B0ECC.x(:,nr) = ( A.x(:,nr) + B.x(:,nr) ) / (4*D)/(2*pi);        
        % Actual Gradient in the presence of eddy currents
        Gout_B0ECC.x(:,nr)      = (ktraj_out_B0ECC.x(2:end,nr)-ktraj_out_B0ECC.x(1:end-1,nr))/dt/gammabar*1e3; % Units (mT/m       
    end
    
    
elseif Naxes == 2
    for  nr = 1:Nreps*Npulses        
        angle_kcc.I_x(:,nr)   = unwrap(angle(k_cc.kx_Pgrad_Px0(:,nr)));
        angle_kcc.II_x(:,nr)  = unwrap(angle(k_cc.kx_Ngrad_Px0(:,nr)));
        angle_kcc.III_x(:,nr) = unwrap(angle(k_cc.kx_Pgrad_Nx0(:,nr)));
        angle_kcc.IV_x(:,nr)  = unwrap(angle(k_cc.kx_Ngrad_Nx0(:,nr)));     

        angle_kcc.I_y(:,nr)   = unwrap(angle(k_cc.ky_Pgrad_Px0(:,nr)));
        angle_kcc.II_y(:,nr)  = unwrap(angle(k_cc.ky_Ngrad_Px0(:,nr)));
        angle_kcc.III_y(:,nr) = unwrap(angle(k_cc.ky_Pgrad_Nx0(:,nr)));
        angle_kcc.IV_y(:,nr)  = unwrap(angle(k_cc.ky_Ngrad_Nx0(:,nr)));        
        
        
        % plots test
        figure()
        subplot(221)
        plot(angle_kcc.I_x(:,nr))
        title('I_x')
        subplot(222)
        plot(angle_kcc.II_x(:,nr))
        title('II_x')
        subplot(223)
        plot(angle_kcc.III_x(:,nr))
        title('III_x')
        subplot(224)
        plot(angle_kcc.IV_x(:,nr))
        title('IV_x')                        
        figure()
        subplot(221)
        plot(angle_kcc.I_y(:,nr))
        title('I_y')
        subplot(222)
        plot(angle_kcc.II_y(:,nr))
        title('II_y')
        subplot(223)
        plot(angle_kcc.III_y(:,nr))
        title('III_y')
        subplot(224)
        plot(angle_kcc.IV_y(:,nr))
        title('IV_y')
        
        
        
        A.x(:,nr) = angle_kcc.I_x(:,nr)  - angle_kcc.II_x(:,nr);
        B.x(:,nr) = angle_kcc.IV_x(:,nr) - angle_kcc.III_x(:,nr);
        A.y(:,nr) = angle_kcc.I_y(:,nr)  - angle_kcc.II_y(:,nr);
        B.y(:,nr) = angle_kcc.IV_y(:,nr) - angle_kcc.III_y(:,nr);        

        % plot test
        figure()
        subplot(221)
        plot(A.x(:,nr))
        title('A_x (I-II)')        
        subplot(222)
        plot(A.y(:,nr))
        title('A_y (I-II)')
        subplot(223)
        plot(B.x(:,nr))
        title('B_x (IV-III)')        
        subplot(224)
        plot(B.y(:,nr))
        title('B_y (IV-III)')
        
        % Phase from B0 eddy currents
        phase_B0ECC.x(:,nr)     = ( A.x(:,nr) - B.x(:,nr) ) / 4;        
        phase_B0ECC.y(:,nr)     = ( A.y(:,nr) - B.y(:,nr) ) / 4;   
        
        % plot tests
        figure()
        subplot(121)
        plot(phase_B0ECC.x(:,nr))
        title('phase - B0ECC_x')
        subplot(122)
        plot(phase_B0ECC.y(:,nr))
        title('phase - B0ECC_y')

        
        
        % Actual k-space trajectory in the presence of eddy currents
        ktraj_out_B0ECC.x(:,nr) = ( A.x(:,nr) + B.x(:,nr) ) / (4*D) / (2*pi);
        ktraj_out_B0ECC.y(:,nr) = ( A.y(:,nr) + B.y(:,nr) ) / (4*D) / (2*pi);
        % Actual Gradient in the presence of eddy currents
        Gout_B0ECC.x(:,nr)      = (ktraj_out_B0ECC.x(2:end,nr)-ktraj_out_B0ECC.x(1:end-1,nr))/dt/gammabar*1e3;
        Gout_B0ECC.y(:,nr)      = (ktraj_out_B0ECC.y(2:end,nr)-ktraj_out_B0ECC.y(1:end-1,nr))/dt/gammabar*1e3;
    end     
    
%     figure()
%     plot(Gout_B0ECC.x)
%     title('Gout_x B0ECC')

elseif Naxes == 3
    for  nr = 1:Nreps*Npulses
        angle_kcc.I_x(:,nr)   = unwrap(angle(k_cc.kx_Pgrad_Px0(:,nr)));
        angle_kcc.II_x(:,nr)  = unwrap(angle(k_cc.kx_Ngrad_Px0(:,nr)));
        angle_kcc.III_x(:,nr) = unwrap(angle(k_cc.kx_Pgrad_Nx0(:,nr)));
        angle_kcc.IV_x(:,nr)  = unwrap(angle(k_cc.kx_Ngrad_Nx0(:,nr)));
        
        angle_kcc.I_y(:,nr)   = unwrap(angle(k_cc.ky_Pgrad_Px0(:,nr)));
        angle_kcc.II_y(:,nr)  = unwrap(angle(k_cc.ky_Ngrad_Px0(:,nr)));
        angle_kcc.III_y(:,nr) = unwrap(angle(k_cc.ky_Pgrad_Nx0(:,nr)));
        angle_kcc.IV_y(:,nr)  = unwrap(angle(k_cc.ky_Ngrad_Nx0(:,nr)));
        
        angle_kcc.I_z(:,nr)   = unwrap(angle(k_cc.kz_Pgrad_Px0(:,nr)));
        angle_kcc.II_z(:,nr)  = unwrap(angle(k_cc.kz_Ngrad_Px0(:,nr)));
        angle_kcc.III_z(:,nr) = unwrap(angle(k_cc.kz_Pgrad_Nx0(:,nr)));
        angle_kcc.IV_z(:,nr)  = unwrap(angle(k_cc.kz_Ngrad_Nx0(:,nr)));
        
        A.x(:,nr) = angle_kcc.I_x(:,nr)  - angle_kcc.II_x(:,nr);
        B.x(:,nr) = angle_kcc.IV_x(:,nr) - angle_kcc.III_x(:,nr);
        A.y(:,nr) = angle_kcc.I_y(:,nr)  - angle_kcc.II_y(:,nr);
        B.y(:,nr) = angle_kcc.IV_y(:,nr) - angle_kcc.III_y(:,nr);
        A.z(:,nr) = angle_kcc.I_z(:,nr)  - angle_kcc.II_z(:,nr);
        B.z(:,nr) = angle_kcc.IV_z(:,nr) - angle_kcc.III_z(:,nr);
        
        % Phase from B0 eddy currents
        phase_B0ECC.x(:,nr)     = ( A.x(:,nr) - B.x(:,nr) ) / 4;        
        phase_B0ECC.y(:,nr)     = ( A.y(:,nr) - B.y(:,nr) ) / 4;                
        phase_B0ECC.z(:,nr)     = ( A.z(:,nr) - B.z(:,nr) ) / 4;                
        % Actual k-space trajectory in the presence of eddy currents
        ktraj_out_B0ECC.x(:,nr) = ( A.x(:,nr) + B.x(:,nr) ) / (4*D)/(2*pi);
        ktraj_out_B0ECC.y(:,nr) = ( A.y(:,nr) + B.y(:,nr) ) / (4*D)/(2*pi);        
        ktraj_out_B0ECC.z(:,nr) = ( A.z(:,nr) + B.z(:,nr) ) / (4*D)/(2*pi);        
        % Actual Gradient in the presence of eddy currents
        Gout_B0ECC.x(:,nr)      = (ktraj_out_B0ECC.x(2:end,nr)-ktraj_out_B0ECC.x(1:end-1,nr))/dt/gammabar*1e3;    
        Gout_B0ECC.y(:,nr)      = (ktraj_out_B0ECC.y(2:end,nr)-ktraj_out_B0ECC.y(1:end-1,nr))/dt/gammabar*1e3;
        Gout_B0ECC.z(:,nr)      = (ktraj_out_B0ECC.z(2:end,nr)-ktraj_out_B0ECC.z(1:end-1,nr))/dt/gammabar*1e3;
    end    
end

    
end