%% 
%
% Generates a spiral k-space trajectory with a method adapted from [1]
% (the paper is a bit buggy).
%
% INPUTS:   * FOV    - field of view (meters)
%           * N      - resolution (Nyquist distance)
%           * Nshots - number of interleaves
%           * alpha  - variable density factor (alpha in [1])
%           * Gmax   - maximum amplitude (in T/m)
%           * SRmax  - maximum slew rate (in T/m/s)
%
% OUTPUTS:  * k      - 
%           * tge    - 
%           * lambda - 
%           * Ts2a   - Time slew regime to amplitude regime
%           * Tend   - Total Time of one shot gradient
%           * Tea    - End time of amplitude regime
%           * Tes    - End time of slew regime
%
% %       - generates da gradient (g1 and g2)
%       - k space design (k1 and k2) of the two spiral regims. 
%           - In case there is a single regim k2 and g2 are 0. 
%           - It also generates the time it takes for each regim.
%
% [1]   "Simple Analytic Variable Density Spiral Design"
%       Dong-hyun Kim, Elfar Adalsteinsson, and Daniel M. Spielman
%       Magnetic Resonance in Medicine 50:214-219 (2003)
%
% Copyright, Matthieu Guerquin-Kern, 2012
% Modified by : Pavan Poojar, MIRC
%               Tiago Fernandes, LASEEB, 2019
 
function [k,g,tDT,lambda,Ts2a,Tend,Tea,Tes]= vds2D_vf(FOV,N,Nshots,alpha,Gmax,SRmax,DT,gamma,slew_safety)

res=FOV/N;

%% Generating first interleave
% gamma = 42.576e6; % in Hz
% gamma = 1;
% gamma = 42.576; 

SRmax = SRmax * slew_safety;

lambda = .5/res; % in m^(-1)
n = (1/(1-(1-Nshots/FOV/lambda)^(1/alpha)));  % number of turns
w = 2*pi*n;
Tea = lambda*w/gamma/Gmax/(alpha+1); % in s
%Tea1 = (lambda*w)/((gamma.*gm)*(alpha+1));
Tes = sqrt(lambda*w^2/SRmax/gamma)/(alpha/2+1); % in s
Ts2a = (Tes^((alpha+1)/(alpha/2+1))*(alpha/2+1)/Tea/(alpha+1))^(1+2/alpha); % in s

if Ts2a<Tes
    tautrans = (Ts2a/Tes).^(1/(alpha/2+1));
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Ts2a)+((t-Ts2a)/Tea + tautrans^(alpha+1)).^(1/(alpha+1)).*(t>Ts2a).*(t<=Tea).*(Tes>=Ts2a);
    Tend = Tea;
else
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Tes);
    Tend = Tes;
end
 
k = @(t) lambda*tau(t).^alpha.*exp(1i*w*tau(t));
dt = Tea*1e-4; % in s
Dt = dt/FOV/abs(k(Tea)-k(Tea-dt)); % in s
% % t = 0:Dt:Tend; % in s
t = 0:DT:Tend; % in s
kt = k(t); % in m^-1

%% Interpolation for the selected Delta Time (DT) (in seconds)

% DT_GE=4e-6;
% % % tDT = 0:DT:Tend;
% % %  
% % %     dist = [0,cumsum(abs(kt(2:end)-kt(1:end-1)))];
% % %     kt_x = interp1(t,real(kt),tDT,'spline'); 
% % %     kt_y = interp1(t,imag(kt),tDT,'spline'); 
% % %     kt_new = complex(kt_x,kt_y);
% % %  
% % % k = kt_new;
% % % 

k   = kt;
tDT = t;

%% Obtain Gradient Waveforms

ktn(1,:) = real(k);
ktn(2,:) = imag(k);

% % [Gx, Gy, ~] = k2g_v2(ktn,DT,gamma);
[Gx, Gy, ~] = AuxSpiralTF_k2g(ktn,DT);

g = Gx + 1*i*Gy;

% figure(103); plot(Gx/gamma); hold on; plot(Gy/gamma);
end
                
 
 
 
 
 


