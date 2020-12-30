%% 18-11-2020: Create Gradient transfer function (GIRF)
% usage: [H] = gradient_create_GIRF(Ginput,ff,Goutput,dt,Npad)
% INPUTS:
%      - Ginput can be MxN - i.e. can accept any number of inputs - time
%      domain
%      - ff = frequencies for which H is defined
%      - Ginput can be MxN - i.e. can accept any number of inputs - time
%      domain
%      - dt = gradient sample dur (default 6.4us if not set)
%      - Npad = number of padding samples to add at start and end (default=0)

function [H] = gradient_create_GIRF(Ginput,ff,Goutput,dt,Npad)

if (~exist('dt','var'))||isempty(dt)
    dt=6.4e-6;
end
if ~exist('Npad','var')
    Npad=0;
end    
Ngradients = size(Ginput,2);

%%% 1 - add some zeros on front and end
G    = cat(1,zeros([Npad Ngradients]),Ginput,zeros([Npad Ngradients]));
Gout = cat(1,zeros([Npad Ngradients]),Goutput,zeros([Npad Ngradients]));
M    = length(G);
t    = (0:M-1)*dt;


%%% 2 - check if frequencies sufficient
df_req = 1/t(end);
df_act = ff(2)-ff(1);
if df_req<df_act
    % New freq range
    ffnew = ff(1):df_req:ff(end);
    ffnew=ffnew-min(abs(ffnew));% centre so we have a zero sample
    % centring may have introduced one frequency out of range
    ffnew(1)=[];
    % resample H
    Hnew = interp1(ff,H,ffnew);
       
    %%% swap
    ff=ffnew;
    H=Hnew;
    if size(H,1)==1
        H=H(:);
    end
    fprintf(1,'GIRF: resampled because frequency resolution was insufficient\n');
    df_act = ff(2)-ff(1);
end
    

%%% 3 - construct Fourier matrix
F = dt*exp(2*pi*1i*ff'*t);
% Also define inverse
Fi = df_act*exp(-2*pi*1i*t'*ff);

%%% 4 - Apply transfer function
H_factor1 = zeros([M 1]);
H_factor2 = zeros([M 1]);

%%% 5 - pre-emph!
for ii=1:Ngradients    
    H_factor1 = H_factor1 + ( conj((F*G(:,ii))) .* (F*Gout(:,ii))   ) ;
    H_factor2 = H_factor2 + ( abs(((F*G(:,ii)))).^2  );    
    % %     H = H + (F*Gout(:,ii))  ./ (F*G(:,ii));
end

H = H_factor1 ./ H_factor2;


%%% 5 - Now make k-spaces before removing the padding
gamma_mT = 2*pi*42577.46778; % rad s^-1 mT^-1

%%% 6 - remove padding 2-3-15
H([1:Npad (M-Npad+1):M],:)=[];
M = length(H);
t = (0:M-1)*dt;


end


