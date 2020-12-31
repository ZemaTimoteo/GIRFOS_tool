%% Obtain external gradients in matlab matrix from pypulseq sequence
% by TTFernandes, 2020 IST
clc

%% --- 1 - Load and initialize vectors
% cd('/home/tfernandes/Documents/Projetos/Project_lfStroke/Tests/1_Spiral_test/Test/test79_Spiraltest3_Nreps-1_Nslices-1_TR-2s_TE-56ms_sm-125_RFbw-1_RFoff-0mm_Diff_Ndirs-2_bvalue-500/cfl_hdr_files')
cd(file_folder_cfl_hdr)
seq_gradX = load('GradX.mat');
seq_gradY = load('GradY.mat');
  
tRecon_Plot = 'Fals';   % See plot - 'True'

i_raster_time = round(1/dt);
ntrials       = (tmp.ndirs + 1)*tmp.nreps;
seq_npoints   = round((tmp.TR*ntrials) / dt);

%% --- 2 - fill vectors for XX
seq_Vect_timeX = zeros(1,seq_npoints);
seq_Vect_gradX = zeros(1,seq_npoints);
nspiX = 1; % for flag of Spiral

seqX_nPulses = size(seq_gradX.time,2);
for i=1:seqX_nPulses
    % vector Time
    seqTime = seq_gradX.time{1,i};
    if size(seqTime,2) == 5
        tInitial = round(seqTime(1+1)*i_raster_time);
        tFinal   = round(seqTime(end)*i_raster_time);
        aux_tx   = [tInitial:1:tFinal];
    else
        tInitial = round(seqTime(1+1)*i_raster_time);
        tFinal   = round(seqTime(end-1)*i_raster_time);
        aux_tx   = tInitial:1:tFinal;                    % have to remove the 1st and last point to compensate effect of pypulseq
    end
    seq_Vect_timeX(1,aux_tx)=1;   
      
    % vector Grad
    aux_GradX = seq_gradX.Grad{1,i};
    if size(seqTime,2) == 5
        %up step
        upGrad_aux_tx = round(seqTime(1+1)*i_raster_time):1:round(seqTime(2+1)*i_raster_time);
        tSlope        = size(upGrad_aux_tx,2)*dt;
        m_up          = (aux_GradX(1,3)-aux_GradX(1,2))/tSlope;
        up_GradX      = [0*dt:dt:tSlope-dt]*m_up;
        seq_Vect_gradX(1,upGrad_aux_tx) = up_GradX;
        
        % platou
        platou_aux_tx   = round(seqTime(2+1)*i_raster_time):1:round(seqTime(3+1)*i_raster_time-1);
        seq_Vect_gradX(1,platou_aux_tx) = aux_GradX(1,3);
        
        % down step
        downGrad_aux_tx = round(seqTime(3+1)*i_raster_time):1:round(seqTime(end)*i_raster_time);
        tSlope          = size(downGrad_aux_tx,2)*dt;
        m_down          = (aux_GradX(1,end)-aux_GradX(1,3+1))/tSlope;
        down_GradX      = aux_GradX(1,3)+[0*dt:dt:tSlope-dt]*m_down;
        seq_Vect_gradX(1,downGrad_aux_tx) = down_GradX;        
        
    else
        seq_Vect_gradX(1,aux_tx)=aux_GradX(1,1+1:end-1);    % have to remove the 1st and last point to compensate effect of pypulseq
        if size(seqTime,2) > 2000
            flagSPI_x(1,nspiX) = aux_tx(1,1)-1;
            nspiX = nspiX+1;
        end
            
    end
    clear seqTime aux_tx aux_GradX upGrad_aux_tx downGrad_aux_tx up_GradX down_GradX m_up m_down tSlope
end

%% --- 3 - fill vectors for YY
seq_Vect_timeY = zeros(1,seq_npoints);
seq_Vect_gradY = zeros(1,seq_npoints);
nspiY = 1; % for flag of Spiral

seqY_nPulses = size(seq_gradY.time,2);
for i=1:seqY_nPulses
    % vector Time
    seqTime = seq_gradY.time{1,i};
    if size(seqTime,2) == 5
        tInitial = round(seqTime(1+1)*i_raster_time);
        tFinal   = round(seqTime(end)*i_raster_time);
        aux_ty   = [tInitial:1:tFinal];
    else
        tInitial = round(seqTime(1+1)*i_raster_time);
        tFinal   = round(seqTime(end-1)*i_raster_time);
        aux_ty   = tInitial:1:tFinal;                    % have to remove the 1st and last point to compensate effect of pypulseq
    end
    seq_Vect_timeY(1,aux_ty)=1;   
      
    % vector Grad
    aux_GradY = seq_gradY.Grad{1,i};
    if size(seqTime,2) == 5
        %up step
        upGrad_aux_ty = round(seqTime(1+1)*i_raster_time):1:round(seqTime(2+1)*i_raster_time);
        tSlope        = size(upGrad_aux_ty,2)*dt;
        m_up          = (aux_GradY(1,3)-aux_GradY(1,2))/tSlope;
        up_GradY      = [0*dt:dt:tSlope-dt]*m_up;
        seq_Vect_gradY(1,upGrad_aux_ty) = up_GradY;
        
        % platou
        platou_aux_ty   = round(seqTime(2+1)*i_raster_time):1:round(seqTime(3+1)*i_raster_time-1);
        seq_Vect_gradY(1,platou_aux_ty) = aux_GradY(1,3);
        
        % down step
        downGrad_aux_ty = round(seqTime(3+1)*i_raster_time):1:round(seqTime(end)*i_raster_time);
        tSlope          = size(downGrad_aux_ty,2)*dt;
        m_down          = (aux_GradY(1,end)-aux_GradY(1,3+1))/tSlope;
        down_GradY      = aux_GradY(1,3)+[0*dt:dt:tSlope-dt]*m_down;
        seq_Vect_gradY(1,downGrad_aux_ty) = down_GradY;        
        
    else
        seq_Vect_gradY(1,aux_ty)=aux_GradY(1,1+1:end-1);    % have to remove the 1st and last point to compensate effect of pypulseq
        if size(seqTime,2) > 2000
            flagSPI_y(1,nspiY) = aux_ty(1,1)-1;
            nspiY = nspiY+1;
        end
    end
    clear seqTime aux_tx aux_GradX upGrad_aux_tx downGrad_aux_tx up_GradX down_GradX m_up m_down tSlope
end


%% --- 4 - Plots

if tRecon_Plot == 'True'
    vectTime = dt*0:dt:(size(seq_Vect_timeY,2)-1)*dt;
    
    figure()
    subplot(121)
    plot(vectTime,seq_Vect_gradX)
    subplot(122)
    plot(vectTime,seq_Vect_gradY)
end
%% --- 5 - Organize data
gradForH_x = zeros(ntrials,nP_before+sizeSpiral);
gradForH_y = zeros(ntrials,nP_before+sizeSpiral);


for j=1:ntrials
    gradForH_x(j,:) = seq_Vect_gradX(1,flagSPI_x(j)-nP_before:flagSPI_x(j)+sizeSpiral-1);
    gradForH_y(j,:) = seq_Vect_gradY(1,flagSPI_y(j)-nP_before:flagSPI_y(j)+sizeSpiral-1);
end

if tRecon_Plot == 'True'   
    figure()
    subplot(1,3,1)
    plot(gradForH_x(1,:))
    subplot(1,3,2)
    plot(gradForH_x(2,:))
    subplot(1,3,3)
    plot(gradForH_x(3,:))
end
