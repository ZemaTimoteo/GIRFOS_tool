%% ========
% This script allows to create a diffusion images based on S=S0exp(-bD),
%    of each tissues adding them together to create an image.
%% ========

%% Part 1 - Data generation for simulation test
if testDiffusion == 'test1'
% % % %     % 0 - Initialization ---
% % % %     name_simu_template = 'simu.xml';
% % % %     
% % % %     tissue = {'CF','GM','WM','LS'}; % Select type tissue 'CF' - CSF; 'GM' - Gray Matter; 'WM' - White Matter; 'LS' - Lesion
% % % %     origFOV    = 180;
% % % %     nOfTissues = size(tissue,2);
% % % %     
% % % %     % 1 - Generate Phatoms and directories
% % % %     clear A
% % % %     
% % % %     for ppp=1:nOfTissues
% % % %         typeTissue = tissue{ppp};
% % % %         cd(dir_seq)
% % % %         dir_tissue{ppp} = [dir_seq,'/',typeTissue];
% % % %         mkdir(dir_tissue{ppp})
% % % % 
% % % %         % --- generate mask per tissue ---
% % % % % %         Sample_Phantom_generator_BRAIN_byRegion_3T % for scanners of 1.5T or 3T
% % % % % %         Sample_Phantom_generator_BRAIN_byRegion
% % % % % %         Sample_Phantom_generator_BRAIN_byRegion_vTest   % all tissues are alike in relaxation times
% % % %         Sample_Phantom_generator_BRAIN_EQUALRegion
% % % % % %         Sample_Phantom_generator_BRAIN_EQUALRegion_vTest
% % % %         
% % % %         % --- read sample per tissue ---
% % % %         sample_data = h5read('sample.h5','/sample/data');
% % % %         
% % % %         % --- plot per tissue ---
% % % %         figure();
% % % %         imagesc(reshape(sample_data(2,:,:),origFOV+2*widenSz,origFOV+2*widenSz))
% % % %         title(typeTissue)
% % % %         colormap gray
% % % %         
% % % %         % --- copy files per tissue ---
% % % %         cd(dir_seq)
% % % %         copyfile(name,dir_tissue{ppp})  
% % % %         copyfile(name_simu_template,dir_tissue{ppp})
% % % % 
% % % %         clear sample_data
% % % %     end
    
elseif testDiffusion == 'test2'
    %% Part 2 -  Organizacao dos dados + organizacao por pasta for diffusion
    
    tissue = {'CF','GM','WM','LS'}; % Select type tissue 'CF' - CSF; 'GM' - Gray Matter; 'WM' - White Matter; 'LS' - Lesion
    origFOV    = 180;
    nOfTissues = size(tissue,2);    
    if sequenc == 'SPI'
        name_simu = ['test_',num2str(test),'_Tissues_phantomBrain_','sm_',num2str(sm),'_T_m_s_gm_',num2str(gm), ...
            '_mT_m_FOV_',num2str(FOVX),'_N_',num2str(N),'_alpha_',num2str(alpha),...
            '_NP_',num2str(NP),'_RES_',num2str(res),'_Time1_',num2str(T1),...
            '_Tend_',num2str(Tend),'_Interleave_',num2str(ninterleaves),'_TR_',num2str(TR)...
            'teste_extgrad'];
    elseif sequenc == 'EPI'
        name_simu = ['EPItest_',num2str(EPItest),'_Tissues_phantomBrain_','sm_',num2str(sm),'_T_m_s_gm_',num2str(gm), ...
            '_T_m_FOV_',num2str(FOVX),'_N_',num2str(N),'_RES_',num2str(res),'_B0_',num2str(B0),'_TE_',num2str(TE),'_TR_',num2str(TR)];
    end
    
    name_simu = strrep(name_simu,'.','_');
    
    % ... Check if there is directory and create if not ...
    dir_IR = [directory_test,'/ImageRecon_files']; % run in asia    
    mkdir(dir_IR)
    cd(dir_IR)
    new_dir = [dir_IR ,'/', name_simu];
    isequal(exist(name_simu, 'dir'),7)
    if ans == 0
        mkdir(name_simu)
        cd(dir_seq)
    end
    clear ans signals
    
    % ... Copy signals data from the tissues (CSF, GM, WM and Lesion) ...
    name_signal_results = 'signals.h5';
    for i=1:size(tissue,2)
        typeTissue = tissue{i};
        % ... copy file ...
        dir_seq_tissue = [dir_seq,'/',typeTissue];
        cd(dir_seq_tissue)
        copyfile(name_signal_results,new_dir)
        % ... select directory ...
        cd(new_dir)
        name_simu = ['signals_',typeTissue,'.mat']; % name for '.mat' file
        % ... load signals ...
        read_h5signals(name_simu,'signals.h5')
        signals.(tissue{i})=load(name_simu);
        delete('signals.h5')
    end
    
elseif testDiffusion == 'test3'
    %% Part 3 - Test Diffusion after simulation
%     if typSamp   == 'DIFF'
        
        % ... resize das máscaras ...
        % % %     scale      = N/(FOVX*1e3);
        % % %     maskRS_CSF = imresize(mask_CSF,scale);
        % % %     maskRS_GM  = imresize(mask_GM,scale);
        % % %     maskRS_WM  = imresize(mask_WM,scale);
        % % %     maskRS_LS  = imresize(mask_LS,scale);
        % % %
        % ... factor due to simulate different samples in JEMRIS
        % % %         factorJemris(1) = mean(abs(complex(signals.CF.M(:,1),signals.CF.M(:,2)))); % data obtained with simulation
        % % %         factorJemris(2) = mean(abs(complex(signals.GM.M(:,1),signals.GM.M(:,2)))); % data obtained with simulation
        % % %         factorJemris(3) = mean(abs(complex(signals.WM.M(:,1),signals.WM.M(:,2)))); % data obtained with simulation
        % % %         factorJemris(4) = mean(abs(complex(signals.LS.M(:,1),signals.LS.M(:,2)))); % data obtained with simulation
        
        factorCSF = 1;%factorJemris(1)/min(factorJemris);
        factorGM = 1;%factorJemris(2)/min(factorJemris);
        factorWM = 1;%factorJemris(3)/min(factorJemris);
        factorLS = 1;%factorJemris(4)/min(factorJemris);
        
        % % %     factorCSF = size(find(maskRS_CSF>0),1);                             % due to simulate different samples in JEMRIS
        D_CSF     = 3.2*1e-3;                                               % Diffusivity of CSF mm^2 * s^-1 - Tese Nunes R.G.
        kdata_CSF = complex(signals.CF.M(:,1),signals.CF.M(:,2))/(factorCSF); % data obtained with simulation
        S_CSF     = kdata_CSF * exp (-b*D_CSF);
        S_CSF_mean  = mean(abs(S_CSF));
        
        % % %     factorGM = size(find(maskRS_GM>0),1);                               % due to simulate different samples in JEMRIS
        D_GM     = 0.8*1e-3;                                                % Diffusivity of GM - Tese Nunes R.G.
        kdata_GM = complex(signals.GM.M(:,1),signals.GM.M(:,2))/(factorGM);   % data obtained with simulation
        S_GM     = kdata_GM * exp (-b*D_GM);
        S_GM_mean  = mean(abs(S_GM));
        
        % % %     factorWM = size(find(maskRS_WM>0),1);                               % due to simulate different samples in JEMRIS
        D_WM     = 0.7*1e-3;                                                % Diffusivity of WM - Tese Nunes R.G.
        kdata_WM = complex(signals.WM.M(:,1),signals.WM.M(:,2))/(factorWM);   % data obtained with simulation
        S_WM     = kdata_WM * exp (-b*D_WM);
        S_WM_mean  = mean(abs(S_WM));
        
        % % %     factorLS = size(find(maskRS_LS>0),1);                               % due to simulate different samples in JEMRIS
        D_LS     = 0.55*1e-3;                                               % Diffusivity of LS -Knight, 2019
        kdata_LS = complex(signals.LS.M(:,1),signals.LS.M(:,2))/(factorLS);   % data obtained with simulation
        S_LS     = kdata_LS * exp (-b*D_LS);
        S_LS_mean  = mean(abs(S_LS));
        
        kdata    = kdata_CSF + kdata_GM + kdata_WM + kdata_LS;  % S0 data      - no diffusion
        kdata_DW = S_CSF + S_GM + S_WM + S_LS;       % S_Total data - diffusion signal from CSF + WM + GM + LS
        
        figure()
        plot(abs(kdata_CSF))
        hold on
        title('T2-weigthed Signal')
        plot(abs(kdata_GM))
        plot(abs(kdata_WM))
        plot(abs(kdata_LS))
        
        figure()
        plot(abs(S_CSF))
        hold on
        title('Diffusion Signal')
        plot(abs(S_GM))
        plot(abs(S_WM))
        plot(abs(S_LS))
               
% %         figure();
% %         plot(abs(kdata_DW))
% %         hold on
% %         plot(abs(kdata))

elseif testDiffusion == 'test4'
    %% Part 4 - Add different signals after fft
    clear m_postcomp m_postcomp_Tissues m_postcompTest imTest_crop imTest imTestTissue_crop
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

    % ... pre-compensation weights (independent of simulated data) ...
    [acf,flag,res] = makePrecompWeights_2D_VORONOI( kTraj, [gridsize/2 gridsize/2] );
    
    
    for i=1:size(tissue,2)
        typeTissue = tissue{i};
        if i == 1
            maskRS_CSF = imresize(mask_CSF,scale);            
            factorCSF  = size(find(maskRS_CSF>0),1);  
            kdata_test{i} = kdata_CSF*factorCSF;
        elseif i ==2
            maskRS_GM = imresize(mask_GM,scale);
            factorGM  = size(find(maskRS_CSF>0),1);  
            kdata_test{i} = kdata_GM*factorGM;
        elseif i ==3
            maskRS_WM = imresize(mask_WM,scale);
            factorWM  = size(find(maskRS_CSF>0),1);  
            kdata_test{i} = kdata_WM*factorWM;
        elseif i ==4
            maskRS_LS  = imresize(mask_LS,scale);
            factorLS  = size(find(maskRS_CSF>0),1);  
            kdata_test{i} = kdata_LS*factorLS;
        end
        
        % ... selecao de parametros ...
        overgridfactor = 3;
        kwidth         = 3; % kwidth = 50;  % kernel width
        kbeta          = pi*sqrt(kwidth^2/overgridfactor^2*(overgridfactor-0.5)^2-0.8 );  % kernel beta
        w              = acf; % w = weigth(kloc);
        
        % ... keiser-bessel gridding + fft ...
        %     trimming       -- 'y' or 'n'
        %     apodize        -- 'y' or 'n'
        %     postcompensate -- 'y' or 'n'
        [m_postcomp_Tissues{i},kwx,kwy] = grid_kb_teste141118(kdata_test{i},kloc,w,gridsize,overgridfactor,kwidth,kbeta,'n','n','y');                                
    end
    
    % ... select tissue ...
    m_postcompTest = m_postcomp_Tissues{1}+m_postcomp_Tissues{2}+m_postcomp_Tissues{3}+m_postcomp_Tissues{4};    
    
    % ... convert in image space ...
    imTest          = fftshift(fft2(fftshift(m_postcompTest)));       
    
    % ... crop image ...
    m_postcompTest_crop = m_postcompTest(gridsize-round(N/2)+1:gridsize+round(N/2),gridsize-round(N/2)+1:gridsize+round(N/2));
    imTest_crop_Tissue_beforeFFT = imTest(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2),round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2));    

    % ... plot figures ...
    figure;
    subplot(121)
    imshow(abs(imTest_crop_Tissue_beforeFFT),[]);    
    
    % ... pre-loop - take mask out of the brain ...
    allmaskRS = maskRS_CSF+maskRS_GM+maskRS_WM+maskRS_LS;
    allmaskRS_noHoles = imfill(allmaskRS,'holes');
    
    if tissueTest == 1
        figure()
        for jj=1:size(m_postcomp_Tissues,2)
            % ... select tissue ...
            m_postcompTestTissue{jj} = m_postcomp_Tissues{jj};                        
            % ... convert in image space ...
            imTestTissue{jj} = fftshift(fft2(fftshift(m_postcompTestTissue{jj}))); 
            % ... crop image ...
            imTestTissue_crop{jj} = imTestTissue{jj}(round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2),round((gridsize*overgridfactor)/2)-round(N/2)+1:round((gridsize*overgridfactor)/2)+round(N/2));
            % ... take mean out of the mask in the brain ...
            imTestTissue_crop{jj}(~allmaskRS_noHoles) = 0;            
            if jj == 1 % getting mask
                mask_Diff4 = maskRS_GM+maskRS_WM+maskRS_LS;
            elseif jj == 2
                mask_Diff4 = maskRS_CSF+maskRS_WM+maskRS_LS;                
            elseif jj == 3
                mask_Diff4 = maskRS_CSF+maskRS_GM+maskRS_LS;
            elseif jj == 4
                mask_Diff4 = maskRS_CSF+maskRS_GM+maskRS_WM;
            end        
            meanMaskTissue(jj) = mean(imTestTissue_crop{jj}(find(mask_Diff4>0)));
            mean_imTestTissue_crop{jj}=imTestTissue_crop{jj}-meanMaskTissue(jj);            
            clear mask_Diff4           
            mean_imTestTissue_crop{jj}(~allmaskRS_noHoles) = 0;                        
            % ... plot das figuras ...
            subplot(2,2,jj)
            imshow(abs(mean_imTestTissue_crop{jj}),[]);
        end     
        
        imTest_crop_Tissue_afterFFT = mean_imTestTissue_crop{1}+mean_imTestTissue_crop{2}+mean_imTestTissue_crop{3}+mean_imTestTissue_crop{4};
        % ... plot figures ...
        figure()
        subplot(122)
        imshow(abs(imTest_crop_Tissue_afterFFT),[]);
    end


    
elseif testDiffusion == 'test5'
    %% Part 5 - Add different DW in the image space
    clear imDW_TestAllTissues imDW_TestTissue_crop
    
    % ... Diffusivity for each Tissue ...
    D{1} = 3.2*1e-3;      % Diffusivity of CSF mm^2 * s^-1 - Tese Nunes R.G.
    D{2} = 0.8*1e-3;      % Diffusivity of GM - Tese Nunes R.G.
    D{3} = 0.7*1e-3;      % Diffusivity of WM - Tese Nunes R.G.
    D{4} = 0.55*1e-3;     % Diffusivity of LS -Knight, 2019
    
    figure()
    for ii=1:size(mean_imTestTissue_crop,2)
        imDW_TestTissue_crop{ii} = mean_imTestTissue_crop{ii}  *  exp (-b*D{ii}); 
        % ... plot das figuras ...
        subplot(2,2,ii)
        imshow(abs(imDW_TestTissue_crop{ii}),[]);
    end    
    % ... sum all tissues ... 
    imDW_TestAllTissues = imDW_TestTissue_crop{1}+imDW_TestTissue_crop{2}+imDW_TestTissue_crop{3}+imDW_TestTissue_crop{4};
    % ... plot figure ...
    figure()
    imshow(abs(imDW_TestAllTissues),[]);

end   

