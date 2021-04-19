%% Compress coils function using BART toolbox
%  -  works for 'getB0_ECC.m' implementation
%  -  outputs the kdata compressed for each of the particular situation
%  see Robison, 2019 paper (10.1002/mrm.27583)
% 
% TTFernandes, IST, Jan 2021

%%

function [k_cc] = compressCoils_BOECC(k,NadcK,Ny,nsl,nc,Npulses,Nreps,Naxes,file_folder)

if Naxes == 1
    
    % ----------------------- x-axis -----------------------------
    kx_ccswp_Pgrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Ngrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Pgrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Ngrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    kx_swp_Pgrad_Px0 = permute(k.kx_Pgrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Ngrad_Px0 = permute(k.kx_Ngrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Pgrad_Nx0 = permute(k.kx_Pgrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Ngrad_Nx0 = permute(k.kx_Ngrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    
    % ... part 5.1 - compression + concatenation ...
    kx_conc_Px0 = [kx_swp_Pgrad_Px0(:,:,:,:,1); kx_swp_Ngrad_Px0(:,:,:,:,1)];
    kx_conc_Nx0 = [kx_swp_Pgrad_Nx0(:,:,:,:,1); kx_swp_Ngrad_Nx0(:,:,:,:,1)];
    
    % ... part 5.1.1 - create ccoef for x-axis
    writecfl('kx_conc_Px0',kx_conc_Px0);
    system(sprintf('bart cc -p %d -M kx_conc_Px0 kx_Px0_cccoef_conc',nc));
    
    writecfl('kx_conc_Nx0',kx_conc_Nx0);
    system(sprintf('bart cc -p %d -M kx_conc_Nx0 kx_Nx0_cccoef_conc',nc));
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        kx_tmp_Pgrad_Px0 = kx_swp_Pgrad_Px0(:,:,:,:,v);
        writecfl('kx_tmp_Pgrad_Px0',kx_tmp_Pgrad_Px0);
        system(sprintf('bart ccapply -p %d kx_tmp_Pgrad_Px0 kx_Px0_cccoef_conc kx_tmpcc_Pgrad_Px0',nc));
        kx_ccswp_Pgrad_Px0(:,:,:,:,v)=readcfl('kx_tmpcc_Pgrad_Px0');
        
        kx_tmp_Ngrad_Px0 = kx_swp_Ngrad_Px0(:,:,:,:,v);
        writecfl('kx_tmp_Ngrad_Px0',kx_tmp_Ngrad_Px0);
        system(sprintf('bart ccapply -p %d kx_tmp_Ngrad_Px0 kx_Px0_cccoef_conc kx_tmpcc_Ngrad_Px0',nc));
        kx_ccswp_Ngrad_Px0(:,:,:,:,v)=readcfl('kx_tmpcc_Ngrad_Px0');
        
        kx_tmp_Pgrad_Nx0=kx_swp_Pgrad_Nx0(:,:,:,:,v);
        writecfl('kx_tmp_Pgrad_Nx0',kx_tmp_Pgrad_Nx0);
        system(sprintf('bart ccapply -p %d kx_tmp_Pgrad_Nx0 kx_Nx0_cccoef_conc kx_tmpcc_Pgrad_Nx0',nc));
        kx_ccswp_Pgrad_Nx0(:,:,:,:,v)=readcfl('kx_tmpcc_Pgrad_Nx0');
        
        kx_tmp_Ngrad_Nx0=kx_swp_Ngrad_Nx0(:,:,:,:,v);
        writecfl('kx_tmp_Ngrad_Nx0',kx_tmp_Ngrad_Nx0);
        system(sprintf('bart ccapply -p %d kx_tmp_Ngrad_Nx0 kx_Nx0_cccoef_conc kx_tmpcc_Ngrad_Nx0',nc));
        kx_ccswp_Ngrad_Nx0(:,:,:,:,v)=readcfl('kx_tmpcc_Ngrad_Nx0');
        
        clear kx_tmpcc_Pgrad_Px0 kx_tmpcc_Ngrad_Px0 kx_tmpcc_Pgrad_Nx0 kx_tmpcc_Ngrad_Nx0
    end
    
    cd(file_folder)
    
    % ... part 5.3 - invert permutation ...
    kx_cc_Pgrad_Px0 = ipermute(kx_ccswp_Pgrad_Px0,[1 3 4 2 5]);
    kx_cc_Ngrad_Px0 = ipermute(kx_ccswp_Ngrad_Px0,[1 3 4 2 5]);
    kx_cc_Pgrad_Nx0 = ipermute(kx_ccswp_Pgrad_Nx0,[1 3 4 2 5]);
    kx_cc_Ngrad_Nx0 = ipermute(kx_ccswp_Ngrad_Nx0,[1 3 4 2 5]);
    
    k_cc.kx_Pgrad_Px0 = reshape(kx_cc_Pgrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kx_Ngrad_Px0 = reshape(kx_cc_Ngrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kx_Pgrad_Nx0 = reshape(kx_cc_Pgrad_Nx0,[NadcK Npulses*Nreps]);
    k_cc.kx_Ngrad_Nx0 = reshape(kx_cc_Ngrad_Nx0,[NadcK Npulses*Nreps]);    
    
    
elseif Naxes == 2    
    % ----------------------- x-axis -----------------------------
    kx_ccswp_Pgrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Ngrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Pgrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Ngrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    kx_swp_Pgrad_Px0 = permute(k.kx_Pgrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Ngrad_Px0 = permute(k.kx_Ngrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Pgrad_Nx0 = permute(k.kx_Pgrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Ngrad_Nx0 = permute(k.kx_Ngrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes        
    
    % ... part 5.1 - compression + concatenation ...
    kx_conc_Px0 = [kx_swp_Pgrad_Px0(:,:,:,:,1); kx_swp_Ngrad_Px0(:,:,:,:,1)];
    kx_conc_Nx0 = [kx_swp_Pgrad_Nx0(:,:,:,:,1); kx_swp_Ngrad_Nx0(:,:,:,:,1)];
    
    % ... part 5.1.1 - create ccoef for x-axis
    writecfl('kx_conc_Px0',kx_conc_Px0);
    system(sprintf('bart cc -p %d -M kx_conc_Px0 kx_Px0_cccoef_conc',nc));

    writecfl('kx_conc_Nx0',kx_conc_Nx0);
    system(sprintf('bart cc -p %d -M kx_conc_Nx0 kx_Nx0_cccoef_conc',nc));        
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        kx_tmp_Pgrad_Px0 = kx_swp_Pgrad_Px0(:,:,:,:,v);
        writecfl('kx_tmp_Pgrad_Px0',kx_tmp_Pgrad_Px0);
        system(sprintf('bart ccapply -p %d kx_tmp_Pgrad_Px0 kx_Px0_cccoef_conc kx_tmpcc_Pgrad_Px0',nc));
        kx_ccswp_Pgrad_Px0(:,:,:,:,v)=readcfl('kx_tmpcc_Pgrad_Px0');
        
        kx_tmp_Ngrad_Px0 = kx_swp_Ngrad_Px0(:,:,:,:,v);
        writecfl('kx_tmp_Ngrad_Px0',kx_tmp_Ngrad_Px0);
        system(sprintf('bart ccapply -p %d kx_tmp_Ngrad_Px0 kx_Px0_cccoef_conc kx_tmpcc_Ngrad_Px0',nc));
        kx_ccswp_Ngrad_Px0(:,:,:,:,v)=readcfl('kx_tmpcc_Ngrad_Px0');
        
        kx_tmp_Pgrad_Nx0=kx_swp_Pgrad_Nx0(:,:,:,:,v);
        writecfl('kx_tmp_Pgrad_Nx0',kx_tmp_Pgrad_Nx0);
        system(sprintf('bart ccapply -p %d kx_tmp_Pgrad_Nx0 kx_Nx0_cccoef_conc kx_tmpcc_Pgrad_Nx0',nc));
        kx_ccswp_Pgrad_Nx0(:,:,:,:,v)=readcfl('kx_tmpcc_Pgrad_Nx0');
        
        kx_tmp_Ngrad_Nx0=kx_swp_Ngrad_Nx0(:,:,:,:,v);
        writecfl('kx_tmp_Ngrad_Nx0',kx_tmp_Ngrad_Nx0);
        system(sprintf('bart ccapply -p %d kx_tmp_Ngrad_Nx0 kx_Nx0_cccoef_conc kx_tmpcc_Ngrad_Nx0',nc));
        kx_ccswp_Ngrad_Nx0(:,:,:,:,v)=readcfl('kx_tmpcc_Ngrad_Nx0');
        
        clear kx_tmpcc_Pgrad_Px0 kx_tmpcc_Ngrad_Px0 kx_tmpcc_Pgrad_Nx0 kx_tmpcc_Ngrad_Nx0
    end
      
    % ... part 5.3 - invert permutation ...
    kx_cc_Pgrad_Px0 = ipermute(kx_ccswp_Pgrad_Px0,[1 3 4 2 5]);
    kx_cc_Ngrad_Px0 = ipermute(kx_ccswp_Ngrad_Px0,[1 3 4 2 5]);
    kx_cc_Pgrad_Nx0 = ipermute(kx_ccswp_Pgrad_Nx0,[1 3 4 2 5]);
    kx_cc_Ngrad_Nx0 = ipermute(kx_ccswp_Ngrad_Nx0,[1 3 4 2 5]);
    
    k_cc.kx_Pgrad_Px0 = reshape(kx_cc_Pgrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kx_Ngrad_Px0 = reshape(kx_cc_Ngrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kx_Pgrad_Nx0 = reshape(kx_cc_Pgrad_Nx0,[NadcK Npulses*Nreps]);
    k_cc.kx_Ngrad_Nx0 = reshape(kx_cc_Ngrad_Nx0,[NadcK Npulses*Nreps]);
     
    
    % ----------------------- y-axis -----------------------------
    ky_ccswp_Pgrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_Ngrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_Pgrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_Ngrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    ky_swp_Pgrad_Px0 = permute(k.ky_Pgrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_Ngrad_Px0 = permute(k.ky_Ngrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_Pgrad_Nx0 = permute(k.ky_Pgrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_Ngrad_Nx0 = permute(k.ky_Ngrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes        
    
    % ... part 5.1 - compression + concatenation ...
    ky_conc_Px0 = [ky_swp_Pgrad_Px0(:,:,:,:,1); ky_swp_Ngrad_Px0(:,:,:,:,1)];
    ky_conc_Nx0 = [ky_swp_Pgrad_Nx0(:,:,:,:,1); ky_swp_Ngrad_Nx0(:,:,:,:,1)];
    
    % ... part 5.1.1 - create ccoef for x-axis
    writecfl('ky_conc_Px0',ky_conc_Px0);
    system(sprintf('bart cc -p %d -M ky_conc_Px0 ky_Px0_cccoef_conc',nc));

    writecfl('ky_conc_Nx0',ky_conc_Nx0);
    system(sprintf('bart cc -p %d -M ky_conc_Nx0 ky_Nx0_cccoef_conc',nc));        
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        ky_tmp_Pgrad_Px0 = ky_swp_Pgrad_Px0(:,:,:,:,v);
        writecfl('ky_tmp_Pgrad_Px0',ky_tmp_Pgrad_Px0);
        system(sprintf('bart ccapply -p %d ky_tmp_Pgrad_Px0 ky_Px0_cccoef_conc ky_tmpcc_Pgrad_Px0',nc));
        ky_ccswp_Pgrad_Px0(:,:,:,:,v)=readcfl('ky_tmpcc_Pgrad_Px0');
        
        ky_tmp_Ngrad_Px0 = ky_swp_Ngrad_Px0(:,:,:,:,v);
        writecfl('ky_tmp_Ngrad_Px0',ky_tmp_Ngrad_Px0);
        system(sprintf('bart ccapply -p %d ky_tmp_Ngrad_Px0 ky_Px0_cccoef_conc ky_tmpcc_Ngrad_Px0',nc));
        ky_ccswp_Ngrad_Px0(:,:,:,:,v)=readcfl('ky_tmpcc_Ngrad_Px0');
        
        ky_tmp_Pgrad_Nx0=ky_swp_Pgrad_Nx0(:,:,:,:,v);
        writecfl('ky_tmp_Pgrad_Nx0',ky_tmp_Pgrad_Nx0);
        system(sprintf('bart ccapply -p %d ky_tmp_Pgrad_Nx0 ky_Nx0_cccoef_conc ky_tmpcc_Pgrad_Nx0',nc));
        ky_ccswp_Pgrad_Nx0(:,:,:,:,v)=readcfl('ky_tmpcc_Pgrad_Nx0');
        
        ky_tmp_Ngrad_Nx0=ky_swp_Ngrad_Nx0(:,:,:,:,v);
        writecfl('ky_tmp_Ngrad_Nx0',ky_tmp_Ngrad_Nx0);
        system(sprintf('bart ccapply -p %d ky_tmp_Ngrad_Nx0 ky_Nx0_cccoef_conc ky_tmpcc_Ngrad_Nx0',nc));
        ky_ccswp_Ngrad_Nx0(:,:,:,:,v)=readcfl('ky_tmpcc_Ngrad_Nx0');
        
        clear ky_tmpcc_Pgrad_Px0 ky_tmpcc_Ngrad_Px0 ky_tmpcc_Pgrad_Nx0 ky_tmpcc_Ngrad_Nx0
    end
      
    % ... part 5.3 - invert permutation ...
    ky_cc_Pgrad_Px0 = ipermute(ky_ccswp_Pgrad_Px0,[1 3 4 2 5]);
    ky_cc_Ngrad_Px0 = ipermute(ky_ccswp_Ngrad_Px0,[1 3 4 2 5]);
    ky_cc_Pgrad_Nx0 = ipermute(ky_ccswp_Pgrad_Nx0,[1 3 4 2 5]);
    ky_cc_Ngrad_Nx0 = ipermute(ky_ccswp_Ngrad_Nx0,[1 3 4 2 5]);
    
    k_cc.ky_Pgrad_Px0 = reshape(ky_cc_Pgrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.ky_Ngrad_Px0 = reshape(ky_cc_Ngrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.ky_Pgrad_Nx0 = reshape(ky_cc_Pgrad_Nx0,[NadcK Npulses*Nreps]);
    k_cc.ky_Ngrad_Nx0 = reshape(ky_cc_Ngrad_Nx0,[NadcK Npulses*Nreps]);
    
    cd(file_folder)

    
    
elseif Naxes == 3
    % ----------------------- x-axis -----------------------------
    kx_ccswp_Pgrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Ngrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Pgrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kx_ccswp_Ngrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    kx_swp_Pgrad_Px0 = permute(k.kx_Pgrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Ngrad_Px0 = permute(k.kx_Ngrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Pgrad_Nx0 = permute(k.kx_Pgrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kx_swp_Ngrad_Nx0 = permute(k.kx_Ngrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes        
    
    % ... part 5.1 - compression + concatenation ...
    kx_conc_Px0 = [kx_swp_Pgrad_Px0(:,:,:,:,1); kx_swp_Ngrad_Px0(:,:,:,:,1)];
    kx_conc_Nx0 = [kx_swp_Pgrad_Nx0(:,:,:,:,1); kx_swp_Ngrad_Nx0(:,:,:,:,1)];
    
    % ... part 5.1.1 - create ccoef for x-axis
    writecfl('kx_conc_Px0',kx_conc_Px0);
    system(sprintf('bart cc -p %d -M kx_conc_Px0 kx_Px0_cccoef_conc',nc));

    writecfl('kx_conc_Nx0',kx_conc_Nx0);
    system(sprintf('bart cc -p %d -M kx_conc_Nx0 kx_Nx0_cccoef_conc',nc));        
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        kx_tmp_Pgrad_Px0 = kx_swp_Pgrad_Px0(:,:,:,:,v);
        writecfl('kx_tmp_Pgrad_Px0',kx_tmp_Pgrad_Px0);
        system(sprintf('bart ccapply -p %d kx_tmp_Pgrad_Px0 kx_Px0_cccoef_conc kx_tmpcc_Pgrad_Px0',nc));
        kx_ccswp_Pgrad_Px0(:,:,:,:,v)=readcfl('kx_tmpcc_Pgrad_Px0');
        
        kx_tmp_Ngrad_Px0 = kx_swp_Ngrad_Px0(:,:,:,:,v);
        writecfl('kx_tmp_Ngrad_Px0',kx_tmp_Ngrad_Px0);
        system(sprintf('bart ccapply -p %d kx_tmp_Ngrad_Px0 kx_Px0_cccoef_conc kx_tmpcc_Ngrad_Px0',nc));
        kx_ccswp_Ngrad_Px0(:,:,:,:,v)=readcfl('kx_tmpcc_Ngrad_Px0');
        
        kx_tmp_Pgrad_Nx0=kx_swp_Pgrad_Nx0(:,:,:,:,v);
        writecfl('kx_tmp_Pgrad_Nx0',kx_tmp_Pgrad_Nx0);
        system(sprintf('bart ccapply -p %d kx_tmp_Pgrad_Nx0 kx_Nx0_cccoef_conc kx_tmpcc_Pgrad_Nx0',nc));
        kx_ccswp_Pgrad_Nx0(:,:,:,:,v)=readcfl('kx_tmpcc_Pgrad_Nx0');
        
        kx_tmp_Ngrad_Nx0=kx_swp_Ngrad_Nx0(:,:,:,:,v);
        writecfl('kx_tmp_Ngrad_Nx0',kx_tmp_Ngrad_Nx0);
        system(sprintf('bart ccapply -p %d kx_tmp_Ngrad_Nx0 kx_Nx0_cccoef_conc kx_tmpcc_Ngrad_Nx0',nc));
        kx_ccswp_Ngrad_Nx0(:,:,:,:,v)=readcfl('kx_tmpcc_Ngrad_Nx0');
        
        clear kx_tmpcc_Pgrad_Px0 kx_tmpcc_Ngrad_Px0 kx_tmpcc_Pgrad_Nx0 kx_tmpcc_Ngrad_Nx0
    end
      
    % ... part 5.3 - invert permutation ...
    kx_cc_Pgrad_Px0 = ipermute(kx_ccswp_Pgrad_Px0,[1 3 4 2 5]);
    kx_cc_Ngrad_Px0 = ipermute(kx_ccswp_Ngrad_Px0,[1 3 4 2 5]);
    kx_cc_Pgrad_Nx0 = ipermute(kx_ccswp_Pgrad_Nx0,[1 3 4 2 5]);
    kx_cc_Ngrad_Nx0 = ipermute(kx_ccswp_Ngrad_Nx0,[1 3 4 2 5]);
    
    k_cc.kx_Pgrad_Px0 = reshape(kx_cc_Pgrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kx_Ngrad_Px0 = reshape(kx_cc_Ngrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kx_Pgrad_Nx0 = reshape(kx_cc_Pgrad_Nx0,[NadcK Npulses*Nreps]);
    k_cc.kx_Ngrad_Nx0 = reshape(kx_cc_Ngrad_Nx0,[NadcK Npulses*Nreps]);
     
    
    % ----------------------- y-axis -----------------------------
    ky_ccswp_Pgrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_Ngrad_Px0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_Pgrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    ky_ccswp_Ngrad_Nx0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    ky_swp_Pgrad_Px0 = permute(k.ky_Pgrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_Ngrad_Px0 = permute(k.ky_Ngrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_Pgrad_Nx0 = permute(k.ky_Pgrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    ky_swp_Ngrad_Nx0 = permute(k.ky_Ngrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes        
    
    % ... part 5.1 - compression + concatenation ...
    ky_conc_Px0 = [ky_swp_Pgrad_Px0(:,:,:,:,1); ky_swp_Ngrad_Px0(:,:,:,:,1)];
    ky_conc_Nx0 = [ky_swp_Pgrad_Nx0(:,:,:,:,1); ky_swp_Ngrad_Nx0(:,:,:,:,1)];
    
    % ... part 5.1.1 - create ccoef for x-axis
    writecfl('ky_conc_Px0',ky_conc_Px0);
    system(sprintf('bart cc -p %d -M ky_conc_Px0 ky_Px0_cccoef_conc',nc));

    writecfl('ky_conc_Nx0',ky_conc_Nx0);
    system(sprintf('bart cc -p %d -M ky_conc_Nx0 ky_Nx0_cccoef_conc',nc));        
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        ky_tmp_Pgrad_Px0 = ky_swp_Pgrad_Px0(:,:,:,:,v);
        writecfl('ky_tmp_Pgrad_Px0',ky_tmp_Pgrad_Px0);
        system(sprintf('bart ccapply -p %d ky_tmp_Pgrad_Px0 ky_Px0_cccoef_conc ky_tmpcc_Pgrad_Px0',nc));
        ky_ccswp_Pgrad_Px0(:,:,:,:,v)=readcfl('ky_tmpcc_Pgrad_Px0');
        
        ky_tmp_Ngrad_Px0 = ky_swp_Ngrad_Px0(:,:,:,:,v);
        writecfl('ky_tmp_Ngrad_Px0',ky_tmp_Ngrad_Px0);
        system(sprintf('bart ccapply -p %d ky_tmp_Ngrad_Px0 ky_Px0_cccoef_conc ky_tmpcc_Ngrad_Px0',nc));
        ky_ccswp_Ngrad_Px0(:,:,:,:,v)=readcfl('ky_tmpcc_Ngrad_Px0');
        
        ky_tmp_Pgrad_Nx0=ky_swp_Pgrad_Nx0(:,:,:,:,v);
        writecfl('ky_tmp_Pgrad_Nx0',ky_tmp_Pgrad_Nx0);
        system(sprintf('bart ccapply -p %d ky_tmp_Pgrad_Nx0 ky_Nx0_cccoef_conc ky_tmpcc_Pgrad_Nx0',nc));
        ky_ccswp_Pgrad_Nx0(:,:,:,:,v)=readcfl('ky_tmpcc_Pgrad_Nx0');
        
        ky_tmp_Ngrad_Nx0=ky_swp_Ngrad_Nx0(:,:,:,:,v);
        writecfl('ky_tmp_Ngrad_Nx0',ky_tmp_Ngrad_Nx0);
        system(sprintf('bart ccapply -p %d ky_tmp_Ngrad_Nx0 ky_Nx0_cccoef_conc ky_tmpcc_Ngrad_Nx0',nc));
        ky_ccswp_Ngrad_Nx0(:,:,:,:,v)=readcfl('ky_tmpcc_Ngrad_Nx0');
        
        clear ky_tmpcc_Pgrad_Px0 ky_tmpcc_Ngrad_Px0 ky_tmpcc_Pgrad_Nx0 ky_tmpcc_Ngrad_Nx0
    end
      
    % ... part 5.3 - invert permutation ...
    ky_cc_Pgrad_Px0 = ipermute(ky_ccswp_Pgrad_Px0,[1 3 4 2 5]);
    ky_cc_Ngrad_Px0 = ipermute(ky_ccswp_Ngrad_Px0,[1 3 4 2 5]);
    ky_cc_Pgrad_Nx0 = ipermute(ky_ccswp_Pgrad_Nx0,[1 3 4 2 5]);
    ky_cc_Ngrad_Nx0 = ipermute(ky_ccswp_Ngrad_Nx0,[1 3 4 2 5]);
    
    k_cc.ky_Pgrad_Px0 = reshape(ky_cc_Pgrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.ky_Ngrad_Px0 = reshape(ky_cc_Ngrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.ky_Pgrad_Nx0 = reshape(ky_cc_Pgrad_Nx0,[NadcK Npulses*Nreps]);
    k_cc.ky_Ngrad_Nx0 = reshape(ky_cc_Ngrad_Nx0,[NadcK Npulses*Nreps]);
           
    
    % ----------------------- z-axis -----------------------------
    kz_ccswp_Pgrad_Py0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kz_ccswp_Ngrad_Py0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kz_ccswp_Pgrad_Ny0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    kz_ccswp_Ngrad_Ny0 = zeros(NadcK,Ny, nsl, nc, Npulses*Nreps);
    
    % ... part 5.0.1 - permute for bart ...
    kz_swp_Pgrad_Px0 = permute(k.kz_Pgrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kz_swp_Ngrad_Px0 = permute(k.kz_Ngrad_Px0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kz_swp_Pgrad_Nx0 = permute(k.kz_Pgrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes
    kz_swp_Ngrad_Nx0 = permute(k.kz_Ngrad_Nx0,[1 3 4 2 5]);   %Nx Ny Nslice Ncoils Nvolumes        
    
    % ... part 5.1 - compression + concatenation ...
    kz_conc_Px0 = [kz_ccswp_Pgrad_Px0(:,:,:,:,1); kz_swp_Ngrad_Px0(:,:,:,:,1)];
    kz_conc_Nx0 = [kz_ccswp_Pgrad_Nx0(:,:,:,:,1); kz_swp_Ngrad_Nx0(:,:,:,:,1)];
    
    % ... part 5.1.1 - create ccoef for x-axis
    writecfl('kz_conc_Px0',kz_conc_Px0);
    system(sprintf('bart cc -p %d -M kz_conc_Px0 kz_Px0_cccoef_conc',nc));

    writecfl('kz_conc_Nx0',kz_conc_Nx0);
    system(sprintf('bart cc -p %d -M kz_conc_Nx0 kz_Nx0_cccoef_conc',nc));        
    
    % ... part 5.2 - apply compression ...
    for v=1:Npulses*Nreps
        kz_tmp_Pgrad_Px0 = kz_swp_Pgrad_Px0(:,:,:,:,v);
        writecfl('kz_tmp_Pgrad_Px0',kz_tmp_Pgrad_Px0);
        system(sprintf('bart ccapply -p %d kz_tmp_Pgrad_Px0 kz_Px0_cccoef_conc kz_tmpcc_Pgrad_Px0',nc));
        kz_ccswp_Pgrad_Px0(:,:,:,:,v)=readcfl('kz_tmpcc_Pgrad_Px0');
        
        kz_tmp_Ngrad_Px0 = kz_swp_Ngrad_Px0(:,:,:,:,v);
        writecfl('kz_tmp_Ngrad_Px0',kz_tmp_Ngrad_Px0);
        system(sprintf('bart ccapply -p %d kz_tmp_Ngrad_Px0 kz_Px0_cccoef_conc kz_tmpcc_Ngrad_Px0',nc));
        kz_ccswp_Ngrad_Px0(:,:,:,:,v)=readcfl('kz_tmpcc_Ngrad_Px0');
        
        kz_tmp_Pgrad_Nx0 = kz_swp_Pgrad_Nx0(:,:,:,:,v);
        writecfl('kz_tmp_Pgrad_Nx0',kz_tmp_Pgrad_Nx0);
        system(sprintf('bart ccapply -p %d kz_tmp_Pgrad_Nx0 kz_Nx0_cccoef_conc kz_tmpcc_Pgrad_Nx0',nc));
        kz_ccswp_Pgrad_Nx0(:,:,:,:,v)=readcfl('kz_tmpcc_Pgrad_Nx0');
        
        kz_tmp_Ngrad_Nx0 = kz_swp_Ngrad_Nx0(:,:,:,:,v);
        writecfl('kz_tmp_Ngrad_Nx0',kz_tmp_Ngrad_Nx0);
        system(sprintf('bart ccapply -p %d kz_tmp_Ngrad_Nx0 kz_Nx0_cccoef_conc kz_tmpcc_Ngrad_Nx0',nc));
        kz_ccswp_Ngrad_Nx0(:,:,:,:,v)=readcfl('kz_tmpcc_Ngrad_Nx0');
        
        clear kz_tmpcc_Pgrad_Px0 kz_tmpcc_Ngrad_Px0 kz_tmpcc_Pgrad_Nx0 kz_tmpcc_Ngrad_Nx0
    end
        
    % ... part 5.3 - invert permutation ...
    kz_cc_Pgrad_Px0 = ipermute(kz_ccswp_Pgrad_Px0,[1 3 4 2 5]);
    kz_cc_Ngrad_Px0 = ipermute(kz_ccswp_Ngrad_Px0,[1 3 4 2 5]);
    kz_cc_Pgrad_Nx0 = ipermute(kz_ccswp_Pgrad_Nx0,[1 3 4 2 5]);
    kz_cc_Ngrad_Nx0 = ipermute(kz_ccswp_Ngrad_Nx0,[1 3 4 2 5]);
    
    k_cc.kz_Pgrad_Px0 = reshape(kz_cc_Pgrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kz_Ngrad_Px0 = reshape(kz_cc_Ngrad_Px0,[NadcK Npulses*Nreps]);
    k_cc.kz_Pgrad_Nx0 = reshape(kz_cc_Pgrad_Nx0,[NadcK Npulses*Nreps]);
    k_cc.kz_Ngrad_Nx0 = reshape(kz_cc_Ngrad_Nx0,[NadcK Npulses*Nreps]);    

    cd(file_folder)    
end
end