function generate_h5_grad( gan,tt, NP, co , dir, seq, type)

%generate_h5_grad Saves the k sapce files in the format h5
%  
name2=num2str(co);
cd(dir)
% --- apaga ficheiros existitentes ---
exist strcat('mGY',name2,'.','h5') dir;
% if ans>0
    delete (strcat(pwd,'/',strcat('mGY',name2,'.','h5')))
% end
exist strcat('mGX',name2,'.','h5') dir;
% if ans>0
    delete (strcat(pwd,'/',strcat('mGX',name2,'.','h5')))
% end

% --- cria ficheiros .h5
if isequal(seq,'SPI') && isequal(type,'SIMU') % --- spiral test ---
    MGY=[tt(1+(co-1)*NP:NP*co)' (imag(gan))'*((2)*pi)];
    MGX=[tt(1+(co-1)*NP:NP*co)' (real(gan))'*((2)*pi)];
    
elseif isequal(seq,'SPI') && isequal(type,'SCAN')% --- spiral test ---
    MGY=[tt(1+(co-1)*NP:NP*co)' (imag(gan))'];
    MGX=[tt(1+(co-1)*NP:NP*co)' (real(gan))'];
    
elseif isequal(seq,'EPI') % --- epi test ---
    MGY=[tt(1+(co-1)*NP:NP*co)' (gan(2,:))'];
    MGX=[tt(1+(co-1)*NP:NP*co)' (gan(1,:))'];
end


%% test saving
% % % name_seq_GY = strcat('mGY',name2,'.','h5');
% % % name_seq_GX = strcat('mGX',name2,'.','h5');
% % % 
% % % GX_afterread = h5read(name_seq_GX,'/extpulse');
% % % GY_afterread = h5read(name_seq_GY,'/extpulse');
% % % 
% % % figure()
% % % plot(MGX(:,1),MGX(:,2),'b')
% % % hold on
% % % plot(MGY(:,1),MGY(:,2),'r')
% % % hold on
% % % plot(GX_afterread(:,1),GX_afterread(:,2),'b')
% % % hold on
% % % plot(GY_afterread(:,1),GY_afterread(:,2),'r')

%% save .h5 files
cd(dir);  % change directory where is going to save file
h5create(strcat('mGY',name2,'.','h5'),'/extpulse',[NP 2])
h5write(strcat('mGY',name2,'.','h5'), '/extpulse', MGY)

h5create(strcat('mGX',name2,'.','h5'),'/extpulse',[NP 2])
h5write(strcat('mGX',name2,'.','h5'), '/extpulse', MGX)

%% Plot - check if data is save properly
% % % name_seq_GY = strcat('mGY',name2,'.','h5');
% % % name_seq_GX = strcat('mGX',name2,'.','h5');
% % % 
% % % GX_afterread = h5read(name_seq_GX,'/extpulse');
% % % GY_afterread = h5read(name_seq_GY,'/extpulse');
% % % 
% % % figure()
% % % plot(MGX(:,1),MGX(:,2),'b')
% % % hold on
% % % plot(MGY(:,1),MGY(:,2),'r')
% % % hold on
% % % plot(GX_afterread(:,1),GX_afterread(:,2),'b')
% % % hold on
% % % plot(GY_afterread(:,1),GY_afterread(:,2),'r')

end