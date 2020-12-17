%% Save sample with two spheres

%% Read signal
cd(dir_Orig)
name_seq_template = 'sample_2Spheres.h5';


sample_data = h5read('sample_2Spheres.h5','/sample/data');
resolutioN = h5read('sample_2Spheres.h5','/sample/resolution');
offsetS = h5read('sample_2Spheres.h5','/sample/offset');


%%
[layers,Nx,Ny] = size(sample_data);
A_wide         = zeros(layers,Nx+2*widenSz,Ny+2*widenSz); % create widen brain
fig1 = figure(33);
for i=1:layers
    A_wide(i,widenSz+1:end-widenSz,widenSz+1:end-widenSz) = sample_data(i,:,:);
    subplot(1,5,i)
    imagesc(reshape(A_wide(i,:,:),Nx+2*widenSz,Ny+2*widenSz))
end
close(fig1)

clear A
A = A_wide;


r = 1;
o = 0;
OFFSET     = [o,o,o];
RESOLUTION = [r,r,r];

%% save sample_2Spheres


cd(cd_dir)
nameSample=['sample']
delete('sample.h5')
SAMPH5=[pwd,strcat('/',nameSample,'.h5')];
h5create(SAMPH5,strcat('/',nameSample,'/data'),size(A));
h5create(SAMPH5,strcat('/',nameSample,'/resolution'),[1 3]);
h5create(SAMPH5,strcat('/',nameSample,'/offset'),[1 3]);
h5write(SAMPH5,strcat('/',nameSample,'/data'), A);
h5write(SAMPH5,strcat('/',nameSample,'/resolution'),RESOLUTION);
h5write(SAMPH5,strcat('/',nameSample,'/offset'),OFFSET);