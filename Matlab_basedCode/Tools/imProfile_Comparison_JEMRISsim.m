%% Comparison from JEMRIS Simulation IRF with Theor
clear all
clc

%% Initialization
Ntest = 3;
testNumber1 = 26;
testNumber2 = 27;


%% Load signals
myCD1 = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/SPI_test_',num2str(testNumber1)];
myCD2 = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/SPI_test_',num2str(testNumber2)];

% comparison test
cd(myCD1)
test2 = load(['test_',num2str(testNumber1),'_imProfile.mat']);
cd(myCD2)
test3 = load(['test_',num2str(testNumber2),'_imProfile.mat']);
myCD3 = ['/home/tfernandes/Documents/Projetos/Project_lfStroke/Code/0_matlabCode/4_Simulation_IRF/Input_sequences/SPI_test_',num2str(test2.testSPI)];
cd(myCD3)
test1 = load(['test_',num2str(test2.testSPI),'_imProfile.mat']);

%% Plots

x_im(1,:)     = [size(test1.Norm_im_crop,1)/10*3  size(test1.Norm_im_crop,1)/10*3];
x_im(2,:)     = [size(test1.Norm_im_crop,1)/10*3  size(test1.Norm_im_crop,1)/10*7];
x_im(3,:)     = [size(test1.Norm_im_crop,1)/10*0  size(test1.Norm_im_crop,1)/10*10];

y_im(1,:)     = [size(test1.Norm_im_crop,1)/10*0 size(test1.Norm_im_crop,1)/10*10];
y_im(2,:)     = [size(test1.Norm_im_crop,1)/10*0 size(test1.Norm_im_crop,1)/10*10];
y_im(3,:)     = [size(test1.Norm_im_crop,1)/10*6 size(test1.Norm_im_crop,1)/10*6];

x_y_line = 2;

figure()
subplot(231)
plot(test1.c_im(:,x_y_line),'r')
hold on
title(['Im-Profile Theor w/ Theor'])
subplot(232)
plot(test2.c_im_noH(:,x_y_line),'b')
hold on
title(['Im-Profile IRF w/ Theor'])
subplot(233)
plot(test2.c_im(:,x_y_line),'g')
hold on
title(['Im-Profile IRF w/ IRF'])
% % subplot(254)
% % plot(test3.c_im_noH(:,x_y_line),'g')
% % hold on
% % title(['Im-Profile IRF w/ Theor - Diff'])
% % subplot(255)
% % plot(test3.c_im(:,x_y_line),'g')
% % hold on
% % title(['Im-Profile IRF w/ IRF - Diff'])

subplot(234)
imshow(test1.Norm_im_noH_crop,[]);
hold on
line(x_im(x_y_line,:),y_im(x_y_line,:))
title('Image Recon Theor w/ Theor')
subplot(235)
imshow(test2.Norm_im_noH_crop,[]);
hold on
line(x_im(x_y_line,:),y_im(x_y_line,:))
title('Image Recon IRF w/ Theor')
subplot(236)
imshow(test2.Norm_im_crop,[]);
hold on
line(x_im(x_y_line,:),y_im(x_y_line,:))
title('Image Recon IRF w/ IRF')
% % subplot(259)
% % imshow(test3.Norm_im_noH_crop,[]);
% % hold on
% % line(x_im(x_y_line,:),y_im(x_y_line,:))
% % title('Image Recon IRF w/ Theor - Diff')
% % subplot(2,5,10)
% % imshow(test3.Norm_im_crop,[]);
% % hold on
% % line(x_im(x_y_line,:),y_im(x_y_line,:))
% % title('Image Recon IRF w/ IRF - Diff')


% % 
% % figure()
% % % subplot(231)
% % plot(test1.c_im(:,x_y_line),'r')
% % hold on
% % % title(['Im-Profile Theor w/ Theor'])
% % % subplot(232)
% % plot(test2.c_im_noH(:,x_y_line),'b')
% % % hold on
% % % title(['Im-Profile IRF w/ Theor'])
% % % subplot(233)
% % plot(test2.c_im(:,x_y_line),'g')
% % % hold on
% % title(['Im-Profile'])