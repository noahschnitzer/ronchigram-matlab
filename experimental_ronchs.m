clear; close all; clc;
%%
fdir = 'data/';
fname_S6 = '6S_a_60_cl2_1_2Mx_300kV_spot7_8cm_m350_350_70im_1us.tif';
fname_S5 = '5S_a_60_cl2_1_2Mx_300kV_spot7_8cm_m350_350_70im_1us.tif';

len = 70;
dim = 2048;
S5 = zeros(dim,dim,len);
S6 = zeros(dim,dim,len);
for it = 1:len
    it
    S5(:,:,it) = imread([fdir fname_S5],it);
    S6(:,:,it) = imread([fdir fname_S6],it);
end


%%
defocuses = -350:10:349;
posx = 380;
posy = 320;
wh = 1200;
crop_x = posx:posx+wh; 
crop_y = posy:posy+wh;
%figure; imagesc(A(crop_y,crop_x,1)); colormap gray; axis image;
%%
% cropped_A = A(crop_y,crop_x,:);
% %figure;imagesc(cropped_A(:,:,1));
% %figure;imagesc(masked(:,:,1));
% masked = A;
% masked(crop_y,crop_x) = 0;
preds_S5 = [];
preds_S6 = [];
load('data/net-26.mat');

for it = 1:len
    it
   im =  S5(crop_y,crop_x,it);
   im = imresize(im,[512 512]);
   im = im - min(im(:));
   im = uint8(im.*255/max(im(:)));
   preds_S5(it) = predict(trainedNet,im);
   
   im =  S6(crop_y,crop_x,it);
   im = imresize(im,[512 512]);
   im = im - min(im(:));
   im = uint8(im.*255/max(im(:)));
   preds_S6(it) = predict(trainedNet,im);

end
%%
caf = 55;
chosen_indices = [21:51 36];
it = 36;
colordef;
c_net_stig = [1 0.5 .5];
figure; 
%plot(defocuses(chosen_indices(1:end-1)),caf.*preds_S5(chosen_indices(1:end-1)),'LineWidth',5,'Color',c_net_stig);
hold on;
plot(defocuses(chosen_indices(1:end-1)),caf.*preds_S6(chosen_indices(1:end-1)),'LineWidth',5,'Color',c_net);

xlim([defocuses(chosen_indices(1)) defocuses(chosen_indices(end-1))]);
ylim([.16 .35].*caf);
set(gca,'FontSize',25);
xlabel('Defocus (nm)','FontSize',35);
ylabel('Prediction (mrad)','FontSize',35);
set(gca,'FontName','Helvetica Neue');

f = figure('Position',[10 10 500 500]);
imagesc(S5(crop_y,crop_x,it)); colormap gray; axis image off;
c = floor(length(crop_x)/2);
viscircles([c c],preds_S5(it).*c,'Color',c_net_stig,'EnhanceVisibility',false,'LineWidth',6);
set(gca,'Position',[0 0 1 1]);

f = figure('Position',[10 10 500 500]);
imagesc(S6(crop_y,crop_x,it)); colormap gray; axis image off;
c = floor(length(crop_x)/2);
viscircles([c c],preds_S6(it).*c,'Color',c_net,'EnhanceVisibility',false,'LineWidth',6);
set(gca,'Position',[0 0 1 1]);

%% Comparing with expected strehl
ab = aberration_generator(1);
ab.mag = ab.mag.*1;
defs = [-150:1:150].*10;
thru_foc_strehl = aberration_series(ab,1024,180,128,1,defs,.8,0);
%% and plotting
figure;
plot(defs,thru_foc_strehl);
%xlim([-150 150]);
%% vids etc below
colordef;
caf = 55;
chosen_indices = [21:51 36];

v = VideoWriter('output/expt_ab_series_plots.avi');
v.FrameRate = 5;
open(v);
for it = chosen_indices(1:end-1)
    it
    colordef;
    c_net_stig = [1 0.5 .5];
    f = figure; 
    %plot(defocuses(chosen_indices(1):it),caf.*preds_S5(chosen_indices(1):it),'LineWidth',5,'Color',c_net_stig);
    hold on;
    plot(defocuses(chosen_indices(1):it),caf.*preds_S6(chosen_indices(1):it),'LineWidth',5,'Color',c_net);

    xlim([defocuses(chosen_indices(1)) defocuses(chosen_indices(end-1))]);
    ylim([.16 .35].*caf);
    set(gca,'FontSize',20);
    xlabel('Defocus (nm)','FontSize',30);
    ylabel('Prediction (mrad)','FontSize',30);
    set(gca,'FontName','Helvetica Neue');
    drawnow;
    print('output/temp.png','-dpng');
    writeVideo(v,imread('output/temp.png'));
    close(f);
end
close(v);
%%
% S5
v = VideoWriter('output/5S_expt_ab_series.avi');
v.FrameRate = 5;
open(v);

for it = chosen_indices
    it
    f = figure('Position',[10 10 500 500]);
    imagesc(S5(crop_y,crop_x,it)); colormap gray; axis image off;
    c = floor(length(crop_x)/2);
    viscircles([c c],preds_S5(it).*c,'Color',c_net_stig,'EnhanceVisibility',false,'LineWidth',6);
    set(gca,'Position',[0 0 1 1]);
    drawnow;
    print('output/temp.png','-dpng');
    writeVideo(v,imread('output/temp.png'));
    close(f);
end
close(v);
%%
% S6
colordef;
v = VideoWriter('output/6S_expt_ab_series.avi');
v.FrameRate = 5;
open(v);

for it = chosen_indices
    it
    f = figure('Position',[10 10 500 500]);
    imagesc(S6(crop_y,crop_x,it)); colormap gray; axis image off;
    c = floor(length(crop_x)/2);
    viscircles([c c],preds_S6(it).*c,'Color',c_net,'EnhanceVisibility',false,'LineWidth',6);
    set(gca,'Position',[0 0 1 1]);

    drawnow;
    print('output/temp.png','-dpng');
    writeVideo(v,imread('output/temp.png'));
    close(f);
end
close(v);
