%%
load('data/metrics_state.mat');
colordef;
ap_size = 128; imdim = 1024; simdim = 180;
%% ronchigram, human ap, network ap, probe for each
idx = 100;
n_choice = net_ap(idx);
h_choice = human_ap(1,idx);
%figure; plot(probe_sizes(idx,:));
ab = abs(idx);
ronch = shifted_ronchigram(ab,[0 0],ap_size,imdim,simdim);
box_dim = 512;%round(sqrt((imdim/simdim*aperture_size)^2/2));
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
ronch = ronch(crop_idx,crop_idx);
figure('Position',[10 10 500 500]); imagesc(ronch); colormap gray; axis image off;
set(gca,'Position',[0 0 1 1]);
c = 256;
viscircles([c c],(n_choice).*imdim/(2*simdim),'Color',c_net,'EnhanceVisibility',false,'LineWidth',6);
viscircles([c c],(h_choice).*imdim/(2*simdim),'Color',c_human,'EnhanceVisibility',false,'LineWidth',6);

[~,n_probe] = calculate_probe(calculate_aberration_function(ab,imdim,simdim),imdim,simdim,n_choice,[0 0]);
[~,h_probe] = calculate_probe(calculate_aberration_function(ab,imdim,simdim),imdim,simdim,h_choice,[0 0]);
delta = 30;extent = imdim/2-delta:imdim/2+delta;
figure; surf(n_probe(extent,extent),'FaceColor','interp'); axis off;
colormap plasma;

figure; surf(h_probe(extent,extent),'FaceColor','interp'); axis off;
colormap plasma;
