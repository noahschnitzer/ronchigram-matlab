%% from strehl_behavior...
load('data/infocus_distribution_w_lims.mat');
colordef;
chosen_ab = abs(212);
scherz_spher = chosen_ab;%abs(300);
scherz_spher.mag = zeros(1,14);
scherz_spher.mag(5) = 1e3;
scherz_spher.mag(1) = 1*calculate_scherzer_focus(scherz_spher);

imdim = 1024;
simdim = 180;
S = strehl_calculator(chosen_ab,imdim,simdim,.8,0);
min_p4 = pi4_calculator(chosen_ab,imdim,simdim);
indiv_p4 = indiv_p4_calculator(chosen_ab,imdim,simdim);

%% spherical only...
spher_domain = [-2 2];
plot_phase_shift(scherz_spher,imdim,simdim/2,0:60,[1,5]);
plot_probe_comparison(figure,scherz_spher,imdim,simdim,pi4_calculator(scherz_spher,imdim,simdim),spher_domain,'red','blue',.5);
%% complex
plot_phase_shift(chosen_ab,imdim,simdim,0:200,[2,3,5,11]);
complex_domain = [-1 1];

plot_probe_comparison(figure,chosen_ab,imdim,simdim,min_p4,complex_domain,'red',c_mw_p4,0.5);
plot_probe_comparison(figure,chosen_ab,imdim,simdim,indiv_p4,complex_domain,'red',c_indiv_p4,0.5);
plot_probe_comparison(figure,chosen_ab,imdim,simdim,S,complex_domain,'red',c_strehl8,0.5);

%% specific probe sizes
ps_s = probe_sizer(chosen_ab,imdim*4,simdim*4,[min_p4,indiv_p4,S]);
px_to_ang(simdim*4,300)*ps_s*2

%% all probe sizes
conv_angles = 1:1:105;
ps_simdim_f = 4;
all_probe_sizes = probe_sizer(chosen_ab,imdim*ps_simdim_f,simdim*ps_simdim_f,conv_angles);
%% probe size vs CA plot
figure; 
plot(conv_angles(:),2*px_to_ang(simdim*ps_simdim_f,300)*all_probe_sizes(:),'LineWidth',4,'Color','black');
%line([S S],[px_to_ang*probe_sizes(idx,strehl_ap(idx)) px_to_ang*max(probe_sizes(idx,:))],'Color',c_strehl,'LineWidth',4);
line([S,S],[0 100],'Color',c_strehl,'LineWidth',4);
line([min_p4,min_p4],[0 100],'Color',c_mw_p4,'LineWidth',4);
line([indiv_p4,indiv_p4],[0 100],'Color',c_indiv_p4,'LineWidth',4);

xlim([3 70]);
ylim([0 3]);
xlabel('Convergence Angle (mrad)');
ylabel(['Probe Size (' char(197) ')']);
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica Neue');
%ylim([0 2*px_to_ang(simdim*ps_simdim_f,300)*all_probe_sizes(plot_domain(1))]);
%% Probe sizes across apertures plots from metrics
load('data/metrics_state.mat','reduced_min_probe_size_ap','reduced_strehl_ap_err','reduced_indiv_p4_ap_err','reduced_mw_p4_ap_err');
% reduced probe sizes are radius'.. hence 2* factor for probe size err
figure; hold on;
min_ap = 5; % min 3
max_ap = 70; % max 95
%px_to_ang = 1;
domain = find(reduced_min_probe_size_ap == min_ap):find(reduced_min_probe_size_ap == max_ap); %3,95
plot(reduced_min_probe_size_ap(domain),2*px_to_ang(simdim,300)*reduced_strehl_ap_err(domain),'LineWidth',4,'Color',c_strehl);
plot(reduced_min_probe_size_ap(domain),2*px_to_ang(simdim,300)*reduced_indiv_p4_ap_err(domain),'LineWidth',4,'Color',c_indiv_p4);
plot(reduced_min_probe_size_ap(domain),2*px_to_ang(simdim,300)*reduced_mw_p4_ap_err(domain),'LineWidth',4,'Color',c_mw_p4);
%plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_net_ap_err(domain),'LineWidth',4,'Color',c_net);
%plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_human_ap_err(domain),'LineWidth',4,'Color',c_human);
plot(nan,'Color',c_min,'LineWidth',4);
%plot(1:70,1:70,'Color','Black','LineWidth',4);

%legend('0.8 Strehl Cutoff','Individual \pi/4','Moving Window \pi/4','Network Prediction','Human Guess','Minimum Probe Size');
xlabel('Convergence Angle (mrad)');
ylabel(['Probe Size Error (' char(197) ')']);
xlim([min_ap max_ap]);
ylim([0 .4]);
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica Neue');

%% aberration fn
chi_map = calculate_aberration_function(chosen_ab,imdim,simdim);%shifted_ronchigram(chosen_ab,[],128,);
%figure; imagesc((chi_map)); colormap(labinterp);
inscribed = 256:256+512;
%figure; imagesc(phase2color(chi_map(inscribed,inscribed))); axis image off;
imwrite(phase2color(chi_map(inscribed,inscribed)),'output/phase_map.jpg');

%% 3d probes
[~, probe_lap] = calculate_probe(chi_map, imdim, simdim, 128,  [0 0]);
[~, probe_sap] = calculate_probe(chi_map, imdim, simdim, 31,  [0 0]);
extent = imdim/2-15:imdim/2+15;
figure; surf(probe_lap(extent,extent),'FaceColor','interp'); axis off;
colormap plasma;
figure; surf(probe_sap(extent,extent),'FaceColor','interp'); axis off;
colormap plasma;

%% %% %% %% %% %%

%% Plotting probes
%plot_probe_comparison(chosen_ab,imdim,simdim,S);
%plot_probe_comparison(chosen_ab,imdim,simdim,min_p4);
plot_probe_comparison(scherz_spher,imdim,simdim,pi4_calculator(scherz_spher,imdim,simdim));
%% Plotting rayleighs
%generic_ab = abs(200);
plot_phase_shift(chosen_ab,imdim,simdim,0:200,[2,3,5,11]);
plot_phase_shift(scherz_spher,imdim,simdim/2,0:60,[1,5]);

%% 3D probe plots
extent = 15;
[X,Y] = meshgrid(-extent:extent,-extent:extent);
X = X.*px_to_ang;
Y = Y.*px_to_ang;
%generic_ab = abs(300);
zero_ab = chosen_ab;
zero_ab.mag = zeros(size(zero_ab.mag));

imdim = 1024;
aperture_size = 128;
simdim = 180;
%magma, inferno, plasma, vega10 20 b c viridis labinterp
[~,~,~,probe,~] = shifted_ronchigram(zero_ab,[0 0], aperture_size, imdim, simdim);
figure; 
subplot(221);
surf(X,Y,ns_norm(probe(512-extent:512+extent,512-extent:512+extent)),'FaceColor','interp');
colormap plasma;
title('lap 0ab');

[~,~,~,probe,~] = shifted_ronchigram(chosen_ab,[0 0], aperture_size, imdim, simdim);
subplot(222);
surf(X,Y,ns_norm(probe(512-extent:512+extent,512-extent:512+extent)),'FaceColor','interp');
colormap plasma;
title('lap ab');

s_ap = 31;
[~,~,~,probe,~] = shifted_ronchigram(zero_ab,[0 0], s_ap, imdim, simdim);
subplot(223); 
surf(X,Y,ns_norm(probe(512-extent:512+extent,512-extent:512+extent)),'FaceColor','interp');
colormap plasma;
title('sap 0ab');

[~,~,~,probe,~] = shifted_ronchigram(chosen_ab,[0 0], s_ap, imdim, simdim);
subplot(224);
surf(X,Y,ns_norm(probe(512-extent:512+extent,512-extent:512+extent)),'FaceColor','interp');
colormap plasma;
title('sap ab');
