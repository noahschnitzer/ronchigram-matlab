%%
clear; close all; clc;
%%
load('data/infocus_distribution_w_lims.mat');
colordef;
chosen_ab = abs(212);
imdim = 1024;
simdim = 180;
S = strehl_calculator(chosen_ab,imdim,simdim,.8,0);
[~,Ss] = strehl_calculator(chosen_ab,imdim,simdim,.8,1);
min_p4 = pi4_calculator(chosen_ab,imdim,simdim);
indiv_p4 = indiv_p4_calculator(chosen_ab,imdim,simdim);
conv_angles = 1:1:105;
ps_simdim_f = 4;
all_probe_sizes = probe_sizer(chosen_ab,imdim*ps_simdim_f,simdim*ps_simdim_f,conv_angles);


%% phase shift
plot_phase_shift(chosen_ab,imdim,simdim,0:200,[2,3,5,11]);

%% probe
domain = [-1 1];

plot_probe_comparison(figure,chosen_ab,imdim,simdim,S,domain,'red',c_strehl8,0.5);

%% Strehl ratio vs CA
figure('Position',[10 10 400 500]); plot(Ss,'LineWidth',4,'Color',c_strehl8); 
xlim([0 80]);
xlabel('Convergence Angle (mrad)');
ylabel('Strehl Ratio');
set(gca,'FontSize',24);
set(gca,'FontName','Helvetica Neue');
%% probe size plot vs CA
figure; 
plot(conv_angles(:),2*px_to_ang(simdim*ps_simdim_f,300)*all_probe_sizes(:),'LineWidth',4,'Color','black');
line([S,S],[0 100],'Color',c_strehl,'LineWidth',4);
line([min_p4,min_p4],[0 100],'Color',c_mw_p4,'LineWidth',4);
line([indiv_p4,indiv_p4],[0 100],'Color',c_indiv_p4,'LineWidth',4);
xlim([3 70]);
ylim([0 3]);
xlabel('Convergence Angle (mrad)');
ylabel(['Probe Size (' char(197) ')']);
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica Neue');
