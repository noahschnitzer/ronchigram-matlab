%% only C3, scherzer focus
% % ab = aberration_generator(1);
% % ab.mag([1:4, 6:end]) = 0; 
% % ronch = shifted_ronchigram(ab, [0 0], aperture_size, imdim, simdim);
% % figure; imagesc(ronch); axis equal off; colormap gray;
% % title('nodef');
% % scherzer = calculate_scherzer_focus(ab);
% % ab.mag(1) = 0;
% % ab.mag(5) = 3;
% % ab.mag(11) = 2;
% % ab.mag 
% % ronch = shifted_ronchigram(ab, [0 0], aperture_size, imdim, simdim);
% % %figure; imagesc(ronch); axis equal off; colormap gray;

%% simple probe
% ab = aberration_generator(1);
% ab.mag(1:4) = 0; ab.mag(6:10) = 0; ab.mag(12:end)=0;
% %scherzer = calculate_scherzer_focus(ab);
% ab.mag(1) = 0;%-scherzer;
% ab.mag 

%% complex probe
ab = aberration_generator(1);
ab.mag = 1*ab.mag; %.4 for nice

%%
load('infocus_distribution_w_lims.mat');
ab = abs(300);

%% params
num_gen = 5000;
imdim = 1024;
aperture_size = 128;
simdim = 180;
full_range = 0;
threshold = .9;

dfs = [-300:10:300];%[-50:2:50];
zfidx = 26;
%cuts = [.5,.6,.7,.8, .9, .99];
cuts = .8;
%%
cuts = [.6,.8,.9];
dfs = [0:5:100];
[results_strehl,results_pi4,results_indiv_p4] = aberration_series(ab,imdim,simdim,aperture_size,1,dfs,cuts,[]);

%%
lims = lims.* abs(1).unit;
ab_range = -.1:.01:.1;
x_ab_res = zeros(length(ab_range),length(ab.mag));
parfor it = 1:length(ab.mag)
    t_ab = ab;
    t_ab.mag(:) = 0;
   x_ab_res(:,it) = aberration_series(t_ab, imdim, simdim, aperture_size, it, ab_range.*lims(it)/ab.unit(it), cuts,[]); 
   it
end

%%
figure;hold on;
entries = {};
for it = 1:length(ab.mag)
    plot(ab_range,x_ab_res(:,it),'LineWidth',4);
    entries{end+1} = ['C' num2str(abs(1).n(it)) num2str(abs(1).m(it)) ': ' num2str(lims(it))];
end
legend(entries);
%xlim([-1 1]);
xlabel('Aberration Magnitude (C/T)');
ylabel('0.8 Strehl Aperture Size (mrad)');
set(gca,'FontSize',20);
%% finding best focus per .8 strehl
idx8 = 2;
[maxv, idx] = max(results_strehl(:,idx8));
maxvs = find(results_strehl(:,idx8)==maxv);
%% plot figs
%line plot
%%%%figure; plot(dfs,results_strehl,'LineWidth',4);
colordef;
figure;hold on;
plot(dfs,results_strehl(:,1),'LineWidth',4,'Color',c_strehl6);
plot(dfs,results_strehl(:,2),'LineWidth',4,'Color',c_strehl8);
plot(dfs,results_strehl(:,3),'LineWidth',4,'Color',c_strehl9);

hold on
plot(dfs,results_pi4,'--','LineWidth',4,'Color',c_mw_p4);
plot(dfs,results_indiv_p4,'--','LineWidth',4,'Color',c_indiv_p4);
labels = {};
for it = 1:length(cuts)
    labels{it} = [ num2str(cuts(it)) ' Strehl Cutoff'  ];
end
labels{end+1} = 'Moving Window \pi/4';
labels{end+1} = 'Individual \pi/4';

legend(labels);
xlabel('Defocus (Å)');
ylabel('Aperture Semi-angle (mrad)');
set(gca,'FontSize',20);
ylim([0 50]);
xlim([0 100])
%image
% % figure; imagesc(results_strehl);
% % xlabel('strehl cutoff');
% % ylabel('defocus');

%% Calculating res
ab.mag(1) = 0;
ress = [];
aperture_sizes = 5:1:60;
for aperture_size = aperture_sizes
    [~,~,~,probe,~] = shifted_ronchigram(ab,[0 0], aperture_size, imdim, simdim);
    lambda = 1.97e-12;
    px_size = lambda/(2*simdim*1e-3);
    res = resolution_test(probe,'probe');
    %title([num2str(aperture_size) 'mrad, ' num2str(res) 'px']);
    ress(end+1) = res;
end
%% Plotting resolution
figure; plot(aperture_sizes,ress,'LineWidth',4); xlabel('Aperture Size (mrad)'); ylabel('Probe Size (px)'); set(gca, 'FontSize',20);
c = get(gca,'colororder');

% for it = 1:length(cuts)
%    ci = c(mod(it-1,7)+1,:);
%    line([results_strehl(zfidx,it) results_strehl(zfidx,it)],[min(ress) max(ress)],'Color',ci); 
%    text(results_strehl(zfidx,it),max(ress),num2str(cuts(it)));
% end
sum_pi4 = pi4_calculator(ab,imdim,simdim);
indiv_pi4 = indiv_p4_calculator(ab, imdim, simdim);

% it = length(cuts)+1;
% line([sum_pi4 sum_pi4], [min(ress) max(ress)],'Color',c(mod(it-1,7)+1,:));
% text(sum_pi4,max(ress),'sum pi/4');
% it = length(cuts)+2;
% line([indiv_pi4 indiv_pi4], [min(ress) max(ress)],'Color',c(mod(it-1,7)+1,:));
% text(indiv_pi4,max(ress),'indiv pi/4');


%% Plotting Ronch
ab.mag(1) = -28;
ronch = shifted_ronchigram(ab, [0 0], 128,imdim, simdim);
box_dim = 512;
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;

figure; imagesc(ronch(crop_idx, crop_idx)); colormap gray; axis equal off;

center = [256 256];
c = get(gca,'colororder');
vals = results_strehl(zfidx,[3 6]);
                %c = c(mod(it-1,7)+1,:);
a = viscircles(center, vals(1)*imdim/(2*simdim),'Color',c(3,:));
%a.Children(1).Color(4) = 0;
viscircles(center, vals(2)*imdim/(2*simdim),'Color',c(6,:));

hold on; plot(nan,'Color',c(3,:),'LineWidth',4); plot(nan,'Color',c(6,:),'LineWidth',4); legend('0.7 Strehl Ratio', '0.99 Strehl Ratio');

set(gca, 'FontSize',20);

%% below is OLD. see criteria_viz for updated scripts

%% Plotting probe
aperture_size = 35;
ab = abs(300);
z_ab = ab;
z_ab.mag(:) = 0;
[~,~,~,probe,~] = shifted_ronchigram(ab,[0 0], aperture_size, imdim, simdim);
lambda = 1.97e-12;
px_to_ang = 2/((2*simdim/1000)/lambda)*10^10;

xd = (0:imdim/2)*px_to_ang;
yd = radial_average(probe)/(pi*(aperture_size*imdim/(2*simdim)).^2).^2;%normalize_data(radial_average(probe),'total');
figure; area(xd,yd,'LineWidth',2, 'FaceColor','blue','FaceAlpha',.5);
%plot(xd,yd,'LineWidth',4);

hold on;
[~,~,~,probe,~] = shifted_ronchigram(z_ab,[0 0], aperture_size, imdim, simdim);
yd = radial_average(probe)/(pi*(aperture_size*imdim/(2*simdim)).^2).^2;%normalize_data(radial_average(probe),'total');
area(-xd,yd,'LineWidth',2, 'FaceColor','red','FaceAlpha',.5);
%plot(xd,yd,'--red','LineWidth',3);

xlabel('Radius (Å)');
ylabel('Normalized Intensity');
xlim([-2 2]);
set(gca,'FontSize',20);


%% Plotting rayleighs
[~,chi0,minp4,~,~] = shifted_ronchigram(ab,[0 0], aperture_size, imdim, simdim);
range = 0:200;
colordef;
figure; 
hold on;
spherical = get_aberration(ab,imdim,simdim,5);
sp5 = get_aberration(ab,imdim,simdim,11);
twofold = get_aberration(ab,imdim,simdim,2);
plot((range)*2*simdim/imdim,spherical(513,512+range),'Color',[200 200 200]/255,'LineWidth',3);
plot((range)*2*simdim/imdim,sp5(513,512+range),'Color',[200 200 200]/255,'LineWidth',3);
plot((range)*2*simdim/imdim,twofold(513,512+range),'Color',c_indiv_p4,'LineWidth',4);

plot((range)*2*simdim/imdim,-chi0(513,512+range),'Color','k','LineWidth',4);
min_v = min(-chi0(513,512+range));
line([0 range(end)*2*simdim/imdim],[min_v+pi/2 min_v+pi/2],'LineStyle','--','Color',c_mw_p4);
line([0 range(end)*2*simdim/imdim],[min_v min_v],'LineStyle','--','Color',c_mw_p4);

line([0 range(end)*2*simdim/imdim],[pi/4 pi/4],'LineStyle','--','Color',c_indiv_p4,'LineWidth',3);
line([0 range(end)*2*simdim/imdim],[-pi/4 -pi/4],'LineStyle','--','Color',c_indiv_p4);
%patch([0 0 range(end)*2*simdim/imdim range(end)*2*simdim/imdim], [-pi/4 pi/4 pi/4 -pi/4],c_indiv_p4,'FaceAlpha',.3,'EdgeColor','none');

patch([0 0 range(end)*2*simdim/imdim range(end)*2*simdim/imdim], [min_v min_v+pi/2 min_v+pi/2 min_v],c_mw_p4,'FaceAlpha',.5,'EdgeColor','none');

%r=rectangle('Position',[0 min_v range(end)*2*simdim/imdim pi/2],'FaceColor',c_mw_p4);
%alpha(r,0.5)
yticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});

xlabel('Convergence Angle (mrad)');
ylabel('Phase Shift \phi');
set(gca,'FontSize',30);


