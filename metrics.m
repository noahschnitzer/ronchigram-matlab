% Metrics:
%   Strehl cutoff of 0.8
%   Individual pi/4
%   Moving window pi/4
%   Human intuition
%   Minimum probe size *(the ground truth)


%% params
n_ab = 1000;
imdim = 1024; simdim=180; aperture_size = 128;

%% Building set of aberrations
%abs  = aberration_generator(n_ab);

%% loading set of aberrations
load('ab_set.mat');
n_ab = length(abs);

%% Identifying apertures for each metric
%% Strehls
threshold = 0.8;
tic
strehl_ap = par_strehl_calculator(abs,imdim,simdim,threshold);
toc
disp('Strehl Calculated');

%% individual pi/4
indiv_p4_ap = zeros(1,n_ab);
for it = 1:n_ab
   indiv_p4_ap(it) = indiv_p4_calculator(abs(it),imdim,simdim); 
end
disp('Indiv \pi/4 Calculated');
%% moving window pi/4
mw_p4_ap = zeros(1,n_ab);
for it = 1:n_ab
   mw_p4_ap(it) = pi4_calculator(abs(it),imdim,simdim); 
end
disp('MW \pi/4 Calculated');
%% human data
n_hum = 2;
load('noah_game_data.mat');
human_ap(1,:) = scaled_res;
load('ssh_game_data.mat');
human_ap(2,:) = scaled_res;
human_ap(human_ap > 101) = 101;
human_ap(human_ap < 1) = 1;
%% network preds
load('net-26.mat');
net_ap = zeros(1,n_ab);
for it = 1:n_ab
   im = shifted_ronchigram(abs(it),[0 0],aperture_size,imdim,simdim); 
   box_dim = 512;%round(sqrt((imdim/simdim*aperture_size)^2/2));
   crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
   im = im(crop_idx,crop_idx);
   imwrite(im, 'temp.png');
   im = imread('temp.png');
   net_ap(it) = predict(trainedNet,im);
   it
end
net_ap = net_ap*aperture_size;
net_ap(net_ap > 101) = 101;
%% probe sizes
tic
probe_sizes = probe_sizer(abs,imdim*2,simdim,[1:1:105]); %useful later!
toc

min_probe_sizes = min(probe_sizes,[],2);
min_probe_size_aps = cell(1,n_ab);
min_probe_size_ap = zeros(1,n_ab);
for it = 1:n_ab
    min_probe_size_aps{it} = find(min_probe_sizes(it)==probe_sizes(it,:));
    min_probe_size_ap(it) = median(min_probe_size_aps{it});
end
disp('Probe Sizes Calculated');
%% Calculating probe size errors for each metric vs minimum probe size aperture
strehl_ap_err = zeros(1,n_ab);
indiv_p4_ap_err = zeros(1,n_ab);
mw_p4_ap_err = zeros(1,n_ab);
net_ap_err = zeros(1,n_ab);
human_ap_err = zeros(1,n_ab);

%errs = zeros(5,n_ab);
for it = 1:n_ab
    strehl_probe_size = probe_sizes(it,strehl_ap(it));
    strehl_ap_err(it) = strehl_probe_size - min_probe_sizes(it);
    indiv_p4_probe_size = probe_sizes(it,round(indiv_p4_ap(it)));
    indiv_p4_ap_err(it) = indiv_p4_probe_size - min_probe_sizes(it);
    mw_p4_probe_size = probe_sizes(it,round(mw_p4_ap(it)));
    mw_p4_ap_err(it) = mw_p4_probe_size - min_probe_sizes(it);
    net_probe_size = probe_sizes(it,round(net_ap(it)));
    net_ap_err(it) = net_probe_size - min_probe_sizes(it);
    for jt = 1:n_hum
        human_probe_size = probe_sizes(it, round(human_ap(jt,it)));
        human_ap_err(it) = human_ap_err(it) + human_probe_size - min_probe_sizes(it);
    end
end

human_ap_err = human_ap_err./n_hum;
%% Averaging down columns (ap sizes vs err)
[reduced_min_probe_size_ap, reduced_strehl_ap_err] = column_means(min_probe_size_ap,strehl_ap_err);
[~,reduced_indiv_p4_ap_err] = column_means(min_probe_size_ap,indiv_p4_ap_err); 
[~, reduced_mw_p4_ap_err] = column_means(min_probe_size_ap,mw_p4_ap_err);
[~, reduced_net_ap_err] = column_means(min_probe_size_ap, net_ap_err);
[~, reduced_human_ap_err] = column_means(min_probe_size_ap, human_ap_err);
%% Averaging down columns (ap sizes vs ap sizes)
[reduced_min_probe_size_ap, reduced_strehl_ap_err] = column_means(min_probe_size_ap,strehl_ap);
[~,reduced_indiv_p4_ap_err] = column_means(min_probe_size_ap,indiv_p4_ap); 
[~, reduced_mw_p4_ap_err] = column_means(min_probe_size_ap,mw_p4_ap);
[~, reduced_net_ap_err] = column_means(min_probe_size_ap, net_ap);
[~, reduced_human_ap_err] = column_means(min_probe_size_ap, human_ap);
%% Probe size unit
lambda = 1.97e-12;
% % % s = imdim/(2*simdim/1000); %px/mrad
% % % f = 2*lambda/s; %mrad m / px, ->diameter not radius
% % % f*6*1e10 %Angstrom
px_to_ang = 2/((2*simdim/1000)/lambda)*10^10;

%% Defining Colors
c = get(gca,'ColorOrder');
c_strehl = c(1,:);
c_indiv_p4 = c(2,:);
c_mw_p4 = c(3,:);
c_net = c(4,:);
c_human = c(5,:);
c_min = 'red';


%% Plotting error
figure; hold on;
min_ap = 10; % min 3
max_ap = 70; % max 95
%px_to_ang = 1;
domain = find(reduced_min_probe_size_ap == min_ap):find(reduced_min_probe_size_ap == max_ap); %3,95
plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_strehl_ap_err(domain),'LineWidth',4,'Color',c_strehl);
plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_indiv_p4_ap_err(domain),'LineWidth',4,'Color',c_indiv_p4);
plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_mw_p4_ap_err(domain),'LineWidth',4,'Color',c_mw_p4);
%plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_net_ap_err(domain),'LineWidth',4,'Color',c_net);
%plot(reduced_min_probe_size_ap(domain),px_to_ang*reduced_human_ap_err(domain),'LineWidth',4,'Color',c_human);
plot(nan,'Color',c_min,'LineWidth',4);
%plot(1:70,1:70,'Color','Black','LineWidth',4);

%legend('0.8 Strehl Cutoff','Individual \pi/4','Moving Window \pi/4','Network Prediction','Human Guess','Minimum Probe Size');
xlabel('Best Aperture (mrad)');
ylabel('Mean Probe Size Error (Å)');
ylim([0 .6]);
set(gca,'FontSize',20);

%% Plotting distribution
figure;
histogram(min_probe_size_ap(min_probe_size_ap>=min_ap & min_probe_size_ap <=max_ap),[min_ap:max_ap] );
xlabel('Best Aperture (mrad)');
ylabel('Number of Examples');
xlim([10 70]);
set(gca,'FontSize',25);
N = sum(min_probe_size_ap>=min_ap & min_probe_size_ap <=max_ap)
%% Individual Aberration
%% Plotting Ronchigram
idx = 100;%100,176,
ronch = shifted_ronchigram(abs(idx),[0 0], aperture_size,imdim,simdim);
bclim = [0.02 .6];
figure; imagesc(ronch,bclim); colormap(gray(2^8)); 
viscircles([imdim/2 imdim/2],min_probe_size_ap(idx)*imdim/(2*simdim),'Color',c_min,'EnhanceVisibility',false,'LineWidth',3);
viscircles([imdim/2 imdim/2],human_ap1(idx)*imdim/(2*simdim),'Color',c_human,'EnhanceVisibility',false,'LineWidth',3);
viscircles([imdim/2 imdim/2],strehl_ap(idx)*imdim/(2*simdim),'Color',c_strehl,'EnhanceVisibility',false,'LineWidth',3);
viscircles([imdim/2 imdim/2],net_ap(idx)*imdim/(2*simdim),'Color',c_net,'EnhanceVisibility',false,'LineWidth',3);
axis equal off;
hold on;
ws = 150;
c = 512;
xlim([c-ws c+ws]);
ylim([c-ws c+ws]);
%% Plotting probe size vs conv angle
figure; 
plot(px_to_ang*probe_sizes(idx,:),'LineWidth',4,'Color','black');
line([strehl_ap(idx) strehl_ap(idx)],[px_to_ang*probe_sizes(idx,strehl_ap(idx)) px_to_ang*max(probe_sizes(idx,:))],'Color',c_strehl,'LineWidth',4);
line([human_ap1(idx) human_ap1(idx)],[px_to_ang*probe_sizes(idx,round(human_ap1(idx))) px_to_ang*max(probe_sizes(idx,:))],'Color',c_human,'LineWidth',4);
line([min_probe_size_ap(idx) min_probe_size_ap(idx)],[px_to_ang*probe_sizes(idx,round(min_probe_size_ap(idx))) px_to_ang*max(probe_sizes(idx,:))],'Color',c_min,'LineWidth',4);
line([net_ap(idx) net_ap(idx)],[px_to_ang*probe_sizes(idx,round(net_ap(idx))) px_to_ang*max(probe_sizes(idx,:))],'Color',c_net,'LineWidth',4);
ylim([0 7]);
xlabel('Convergence Angle (mrad)');
ylabel('Probe Size (Å)');
set(gca,'FontSize',20);
set(gca,'color','none')

%% metric aggregate errors
human_aggregate_err = mean(human_ap_err(min_probe_size_ap>10 & min_probe_size_ap < 70))*px_to_ang*100
net_aggregate_err = mean(net_ap_err(min_probe_size_ap>10 & min_probe_size_ap < 70))*px_to_ang*100
strehl_aggregate_err = mean(strehl_ap_err(min_probe_size_ap>10 & min_probe_size_ap < 70))*px_to_ang*100
mw_p4_aggregate_err = mean(mw_p4_ap_err(min_probe_size_ap>10 & min_probe_size_ap < 70))*px_to_ang*100
indiv_p4_aggregate_err = mean(indiv_p4_ap_err(min_probe_size_ap>10 & min_probe_size_ap < 70))*px_to_ang*100

true_avg_size = mean(min_probe_sizes)*px_to_ang*100;
pct_improved = (net_aggregate_err-human_aggregate_err)/true_avg_size*100

