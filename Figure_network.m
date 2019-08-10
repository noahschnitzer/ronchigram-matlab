%%
load('data/dataset_24_vals.mat');

%%
figure;
histogram(val_table.Strehl,[0:max(val_table.Strehl)]);
xlabel('0.8 Strehl Aperture (mrad)');
ylabel('Number of Examples');
xlim([0 max(val_table.Strehl)]);
set(gca,'FontSize',25);
%N = sum(min_probe_size_ap>=min_ap & min_probe_size_ap <=max_ap)
