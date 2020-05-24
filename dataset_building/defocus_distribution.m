%%
load('dataset_24_vals_2.mat');
datasets_path = 'datasets/';
out_dir = 'dataset_26';
num_gen = size(val_table,1);
filepaths = cell(num_gen, 1);
defocii = cell(size(val_table,1),1);
strehl_apertures = zeros(size(val_table,1),1);

aberrations = meta{8};

imdim = 1024;
aperture_size = 128;
simdim = 180;
threshold = .8;

mkdir([datasets_path out_dir]);
for it = 0:1000:99000
  mkdir([datasets_path out_dir '/' num2str(it)]); 
end
%%
parfor it = 1:size(val_table,1)
    it
    tic

    ab = aberrations(it);
    dfs = -30:2:30;
    zdx = 16;
    [results_strehl,~,~] = aberration_series(ab, imdim, simdim, aperture_size, 1,dfs, [.8], []);
    [maxv, idx] = max(results_strehl);
    maxvs = find(results_strehl==maxv);
    strehl_apertures(it) = maxv;
    defocii{it} = results_strehl;
    if ~any(maxvs==zdx)
        ab.mag(1) = dfs(ceil(median(maxvs)));
    else
        ab.mag(1) = 0;
    end
    %disp(['it: ' num2str(it) ', def:' num2str(ab.mag(1))]);

    aberrations(it) = ab;
    
    shifts = round(30.*(rand(1,2)-.5));
    cur_sub = num2str(floor(it/1000)*1000);
    file_str = [datasets_path out_dir '/' cur_sub '/' num2str(it) '.png'];
    filepaths{it} = file_str;
    [im,~,~,~] = shifted_ronchigram(ab, shifts, aperture_size, imdim, simdim);
    box_dim = 512;%round(sqrt((imdim/simdim*aperture_size)^2/2));
    crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
    im = im(crop_idx,crop_idx);
    imwrite(im, file_str); 
    toc
end



val_table = table(filepaths,  strehl_apertures);

val_table.Properties.VariableNames = {'Predictors', 'Strehl'};
meta = {num_gen,imdim,aperture_size,simdim,threshold,datasets_path,out_dir,aberrations};
save([out_dir '_vals.mat'],'val_table','meta');


