%%
num_gen = 99e3;
imdim = 1024;
aperture_size = 128;
simdim = 180;
full_range = 0;
threshold = .8;%threshold = .9;
datasets_path = 'datasets/';
out_dir = 'dataset_99';

%%
%[aberrations] = aberration_generator(num_gen);
[aberrations] = distribution_generator(num_gen,[1 90],1024,180);
save('wip_ab_dist.mat','aberrations');
disp('Distribution Built');
%%
[strehl_apertures,~] = strehl_calculator(aberrations, 256, simdim, threshold, full_range);
strehl_apertures = strehl_apertures / aperture_size;
save('wip_strehl.mat','strehl_apertures');
disp('Strehl Apertures Calculated');

%%
load('wip_strehl.mat')
load('wip_ab_dist.mat')
strehl_apertures = strehl_apertures(1:num_gen);
filepaths = cell(num_gen, 1);

mkdir([datasets_path out_dir]);
mkdir([datasets_path out_dir '/0']);
cur_sub = '0';
for it = 1:num_gen
    if mod(it,1000) == 0
        mkdir([datasets_path out_dir '/' num2str(it)]);
        cur_sub = num2str(it);
        disp([num2str(it) '/' num2str(num_gen)]); 
    end

    
    shifts = round(30.*(rand(1,2)-.5));
    file_str = [datasets_path out_dir '/' cur_sub '/' num2str(it) '.png'];
    filepaths{it} = file_str;
    [im, ~, ~, ~] = shifted_ronchigram(aberrations(it), shifts, aperture_size, imdim, simdim);
    box_dim = 512;%round(sqrt((imdim/simdim*aperture_size)^2/2));
    crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
    im = im(crop_idx,crop_idx);
    im0 = imresize(im,[227 227]);
    im = cat(3,im0,im0);
    im = cat(3,im,im0);
    imwrite(im, file_str); 
% %     try 
% %         imread(file_str);
% %     catch
% %         display(['Unable to open' num2str(it)]);
% %         imwrite(im,file_str);
% %     end
end


val_table = table(filepaths,  strehl_apertures');

val_table.Properties.VariableNames = {'Predictors', 'Strehl'};
meta = {num_gen,imdim,aperture_size,simdim,threshold,datasets_path,out_dir,aberrations};
save([out_dir '_vals.mat'],'val_table','meta');


