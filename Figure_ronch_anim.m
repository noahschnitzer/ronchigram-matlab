%% set up
load('data/infocus_distribution_w_lims.mat');
colordef;
chosen_ab = abs(212);
anim_abs(1) = chosen_ab;
aperture_size = 128;
simdim = 180;
imdim = 1024;

%% ab series
series_incl = zeros(14,14);
series_incl(1,1:end) = 1;
series_incl(5,2:end) = 1;
series_incl(2,3:end) = 1;
series_incl(3,4:end) = 1;
series_incl(4,5:end) = 1;
for it = 6:14
   series_incl(it,it:end) = 1; 
end
for it = 1:14
   anim_abs(it) = chosen_ab;
   anim_abs(it).mag(series_incl(:,it)==0) = 0;
end
%%
v = VideoWriter('output/abseries_out.avi');
v.FrameRate = 1;
open(v);
fig = figure;

for ab = anim_abs
   ronch_im = shifted_ronchigram(ab,[0 0],aperture_size,imdim,simdim); 
   imagesc(ronch_im);axis image; colormap gray;
   drawnow;
   writeVideo(v,ronch_im);
end
close(v);