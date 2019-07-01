figure; imagesc(shifted_ronchigram_o(abs(100),[0 0],128, 1024, 180));colormap gray;
title('original');
figure; imagesc(shifted_ronchigram(abs(100),[0 0],128, 1024, 180));colormap gray;
title('new');

%%
chi0s = calculate_aberration_function(abs(1:5),1024,180);
probes = calculate_probe(chi0s, 1024, 180, 128, [0 0]);

%%
hold3 = probe_sizer(abs(1:100:1000),1024,180,[1:128]);
strehl_calculator(abs(1:100:1000),1024,180,.8,0);

%%

ab = aberration_generator(1);
ab.mag = ab.mag*10;
%%
ronch = zeros(256);
for it = 1:1
   ronch_0 = shifted_ronchigram(ab,[0 0],50,256,70);
   ronch = ronch+ronch_0;
end
%%
figure; imagesc(ronch_0); title('single'); colormap gray;
figure; imagesc(ronch); title('multi'); colormap gray;



%%
idx = 100;%100,176,
im = shifted_ronchigram(abs(idx),[0 0], aperture_size,imdim,simdim);
box_dim = 512;%round(sqrt((imdim/simdim*aperture_size)^2/2));
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
im = im(crop_idx,crop_idx);
figure; imagesc(im);
axis image off; colormap gray;
