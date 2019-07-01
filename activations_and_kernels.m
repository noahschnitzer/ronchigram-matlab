%% Load network
load('net-26.mat');

%% Trial ronchigram
%load('through_focal_basis.mat');
new_ab = abs(100);
%new_ab.mag(:) = 0;
ronch = shifted_ronchigram(new_ab,[0 0], 128,1024, 180);
ronch = uint8(ronch*255);
box_dim = 512;
imdim = 1024;
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
ronch = ronch(crop_idx,crop_idx);
%ronch = imread('/home/hlab/aberration_correction/datasets/dataset_23/0/1.png');
pred = predict(trainedNet,ronch)*128;
figure; imagesc(ronch); colormap gray; axis equal off;
% gausian
% % ronch  = fspecial('gaussian', [512 512], 40);
% % ronch = uint8(ronch*255/max(ronch(:)));

% 2 freq random noise
% n1 = rand([16 16]);
% n2 = rand([32 32]);
% n1 = uint8(imresize(n1,512/16)*255);
% n2 = uint8(255*imresize(n2,512/32));
% ronch(1:256,1:512) = n1(1:256,1:512);
% ronch(257:512,1:512) = n2(257:512,1:512);
% figure; imagesc(ronch); colormap gray; axis equal off;

%% Generate activations
activation_layers = [4 7 11 14 18 21 25 28];
for act_idx = [4 28]%activation_layers
    
    act = activations(trainedNet,ronch,act_idx);
    sz = size(act); % insx, insy, nch
    acts = num2cell(act,[1 2]);
    figure;
    imagesc(imtile((acts)));colormap gray; axis equal;
end
%%
kernel_layers = [2 5 9 12 16 19 23 26];

for kernel_idx = 2
    kernel = trainedNet.Layers(kernel_idx).Weights;
    kernels = num2cell(kernel,[1 2]);
    figure;
    imagesc(imtile(kernels)); axis equal off;
    
    sz = size(kernel)
    c_sz = sz(3)*sz(4);
    for n = 0:sqrt(c_sz)-1
        for m = 0:sqrt(c_sz)-1
        line([n*sz(1)+.5 (n+1)*sz(1)+.5], [m*sz(2)+.5 m*sz(2)+.5]);
        line([n*sz(1)+.5 (n)*sz(1)+.5], [m*sz(2)+.5 (m+1)*sz(2)+.5]);
        end
    end
    
    
    %sz = size(kernel); % filtersz filtersz nchin numkernelout
    
end


%act4 = activations(trainedNet, ronch,4);
%sz = size(act4);
%act4 = reshape(act4, [sz(1) sz(2) 1 sz(3)]);
%imshow(imtile(mat2gray(act4),'GridSize',[4 4]));