%%
fdir = './';
fname = 'foc_stack.tif';
len = 70;
dim = 2048;
A = zeros(dim,dim,len);
for it = 1:len
    A(:,:,it) = imread([fdir fname],it);
end

%%
load('net-26.mat');
defocuses = -350:10:349;
posx = 360;
posy = 300;
wh = 1300;
crop_x = posx:posx+wh; 
crop_y = posy:posy+wh;
% cropped_A = A(crop_y,crop_x,:);
% %figure;imagesc(cropped_A(:,:,1));
% %figure;imagesc(masked(:,:,1));
% masked = A;
% masked(crop_y,crop_x) = 0;
preds = [];
for it = 1:len
   im =  A(crop_y,crop_x,it);
   im = imresize(im,[512 512]);
   im = im - min(im(:));
   im = uint8(im.*255/max(im(:)));
   preds(it) = predict(trainedNet,im);
end
figure; plot(preds);