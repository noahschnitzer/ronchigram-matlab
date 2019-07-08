%% for protoyping... continued in sawada ab meas
clear all; close all; clc;
%% ab
ab = aberration_generator(1);
aperture_size = 128;
imdim = 1024;
simdim = 180;
box_dim = 512;
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
%% ronch
[ronch, ~,~,~,~] = shifted_ronchigram(ab,[0 0], aperture_size, imdim, simdim);
figure; imagesc(ronch(crop_idx,crop_idx));
crop_ronch = ronch(crop_idx,crop_idx);
%% sram
n = 8;
ac_dim = box_dim/n;
sram = [];
%acs = zeros(box_dim/n,box_dim/n,n,n);
for it = 1:n
    for jt = 1:n
        offset_x = (it-1)*ac_dim+1;
        offset_y = (jt-1)*ac_dim+1;
        %crop_ronch(offset_x:offset_x+ac_dim,offset_y:offset_y+ac_dim) = it+jt;
        %boxes(:,:,it,jt) = crop_ronch(
        sram(:,:,sub2ind([n n],it,jt)) = fftshift(ns_norm(ifft2(abs(fft2(crop_ronch(offset_x:offset_x+ac_dim-1,offset_y:offset_y+ac_dim-1))).^2)));
        %acs(:,:,sub2ind([n n],it,jt)) = (xcorr2(crop_ronch(offset_x:offset_x+ac_dim-1,offset_y:offset_y+ac_dim-1)));
    end
end
%Ge=  datacube(acs);
%Ge.vis4D;
figure;montage((sram),'Size',[n n]);
%% fit all
x = [];
fit_contour = [];
for it = sub2ind([8 8],4,5) %1:size(sram,3)
    clb = .5; cub = .7;
    ell_iso = bwmorph(sram(:,:,it) < cub & sram(:,:,it) > clb, 'thin');
    [r,c] = find(ell_iso);
    if ~isempty(r)
        dr = max(r) - min(r); dc = max(c) - min(c);
        r_0 = min(r)+dr/2; c_0 = min(c)+dc/2;
        if dr > dc
            sf = 1/dr;
        else
            sf = 1/dc;
        end
        xdata = cat(2,(r-r_0).*sf,(c-c_0).*sf); %figure; scatter(xdata(:,1),xdata(:,2))
        ydata = .4.*ones(length(r),1);
        x0 = [1.5,-.5,1.5];
        lb = [-inf, -inf, -inf];
        ub = [inf, inf, inf];
        ell_eq = @(param,indep)  (param(1).*indep(:,1)+param(2).*indep(:,2)).^2+(param(2).*indep(:,1)+param(3).*indep(:,2)).^2;
        tol = 1e-16;
        maxFunEvals = 1e6;
        maxIter = 1e3;
        opt = optimset('Display','Iter','TolFun',tol,'TolX',tol,'MaxFunEvals',maxFunEvals,'MaxIter',maxIter);
        x(it,:) = lsqcurvefit(ell_eq,x0,xdata,ydata,lb,ub,opt);
    else
        x(it,:) = [0 0 0];
    end
    
    %%%
    domain = linspace(-.5,.5,size(sram,1));
    [X,Y] = meshgrid(domain,domain);
    indep_mesh = cat(3,Y,X);
    fit_ell = (x(it,1).*indep_mesh(:,:,1)+x(it,2).*indep_mesh(:,:,2)).^2+(x(it,2).*indep_mesh(:,:,1)+x(it,3).*indep_mesh(:,:,2)).^2;
    fit_ell_exp = exp(-fit_ell);


    figure;
    subplot(2,2,1); imagesc(sram(:,:,it)); axis image; title(num2str(idx));
    subplot(2,2,2); imagesc(bwmorph(sram(:,:,it) < cub & sram(:,:,idx) > clb, 'thin')); axis image;
    subplot(2,2,3); imagesc(fit_ell_exp); axis image;

    fit_contour(:,:,it) = fit_ell_exp > clb & fit_ell_exp < cub;%fit_ell_exp > min(exp(-ell_eq(x(it,:),xdata))) & fit_ell_exp < max(exp(-ell_eq(x(it,:),xdata)));
    subplot(2,2,4); imagesc(fit_contour(:,:,it)); axis image;
end
%% checking fit (plotting ellipses)
idx = sub2ind([8 8],4,4);
domain = linspace(-.5,.5,size(sram,1))./(2.*sf);
[X,Y] = meshgrid(domain,domain);
indep_mesh = cat(3,X,Y);
fit_ell = (x(1).*indep_mesh(:,:,1)+x(2).*indep_mesh(:,:,2)).^2+(x(2).*indep_mesh(:,:,1)+x(3).*indep_mesh(:,:,2)).^2;
fit_ell_exp = exp(-fit_ell);


figure;
subplot(2,2,1); imagesc(sram(:,:,idx)); axis image; title(num2str(idx));
subplot(2,2,2); imagesc(bwmorph(sram(:,:,idx) < cub & sram(:,:,idx) > clb, 'thin')); axis image;
subplot(2,2,3); imagesc(fit_ell_exp); axis image;

%% fit_z
clb = .5;
cub = .7; idx = 4;
ell_iso = bwmorph(sram(:,:,idx) < cub & sram(:,:,idx) > clb, 'thin');
[r,c] = find(ell_iso);
dr = max(r) - min(r); dc = max(c) - min(c);
r_0 = min(r)+dr/2; c_0 = min(c)+dc/2;
if dr > dc
    sf = 1/dr;
else
    sf = 1/dc;
end
xdata = cat(2,(r-r_0).*sf,(c-c_0).*sf); %figure; scatter(xdata(:,1),xdata(:,2))
ydata = .4.*ones(length(r),1);
x0 = [1.5,-.5,1.5];
lb = [-inf, -inf, -inf];
ub = [inf, inf, inf];
ell_eq = @(param,indep)  (param(1).*indep(:,1)+param(2).*indep(:,2)).^2+(param(2).*indep(:,1)+param(3).*indep(:,2)).^2;
tol = 1e-16;
maxFunEvals = 1e6;
maxIter = 1e3;
opt = optimset('Display','Iter','TolFun',tol,'TolX',tol,'MaxFunEvals',maxFunEvals,'MaxIter',maxIter);
x = lsqcurvefit(ell_eq,x0,xdata,ydata,lb,ub,opt);

domain = linspace(-.5,.5,size(sram,1))./(2.*sf);
[X,Y] = meshgrid(domain,domain);
indep_mesh = cat(3,X,Y);
fit_ell = (x(1).*indep_mesh(:,:,1)+x(2).*indep_mesh(:,:,2)).^2+(x(2).*indep_mesh(:,:,1)+x(3).*indep_mesh(:,:,2)).^2;
fit_ell_exp = exp(-fit_ell);
figure; 
subplot(2,2,1); imagesc(sram(:,:,idx)); axis image; title(num2str(idx));
subplot(2,2,2); imagesc(ell_iso); axis image;
subplot(2,2,3); imagesc(fit_ell_exp); axis image;
subplot(2,2,4); imagesc(fit_ell_exp > min(exp(-ell_eq(x,xdata))) & fit_ell_exp < max(exp(-ell_eq(x,xdata))));

%% fit_y
clb = .5;
cub = .7;
idx = 44;
ell = bwmorph(sram(:,:,idx) < cub & sram(:,:,idx) > clb, 'thin');
[r,c]=find(ell);
xdata = ns_norm(cat(2,r,c))-.5; %%%% this is an issue! not properly placing on grid
% should calculate center, and scale axis equally
% then feed scaling factors back in for display
ydata = .4.*ones(length(r),1);
x0 = [1.5,-.5,1.5];
lb = [-inf, -inf, -inf];
ub = [inf, inf, inf];
ell_eq = @(param,indep)  (param(1).*indep(:,1)+param(2).*indep(:,2)).^2+(param(2).*indep(:,1)+param(3).*indep(:,2)).^2;
%fun = @(param,indep) [(diff(ell_eq(param,indep))); 0];

tol = 1e-16;
maxFunEvals = 1e6;
maxIter = 1e3;
opt = optimset('Display','Iter','TolFun',tol,'TolX',tol,'MaxFunEvals',maxFunEvals,'MaxIter',maxIter);
x = lsqcurvefit(ell_eq,x0,xdata,ydata,lb,ub,opt);

%figure;scatter3(xdata(:,1),xdata(:,2),ell_eq(x,xdata));
figure; scatter3(xdata(:,1),xdata(:,2),exp(-ell_eq(x,xdata)));
domain = linspace(-.5,.5,size(sram,1));
[X,Y] = meshgrid(domain,domain);
indep_mesh = cat(3,X,Y);
fit_ell = (x(1).*indep_mesh(:,:,1)+x(2).*indep_mesh(:,:,2)).^2+(x(2).*indep_mesh(:,:,1)+x(3).*indep_mesh(:,:,2)).^2;
fit_ell_exp = exp(-fit_ell);
figure; 
subplot(2,2,1); imagesc(sram(:,:,idx)); axis image; title(num2str(idx));
subplot(2,2,2); imagesc(ell); axis image;
subplot(2,2,3); imagesc(fit_ell_exp); axis image;
subplot(2,2,4); imagesc(fit_ell_exp > min(exp(-ell_eq(x,xdata))) & fit_ell_exp < max(exp(-ell_eq(x,xdata))));

%%
x = [.3 -.2 .3];
fit_ell = (x(1).*indep_mesh(:,:,1)+x(2).*indep_mesh(:,:,2)).^2+(x(2).*indep_mesh(:,:,1)+x(3).*indep_mesh(:,:,2)).^2;
figure; 
subplot(1,3,1); imagesc(sram(:,:,idx)); axis image; title(num2str(idx));
subplot(1,3,2); imagesc(ell); axis image;
subplot(1,3,3); imagesc(exp(-fit_ell)); axis image;

%% fit_x
domain = linspace(-1,1,size(sram,1));
[X,Y] = meshgrid(domain,domain);
xdata = cat(3,X,Y);
clb = .7;
cub = .8;
idx = 5;
ydata = double(bwmorph(sram(:,:,idx) < cub & sram(:,:,idx) > clb, 'thin')); %skel or thin..
x0 = [2, 1, 1];
lb = [-inf, -inf, -inf];
ub = [inf, inf, inf];

fun = @(param,indep) double(.1>abs(.75- (param(1).*indep(:,:,1)+param(2).*indep(:,:,2)).^2+(param(2).*indep(:,:,1)+param(3).*indep(:,:,2)).^2));
tol = 1e-10;
maxFunEvals = 1e6;
maxIter = 1e3;
opt = optimset('Display','Iter','TolFun',tol,'TolX',tol,'MaxFunEvals',maxFunEvals,'MaxIter',maxIter);
x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,opt);

figure; imagesc(ydata);
figure; imagesc(fun(x,xdata));