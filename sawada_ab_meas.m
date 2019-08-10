%%
clear; close all; clc;
%%
ab = aberration_generator(1);
%%
%ab.mag = zeros(size(ab.mag));
aperture_size = 128;
imdim = 1024;
simdim = 180;
box_dim = 512;
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;

over_ab = ab; over_ab.mag(1) = 50; %over_ab.mag(2) = -2;
under_ab = ab;under_ab.mag(1) = -50;
[over_ronch, ~,~,~,~] = shifted_ronchigram(over_ab,[0 0], aperture_size, imdim, simdim);
[under_ronch, ~,~,~,~] = shifted_ronchigram(under_ab,[0 0], aperture_size, imdim, simdim);
over_ronch = over_ronch(crop_idx,crop_idx);
under_ronch = under_ronch(crop_idx,crop_idx);

%figure;montage((over_sram),'Size',[n n]);
%figure;montage((under_sram),'Size',[n n]);
n = 8;
over_sram = get_sram(over_ronch,n);
under_sram = get_sram(under_ronch,n);
idx = 1:size(over_sram,3);%sub2ind([n,n],4,4);
[over_xs, over_contours,over_fit_ell_exps] = fit_sram(over_sram,idx,.6,.8,0);
[under_xs, under_contours,under_fit_ell_exps] = fit_sram(under_sram,idx,.6,.8,0);

figure; montage(under_sram); 
figure; montage(under_fit_ell_exps); 

figure; montage(over_sram); 
figure; montage(over_fit_ell_exps); 

%%
chi = get_ab_fun(under_xs,over_xs);
figure; imagesc(chi);
%% extracting abs
good_segments = ones(n^2,1);
for it = 1:n^2
    good_segments(it) = ~isequal(under_xs(it,:),[0 0 0]) & ~isequal(over_xs(it,:),[0 0 0]);
end

[abs] = fit_aberrations(over_xs, under_xs,n);

%% fns

function chi = get_ab_fun(under_fits,over_fits,param)
    for it = 1:8
        for jt = 1:8
            fit = under_fits(sub2ind([8 8],it,jt),:);
            chi(it,jt) = fit(2)+fit(3);
            
        end
    end
end

function [sram] = get_sram(ronch, n)
    ac_dim = size(ronch,1) / n;
    for it = 1:n
        for jt = 1:n
            offset_x = (it-1)*ac_dim+1;
            offset_y = (jt-1)*ac_dim+1;
            %crop_ronch(offset_x:offset_x+ac_dim,offset_y:offset_y+ac_dim) = it+jt;
            %boxes(:,:,it,jt) = crop_ronch(
            sram(:,:,sub2ind([n n],it,jt)) = fftshift(ns_norm(ifft2(abs(fft2(ronch(offset_x:offset_x+ac_dim-1,offset_y:offset_y+ac_dim-1))).^2)));
            %acs(:,:,sub2ind([n n],it,jt)) = (xcorr2(crop_ronch(offset_x:offset_x+ac_dim-1,offset_y:offset_y+ac_dim-1)));
        end
    end


end



function [xs, fit_contour,fit_ell_exps] = fit_sram(sram, idxs, clb,cub,doplot)
    xs = [];
    fit_contour = [];
    for it = idxs %1:size(sram,3)
        ell_iso = bwmorph(sram(:,:,it) < cub & sram(:,:,it) > clb, 'thin');
        [r,c] = find(ell_iso);
        sf = 0;
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
            xs(it,:) = lsqcurvefit(ell_eq,x0,xdata,ydata,lb,ub,opt);
        else
            xs(it,:) = [0 0 0];
        end

        %%%
        domain = linspace(-.5,.5,size(sram,1)).*(sf.*size(sram,1));
        [X,Y] = meshgrid(domain,domain);
        indep_mesh = cat(3,Y,X);
        fit_ell = (xs(it,1).*indep_mesh(:,:,1)+xs(it,2).*indep_mesh(:,:,2)).^2+(xs(it,2).*indep_mesh(:,:,1)+xs(it,3).*indep_mesh(:,:,2)).^2;
        fit_ell_exp = exp(-fit_ell);
        fit_ell_exps(:,:,it) = ns_norm(fit_ell_exp);
        fit_contour(:,:,it) = fit_ell_exp > clb & fit_ell_exp < cub;%fit_ell_exp > min(exp(-ell_eq(x(it,:),xdata))) & fit_ell_exp < max(exp(-ell_eq(x(it,:),xdata)));
        if doplot
            figure;
            subplot(2,2,1); imagesc(sram(:,:,it)); axis image; title(num2str(it));
            subplot(2,2,2); imagesc(bwmorph(sram(:,:,it) < cub & sram(:,:,it) > clb, 'thin')); axis image;
            subplot(2,2,3); imagesc(fit_ell_exp); axis image;

            subplot(2,2,4); imagesc(fit_contour(:,:,it)); axis image;
        end
    end

end

function [abs] = fit_aberrations(over_coefs, under_coefs,n)
    %xdata: alpha, theta
    %ydata: A,B,Cs
    %soln: abs
    
    % assembling xdata
    al_max = 180 * 10^-3; 
    al_vec = (linspace(-al_max,al_max,n));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    alpha = sqrt(alxx.^2 + alyy.^2);
    theta = atan2(alyy,alxx);
    xdata = cat(2,alpha(:),theta(:));
    ydata = cat(3,over_coefs,under_coefs);
    
    %ydata = over_coefs;%cat(3,over_coefs,under_coefs);
    guess = zeros(2,12);
    lb = ones(2,12).*-inf;
    ub = ones(2,12).*inf;
    
    tol = 1e-16;
    maxFunEvals = 1e6;
    maxIter = 1e3;
    opt = optimset('Display','Iter','TolFun',tol,'TolX',tol,'MaxFunEvals',maxFunEvals,'MaxIter',maxIter);
    fun = @(param, indep) calc_coeffs(param, indep);
    abs = lsqcurvefit(fun,guess,xdata,ydata,lb,ub,opt);
    %abs_over = lsqcurvefit(fun,guess,xdata,over_coefs,lb,ub,opt);
    %abs_under = lsqcurvefit(fun,guess,xdata,under_coefs,lb,ub,opt);
    %abs = cat(3,abs_over,abs_under);
end


function [abs] = fit_aberrations_v1(over_coefs, under_coefs,n)
    %xdata: alpha, theta
    %ydata: A,B,Cs
    %soln: abs
    
    % assembling xdata
    al_max = 180 * 10^-3; 
    al_vec = (linspace(-al_max,al_max,n));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    alpha = sqrt(alxx.^2 + alyy.^2);
    theta = atan2(alyy,alxx);
    xdata = cat(2,alpha(:),theta(:));
    
    
    %ydata = over_coefs;%cat(3,over_coefs,under_coefs);
    guess = zeros(2,12);
    lb = ones(2,12).*-inf;
    ub = ones(2,12).*inf;
    
    tol = 1e-16;
    maxFunEvals = 1e6;
    maxIter = 1e3;
    opt = optimset('Display','Iter','TolFun',tol,'TolX',tol,'MaxFunEvals',maxFunEvals,'MaxIter',maxIter);
    fun = @(param, indep) calc_coeffs(param, indep);
    abs_over = lsqcurvefit(fun,guess,xdata,over_coefs,lb,ub,opt);
    abs_under = lsqcurvefit(fun,guess,xdata,under_coefs,lb,ub,opt);
    abs = cat(3,abs_over,abs_under);
end

function [coeffs] = calc_coeffs(param,indep)
   mags = param(1,:);
   angles = param(2,:);
   alpha = indep(:,1);
   theta = indep(:,2);
   coeffs = transpose([calc_A(mags,angles,alpha,theta); calc_B(mags,angles,alpha,theta); calc_C(mags,angles,alpha,theta)]);

end

function [out] = calc_A(mags,angles,alpha,theta)
    t = zeros(12,length(alpha));
    %O2
    t(1,:) = mags(1); 
    %A2
    t(2,:) = mags(2).*cos(2.*angles(2)); 
    %P3
    t(3,:) = 4/3.*mags(3).*alpha.*cos(theta-angles(3))+2/3.*mags(3).*alpha.*cos(theta+angles(3));
    %A3
    t(4,:) = 2.*mags(4).*alpha.*cos(theta-3.*angles(4)); 
    %O4
    t(5,:) = 2.*mags(5).*alpha.^2+mags(5).*alpha.^2.*cos(2.*theta);
    %Q4
    t(6,:) = 3/2.*mags(6).*alpha.^2.*cos(angles(6)-2.*theta)+3/2.*mags(6).*alpha.^2;
    %A4
    t(7,:) = 3.*mags(7).*alpha.^2.*cos(2.*theta-4.*angles(7));
    %P5
    t(8,:) = 12/5.*mags(8).*alpha.^3.*cos(theta-angles(8))+6/5.*mags(8).*alpha.^3.*cos(theta+angles(8))+2/5.*mags(8).*alpha.^3.*cos(3.*theta-angles(8));
    %R5
    t(9,:) = 8/5.*mags(9).*alpha.^3.*cos(3.*theta-angles(9))+12/5.*mags(9).*alpha.^3.*cos(theta-angles(9));
    %A5
    t(10,:) = 4.*mags(10).*alpha.^3.*cos(3.*theta-5.*angles(10));
    %O6
    t(11,:) = 3.*mags(11).*alpha.^4+2.*alpha.^4.*mags(11).*cos(2.*theta);
    %A6
    t(12,:) = 5.*mags(12).*alpha.^4.*cos(4.*theta-6.*angles(12));

    out = sum(t,1);
end

function [out] = calc_B(mags,angles,alpha,theta)
    t = zeros(12,length(alpha));
    %O2
    t(1,:) = 0; 
    %A2
    t(2,:) = mags(2).*sin(2.*angles(2)); 
    %P3
    t(3,:) = 2/3.*mags(3).*alpha.*sin(theta+angles(3));
    %A3
    t(4,:) = -2.*mags(4).*alpha.*sin(theta-3.*angles(4)); 
    %O4
    t(5,:) = mags(5).*alpha.^2.*sin(2.*theta);
    %Q4
    t(6,:) = 0;
    %A4
    t(7,:) = -3.*mags(7).*alpha.^2.*sin(2.*theta-4.*angles(7));
    %P5
    t(8,:) = 6/5.*mags(8).*alpha.^3.*sin(theta+angles(8))+2/5.*mags(8).*alpha.^3.*sin(theta-angles(8));
    %R5
    t(9,:) = -12/5.*mags(9).*alpha.^3.*sin(theta-angles(9));
    %A5
    t(10,:) = -4.*mags(10).*alpha.^3.*sin(3.*theta-5.*angles(10));
    %O6
    t(11,:) = 2.*mags(11).*alpha.^4.*sin(2.*theta);
    %A6
    t(12,:) = -5.*mags(12).*alpha.^4.*sin(4.*theta-6.*angles(12));

    out = sum(t,1);
end

function [out] = calc_C(mags,angles,alpha,theta)
    t = zeros(12,length(alpha));
    %O2
    t(1,:) = mags(1); 
    %A2
    t(2,:) = -mags(2).*cos(2.*angles(2)); 
    %P3
    t(3,:) = 4/3.*mags(3).*alpha.*cos(theta-angles(3))-2/3.*mags(3).*alpha.*cos(theta+angles(3));
    %A3
    t(4,:) = -2.*mags(4).*alpha.*cos(theta-3.*angles(4)); 
    %O4
    t(5,:) = 2.*mags(5).*alpha.^2-mags(5).*alpha.^2.*cos(2.*theta);
    %Q4
    t(6,:) = 3/2.*mags(6).*alpha.^2.*cos(angles(6)-2.*theta)-3/2.*mags(6).*alpha.^2;
    %A4
    t(7,:) = -3.*mags(7).*alpha.^2.*cos(2.*theta-4.*angles(7));
    %P5
    t(8,:) = 12/5.*mags(8).*alpha.^3.*cos(theta-angles(8))-6/5.*mags(8).*alpha.^3.*cos(theta+angles(8))-2/5.*mags(8).*alpha.^3.*cos(3.*theta-angles(8));
    %R5
    t(9,:) = 8/5.*mags(9).*alpha.^3.*cos(3.*theta-angles(9))-12/5.*mags(9).*alpha.^3.*cos(theta-angles(9));
    %A5
    t(10,:) = -4.*mags(10).*alpha.^3.*cos(3.*theta-5.*angles(10));
    %O6
    t(11,:) = 3.*mags(11).*alpha.^4-2.*alpha.^4.*mags(11).*cos(2.*theta);
    %A6
    t(12,:) = -5.*mags(12).*alpha.^4.*cos(4.*theta-6.*angles(12));

    out = sum(t,1);
end

