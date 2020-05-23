%%% Generates ronchigrams using white noise phase grating with eikonal
%%% approximation to model an amorphous sample.
%%%
%%% By Suk Hyun Sung
%%% revised by noah schnitzer
%%%
%%% Parameters:


function [im, chi0, min_p4, S] = shifted_ronchigram(aberrations, shifts, aperture_size, imdim, simdim)


    % units
    ang = 10^-10;
    nm  = 10^-9;
    um  = 10^-6;
    mm  = 10^-3;
    deg = pi/180;
    numPx = 256;
    if nargin > 3
        numPx = imdim;
    end
    
    numAb = length(aberrations.n);
    
    im = zeros(numPx,numPx);
    kev = 300;
    lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10; % E-wavelength in **meter**
    
    
    % constructing probe
    % alpha (in rad) grid
    al_max = 70 * 10^-3; %rad %orig 70
    if nargin > 4
       al_max = simdim * 10^-3; 
    end
    al_vec = (linspace(-al_max,al_max,numPx));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    % Polar grid
    al_rr = sqrt(alxx.^2 + alyy.^2);
    al_pp = atan2(alyy,alxx);
    % Objective Aperture
    obj_ap_r  = 60 * 10^-3; %rad
    if nargin > 2
        obj_ap_r = aperture_size.*10^-3;
    end
    obj_ap    = al_rr<= obj_ap_r;    
    % creating amorphous specimen by generating random noise kernel,
    % and downsampling
    noise_kernel_size = numPx/8; %256 ->32, 512 ~ 32, 1024 -> 128
    resize_factor = numPx./noise_kernel_size;
    noise_fn = randn(noise_kernel_size,noise_kernel_size);
    noise_fun = imresize(noise_fn,resize_factor);
    charge_e = 1.602e-19;
    mass_e = 9.11e-31;
    c = 3e8;
    % Transmission Operator from Random pi/4 shift specimen +
    % interaction parameter

    interaction_param = 2*pi./(lambda.*kev./charge_e.*1000).*(mass_e*c.^2+kev.*1000)./(2*mass_e*c.^2+kev.*1000);
    interaction_param_0 = 1.7042e-12; %300 kV normalization factor
    trans = exp(-1i*pi/4*noise_fun*interaction_param/interaction_param_0);
    

    % Calculate Probe
    chi = zeros(size(al_rr));
    for kt = 1: numAb
        Cnm = aberrations.mag(1, kt) * aberrations.unit(kt) ;
        m = aberrations.m(kt);
        n = aberrations.n(kt);
        phinm = aberrations.angle(1,kt) * deg;
        chi = chi +Cnm*cos(m*(al_pp-phinm)).*(al_rr.^(n+1))/(n+1);
    end
    chi0 = 2*pi/lambda * chi; %aberration function
    max_p4 = 0;
    max_center = 0;
    for lim_center = -pi/4:pi/20:pi/4
        
        lb = lim_center - pi/4;
        ub = lim_center + pi/4;
        chi0_p4 = (chi0 < lb) | (chi0 > ub);
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_p4 = min(al_rr_p4(:))*1000;
        if min_p4 > max_p4
            max_p4 = min_p4;
            max_center = lim_center;
        end
    end
    
    min_p4 = max_p4;
    
    obj_ap = imtranslate(obj_ap,shifts);
    
    
    expchi0 = exp(-1i * chi0).*obj_ap;          % Apply Objectived Ap
    
    
    psi_p = fft2(expchi0);                     % Probe Wavefunction
    max_psi_p = max(abs(psi_p(:)).^2);
    %figure;
    %imagesc(abs(psi_p).^2);
    
    A = double(sum(obj_ap(:)));
    
    S = max_psi_p /A^2;
    % Transmit Probe
    psi_t = trans.*psi_p;                       % Transmitted Wavefnction

    % Form ronchigram, apply obj astandard unp again to kill outside
    
    ronch = ifft2(psi_t) .*obj_ap;
    ronch = imtranslate(ronch, -shifts);
    ronch = abs(ronch).^2;
    ronch = ronch/max(ronch(:));                % Normalize

    im = ronch;
    
    %probe = abs(ifftshift(psi_p)).^2;
end


function [ronch, chi0, min_p4, probe, S] = shifted_ronchigram_plus(aberrations, shifts, aperture_size, imdim, simdim)
    kev = 300;
    lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10; % E-wavelength in **meter**
    
    al_max = simdim * 10^-3; 
    al_vec = (linspace(-al_max,al_max,imdim));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    al_rr = sqrt(alxx.^2 + alyy.^2);
    al_pp = atan2(alyy,alxx);
    obj_ap_r = aperture_size.*10^-3;
    obj_ap    = al_rr<= obj_ap_r;    
    if ~isequal(shifts,[0 0])
        obj_ap = imtranslate(obj_ap,shifts);
    end

    chi0 = calculate_aberration_function(aberrations,imdim,simdim);
    [psi_p, probe] = calculate_probe(chi0, imdim, simdim, aperture_size,  shifts);
    % ronchigram
    ronch_t = zeros(size(al_pp));
    nsim = 1;
    nnoise = 1;
    noisefact = 8;%16;
    for simnum = 1:nsim
        noise_kernel_size = imdim/noisefact; %256 ->32, 512 ~ 32, 1024 -> 128
        resize_factor = imdim./noise_kernel_size;


        noise_fn = zeros(noise_kernel_size,noise_kernel_size);%randn(noise_kernel_size,noise_kernel_size);
        for it = 1:nnoise
           noise_fn = noise_fn + randn(noise_kernel_size,noise_kernel_size);
        end
        noise_fn = noise_fn./nnoise;
        noise_fn = noise_fn + 1;
        noise_fn = noise_fn./2;
        %figure; histogram(noise_fn(:));
        %for finite sized probe
        %noise_fun = imfilter(imresize(noise_fn,resize_factor),fspecial('gaussian',[10 10],1));
        
        noise_fun = (imresize(noise_fn,resize_factor)); %for normal probe
        
        %noise_fun = perlin_noise(imdim/noisefact,imdim/noisefact);
        %noise_fun = imresize(noise_fun,noisefact);
        %figure; imagesc(noise_fun);

        charge_e = 1.602e-19;
        mass_e = 9.11e-31;
        c = 3e8;
        interaction_param = 2*pi./(lambda.*kev./charge_e.*1000).*(mass_e*c.^2+kev.*1000)./(2*mass_e*c.^2+kev.*1000);
        interaction_param_0 = 1.7042e-12; %300 kV normalization factor
        trans = exp(-1i*pi/4*noise_fun*interaction_param/interaction_param_0);
        psi_t = trans.*psi_p;                       % Transmitted Wavefnction    
        ronch = ifft2(psi_t) .*obj_ap;
        if ~isequal(shifts,[0 0])
            ronch = imtranslate(ronch, -shifts);
        end
        ronch_t = ronch_t + abs(ronch).^2;
    end

    ronch_t = ronch_t - min(ronch_t(:));
    ronch_t = ronch_t./max(ronch_t(:));                % Normalize
    ronch = ronch_t;
    % misc    
    max_p4 = 0;
    max_center = 0;
    for lim_center = -pi/4:pi/80:pi/4
        lb = lim_center - pi/4;
        ub = lim_center + pi/4;
        chi0_p4 = (chi0 < lb) | (chi0 > ub);
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_p4 = min(al_rr_p4(:))*1000;
        if min_p4 > max_p4
            max_p4 = min_p4;
            max_center = lim_center;
        end
    end
    min_p4 = max_p4;
    max_psi_p = max(abs(psi_p(:)).^2);
    A = double(sum(obj_ap(:)));
    S = max_psi_p /A^2;
   
end