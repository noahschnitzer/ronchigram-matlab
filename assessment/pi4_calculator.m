%%% Generates ronchigrams using white noise phase grating with eikonal
%%% approximation to model an amorphous sample.
%%%
%%% By Suk Hyun Sung
%%% revised by noah schnitzer
%%%
%%% Parameters:

function [min_p4] = pi4_calculator(aberrations, imdim, simdim)


    % units
    deg = pi/180;
    numPx = 256;
    if nargin > 3
        numPx = imdim;
    end
    
    numAb = length(aberrations.n);
    
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
    for lim_center = -pi/4:pi/20:pi/4
        
        lb = lim_center - pi/4;
        ub = lim_center + pi/4;
        chi0_p4 = (chi0 < lb) | (chi0 > ub);
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_p4 = min(al_rr_p4(:))*1000;
        if min_p4 > max_p4
            max_p4 = min_p4;
        end
    end
    
    min_p4 = max_p4;

end