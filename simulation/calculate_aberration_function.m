function [chi0] = calculate_aberration_function(abs, imdim,simdim)
    
    kev = 300;
    lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10; % E-wavelength in **meter**
    al_max = simdim * 10^-3; 
    al_vec = (linspace(-al_max,al_max,imdim));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    al_rr = sqrt(alxx.^2 + alyy.^2);
    al_pp = atan2(alyy,alxx);
    deg = pi/180;
    numAb = length(abs(1).n);
    numFns = length(abs);
    chi = zeros(imdim,imdim,numFns);
    for it = 1:numFns
        aberrations = abs(it);

        for kt = 1: numAb
            Cnm = aberrations.mag(1, kt) * aberrations.unit(kt) ;
            m = aberrations.m(kt);
            n = aberrations.n(kt);
            phinm = aberrations.angle(1,kt) * deg;
            chi(:,:,it) = chi(:,:,it) +Cnm*cos(m*(al_pp-phinm)).*(al_rr.^(n+1))/(n+1);
        end
    end
    chi0 = 2*pi/lambda * chi; %aberration function
end