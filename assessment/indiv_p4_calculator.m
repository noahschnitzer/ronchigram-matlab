function [global_min_indiv_p4] = indiv_p4_calculator(aberrations, imdim, simdim)
   deg = pi/180;
   numPx = imdim;
   numAb = length(aberrations(1).n);
   kev = 300;
   lambda = 12.3986./sqrt((2*511.0+kev).*kev)* 10^-10; % E-wavelength in **meter**
   al_max = simdim * 10^-3; %rad %orig 70
   al_vec = (linspace(-al_max,al_max,numPx));
   [alxx,alyy] = meshgrid(al_vec,al_vec);
   % Polar grid
   al_rr = sqrt(alxx.^2 + alyy.^2);
   al_pp = atan2(alyy,alxx);
   %sample
   chi = zeros(size(al_rr));
   global_min_indiv_p4 = inf;

   for kt = 1: numAb
       Cnm = aberrations.mag(1, kt) * aberrations.unit(kt) ;
       m = aberrations.m(kt);
       n = aberrations.n(kt);
       phinm = aberrations.angle(1,kt) * deg;
       ab = 2*pi/lambda*Cnm*cos(m*(al_pp-phinm)).*(al_rr.^(n+1))/(n+1);
       chi0_p4 = abs(ab) > pi/4;
       al_rr_p4 = chi0_p4 .* al_rr;
       al_rr_p4( al_rr_p4 == 0 ) = inf;
       min_indiv_p4 = min(al_rr_p4(:))*1000;
       if min_indiv_p4 < global_min_indiv_p4
          global_min_indiv_p4 = min_indiv_p4;
          %kt
          %figure; plot(al_rr_p4);
       end
       %imagesc(chi0_p4);
       %figure;
       %imagesc(ab);

       %chi = chi +Cnm*cos(m*(al_pp-phinm)).*(al_rr.^(n+1))/(n+1);
   end
   chi0 = 2*pi/lambda * chi; %aberration function



end