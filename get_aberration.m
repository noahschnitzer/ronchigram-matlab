function [chi0] = get_aberration(ab,imdim,simdim,idx)

   deg = pi/180;
   numPx = imdim;
   kev = 300;
   lambda = 12.3986./sqrt((2*511.0+kev).*kev)* 10^-10; % E-wavelength in **meter**
   al_max = simdim * 10^-3; %rad %orig 70
   al_vec = (linspace(-al_max,al_max,numPx));
   [alxx,alyy] = meshgrid(al_vec,al_vec);
   % Polar grid
   al_rr = sqrt(alxx.^2 + alyy.^2);
   al_pp = atan2(alyy,alxx);
   %sample
   %chi = zeros(size(al_rr));
   Cnm = ab.mag(1, idx) * ab.unit(idx) ;
   m = ab.m(idx);
   n = ab.n(idx);
   phinm = ab.angle(1,idx) * deg;
   chi0 = 2*pi/lambda*Cnm*cos(m*(al_pp-phinm)).*(al_rr.^(n+1))/(n+1);

end