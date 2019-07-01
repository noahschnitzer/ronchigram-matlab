
function [threshold_apertures, Ss] = strehl_calculator(aberrations, imdim, simdim, threshold, full_range)
    if nargin < 5
        full_range = false;
    end
    mrad_range = simdim:-1:1;
    threshold_apertures = zeros(1,length(aberrations));
    Ss = zeros(length(aberrations),length(mrad_range));

    deg = pi/180;
    numPx = imdim;
    numAb = length(aberrations(1).n);
    kev = 300;
    lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10; % E-wavelength in **meter**
    al_max = simdim * 10^-3; %rad %orig 70
    al_vec = (linspace(-al_max,al_max,numPx));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    % Polar grid
    al_rr = sqrt(alxx.^2 + alyy.^2);
    al_pp = atan2(alyy,alxx);
    %sample

    
    for it = 1:length(threshold_apertures)
        
        if mod(it,1000) == 0
            disp(['Strehls: ' num2str(it) '/' num2str(length(threshold_apertures))]); 
        end

        % Calculate Probe
        chi = zeros(size(al_rr));
        for kt = 1: numAb
            Cnm = aberrations(it).mag(1, kt) * aberrations(it).unit(kt) ;
            m = aberrations(it).m(kt);
            n = aberrations(it).n(kt);
            phinm = aberrations(it).angle(1,kt) * deg;
            chi = chi +Cnm*cos(m*(al_pp-phinm)).*(al_rr.^(n+1))/(n+1);
        end
        chi0 = 2*pi/lambda * chi; %aberration function
        for jt = mrad_range
            obj_ap_r  = jt * 10^-3; %rad

            obj_ap    = al_rr<= obj_ap_r;
            
            expchi0 = exp(-1i * chi0).*obj_ap;          % Apply Objectived Ap
            psi_p = fft2(expchi0);                     % Probe Wavefunction
            max_psi_p = max(abs(psi_p(:)).^2);
            %figure;
            %imagesc(abs(psi_p).^2);

            A = double(sum(obj_ap(:)));

            S = max_psi_p /A^2;
            Ss(it,jt) = S;
            if S > threshold
                threshold_apertures(it) = jt;
                if ~full_range
                    break
                end
            end
        end
    end


end