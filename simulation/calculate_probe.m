function [psi_p, probe] = calculate_probe(chi0,imdim, simdim, aperture_size, shifts)
    al_max = simdim * 10^-3; 
    al_vec = (linspace(-al_max,al_max,imdim));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    al_rr = sqrt(alxx.^2 + alyy.^2);


    obj_ap_r = aperture_size.*10^-3;
    obj_ap    = al_rr<= obj_ap_r;    
    if ~isequal(shifts,[0 0])
        obj_ap = imtranslate(obj_ap,shifts);
    end
    expchi0 = exp(-1i * chi0).*obj_ap;          % Apply Objectived Ap
    psi_p = fft2(expchi0);                     % Probe Wavefunction
    probe = abs(ifftshift(psi_p)).^2;
end