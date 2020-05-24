function [mask] = aperture_mask(imdim,simdim,ap_size)
    al_max = simdim; 
    al_vec = (linspace(-al_max,al_max,imdim));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    al_rr = sqrt(alxx.^2 + alyy.^2);
    mask = al_rr<= ap_size;
end