function [probe_sizes] = probe_sizer(abs,imdim,simdim,aperture_sizes)
    n_ab = length(abs);
    probe_sizes = zeros(n_ab,length(aperture_sizes));
    for it = 1:n_ab
        ab = abs(it);
        ab_probe_sizes = zeros(1,length(aperture_sizes));
        chi0 = calculate_aberration_function(ab,imdim,simdim);
        for jt = 1:length(aperture_sizes)
            aperture_size = aperture_sizes(jt);
            %[~,~,~,probe,~] = shifted_ronchigram_o(ab,[0 0], aperture_size, imdim, simdim);
            [~,probe] = calculate_probe(chi0,imdim,simdim,aperture_size,[0 0]);
            ab_probe_sizes(jt) = resolution_test(probe,'effprobe');
        end
        %probe_sizes(it,jt) = res;
        probe_sizes(it,:) = ab_probe_sizes;
    end
end