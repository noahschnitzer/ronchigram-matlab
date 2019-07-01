function def = optimize_defocus(ab)
    %heuristics:
        %- max strehl
        %- max pi/4
        %- probe...
        %- cumulative aberration fn, pos. of zeros?

    %Kirkland pg 33: delta_f = sqrt((2n_D - .5) C_s lambda) for n_D = 1 ->
    %scherzer
    
    %C3:5, C5:11
    C_s = ab.mag(5).*ab.unit(5) + ab.mag(11).*ab.unit(11);
    kev = 300;
    lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10; % E-wavelength in **meter**
    def = -sign(C_s)*sqrt((2*1-.5)*abs(C_s)*lambda)/ab.unit(1);
    %def = -C_s/ab.unit(1);
end