function def = calculate_scherzer_focus(ab)

    %Kirkland pg 33: delta_f = sqrt((2n_D - .5) C_s lambda) for n_D = 1 ->
    %scherzer
    
    %C3:5, C5:11
    C_s = ab.mag(5).*ab.unit(5) + ab.mag(11).*ab.unit(11);
    kev = 300;
    lambda = 12.3986./sqrt((2*511.0+kev).*kev) * 10^-10; % E-wavelength in **meter**
    %def = -sign(C_s)*sqrt((2*1-.5)*abs(C_s)*lambda)/ab.unit(1); %kirkland
    %def = -sign(C_s)*1.2*sqrt(abs(C_s)*lambda)/ab.unit(1); %wiki
    def = -sign(C_s)*sqrt(lambda*abs(C_s))/ab.unit(1); %weyland, muller
    %def = -C_s/ab.unit(1);
end