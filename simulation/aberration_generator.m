function [aberrations] = aberration_generator(num_gen)
    aberrations = struct;
    for abit = 1:num_gen
        %num_gen = 1;
        ang = 1e-10;
        nm = 1e-9;
        um = 1e-6;
        mm = 1e-3;

        aberrations(abit).n =    [  1, 1,  2,   2,   3,   3,   3,   4,   4,   4,   5, 5, 5, 5];
        aberrations(abit).m =    [  0, 2,  1,   3,   0,   2,   4,   1,   3,   5,   0, 2, 4, 6];
        aberrations(abit).angle = 180*(2*rand(1,length(aberrations(abit).n)))./aberrations(abit).n; %zeros(1,length(aberrations(abit).n));%
        aberrations(abit).unit = [ang,  nm,nm,  nm,  um,  um,  um,  mm,  mm,  mm,  mm,mm,mm,mm];
        %                    x , 157,95.5,10.4,10.4,5.22, 0.1,  xx,  xx,  10,10,10,10
        lims = [50, 50, 157,95.5, 10.4, 10.4,5.22,0.1,0.5,0.5,10,10,10,10]; % Kirkland Ultramicroscopy 2011
        scaling = .5*rand();
        for it = 1:length(lims)
            ab_val = scaling*lims(it)*(rand()-.5);
            aberrations(abit).mag(it) = ab_val;
        end
    end
    
end