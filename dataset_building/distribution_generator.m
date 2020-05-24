function [aberrations] = distribution_generator(num_abs, range, imdim, simdim)
    aberrations = aberration_generator(1);
    bins = range(1):5:range(2);
    gen_per_iteration = length(bins)*100;
    for iteration = 1:num_abs/gen_per_iteration
        for bin = bins
            aberrations = cat(2,aberrations, sub_distribution_generator(100, [bin bin+4], imdim, simdim));
            display(['distribution: ' num2str(length(aberrations)) '/' num2str(num_abs)]);
        end
    end
    aberrations(1) = [];

end



%generates a somewhat uniform/normal distribution of pi/4 ratios for
%at least num_abs ronchigrams across range mrads. Calls aberration_generator to
%generate basis aberrations. Range is inclusive inclusive? 
%[range_0, range_end)]
function [aberrations] = sub_distribution_generator(num_abs, range, imdim, simdim)
    %terrible strategy: for each bin, generate aberrations until I have
    %enough to fill it!
    
    %better: generate aberrations -- if useful, add to set. If not, throw
    %away. Continue until no bins need to be filled!
    % further: incorporate heuristic -- distribution mean? median? -- to
    % apply factor to abs to shift one way or other to efficiently fill all
    % bins to the extrema
    tic
    summed = 0;
    total = 0;
    speed = .01;%0.0001;
    factor = 1;%.01;
    mean_diff = 0;
    
    bins = zeros(1, 500); %reasonable upper bound for generated aberrations and desired ap size
    bins(range(1):range(2)) = ceil(num_abs/length(range(1):range(2)));
    aberrations = aberration_generator(1);
    rem = sum(bins);
    
    while rem ~= 0
        new_ab = aberration_generator(1);
        if factor < speed *1.5
           speed = speed/10;
        elseif factor > 150*speed
            speed = speed*10;
        end
        factor = factor + speed*sign(mean_diff);
        new_ab.mag = new_ab.mag*factor;
        pi4 = floor(pi4_calculator(new_ab, imdim, simdim));
        if pi4 ~= inf
            
        
            total = total+1;
            summed = summed + pi4;
            mean_diff = summed/total-(sum(          bins(range(1):range(2)).* (range(1): range(2))     )/sum(bins(range(1): range(2))) );
        else
            mean_diff = 1;
            total = total+1;
            summed = summed+100;
        end
%         if total > 10 || mean_diff == inf
%            total = 0;
%            summed = 0;
%         end
        
        if pi4 > 0 && pi4 < 500 && bins(pi4) > 0
            aberrations(end+1) = new_ab;
            bins(pi4) = bins(pi4) - 1;
            rem = sum(bins);
            %factor = 1;
            total = 0;
            summed = 0;
        end
        
    end
    aberrations(1) = [];
    toc
    
    %strehls = strehl_calculator(aberrations, 256, simdim, .9, 0);
    %histogram(strehls,[0:100]);
    
    %0
end