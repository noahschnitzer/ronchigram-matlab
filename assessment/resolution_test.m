function res = resolution_test(probe, method)

    if strcmp(method,'rayleigh')
        res = rayleigh_resolution_test(probe);
    elseif strcmp(method,'probe')
        res = probe_size_resolution_test(probe);
    elseif strcmp(method,'effprobe')
        res = eff_probe_size_resolution_test(probe);
    else
        res = -1;
    end
    
end

function res = rayleigh_resolution_test(probe)
    probe_max = max(probe(:));
    threshold = .5;
    for disp = 10:-1:1
        disp
        for rot1 = 0:10:359
            for rot2 = 0:10:359
                %also need to rtoate p1
                p1 = imrotate(probe,rot1,'crop');
                p1 = imtranslate(p1,[0 disp]);
                p2 = imrotate(probe, rot2, 'crop');
                p2 = imtranslate(p2,[0 -disp]);
                super = p1+p2;
                over = super( floor(size(probe,1)/2+1),floor(size(probe,2)/2+1)); %floor or ceil?
                rat = over/probe_max;
                %figure; imagesc(super);

                if rat > threshold
                    figure; imagesc(super);
                   res = disp+1;
                   return;
                end
            end
        end

    end
end

% * RETURNS RADIUS NOT DIAMETER
function res = probe_size_resolution_test(probe)
    %probe size ~ diameter containing half of probe current per kirkland
    %ultramicroscp 2011, pg 1529
    numPx = size(probe,1);
    [x,y] = meshgrid(1:numPx, 1:numPx);
    cx = numPx/2+1; cy = cx;
    total_current = sum(probe(:));
    res = -1;
    for r = 1:numPx/2
       mask = (x-cx).^2+(y-cy).^2 <= r.^2;
       masked_probe = mask.*probe;
       masked_int = sum(masked_probe(:));
       if masked_int > 1/2*total_current
          res = r;
          return;
       end
    end
end

% based off https://www.mathworks.com/matlabcentral/fileexchange/56271-binarysearch-a-n-num
% now considering mid, mid+1, looking to see one on either side of 0.5 int.
function res = eff_probe_size_resolution_test(probe)
    numPx = size(probe,1);
    [x,y] = meshgrid(1:numPx, 1:numPx);
    cx = numPx/2+1; cy = cx;
    total_current = sum(probe(:));
    thrsh = .41*total_current;
    
    r = 1:numPx/2+1;
    left = 1;
    right = numPx/2;
    while(left <= right)
        mid = ceil((left+right)/2);
        mask1 = (x-cx).^2+(y-cy).^2 <= r(mid).^2;
        masked_probe1 = mask1.*probe;
        masked_int1 = sum(masked_probe1(:));
        mask2 = (x-cx).^2+(y-cy).^2 <= r(mid+1).^2;
        masked_probe2 = mask2.*probe;
        masked_int2 = sum(masked_probe2(:));
        %masked_int2 is wider, should always be larger
        if masked_int1 <= thrsh && masked_int2 >= thrsh
            res = r(mid+1);
            return 
        else
            if masked_int1 > thrsh
                right = mid - 1;
            else 
                left = mid + 1;
            end
            
        end
    end
    res = -1;
end