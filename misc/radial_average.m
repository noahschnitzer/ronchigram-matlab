function [outfn] = radial_average(im)
    imdim = size(im,1);
    [mX, mY] = meshgrid(1:imdim,1:imdim);
    c = imdim/2+1;
    rad_im = sqrt((mX-c).^2+(mY-c).^2);
    edges = 0:c;
    outfn = zeros(1,length(edges)-1);
    for it = 2:length(edges)-1
        roi = rad_im >= edges(it-1) & rad_im < edges(it);
       outfn(it-1) = sum(im(   roi     )) / sum(roi(:));
    end

end