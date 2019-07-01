function im = perlin_noise(n,m)

    im = zeros(n,m);
    i = 0;
    w = sqrt(n*m);

    while w > 3
        i = i + 1;
        d = interp2(randn(n, m), i-1, 'spline');
        im = im + i * d(1:n, 1:m);
        w = w - ceil(w/2 - 1);
    end
    
    im = im - min(im(:));
    im = im/max(im(:));
end