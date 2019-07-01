function plot_probe_comparison(ab,imdim,simdim,aperture_size);
    z_ab = ab;
    z_ab.mag(:) = 0;
    [~,~,~,probe,~] = shifted_ronchigram(ab,[0 0], aperture_size, imdim, simdim);
    %lambda = 1.97e-12;
    %px_to_ang = lambda*1/(2*simdim/1000)*10^10; 
    xd = (0:imdim/2)*px_to_ang(simdim,300);
    yd = radial_average(probe)/(pi*(aperture_size*imdim/(2*simdim)).^2).^2;%normalize_data(radial_average(probe),'total');
    figure; area(xd,yd,'LineWidth',2, 'FaceColor','blue','FaceAlpha',.5);
    hold on;
    [~,~,~,probe,~] = shifted_ronchigram(z_ab,[0 0], aperture_size, imdim, simdim);
    yd = radial_average(probe)/(pi*(aperture_size*imdim/(2*simdim)).^2).^2;%normalize_data(radial_average(probe),'total');
    area(-xd,yd,'LineWidth',2, 'FaceColor','red','FaceAlpha',.5);
    xlabel('Radius (Å)');
    ylabel('Normalized Intensity');
    xlim([-2 2]);
    set(gca,'FontSize',20);
    set(gca,'FontName','Helvetica Neue');

    set(gcf,'Position',[1 1 300 300]);


end