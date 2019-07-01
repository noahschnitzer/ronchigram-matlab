% plots phase shift plot for ab. Aberrations length scale set by imdim,
% simdim
% Range set by range
% aberrations to plot set by indices.
function plot_phase_shift(ab, imdim, simdim, range, indices)
    scaled_range = range*2*simdim./(imdim);
    srange_x = imdim/2+1;
    srange_y = imdim/2+range;
    colordef;
    al_max = simdim * 10^-3; 
    al_vec = (linspace(-al_max,al_max,imdim));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    al_rr = sqrt(alxx.^2 + alyy.^2);
    al_rr = al_rr(srange_x,srange_y);
    
    % calculating total aberration fn and pi4 cutoff
    chi0 = calculate_aberration_function(ab,imdim,simdim);
    chi0 = chi0(srange_x,srange_y);
    max_p4 = 0;
    max_center = 0;
    for lim_center = -pi/4:pi/80:pi/4
        lb = lim_center - pi/4;
        ub = lim_center + pi/4;
        chi0_p4 = (chi0 < lb) | (chi0 > ub);
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_p4 = min(al_rr_p4(:))*1000;
        %display([num2str(lim_center) ' : ' num2str(min_p4) ]);
        if min_p4 > max_p4
            max_p4 = min_p4;
            max_center = lim_center;
        end
    end
    min_p4 = max_p4;
    %display(['pick: ' num2str(max_center)]);

    
    % calculating individual aberrations and plotting
    %ab_maps = zeros(imdim,imdim,length(indices));
    %pi4s = zeros(1,length(indices));
    indiv_chis = squeeze(zeros([size(al_rr) length(indices)]));
    figure; hold on;
    plot(scaled_range,chi0,'Color','k','LineWidth',4);

    
    
    for it = 1:length(indices)
        idx = indices(it);
        tchi0 = get_aberration(ab, imdim, simdim, idx);
        tchi0 = tchi0(srange_x,srange_y);
        chi0_p4 = abs(tchi0) > pi/4;
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_indiv_p4(it) = min(al_rr_p4(:))*1000;
        indiv_chis(:,it) = tchi0;
        %ab_maps(:,:,it) = get_aberration(ab, imdim, simdim, idx);
        %pi4s(it) = max(al_rr.*(ab_maps(:,:,it)>-pi/4 & ab_maps(:,:,it) < pi/4),[],'all');
    end
    
    [mp4val,mp4idx] = min(min_indiv_p4);
    for it = 1:length(indices)
       if it==mp4idx
           plot(scaled_range, indiv_chis(:,it),'Color',c_indiv_p4,'LineWidth',3);
       else
           plot(scaled_range, indiv_chis(:,it),'Color',c_grey,'LineWidth',3);
       end
        
    end
    % total p4 box
    line([0 scaled_range(end)],[max_center+pi/4 max_center+pi/4],'LineStyle','--','Color',c_mw_p4);
    line([0 scaled_range(end)],[max_center-pi/4 max_center-pi/4],'LineStyle','--','Color',c_mw_p4);
    patch([0 0 scaled_range(end) scaled_range(end)], [max_center-pi/4 max_center+pi/4 max_center+pi/4 max_center-pi/4],c_mw_p4,'FaceAlpha',.5,'EdgeColor','none');

    % indiv p4s bounds
    line([0 scaled_range(end)],[pi/4 pi/4],'LineStyle','--','Color',c_indiv_p4,'LineWidth',3);
    line([0 scaled_range(end)],[-pi/4 -pi/4],'LineStyle','--','Color',c_indiv_p4,'LineWidth',3);
    % ticks in pis
    yticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
    ylim([-pi pi]);
    xlim([0 scaled_range(end)]);
    set(gca,'FontSize',20);
    set(gca,'FontName','Helvetica Neue');
    xlabel('Convergence Angle');
    ylabel('Phase Shift');
    set(gcf,'Position',[1 1 500 300]);
end




function plot_phase_shift_v1(ab, imdim, simdim, range, indices)
    scaled_range = range*2*simdim./(imdim);
    srange_x = imdim/2+1;
    srange_y = imdim/2+range;
    colordef;
    al_max = simdim * 10^-3; 
    al_vec = (linspace(-al_max,al_max,imdim));
    [alxx,alyy] = meshgrid(al_vec,al_vec);
    al_rr = sqrt(alxx.^2 + alyy.^2);
    
    % calculating total aberration fn and pi4 cutoff
    chi0 = calculate_aberration_function(ab,imdim,simdim);
    max_p4 = 0;
    max_center = 0;
    for lim_center = -pi/4:pi/20:pi/4
        lb = lim_center - pi/4;
        ub = lim_center + pi/4;
        chi0_p4 = (chi0 < lb) | (chi0 > ub);
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_p4 = min(al_rr_p4(:))*1000;
        if min_p4 > max_p4
            max_p4 = min_p4;
            max_center = lim_center;
        end
    end
    min_p4 = max_p4;

    
    % calculating individual aberrations and plotting
    ab_maps = zeros(imdim,imdim,length(indices));
    pi4s = zeros(1,length(indices));
    figure; hold on;
    plot(scaled_range,chi0(imdim/2+1,imdim/2+range),'Color','k','LineWidth',4);

    
    
    for it = 1:length(indices)
        idx = indices(it);
        chi0_p4 = abs(get_aberration(ab, imdim, simdim, idx)) > pi/4;
        al_rr_p4 = chi0_p4 .* al_rr;
        al_rr_p4( al_rr_p4 == 0 ) = inf;
        min_indiv_p4(it) = min(al_rr_p4(:))*1000;
        %ab_maps(:,:,it) = get_aberration(ab, imdim, simdim, idx);
        %pi4s(it) = max(al_rr.*(ab_maps(:,:,it)>-pi/4 & ab_maps(:,:,it) < pi/4),[],'all');
    end
    
    [mp4val,mp4idx] = min(min_indiv_p4);
    for it = 1:length(indices)
       if it==mp4idx
           plot(scaled_range, ab_maps(imdim/2+1,imdim/2+range,it),'Color',c_indiv_p4,'LineWidth',3);
       else
           plot(scaled_range, ab_maps(imdim/2+1,imdim/2+range,it),'Color',c_grey,'LineWidth',3);
       end
        
    end
    % total p4 box
    line([0 scaled_range(end)],[max_center+pi/4 max_center+pi/4],'LineStyle','--','Color',c_mw_p4);
    line([0 scaled_range(end)],[max_center-pi/4 max_center-pi/4],'LineStyle','--','Color',c_mw_p4);
    patch([0 0 scaled_range(end) scaled_range(end)], [max_center-pi/4 max_center+pi/4 max_center+pi/4 max_center-pi/4],c_mw_p4,'FaceAlpha',.5,'EdgeColor','none');

    % indiv p4s bounds
    line([0 scaled_range(end)],[pi/4 pi/4],'LineStyle','--','Color',c_indiv_p4,'LineWidth',3);
    line([0 scaled_range(end)],[-pi/4 -pi/4],'LineStyle','--','Color',c_indiv_p4);

end
%        plot((range)*2*simdim/imdim, ab_maps(imdim/2+1,imdim/2+range,it),'Color',[200 200 200]/255,'LineWidth',3);

function plot_phase_shift_v0(ab, imdim, simdim, range)


    [~,chi0,minp4,~,~] = shifted_ronchigram(ab,[0 0], 128, imdim, simdim);
    %range = 0:200;
    colordef;
    figure; 
    hold on;
    spherical = get_aberration(ab,imdim,simdim,5);
    sp5 = get_aberration(ab,imdim,simdim,11);
    twofold = get_aberration(ab,imdim,simdim,2);
    plot((range)*2*simdim/imdim,spherical(513,512+range),'Color',[200 200 200]/255,'LineWidth',3);
    plot((range)*2*simdim/imdim,sp5(513,512+range),'Color',[200 200 200]/255,'LineWidth',3);
    plot((range)*2*simdim/imdim,twofold(513,512+range),'Color',c_indiv_p4,'LineWidth',4);

    plot((range)*2*simdim/imdim,-chi0(513,512+range),'Color','k','LineWidth',4);
    min_v = min(-chi0(513,512+range));
    line([0 range(end)*2*simdim/imdim],[min_v+pi/2 min_v+pi/2],'LineStyle','--','Color',c_mw_p4);
    line([0 range(end)*2*simdim/imdim],[min_v min_v],'LineStyle','--','Color',c_mw_p4);

    line([0 range(end)*2*simdim/imdim],[pi/4 pi/4],'LineStyle','--','Color',c_indiv_p4,'LineWidth',3);
    line([0 range(end)*2*simdim/imdim],[-pi/4 -pi/4],'LineStyle','--','Color',c_indiv_p4);
    %patch([0 0 range(end)*2*simdim/imdim range(end)*2*simdim/imdim], [-pi/4 pi/4 pi/4 -pi/4],c_indiv_p4,'FaceAlpha',.3,'EdgeColor','none');

    patch([0 0 range(end)*2*simdim/imdim range(end)*2*simdim/imdim], [min_v min_v+pi/2 min_v+pi/2 min_v],c_mw_p4,'FaceAlpha',.5,'EdgeColor','none');

    %r=rectangle('Position',[0 min_v range(end)*2*simdim/imdim pi/2],'FaceColor',c_mw_p4);
    %alpha(r,0.5)
    yticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});

    xlabel('Convergence Angle (mrad)');
    ylabel('Phase Shift \phi');
    set(gca,'FontSize',30);

end