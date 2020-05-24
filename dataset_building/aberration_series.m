function [results_strehl,results_pi4,results_indiv_p4] = aberration_series(ab, imdim, simdim, aperture_size, aberration_index, aberration_range, cuts, plot_ronch_indices)
   dfs = aberration_range;
   results_strehl = zeros(length(dfs),length(cuts));
   results_pi4 = zeros(length(dfs),1);
   results_indiv_p4 = zeros(length(dfs),1);    %%
   for jt = 1:length(dfs)
       df = dfs(jt);
       %% calculations
       ab.mag(aberration_index) = df;
       ronch = shifted_ronchigram(ab, [0 0], aperture_size, imdim, simdim);
       [strehls, Ss] = strehl_calculator(ab, 128, simdim, .9, 1);
       aps = zeros(size(cuts));
       labels = {};
       for it = 1:length(cuts)
           t = find(Ss > cuts(it));
           if length(t) > 0
           aps(it) = t(end);
           else
              aps(it) = 1;
           end
           results_strehl(jt, it) = aps(it);
           labels{it} = [num2str(cuts(it)) ' Strehl Ratio:' num2str(aps(it)) 'mrad'];
       end
       results_pi4(jt) = pi4_calculator(ab,imdim,simdim);
       results_indiv_p4(jt) = 0;%indiv_p4_calculator(ab, imdim, simdim);        %% plot ronch
       if any(plot_ronch_indices == jt)
   % %     if(abs(df) < 20)
           box_dim = 512;
           crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;
           figure; imagesc(ronch(crop_idx,crop_idx)); axis equal off; colormap gray;
           hold on;
           %figure; plot(Ss);
           center = [256 256];
           for it = 1:length(aps)
               c = get(gca,'colororder');
               c = c(mod(it-1,7)+1,:);
               viscircles(center, aps(it)*imdim/(2*simdim),'Color',c);
               plot(nan);
               %viscircles(center, 128*imdim/(2*simdim));
           end
           legend(labels);
           title([num2str(df) ' Å Defocus']);
           set(gca,'FontSize',12);    % %    end
       end
       %% plot probe    
   end
end
