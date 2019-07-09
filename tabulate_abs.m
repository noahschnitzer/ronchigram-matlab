function tabulate_abs(fig,ab,include_idxs)
    figure(fig);
    axis equal off;
    ab_names = {};
    xlim([0 1]);
    ylim([0 1]);
    for it = 1:length(ab.mag)
        %ab_names{it,1} = ['C_' num2str(ab.n(it)) '_' num2str(ab.m(it))];
        %ab_names{it,2} = num2str(ab.mag(it).*ab.unit(it));
        ab_name = ['C_' num2str(ab.n(it)) '_' num2str(ab.m(it))];
        ab_val = ab.mag(it).*ab.unit(it);
        ab_unit = '';
        if abs(ab_val) < 1e-9
           ab_val = ab_val.*1e10;
           ab_unit = char(197);
        elseif abs(ab_val) < 1e-6
            ab_val = ab_val.*1e9;
            ab_unit = 'nm';
        elseif abs(ab_val) < 1e-3
            ab_val = ab_val.*1e6;
            ab_unit = [char(181) 'm'];
        else
            ab_val = ab_val.*1e3;
            ab_unit = 'mm';
        end
            text(0.0,1-it/length(ab.mag),ab_name);
        if(include_idxs(it))
            text(0.5,1-it/length(ab.mag),[num2str(ab_val,3) ' ' ab_unit],'HorizontalAlignment', 'right');
        end
        
    end
    %set(gcf,'FontName','Helvetica Neue');
    set(findall(gcf,'-property','FontName'),'FontName','Helvetica Neue');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);

    %uit = uitable(fig,'Data',ab_names);
end