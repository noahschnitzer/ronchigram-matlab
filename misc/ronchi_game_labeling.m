function ronchi_game
rng('shuffle')
vert_dim = 550;
horz_dim = 700;

f = figure('Visible','off','Position',[50,50,horz_dim,vert_dim]);


movegui(f,'center')

ha = axes('Units','pixels','Position',[25,25,500,500],'Units','normalized');
hstart    = uicontrol('Style','pushbutton',...
             'String','Start','Position',[530,475,150,50],...
             'Callback',@startbutton_Callback,'FontSize',15,'Units','normalized');
hnext    = uicontrol('Style','pushbutton',...
             'String','Next','Position',[530,475,150,50],...
             'Callback',@nextbutton_Callback,'FontSize',15,'Units','normalized','Visible','off');

it = 0;
user_sel = [];
f.Visible = 'on';
load ab_set.mat abs
abs = abs(1:end);
res = zeros(1,length(abs));
perms = randperm(length(abs));
%abs_permd = abs(perms);

imdim = 1024;
simdim = 180;
aperture_size = 128;

    function startbutton_Callback(source,eventdata)
        it = 1;
        hstart.Visible = 'off';
        interaction()
        hnext.Visible = 'on';
    end

    function nextbutton_Callback(source,eventdata)
       res(perms(it)) = user_sel.Radius;
       it = it + 1;
       if it > length(abs)
          res
          scaled_res = res*2*simdim/imdim
          save('game_results/res_1.mat','res','scaled_res');
       end
       hnext.Visible = 'off';
       interaction();
       hnext.Visible = 'on';
    end

    function interaction()
        box_dim = 512;
        shifts = [0 0];
        im = shifted_ronchigram(abs(perms(it)),shifts,aperture_size,imdim,simdim);
        crop_idx = imdim/2-box_dim/2+1:imdim/2+box_dim/2;
        %im = im(crop_idx,crop_idx);
        imagesc(im);
        colormap gray;
        axis equal off;
        title([num2str(it)]);

        %title([num2str(it) ' (' num2str(perms(it)) ')'  ]);
        user_sel = drawcircle('FaceAlpha',0,'Color','red');
    end



end
