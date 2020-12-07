function ronchi_game( varargin )
rng('shuffle')
vert_dim = 550;
horz_dim = 700;

it = 0;
user_sel = [];
mode = '';

if nargin == 1
    abers = varargin{1};
    numMetric = 4;
elseif nargin == 2
    abers = varargin{1};
    cnn   = varargin{2};
    numMetric = 5;
elseif nargin == 3
    abers = varargin{1};
    cnn = varargin{2};
    numMetric = 5;
    mode = varargin{3};
end
    

numAbs = length(abers);
res = zeros(1,numAbs);
perms = randperm(numAbs);

imdim = 1024;
simdim = 180;
aperture_size = 128;
center = imdim/2;
box_dim = 512;%round(sqrt((imdim/simdim*aperture_size)^2/2));
crop_idx = imdim/2-box_dim/2 +1: imdim/2+ box_dim/2;

fig.f = figure('Visible','off','Position',[50,50,horz_dim,vert_dim]);
movegui(fig.f,'center')
fig.f.Visible = 'on';

fig.ax = axes('Units','pixels','Position',[25,25,500,500],'Units','normalized');

% Buttons
fig.button.start    = uicontrol('Style','pushbutton',...
             'String','Start','Position',[530,475,150,50],...
             'Callback',@startbutton_Callback,'FontSize',15,'Units','normalized'); 
fig.button.eval    = uicontrol('Style','pushbutton',...
             'String','Evaluate','Position',[530,475,150,50],...
             'Callback',@evalbutton_Callback,'FontSize',15,'Units','normalized','Visible','off');
fig.button.next    = uicontrol('Style','pushbutton',...
             'String','Next','Position',[530,475,150,50],...
             'Callback',@nextbutton_Callback,'FontSize',15,'Units','normalized','Visible','off');
% Checkbox
fig.feedBack = uicontrol('Style','checkbox','String', 'Show Feedback', ...
            'Value',0,'Position',[530 450 150 25],'FontSize',15);

% Image
im = zeros(imdim);
fig.im = imagesc( fig.ax, im );
axis(fig.ax,[0,imdim,0,imdim],'off')
colormap gray

% Title
fig.title = title( sprintf( '0 / %d',numAbs) );


labels = {'Guess', '0.8 Strehl',  'Pi/4', 'Individual Pi/4', 'CNN' };
colors = {'red', 'blue', 'green' 'cyan', 'magenta'};
radii  = [0,0,0,0];
hs = [];
circlesInit();

    function startbutton_Callback(source,eventdata)
        it = 1;
        fig.button.start.Visible = 'off';
        interaction()
        fig.button.eval.Visible = 'on';
    end

    function evalbutton_Callback(source,eventdata)
        fig.button.eval.Visible = 'off';
        res(perms(it)) = user_sel.Radius;
        user_sel.Visible = 'off';
        if fig.feedBack.Value == 1
            radii(1) = user_sel.Radius/(1/2*imdim/simdim);
            if hs(2).hCk.Value == 1
                radii(2) = strehl_calculator(abers(perms(it)), imdim, simdim, 0.8,0);
            end
            if hs(4).hCk.Value == 1
                radii(4) = indiv_p4_calculator(abers(perms(it)), imdim, simdim);
            end
            if hs(5).hCk.Value == 1
                radii(5) = cnn.predict(uint8(im(crop_idx,crop_idx)*255))*128;
            end
            circlesDraw( );

            if it < numAbs
                fig.button.next.Visible = 'on';
            elseif it == numAbs
                scaled_res = res*2*simdim/imdim;
                save('results/res_1.mat','res','scaled_res');
                disp('Saved results to results/res_1.mat');
                fig.title.String = 'Finished!';
            end
        else
            if it < numAbs
               it = it + 1;
               circlesOff()
               fig.button.next.Visible = 'off';
               interaction();
               fig.button.eval.Visible = 'on';
            elseif it == numAbs
                scaled_res = res*2*simdim/imdim;
                save('results/res_1.mat','res','scaled_res');
                disp('Saved results to results/res_1.mat');
                fig.title.String = 'Finished!';
            end
        end

    end
    function nextbutton_Callback(source,eventdata)
       if it <= numAbs
           it = it + 1;
           circlesOff()
           fig.button.next.Visible = 'off';
           interaction();
           fig.button.eval.Visible = 'on';
       end
    end

    function interaction()
        shifts = [0 0];
        if strcmp(mode,'legacy')
            [im,~,minp4,~] = shifted_ronchigram_o(abers(perms(it)),shifts,aperture_size,imdim,simdim);
        else
            [im,~,minp4,~] = shifted_ronchigram(abers(perms(it)),shifts,aperture_size,imdim,simdim);
        end
        radii(3) = minp4;
        fig.im.CData = im;
        fig.title.String=(sprintf('%d / %d',it, numAbs));

        %title([num2str(it) ' (' num2str(perms(it)) ')'  ]);
        user_sel = drawcircle('FaceAlpha',0,'Color','red');
    end
        
    function circlesInit()
        for ind = 1:numMetric
            hs(ind).label = labels{ind};
            hs(ind).hC = ...
                drawcircle('Center',[center, center], 'Radius', 50,...
                'Visible','off','InteractionsAllowed','none',...
                'FaceAlpha',0,'Color',colors{ind});
            hs(ind).str = uicontrol('Style','text','String','',...
                'Position',[530,50+50*ind,150,50],'FontSize',15,...
                'ForegroundColor',colors{ind});
            
            hs(ind).hCk = ...
                uicontrol('Style','checkbox',...
                'String', ['Show ', hs(ind).label], ...
                'Value',0,'Position',[530, 300+25*ind, 200, 25],'FontSize',15);
                
        end
        hs(1).hCk.Value = 1;
        hs(1).hCk.Visible = 'off';
    end
    function circlesOff()
        for ind = 1:numMetric
            hs(ind).hC.Visible = 'off';
        end
    end
    function circlesDraw( )
        for ind = 1:numMetric
            if hs(ind).hCk.Value == 1
                hs(ind).hC.Radius = radii(ind)*1/2*imdim/simdim;
                hs(ind).hC.Visible = 'on';
                hs(ind).str.String = sprintf('%s: %.1f mrad',hs(ind).label,radii(ind));
            else
                hs(ind).str.String = sprintf('%s: N/A',hs(ind).label);
            end
        end
    end

end
