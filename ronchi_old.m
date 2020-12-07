function gui_driver
rng('shuffle')
vert_dim = 550;
horz_dim = 700;
f = figure('Visible','off','Position',[50,50,horz_dim,vert_dim]);
%buttons
hnext    = uicontrol('Style','pushbutton',...
             'String','Start Playing','Position',[530,475,150,50],...
             'Callback',@nextbutton_Callback,'FontSize',15);
hfinish    = uicontrol('Style','pushbutton',...
             'String','Finish Playing','Position',[530,0,150,25],...
             'Callback',@finishbutton_Callback,'FontSize',15);
% checkboxes
hcheckboxfeedback = uicontrol('Style','checkbox','String', 'Show Feedback', ...
            'Value',0,'Position',[530 450 200 25],'FontSize',15);
hcheckboxindivP4 = uicontrol('Style','checkbox','String', 'Show individual pi/4', ...
            'Value',0,'Position',[530 425 200 25],'FontSize',15);
hcheckboxStrehl = uicontrol('Style','checkbox','String', 'Show Strehl',...
            'Value',0,'Position',[530 400 150 25],'FontSize',15);
hcheckboxPlotProbe = uicontrol('Style','checkbox','String', 'Plot Probe',...
            'Value',0,'Position',[530 375 150 25],'FontSize',15);
% texts
htextStrehl  = uicontrol('Style','text','String','',...
           'Position',[530,250,150,50],'FontSize',15,'ForegroundColor','blue');
htextP4 = uicontrol('Style','text','String','',...
           'Position',[530,200,150,50],'FontSize',15,'ForegroundColor','cyan');
htextGuess = uicontrol('Style','text','String','',...
           'Position',[530,150,150,50],'FontSize',15,'ForegroundColor','red');
htextindivP4 = uicontrol('Style','text','String','',...
           'Position',[530,100,150,50],'FontSize',15,'ForegroundColor','green');
htextStatus = uicontrol('Style', 'text','String', '',...
            'Position', [530, 300, 150, 50], 'FontSize',15);
htextTag = uicontrol('Style','text','String', 'Data Tag',...
            'Position', [530,60,150,30],'FontSize',15);
htextCredit = uicontrol('Style','text','String','Noah Schnitzer and Suk Hyun Sung @ Hovden Lab',...
            'Position', [10,0,300,20],'FontSize',12,'ForegroundColor',[1 1 1]*100/255);
ha = axes('Units','pixels','Position',[25,25,500,500]);
%input
htag = uicontrol('Style','edit', 'Position', [530, 30, 150, 30]);
%%% For NN
%%% load net_12.mat net
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
ha.Units = 'normalized';
hnext.Units = 'normalized';
hfinish.Units = 'normalized';
hcheckboxfeedback.Units = 'normalized';
hcheckboxindivP4.Units = 'normalized';
hcheckboxStrehl.Units = 'normalized';
hcheckboxPlotProbe.Units = 'normalized';
htextStrehl.Units = 'normalized';
htextP4.Units = 'normalized';
htextGuess.Units = 'normalized';
htextindivP4.Units = 'normalized';
htextStatus.Units = 'normalized';
htextTag.Units = 'normalized';
htag.Units = 'normalized';
htextCredit.Units = 'normalized';
% Assign the a name to appear in the window title.
f.Name = 'Ronchigram Game';
% Move the window to the center of the screen.
movegui(f,'center')
imshow(zeros(256,256));
it = 0;
% Make the window visible.
f.Visible = 'on';
results = [];
aberrations = aberration_generator(1);
user_sel = [];
imdim = 256;
simdim = 110; 
ap_size = 90;
shifts = [];
min_p4 = 0;
  function nextbutton_Callback(source,eventdata) 
    htextStrehl.String = '';
    htextGuess.String = '';
    htextindivP4.String = '';
    htextP4.String = '';
    htextStatus.String = '';
    % Prompt user to draw circle
    if strcmp(hnext.String,'Next Ronchigram') | strcmp(hnext.String,'Start Playing')
        it = it + 1;
        shifts = round(40.*(rand(1,2)-.5));
        new_ab = aberration_generator(1);
        aberrations(it) = new_ab(1);
        [im, ~, min_p4,~, ~] = shifted_ronchigram(aberrations(end),shifts,ap_size,imdim,simdim);
        %aberrations(end).mag.*aberrations(end).unit
        %aberrations(end).angle
        imagesc(im);
        colormap gray;
        axis equal off;
        title(['Trial ' num2str(it)]);
        hnext.Visible = false;
        user_sel = drawcircle('FaceAlpha',0,'Color','red');
        %thread halts until circle drawn
        hnext.Visible = true;
        hnext.String = 'Done';
        % now wait for user to press hnext
    else 
        % User has drawn circle, pressed hnext to indicate they are finished
        guessed_radius = user_sel.Radius.*2*simdim/imdim; 
        pred = 0;
        strehl_ap = 0;
        %%% for NN
        %%%     if hcheckboxindivP4.Value
        %%%         pred = double(net.predict(im.*256));
        %%%     end
        %uncalculated values -> 0 unless user chooses to calculate feedback now
        results(end+1,:) = [strehl_ap, 0, guessed_radius, pred];
        if hcheckboxfeedback.Value
            htextStatus.String = 'Radii:';
            cutoff1 = .8;
            [strehl_ap, ~] = strehl_calculator(aberrations(end), imdim, simdim, cutoff1,0);
            [indiv_p4_ap] = indiv_p4_calculator(aberrations(end), imdim, simdim);
            results(end,:) = [min_p4,strehl_ap, guessed_radius, indiv_p4_ap];
            htextP4.String = [ ' Pi/4: ' num2str(min_p4) 'mrad'];
            htextGuess.String = ['Guess: ' num2str(round(guessed_radius,1)) 'mrad'];
            viscircles([imdim/2-shifts(1) imdim/2-shifts(2)], min_p4.*1/2*imdim/simdim,'Color','cyan');
            vis_aps = [guessed_radius,min_p4];
            if hcheckboxStrehl.Value ~= 0
                htextStrehl.String = [num2str(cutoff1) 'Strehl: ' num2str(round(strehl_ap,1)) 'mrad'];
                viscircles([imdim/2-shifts(1) imdim/2-shifts(2)], strehl_ap.*1/2*imdim/simdim,'Color','blue');
                vis_aps(end+1) = strehl_ap;
            else
                htextStrehl.String = '';
            end
            if hcheckboxindivP4.Value ~= 0
                htextindivP4.String = ['Individual pi/4: ' num2str(indiv_p4_ap) 'mrad'];
                viscircles([imdim/2-shifts(1) imdim/2-shifts(2)], indiv_p4_ap.*1/2*imdim/simdim,'Color','green');
                vis_aps(end+1) = indiv_p4_ap;
            else
                htextindivP4.String = '';
            end
            if hcheckboxPlotProbe.Value 
                probe_vis(figure,aberrations(end), vis_aps,simdim, imdim);
            end
            hnext.String = 'Next Ronchigram';
        else %feedback
            hnext.String = 'Next Ronchigram';
            nextbutton_Callback;
        end
    end
  end
  function finishbutton_Callback(source,eventdata) 
  % Display mesh plot of the currently selected data.
        cla
        err = results(:,1)-results(:,3);
        RMSE = sqrt(mean(err.^2));
        htag.Visible = false;
        htextTag.String = ['Tag: ' htag.String];
        if hcheckboxfeedback.Value
            htextStatus.String = ['RMSE: ' num2str(RMSE) ' mrad'];
            plot(results(:,1),'o','MarkerSize',15,'MarkerEdgeColor','blue');
            hold on;
            plot(results(:,2),'x','MarkerSize',15,'MarkerEdgeColor','cyan');
            plot(results(:,3),'^','MarkerSize',15,'MarkerEdgeColor','red');
            plot(results(:,4),'*','MarkerSize',15,'MarkerEdgeColor','green');
            legend('Sum Pi/4','Strehl','Guess','Individual Pi/4');
            set(gca,'FontSize',15);
            xlabel('Trial');
            ylabel('Aperture Size');
        end
        hnext.Visible = false;
        hfinish.Visible = false;
        hcheckboxfeedback.Visible = false;
        hcheckboxindivP4.Visible = false;
        hcheckboxStrehl.Visible = false;
        hcheckboxPlotProbe.Visible = false;
        tag = htag.String;
        t=datetime('now');
        format = 'yyyymmdd_HHMMSS';
        fname = ['game_results/results_' datestr(t,format) '_tag_' tag];
        aberrations = aberrations(2:end);
        if ~exist('game_results','dir')
            mkdir('game_results');
        end
        save(fname,'results','aberrations');
  end
end