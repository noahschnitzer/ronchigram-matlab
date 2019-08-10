%% setup
figure('Position',[1 1 800 700]); 
hold on;
cone_zfactor = 4;
n_at = 100;
spec_sf_x  = .33;
spec_sf_z = .2;
cone_samp = 20;
ap_samp = 100;
ap_extent = .4;
ronch_samp = 512;
ronch_extent = .35;

samp_color = [99 129  39]/255;
cone_color = [0 134 192]/255;
cone_alpha = .25;
ap_color = [119 108 193]/255;
%% prespecimen cone
r = linspace(0,1/cone_zfactor,cone_samp) ;
th = linspace(0,2*pi,cone_samp) ;
[R,T] = meshgrid(r,th) ;
X = R.*cos(T) ;
Y = R.*sin(T) ;
Z = cone_zfactor.*R ;
surf(X,Y,Z,'FaceColor',cone_color,'FaceAlpha',cone_alpha);
%% aperture
[X,Y] = meshgrid(linspace(-1,1,ap_samp),linspace(-1,1,ap_samp));
X = X.* ap_extent; Y=  Y.* ap_extent;
Z = .5.*ones(size(X));  %placing aperture halfway up cone
Z(X.^2+Y.^2 < (1/cone_zfactor./2).^2) = nan;
surf(X,Y,Z,'FaceColor',ap_color,'FaceAlpha',.8,'EdgeColor','none');
[X,Y] = meshgrid(linspace(-1,1,ap_samp/6),linspace(-1,1,ap_samp/6));
X = X.* ap_extent; Y=  Y.* ap_extent;
Z = .5.*ones(size(X));  %placing aperture halfway up cone
Z(X.^2+Y.^2 < (1/cone_zfactor./2).^2) = nan;
surf(X,Y,Z,'FaceColor','none','FaceAlpha',.8);
%% specimen
rand_pos = rand(3,n_at);
rand_pos(1:2,:) = rand_pos(1:2,:).*spec_sf_x-spec_sf_x/2;
rand_pos(3,:) = rand_pos(3,:).*spec_sf_z - spec_sf_z/2;
scatter3(rand_pos(1,:),rand_pos(2,:),rand_pos(3,:),15,'filled','MarkerFaceColor',samp_color,'MarkerEdgeColor','k');

%% post specimen cone
X = R.*cos(T) ;
Y = R.*sin(T) ;
Z = -cone_zfactor.*R ;
surf(X,Y,Z,'FaceColor',cone_color,'FaceAlpha',cone_alpha);

%% ronch
[X,Y] = meshgrid(linspace(-ronch_extent,ronch_extent,ronch_samp),linspace(-ronch_extent,ronch_extent,ronch_samp));
Z = ones(size(X)).*-1;

ab = aberration_generator(1);
ronch = shifted_ronchigram(ab,[0 0],128,ronch_samp,180);

a = surf(X,Y,Z,...    % Plot the surface
     'CData',(ronch-.5).*2,...
     'FaceColor','texturemap', 'EdgeColor','none'); colormap gray;


%% final
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
view(3);
axis off;

