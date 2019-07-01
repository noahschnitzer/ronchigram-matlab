%% load dataset
load('dataset_24_vals.mat');

%% load dists
dist = meta{8};

%% selection rule
%mult = 37;
%selection = 1:mult:1000*mult-1;%
%selection = 1:length(dist);

% 1000 aberrations s.t.:
% they comrpise a uniform distribution
% they range from 2 to 101 mrad (100 bins)
% hence, 10 at each
selection = [];
for it = 2:101
   indices = find(val_table.Strehl == it);
   selection = [selection; indices(1:10)];
end
%% plotting
figure;histogram(val_table(selection,:).Strehl,0:110);
abs = dist(selection);

%% distribution of each aberrations
figure; 
for it = 1:length(abs(1).mag)
   l = [];
   order = abs(1).n(it);
   degree = abs(1).m(it);
   T = (order+1)*1.969e-12/(8*(50e-3)^(order+1));
   for jt = 1:length(abs)
      l(end+1) = abs(jt).mag(it).*abs(jt).unit(it);
   end
   stdev = std(l);
   subplot(4,4,it); histogram(l/T,[-10:.01:10]);
   title(['C' num2str(order) num2str(degree) ', T:' num2str(T) ', std:' num2str(stdev)]);
end
