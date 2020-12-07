addpath('../simulation')
addpath('../assessment')
clear
%%
load net-26.mat trainedNet
load ab_set.mat abers
abers = abers(1:end);

ronchi_game(abers, trainedNet)

% For optimal performance with the net-26.mat weights, enable  
% the legacy ronchigram simulation mode to match training 
% simulation parameters:
% ronchi_game(abers, trainedNet, 'legacy');