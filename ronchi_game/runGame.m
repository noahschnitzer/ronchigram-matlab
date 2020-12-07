addpath('../simulation')
addpath('../assessment')
clear
%%
load net-26.mat trainedNet
load ab_set.mat abers
abers = abers(1:end);

ronchi_game(abers, trainedNet)