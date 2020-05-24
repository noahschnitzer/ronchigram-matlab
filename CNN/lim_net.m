%% Loading data
%load from generated dataset. Easiest to format training and validation
%data as a table with a column 'Predictors' with filenames for Ronchigrams
%and a column with the label e.g. 'Strehl' with the 0.8 Strehl ratio
%aperture size
%For code below, assume training data in a table 'train_table', validation
%data in 'validation_table'
%e.g.:
%   train_table = table(filepaths,  strehl_apertures);
%   train_table.Properties.VariableNames = {'Predictors', 'Strehl'};

%% Training

layers = [ ...
    imageInputLayer([512 512 1],'Normalization','none')
    %filterSize, numFilters, 'Stride', stride, 'Padding', padding
    
    convolution2dLayer(11,32,'Stride',1,'Padding','same') %5,8 %11,16
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(5,32,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2, 'Stride',2)

    
    convolution2dLayer(3,64,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer

    averagePooling2dLayer(2, 'Stride',2)
    
    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer

    averagePooling2dLayer(2, 'Stride',2)
    
    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer

    convolution2dLayer(3,128,'Stride',1,'Padding','same')
    batchNormalizationLayer
    reluLayer

    averagePooling2dLayer(2, 'Stride',2)
    
    %fullyConnectedLayer(32)
    %reluLayer
    
    fullyConnectedLayer(16)
    reluLayer
    
    fullyConnectedLayer(8)
    reluLayer
    
    fullyConnectedLayer(1)
    regressionLayer];

options = trainingOptions('sgdm',...
    'InitialLearnRate',5e-5,...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.1,...
    'LearnRateDropPeriod',50,...
    'MaxEpochs',10,...
    'Shuffle','every-epoch',...
    'ValidationData', validation_table,...
    'ValidationPatience',inf,...
    'ValidationFrequency',1000,...
    'Verbose',true,...
    'Plots','training-progress',...
    'MiniBatchSize', 32,...
    'CheckpointPath' , 'checkpoints'); %64 ->16 for 1024 imgs

[trainedNet,traininfo] = trainNetwork(train_table,layers,options);
save('trained_net.mat',trainedNet,traininfo,options,meta);