% Load and preprocess the data
imds = imageDatastore('path/to/images');
imds.ReadFcn = @(filename)readAndPreprocessImage(filename);

% Split the data into training and testing sets
[imdsTrain, imdsTest] = splitEachLabel(imds, 0.7, 'randomized');

% Define the CNN architecture
layers = [
    imageInputLayer([32 32 3])
    
    convolution2dLayer(3, 16, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2, 'Stride', 2)
    
    convolution2dLayer(3, 32, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2, 'Stride', 2)
    
    convolution2dLayer(3, 64, 'Padding', 'same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer
];

% Train the CNN
options = trainingOptions('sgdm', ...
    'InitialLearnRate', 0.01, ...
    'MaxEpochs', 10, ...
    'MiniBatchSize', 128, ...
    'Plots', 'training-progress');
net = trainNetwork(imdsTrain, layers, options);

% Evaluate the CNN
YPred = classify(net, imdsTest);
YTest = imdsTest.Labels;
accuracy = sum(YPred == YTest)/numel(YTest);
disp("Accuracy: " + accuracy);
