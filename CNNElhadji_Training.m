imds = imageDatastore('FacesForCNN_PIR',...
    'IncludeSubfolders',true,...
    'LabelSource','foldernames');

labelCount = countEachLabel(imds)

imds.ReadFcn =@(loc)imresize(imread(loc),[224,224]);

[imTrain,imValidation] = splitEachLabel(imds,0.7,'randomized');
numTrainImages = numel(imTrain.Labels);


layers = [
    imageInputLayer([224 224 3])
    
    convolution2dLayer(3,9,'Padding','same') %la couche convolution avec 8 filtres de taille 3x3 le padding=same permet d'avoir la même taille d'image
                                             %en entrée qu'en sortie.
                                             
    batchNormalizationLayer
    %tanhLayer
    leakyReluLayer
    %reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,18,'Padding','same')
    batchNormalizationLayer
    %tanhLayer
    leakyReluLayer
    %reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,36,'Padding','same')
    batchNormalizationLayer
    %tanhLayer
    leakyReluLayer
    %reluLayer
    
%     fullyConnectedLayer(10)
    fullyConnectedLayer(5)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.0001, ...
    'MaxEpochs',15, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imValidation, ...
    'ValidationFrequency',10, ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(imTrain,layers,options);
YPred = classify(net,imValidation);
YValidation = imValidation.Labels;
accuracy = sum(YPred == YValidation)/numel(YValidation)
save('net_training', 'net');