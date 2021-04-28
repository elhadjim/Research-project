%Training Data 
imds = imageDatastore('FacesForCNN_PIR', 'IncludeSubfolders', 1, 'LabelSource', 'foldernames');
imds.ReadFcn = @(loc)imresize(imread(loc),[227,227]); %Fonction qui sert � redimensionner les images (227x227 impos� par Alexnet)
imds.Labels; %V�rifier que imds contient effectivement des donn�es

[trainingSet, testSet] = splitEachLabel(imds, 0.6, 'randomize'); 

%Alexnet : CNN
net = alexnet;

layers = net.Layers; % Take a look at the layers

layers; % notice the 1000 in the last fully connected layer. This is for the 1000 categories AlexNet knows.

numTrainImages = numel(trainingSet.Labels)
idx = randperm(numTrainImages,16)
figure
for i = 1:16
    subplot(4,4,i)
    I = readimage(trainingSet,idx(i));
    imshow(I)
end

inputSize = net.Layers(1).InputSize;

layersTransfer = net.Layers(1:end-3);
numClasses = numel(categories(trainingSet.Labels))

layers = [
    layersTransfer
    fullyConnectedLayer(numClasses,'WeightLearnRateFactor',20,'BiasLearnRateFactor',20)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'MiniBatchSize',10, ...
    'MaxEpochs',6, ...
    'InitialLearnRate',1e-4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',testSet, ...
    'ValidationFrequency',3, ...
    'Verbose',false, ...
    'Plots','training-progress');

netTransfer = trainNetwork(trainingSet,layers,options);

[YPred,scores] = classify(netTransfer,testSet);

idx = randperm(numel(testSet.Files),4);
figure
for i = 1:4
    subplot(2,2,i)
    I = readimage(imds,idx(i));
    imshow(I)
    label = YPred(idx(i));
    title(string(label));
end

save('alexnet_training', 'netTransfer');