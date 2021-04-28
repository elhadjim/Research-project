function newBoxPolygon = SUFR_GUI(guiimage1, guiimage2)
disp("SURF is launched");
if(size(guiimage1, 3)==3)
    boxImage = rgb2gray(guiimage1);
else
    boxImage = guiimage1;
end
% boxImage = rgb2gray(imread('waldo.png'));
% figure;
% imshow(boxImage);
% title('Image of a face');

if(size(guiimage2, 3)==3)
    sceneImage = rgb2gray(guiimage2);
else
    sceneImage = guiimage2;
end
% sceneImage = rgb2gray(imread('fwaldo.png'));
%sceneImage = imresize(imrotate(boxImage,-20),1.2);
% figure;
% imshow(sceneImage);
% title('Image of another face');

boxPoints = detectSURFFeatures(boxImage);
scenePoints = detectSURFFeatures(sceneImage);

% figure;
% imshow(boxImage);
% title('100 Strongest Feature Points from Box Image');
% hold on;
% plot(selectStrongest(boxPoints, 10000));

% figure;
% imshow(sceneImage);
% title('300 Strongest Feature Points from Scene Image');
% hold on;
% plot(selectStrongest(scenePoints, 500));

[boxFeatures, boxPoints] = extractFeatures(boxImage, boxPoints);
[sceneFeatures, scenePoints] = extractFeatures(sceneImage, scenePoints);

boxPairs = matchFeatures(boxFeatures, sceneFeatures);

matchedBoxPoints = boxPoints(boxPairs(:, 1), :);
matchedScenePoints = scenePoints(boxPairs(:, 2), :);
% figure;
% showMatchedFeatures(boxImage, sceneImage, matchedBoxPoints, ...
%     matchedScenePoints, 'montage');
% title('Putatively Matched Points (Including Outliers)');

[tform, inlierBoxPoints, inlierScenePoints] = ...
    estimateGeometricTransform(matchedBoxPoints, matchedScenePoints, 'affine');

% figure;
% showMatchedFeatures(boxImage, sceneImage, inlierBoxPoints, ...
%     inlierScenePoints, 'montage');
% title('Matched Points (Inliers Only)');

boxPolygon = [1, 1;...                           % top-left
        size(boxImage, 2), 1;...                 % top-right
        size(boxImage, 2), size(boxImage, 1);... % bottom-right
        1, size(boxImage, 1);...                 % bottom-left
        1, 1];                   % top-left again to close the polygon
    
    newBoxPolygon = transformPointsForward(tform, boxPolygon);
    
%     figure;
% imshow(sceneImage);
% hold on;
% line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'y');
% title('Detected Box');

disp("SURF is done");