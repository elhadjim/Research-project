close all;

faces_trainingset_construction = true;
nonfaces_trainingset_construction = true;

%% Settings for the construction of the LBP images and their global (spatially enhanced) histograms
windows = 7;   % number of regions per line and per columns (total number of regions = windows²)
m = windows*windows;    % number of regions
nbins = 257;    % number of bins in a histogram
uniform_patterns = false; % considering uniform patterns when computing histograms
weighted = true; % using weighted chi square distance
firstBar = 0; % number of the first histogram bin to consider when computing the distance (from 0 to 255)
lastBar = 255; % number of the last hisotgram bin to consider when computing the distance (from 0 to 255, or 0 to 255 when using uniform patterns, since non-uniform ones are stored in the 256th bin)


%% Setting up
facebasedirectory = 'LBP_TrainingSet/Faces/';
nonfacebasedirectory = 'LBP_TrainingSet/NonFaces/';
imagetype = '**/*.gif';
faceimagefiles = dir(fullfile(facebasedirectory, imagetype));
nonfaceimagefiles = dir(fullfile(nonfacebasedirectory, imagetype));
l_face = length(faceimagefiles); % number of face images
l_nonface = length(nonfaceimagefiles); % number of non-face images
sz = size(imread([faceimagefiles(1).folder '\' faceimagefiles(1).name]));   % dimensions of an image
face_images = zeros(sz(1), sz(2), l_face);    % (3D matrix) - "list" of l images (as matrices)    
nonface_images = zeros(sz(1), sz(2), l_nonface);

% extraction of the images as 2D matrices of intensity pixels
for i=1:l_face
    currentimage = imread([faceimagefiles(i).folder '\' faceimagefiles(i).name]);
    if(size(currentimage, 3)==3)
        graycurrentimage = rgb2gray(currentimage);  % converting in gray levels
%         face_images(:,:,i) = im2double(graycurrentimage);
        face_images(:,:,i) = (graycurrentimage);
    else
%         face_images(:,:,i) = im2double(currentimage);
        face_images(:,:,i) = (currentimage);
    end
end

for i=1:l_nonface
    currentimage = imread([nonfaceimagefiles(i).folder '\' nonfaceimagefiles(i).name]);
    if(size(currentimage, 3)==3)
        graycurrentimage = rgb2gray(currentimage);
%         nonface_images(:,:,i) = im2double(graycurrentimage);
        nonface_images(:,:,i) = (graycurrentimage);
    else
%         nonface_images(:,:,i) = im2double(currentimage);
        nonface_images(:,:,i) = (currentimage);
    end
end

pixels_window_h = fix(sz(2)/windows);   % number of pixels per region in a line
pixels_window_v = fix(sz(1)/windows);   % number of pixels per region in a column
neighborhood = zeros(3,3);  % temporary matrix used to compute the new value of each pixel using its 3-by-3 neighborhood


%% Faces
% Construction of the histograms of faces
if(faces_trainingset_construction)
    LBP_faceimages = uint8(zeros(sz(1),sz(2),l_face));  % "LBP" images of faces
    facehistograms = zeros(nbins, m, l_face);  % "list" of spatially enhanced histograms of all faces (as matrices)
    for t=1:l_face  % for every image
        global_histogram = zeros(nbins, m);   % list of vectors (one vector representing a local (regional) histogram)
        % Construction of the LBP image
        for i=1:sz(1) % parcours de tous les pixels de la matrice
            for j=1:sz(2)
                for k=-1:1  % parcours des voisins du pixel (i,j)
                    for z=-1:1
                        k2 = k;
                        z2 = z;
                        if(i+k<1); k2 = k+1; end
                        if(i+k>sz(1)); k2 = k-1; end
                        if(j+z<1); z2 = z+1; end
                        if(j+z>sz(2)); z2 = z-1; end                                          
                        
                        if(face_images(i+k2,j+z2,t)>=face_images(i,j,t))    % tresholding the neighbors' values with the pixel's one
                            neighborhood(2+k, 2+z) = 1;
                        else
                            neighborhood(2+k, 2+z) = 0;
                        end
                    end
                end
                binary = [neighborhood(2,1) neighborhood(3,1) neighborhood(3,2) neighborhood(3,3) neighborhood(2,3) neighborhood(1,3) neighborhood(1,2) neighborhood(1,1)]; % combining all the neighbors' values into a binary number
%                 binary = [neighborhood(2,1) 1 neighborhood(3,2) 1 neighborhood(2,3) 1 neighborhood(1,2) 1];
                LBP_faceimages(i,j,t) = bi2de(binary, 'right-msb');  % converting the binary value into a decimal one
            end
        end

%         figure;
%         imshow(face_images(:,:,t));
%         figure;
%         imshow(LBP_faceimages(:,:,t));

        % Construction of the histograms
%         imwrite(LBP_faceimages(:,:,t), ("LBP_TrainingSet/LBP/" + num2str(t) + ".gif"));
        for i=1:windows    % parcours de toutes les régions de l'image
            for j=1:windows
                if(i==windows && j==windows)
                    yrange = (i-1)*pixels_window_v+1:length(LBP_faceimages(:,1,t));
                    xrange = (j-1)*pixels_window_h+1:length(LBP_faceimages(1,:,t));
                    width = length(LBP_faceimages(1,:,t))-((j-1)*pixels_window_h+1);
                    height = length(LBP_faceimages(:,1,t))-((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                elseif(i==windows)
                    yrange = (i-1)*pixels_window_v+1:length(LBP_faceimages(:,1,t));
                    xrange = (j-1)*pixels_window_h+1:j*pixels_window_h;
                    width = j*pixels_window_h - ((j-1)*pixels_window_h+1);
                    height = length(LBP_faceimages(:,1,t)) - ((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                elseif(j==windows)
                    yrange = (i-1)*pixels_window_v+1:i*pixels_window_v;
                    xrange = (j-1)*pixels_window_h+1:length(LBP_faceimages(1,:,t));
                    width = length(LBP_faceimages(1,:,t)) - ((j-1)*pixels_window_h+1);
                    height = i*pixels_window_v - ((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                else
                    yrange = (i-1)*pixels_window_v+1:i*pixels_window_v;
                    xrange = (j-1)*pixels_window_h+1:j*pixels_window_h;
                    width = j*pixels_window_h - ((j-1)*pixels_window_h+1);
                    height = i*pixels_window_v - ((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                end
                LBP_faceimages = uint8(LBP_faceimages);
                c = histo(LBP_faceimages(yrange, xrange, t), 257, uniform_patterns);
%                 [c,~] = imhist(LBP_faceimages(yrange, xrange, t), 256); % c is a (column) vector containing the histogram counts for one region (the local histogram)
                global_histogram(:,(i-1)*windows+j) = c;    % adding the local histogram in the global one (for the whole image t)
%                 rectangle('Position', [x y width height], 'EdgeColor', 'blue', 'LineWidth', 1);
            end
        end
        facehistograms(:,:,t) = global_histogram;   % addind the global histogram of image t in the list of histograms
    end
    
    % Calculation of the histogram model for faces
    Mface = zeros(nbins, m); % mean of all histograms for faces (dimension of a spatially enhanced histogram)
    for t=1:l_face
        Mface = Mface + facehistograms(:,:,t);
    end
    Mface = Mface / l_face;

% plot_hist(Mface, windows, nbins, uniform_patterns, firstBar, lastBar);   % displaying the mean spatially enhanced histogram Mface
    
end

%% Non-faces
% Construction of the histograms for non-faces
if(nonfaces_trainingset_construction)
    LBP_nonfaceimages = uint8(zeros(sz(1),sz(2),l_nonface));  % "LBP" images of non-faces
    nonfacehistograms = zeros(nbins, m, l_nonface);  % "list" of spatially enhanced histograms of non-faces (as matrices)
    for t=1:l_nonface
        global_histogram = zeros(nbins, m);   % list of vectors (one vector representing a local (regional) histogram)
        for i=1:sz(1) % parcours de tous les pixels de la matrice
            for j=1:sz(2)
                for k=-1:1  % parcours des voisins du pixel (i,j)
                    for z=-1:1
                        k2 = k;
                        z2 = z;
                        if(i+k<1); k2 = k+1; end
                        if(i+k>sz(1)); k2 = k-1; end
                        if(j+z<1); z2 = z+1; end
                        if(j+z>sz(2)); z2 = z-1; end
                        
                        if(nonface_images(i+k2,j+z2,t)>=nonface_images(i,j,t))    % tresholding the neighbors' values with the pixel's one
                            neighborhood(2+k, 2+z) = 1;
                        else
                            neighborhood(2+k, 2+z) = 0;
                        end
                    end
                end
                binary = [neighborhood(2,1) neighborhood(3,1) neighborhood(3,2) neighborhood(3,3) neighborhood(2,3) neighborhood(1,3) neighborhood(1,2) neighborhood(1,1)]; % combining all the neighbors' values into a binary number
%                 binary = [neighborhood(2,1) 0 neighborhood(3,2) 0 neighborhood(2,3) 0 neighborhood(1,2) 0];
                LBP_nonfaceimages(i,j,t) = bi2de(binary, 'right-msb');  % converting the binary value into a decimal one
            end
        end

%         figure;
%         imshow(nonface_images(:,:,t));
%         figure;
%         imshow(LBP_nonfaceimages(:,:,t));
        for i=1:windows    % parcours de toutes les régions de l'image
            for j=1:windows
                if(i==windows && j==windows)
                    yrange = (i-1)*pixels_window_v+1:length(LBP_nonfaceimages(:,1,t));
                    xrange = (j-1)*pixels_window_h+1:length(LBP_nonfaceimages(1,:,t));
                    width = length(LBP_nonfaceimages(1,:,t))-((j-1)*pixels_window_h+1);
                    height = length(LBP_nonfaceimages(:,1,t))-((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                elseif(i==windows)
                    yrange = (i-1)*pixels_window_v+1:length(LBP_nonfaceimages(:,1,t));
                    xrange = (j-1)*pixels_window_h+1:j*pixels_window_h;
                    width = j*pixels_window_h - ((j-1)*pixels_window_h+1);
                    height = length(LBP_nonfaceimages(:,1,t)) - ((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                elseif(j==windows)
                    yrange = (i-1)*pixels_window_v+1:i*pixels_window_v;
                    xrange = (j-1)*pixels_window_h+1:length(LBP_nonfaceimages(1,:,t));
                    width = length(LBP_nonfaceimages(1,:,t)) - ((j-1)*pixels_window_h+1);
                    height = i*pixels_window_v - ((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                else
                    yrange = (i-1)*pixels_window_v+1:i*pixels_window_v;
                    xrange = (j-1)*pixels_window_h+1:j*pixels_window_h;
                    width = j*pixels_window_h - ((j-1)*pixels_window_h+1);
                    height = i*pixels_window_v - ((i-1)*pixels_window_v+1);
                    x = (j-1)*pixels_window_h+1;
                    y = (i-1)*pixels_window_v+1;
                end
                LBP_nonfaceimages = uint8(LBP_nonfaceimages);
                c = histo(LBP_nonfaceimages(yrange, xrange, t), 257, uniform_patterns);
%                 [c,~] = imhist(LBP_nonfaceimages(yrange, xrange, t), 256);  % c is a (column) vector containing the histogram counts for one region (the local histogram)
                global_histogram(:,(i-1)*windows+j) = c;    % adding the local histogram in the global one (for the whole image t)
%                 rectangle('Position', [x y width height], 'EdgeColor', 'red', 'LineWidth', 1);
            end
        end
        nonfacehistograms(:,:,t) = global_histogram;   % addind the global histogram of image t in the list of histograms
    end

% Calculation of the histogram model for non-faces
    Mnonface = zeros(nbins, m); % mean of all histograms for faces (dimension of a spatially enhanced histogram)
    for t=1:l_nonface
        Mnonface = Mnonface + nonfacehistograms(:,:,t);
    end
    Mnonface = Mnonface / l_nonface;
    
%     plot_hist(Mnonface, windows, nbins, uniform_patterns, firstBar, lastBar);  % displaying Mnonface
    
end


%% Display some LBP images
% figure;
% for i=2:2:10
%     subplot(2,5,i/2);
%     imshow(mat2gray(LBP_faceimages(:,:,i-1)));
% end
% for i=2:2:10
%     subplot(2,5,5+i/2);
%     imshow(mat2gray(LBP_nonfaceimages(:,:,i-1)));
% end


%% Saving the variables
save('LBP_training', 'windows', 'm', 'nbins', 'uniform_patterns', 'weighted', 'firstBar', 'lastBar', 'Mface', 'Mnonface');




%% Histogram
function c = histo(Image, nbins, uniform_patterns)
    c = zeros(nbins, 1);
    sz = size(Image);
    n = sz(1)*sz(2);
    for i=1:sz(1)
        for j=1:sz(2)
            if(uniform_patterns==true && uniform(Image(i,j))~=true)
                c(257) = c(257) + 1;
            else
                c(Image(i,j)+1) = c(Image(i,j)+1) + 1;
            end
        end
    end
    c = c/n;
end

%% Uniform Pattern?
function bool = uniform(n)
    array = de2bi(n,8);
    sz = size(array);
    nbTransitions = 0;
    for i=1:(sz(2))-1
       if(array(i)~= array(i+1))
           nbTransitions = nbTransitions + 1;
       end
    end
    bool = (nbTransitions <= 2);
end


%% Display a spatially enhanced histogram H
function plot_hist(H2, windows, ~, uniform_patterns, firstBar, lastBar)
    figure;
    axis off;
    for i = 1:windows
        for j=1:windows
           H = H2((firstBar+1:lastBar+1),(i-1)*windows+j);
           
           if(uniform_patterns)
               d = [];
               k = 1;
               for u=1:size(H)
                   if(H(u)~=0)
                       d(k)=H(u);
                       k=k+1;
                   end
               end
               H=d(:);
           end
           
           subplot(windows, windows, (i-1)*windows+j);
           axis off;
%            bar(0:1:nbins-1, H(:,(i-1)*windows+j));
           bar(firstBar:1:size(H)-1+firstBar, H);
           axis off;
        end
    end
end
