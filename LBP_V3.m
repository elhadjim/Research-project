close all;

%% What part of the code do I want to execute?
new_face_detection = false;
faces_trainingset_construction = false;
nonfaces_trainingset_construction = false;
human_treshold_analysis = true;

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
        face_images(:,:,i) = (graycurrentimage);
    else
        face_images(:,:,i) = (currentimage);
    end
end

for i=1:l_nonface
    currentimage = imread([nonfaceimagefiles(i).folder '\' nonfaceimagefiles(i).name]);
    if(size(currentimage, 3)==3)
        graycurrentimage = rgb2gray(currentimage);
        nonface_images(:,:,i) = (graycurrentimage);
    else
        nonface_images(:,:,i) = (currentimage);
    end
end

pixels_window_h = fix(sz(2)/windows);   % number of pixels per region in a line
pixels_window_v = fix(sz(1)/windows);   % number of pixels per region in a column
neighborhood = zeros(3,3);  % temporary matrix used to compute the new value of each pixel using its 3-by-3 neighborhood


%% Faces
% Construction of the histograms of faces
if(faces_trainingset_construction)
    disp("Construction of the histograms of faces");
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
    disp("Construction of the histograms for non-faces");
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


%% Analysis of a new image
mindist=inf;
minx=0; miny=0; minheight=0; minwidth = 0;
if(new_face_detection)
    file = input('Enter image file name: ', 's');
    figure;
    hold on;
    imshow(imread(file));
    imagematrix = imread(file);
    if(size(imagematrix, 3)==3)
        imagematrix = (rgb2gray(imagematrix));    % converting the image as a gray levels image in case it is colored
    end
    original_size = size(imagematrix);

    for scale = 1:5    % every scale
        disp(scale);
        width = original_size(2)/scale;
        height = 4/3*width;
        x = 1;
        while (x+width-1<=original_size(2))    % every location
            y = 1;
            while(y+height-1<=original_size(1))
                sample = imagematrix(floor(y:y+height-1), floor(x:x+width-1));
                sz = size(sample);
                pixels_window_h = fix(sz(2)/windows);   % number of pixels per region in a line
                pixels_window_v = fix(sz(1)/windows);   % number of pixels per region in a column
                LBP_image = zeros(sz(1),sz(2));  % "LBP" image of the sample
                neighborhood = zeros(3,3);  % temporary matrix used to compute the new value of each pixel using its 3-by-3 neighborhood

                % LBP computation for the sample
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
                                
                                if(sample(i+k2,j+z2)>=sample(i,j))    % tresholding the neighbors' values with the pixel's one
                                    neighborhood(2+k, 2+z) = 1;
                                else
                                    neighborhood(2+k, 2+z) = 0;
                                end
                            end
                        end
                        binary = [neighborhood(2,1) neighborhood(3,1) neighborhood(3,2) neighborhood(3,3) neighborhood(2,3) neighborhood(1,3) neighborhood(1,2) neighborhood(1,1)]; % combining all the neighbors' values into a binary number
%                         binary = [neighborhood(2,1) 0 neighborhood(3,2) 0 neighborhood(2,3) 0 neighborhood(1,2) 0];
                        LBP_image(i,j) = bi2de(binary, 'right-msb');  % converting the binary value into a decimal one
                    end
                end
                LBP_image = uint8(LBP_image);
                
                % Computation of the histogram of the sample
                for i=1:windows    % parcours de toutes les régions de l'image
                    for j=1:windows
                        if(i==windows && j==windows)
                            lbpimage = LBP_image((i-1)*pixels_window_v+1:length(LBP_image(:,1)), (j-1)*pixels_window_h+1:length(LBP_image(1,:)));
                            c = histo(lbpimage, nbins, uniform_patterns);                        
%                             [c,b] = imhist(LBP_image((i-1)*pixels_window_v+1:length(LBP_image(:,1)), (j-1)*pixels_window_h+1:length(LBP_image(1,:))), nbins);                        
                        elseif(i==windows)
                            lbpimage = LBP_image((i-1)*pixels_window_v+1:length(LBP_image(:,1)), (j-1)*pixels_window_h+1:j*pixels_window_h);
                            c = histo(lbpimage, nbins, uniform_patterns);
                        elseif(j==windows)
                            lbpimage = LBP_image((i-1)*pixels_window_v+1:i*pixels_window_v, (j-1)*pixels_window_h+1:length(LBP_image(1,:)));
                            c = histo(lbpimage, nbins, uniform_patterns);
                        else
                            lbpimage = LBP_image((i-1)*pixels_window_v+1:i*pixels_window_v, (j-1)*pixels_window_h+1:j*pixels_window_h);
                            c = histo(lbpimage, nbins, uniform_patterns);
                        end
                        global_histogram(:,(i-1)*windows+j) = c;    % adding the local histogram in the global one (for the whole image t)
                    end
                end
                
%                 plot_hist(global_histogram, windows, nbins, uniform_patterns);

                % Comparison of this histogram with the model for faces and non-faces
                dist_face = chisq(global_histogram, Mface, weighted, m, windows, firstBar, lastBar);
                dist_nonface = chisq(global_histogram, Mnonface, weighted, m,  windows, firstBar, lastBar);
                if(dist_face<dist_nonface)
%                     disp(["this is a face: " dist_face ", " dist_nonface]);
                    rectangle('Position', [x y width height], 'EdgeColor', 'red', 'LineWidth', 1);
%                     newface = imcrop(imread(file), [x y width height]);
%                     imwrite(newface, "LBP_TrainingSet/DetectedFaces/face" + num2str(i) + num2str(x) + num2str(y) + ".jpg");
                    if(dist_face<mindist)
                        mindist = dist_face;
                        minx = x; minwidth = width; miny = y; minheight = height;
                    end
                else
%                     disp(["this is not a face: " dist_face ", " dist_nonface]);
                end
                y = y + height/2;
            end
            x = x + width/2;
        end
    end
    rectangle('Position', [minx miny minwidth minheight], 'EdgeColor', 'blue', 'LineWidth', 1);
%     newface = imcrop(imread(file), [minx miny minwidth minheight]);
%     imwrite(newface, "LBP_TrainingSet/DetectedFaces/bestface.jpg");
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Analysis of distance to the face space to define a treshold     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(human_treshold_analysis)
    disp("Analysis of the distances to face and non-face histograms");
    basedirectory3 = 'FacesVsNonFaces/';
    imagetype3 = '**/*.gif';
    imagefiles3 = dir(fullfile(basedirectory3, imagetype3));
    M3 = length(imagefiles3); % number of images in the base
    H3 = 0; %number of pictures of humans
    NH3 = 0; %number of pictures of non humans
    sz = size(imread([imagefiles3(1).folder '\' imagefiles3(1).name]));
    testimages = zeros(sz(1), sz(2), M3);
    facetofacedistances = [];
    nonfacetofacedistances = [];
    facetononfacedistances = [];
    nonfacetononfacedistances = [];
    
    for i=1:M3
        currentimage = imread([imagefiles3(i).folder '\' imagefiles3(i).name]);
        if(size(currentimage, 3)==3)
            graycurrentimage = rgb2gray(currentimage);
            testimages(:,:,i) = graycurrentimage;
        else
            testimages(:,:,i) = currentimage;
        end
        if(contains(imagefiles3(i).name,'human'))
            H3=H3+1;
        else
            NH3=NH3+1;
        end
    end
    
    LBP_testimages = uint8(zeros(sz(1),sz(2),M3));  % "LBP" images
%     testhistograms = zeros(nbins, m, M3);  % "list" of spatially enhanced histograms of all images (as matrices)
    for t=1:M3  % for every image
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

                        if(testimages(i+k2,j+z2,t)>=testimages(i,j,t))    % tresholding the neighbors' values with the pixel's one
                            neighborhood(2+k, 2+z) = 1;
                        else
                            neighborhood(2+k, 2+z) = 0;
                        end
                    end
                end
                binary = [neighborhood(2,1) neighborhood(3,1) neighborhood(3,2) neighborhood(3,3) neighborhood(2,3) neighborhood(1,3) neighborhood(1,2) neighborhood(1,1)]; % combining all the neighbors' values into a binary number
    %                 binary = [neighborhood(2,1) 1 neighborhood(3,2) 1 neighborhood(2,3) 1 neighborhood(1,2) 1];
                LBP_testimages(i,j,t) = bi2de(binary, 'right-msb');  % converting the binary value into a decimal one
            end
        end

        % Construction of the histograms
        for i=1:windows    % parcours de toutes les régions de l'image
            for j=1:windows
                if(i==windows && j==windows)
                    yrange = (i-1)*pixels_window_v+1:length(LBP_testimages(:,1,t));
                    xrange = (j-1)*pixels_window_h+1:length(LBP_testimages(1,:,t));
                elseif(i==windows)
                    yrange = (i-1)*pixels_window_v+1:length(LBP_testimages(:,1,t));
                    xrange = (j-1)*pixels_window_h+1:j*pixels_window_h;
                elseif(j==windows)
                    yrange = (i-1)*pixels_window_v+1:i*pixels_window_v;
                    xrange = (j-1)*pixels_window_h+1:length(LBP_testimages(1,:,t));
                else
                    yrange = (i-1)*pixels_window_v+1:i*pixels_window_v;
                    xrange = (j-1)*pixels_window_h+1:j*pixels_window_h;
                end
                LBP_testimages = uint8(LBP_testimages);
                c = histo(LBP_testimages(yrange, xrange, t), 257, uniform_patterns);
    %                 [c,~] = imhist(LBP_faceimages(yrange, xrange, t), 256); % c is a (column) vector containing the histogram counts for one region (the local histogram)
                global_histogram(:,(i-1)*windows+j) = c;    % adding the local histogram in the global one (for the whole image t)
    %                 rectangle('Position', [x y width height], 'EdgeColor', 'blue', 'LineWidth', 1);
            end
        end
%         testhistograms(:,:,t) = global_histogram;   % addind the global histogram of image t in the list of histograms
        testdist = chisq(global_histogram, Mface, weighted, m, windows, firstBar, lastBar);
        testdist2 = chisq(global_histogram, Mnonface, weighted, m, windows, firstBar, lastBar);
        if(contains(imagefiles3(t).name,'human'))
            facetofacedistances = [facetofacedistances testdist];
            facetononfacedistances = [facetononfacedistances testdist2];
        else
            nonfacetofacedistances = [nonfacetofacedistances testdist];
            nonfacetononfacedistances = [nonfacetononfacedistances testdist2];
        end
    end
    
%     disp(humandistances);
    figure;
    hold on;
    ax3(1) = subplot(1,2,1);
    boxplot(facetofacedistances);
    title('Distance btw hist of a face and face model');
    ax3(2) = subplot(1,2,2);
    boxplot(facetononfacedistances);
    title('Distance btw hist of a face and non-face model');
    linkaxes(ax3,'y');
    figure;
    hold on;
    ax4(1) = subplot(1,2,1);
    boxplot(nonfacetofacedistances);
    title('Distance btw hist of a non-face and face model');
    ax4(2) = subplot(1,2,2);
    boxplot(nonfacetononfacedistances);
    title('Distance btw hist of a non-face and non-face model');
    linkaxes(ax4,'y');
end


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
    for i=1:(sz(2))
       if(array(i)~= array(mod(i,sz(2))+1))
           nbTransitions = nbTransitions + 1;
       end
    end
    bool = (nbTransitions <= 2);
end

%% Chi Square distance calculation
function dist = chisq(S, M, weighted, m, windows, firstBar, lastBar)
    dist = 0;
    
    weights = ones(windows, windows);
    if(weighted)
        weights(1,1)=0; weights(1,windows)=0;
        weights(windows-2,1)=0; weights(windows-1,1)=0; weights(windows,1)=0; weights(windows,2)=0;
        weights(windows-2, windows)=0; weights(windows-1, windows)=0; weights(windows, windows)=0; weights(windows, windows-1)=0;
        weights(3,3)=2; weights(3,4)=2; weights(4,3)=2; weights(4,4)=2;
        weights(3,5)=2; weights(3,6)=2; weights(4,5)=2; weights(4,5)=2;
        weights(6,4)=2;
        weights(4,4)=1.5; weights(5,3)=1.5; weights(5,4)=1.5; weights(5,5)=1.5; weights(6,3)=1.5; weights(6,5)=1.5;
    end
    
    sum=0;
    for j=1:m
        dist = 0;
        for i=firstBar+1:lastBar+1
            if(S(i,j)~=0 && M(i,j)~=0)
                dist = dist + (((S(i,j)-M(i,j))*(S(i,j)-M(i,j)))/abs(S(i,j)+M(i,j)));
            end
        end
        sum = sum+(dist*weights(fix((j-1)/windows)+1, j - (windows*(fix((j-1)/windows)))));
    end
    dist = sum;
end


%% Display a spatially enhanced histogram H
function plot_hist(H2, windows, nbins, uniform_patterns, firstBar, lastBar)
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