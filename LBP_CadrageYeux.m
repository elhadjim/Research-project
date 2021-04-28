function [xs, ys, widths, heights] = LBP_CadrageYeux(GUIimage)
    load('LBP_Eyetraining');
    disp("LBP is launched!");

    mindist=inf;
    minx=0; miny=0; minheight=0; minwidth = 0;
    xs = []; ys = []; widths = []; heights = []; distances = [];
    weights = ones(windows, windows);
    if(weighted)
        weights(1,1)=0.5; weights(1,3)=0.5; weights(3,1)=0.5; weights(3,3)=0.5; weights(2,2)=0.7;
        weights(1,2)=0.3; weights(3,2)=0.3;
        weights(2,1)=1; weights(2,3)=1;
    end
    imagematrix = GUIimage;
    if(size(imagematrix, 3)==3)
        imagematrix = (rgb2gray(imagematrix));    % converting the image as a gray levels image in case it is colored
    end
    original_size = size(imagematrix);
    
    for scale = 1:0.2:2    % every scale
%         disp([int2str(scale) '/' int2str(scaleMax)]);
        width = original_size(2)/scale;
        height = 9/16*width;
        x = 1;
        while (x+width-1<=original_size(2))    % every location
            y = 1;
            while(y+height-1<=original_size(1))
%                 for angle=-40:10:40
                    sample = imagematrix(floor(y:y+height-1), floor(x:x+width-1));
%                     sample = imrotate(sample, angle);
                    sample = padarray(sample, [1 1], 0, 'pre');
                    sample = padarray(sample, [1 1], 0, 'post');
                    sz = size(sample);
                    pixels_window_h = fix(sz(2)/windows);   % number of pixels per region in a line
                    pixels_window_v = fix(sz(1)/windows);   % number of pixels per region in a column
                    LBP_image = zeros(sz(1),sz(2));  % "LBP" image of the sample
                    neighborhood = zeros(3,3);  % temporary matrix used to compute the new value of each pixel using its 3-by-3 neighborhood

            % LBP computation for the sample
    %                 LBP_image = nlfilter(sample,[3 3], LBP);
                    for i=2:sz(1)-1 % parcours de tous les pixels de la matrice
                        for j=2:sz(2)-1
    %                         LBP_image(i,j) = LBP(sample(i-1:i+1, j-1:j+1));
                            for k=-1:1  % parcours des voisins du pixel (i,j)
                                for z=-1:1
    %                                 k2 = k;
    %                                 z2 = z;
    %                                 if(i+k<1); k2 = k+1; end
    %                                 if(i+k>sz(1)); k2 = k-1; end
    %                                 if(j+z<1); z2 = z+1; end
    %                                 if(j+z>sz(2)); z2 = z-1; end
    % %                                 k2 = mod(i+k-1, sz(1))+1;
    % %                                 z2 = mod(j+z-1, sz(2))+1;
    %                                 
                                    neighborhood(2+k, 2+z) = (sample(i+k,j+z)>=sample(i,j)); % tresholding the neighbors' values with the pixel's one
                                end
                            end
                              LBP_image(i,j) = neighborhood(1,1)*2^7 + neighborhood(1,2)*2^6 + neighborhood(1,3)*2^5 + neighborhood(2,3)*2^4 + neighborhood(3,3)*2^3 + neighborhood(3,2)*2^2 + neighborhood(3,1)*2 + neighborhood(2,1)*1; % combining all the neighbors' values into a binary number
    %                           LBP_image(i,j) = neighborhood(2,1)*2^7 + neighborhood(3,1)*2^6 + neighborhood(3,2)*2^5 + neighborhood(3,3)*2^4 + neighborhood(2,3)*2^3 + neighborhood(1,3)*2^2 + neighborhood(1,2)*2 + neighborhood(1,1)*1; % combining all the neighbors' values into a binary number
    %                         binary = [neighborhood(2,1) neighborhood(3,1) neighborhood(3,2) neighborhood(3,3) neighborhood(2,3) neighborhood(1,3) neighborhood(1,2) neighborhood(1,1)]; % combining all the neighbors' values into a binary number
    %                         LBP_image(i,j) = bi2de(binary, 'right-msb');  % converting the binary value into a decimal one
                        end
                    end
                    LBP_image = uint8(LBP_image(2:sz(1)-1, 2:sz(2)-1)); % caster l'image en uint8 pour que ce soit bien entre 0 et 255, et non pas entre 0 et 1

            % Computation of the histogram of the sample
                    global_histogram = zeros(nbins, m);   % list of vectors (one vector representing a local (regional) histogram)
                    for i=1:windows    % parcours de toutes les régions de l'image
                        for j=1:windows
                            if(i==windows && j==windows)
                                lbpimage = LBP_image((i-1)*pixels_window_v+1:length(LBP_image(:,1)), (j-1)*pixels_window_h+1:length(LBP_image(1,:)));
                                c = histo(lbpimage, nbins, uniform_patterns);                        
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

                    % Comparison of this histogram with the model for faces and non-faces
                    dist_face = chisq(global_histogram, Mface, m, windows, firstBar, lastBar, weights);
                    dist_nonface = chisq(global_histogram, Mnonface, m,  windows, firstBar, lastBar, weights);
                    seuil = 10;
                    if(dist_face<dist_nonface && (dist_nonface-dist_face)/(dist_nonface+dist_face)*100>seuil)
                        xs=[xs x]; ys=[ys y]; widths=[widths width]; heights=[heights height]; distances=[distances dist_face];
                        if(dist_face<mindist)
                            mindist = dist_face;
                            minx = x; minwidth = width; miny = y; minheight = height;
                        end
%                         figure;
%                         imshow(sample);
                    end
%                 end
                y = y + height/2;
            end
            x = x + width/3;
        end
    end

    
% Regroupage des rectangles pour n'avoir qu'un seul rectangle qui encadre
% les yeux - on commence le regroupement par le rectangle avec la plus
% petite distance - si un rectangle n'appartient pas à ce groupe, alors il
% ne représente probablement pas des yeux (un seul visage n'a qu'une seule
% paire d'yeux), donc on l'ignore
    div = 5; % 2 rectangles sont mis dans le même groupe s'ils sont distants de moins de width/div
    mult = 1; % pondération des coordonnées d'un rectangle en mult*1/distance
    [distances, ind] = sort(distances); % tri du vecteur contenant les distances (ordre croissant)
    xs = xs(ind);
    ys = ys(ind);
    widths = widths(ind);
    heights = heights(ind);
    % tri de xs, ys, widths et heights selon distances croissantes
    % pour que le parcours de la liste des rectangles pour les regrouper
    % commence par celui qui a la plus petite distance,
    % et ensuite on y a ajoute les autres rectangles, aux distances plus élevées
    if(isempty(xs)==false) % si des visages ont été détectés
        eyesRect = [xs(1) xs(1)+widths(1) ys(1) ys(1)+heights(1)]; % vecteur qui va contenir un rectangle, caractérisé par x, y, width et height,
                                                   % il contient d'abord le rectangle à la distance minimale
        for i=2:length(xs) % pour chaque rectangle contenant un visage détecté
            if(((xs(i)>=(eyesRect(1)-(widths(i)/div)) && xs(i)<=(eyesRect(1)+(widths(i)/div))) || ((xs(i)+widths(i))>=(eyesRect(2)-(widths(i)/div)) && (xs(i)+widths(i))<=(eyesRect(2)+(widths(i)/div)))) && ((ys(i)>=eyesRect(3)-(heights(i)/div) && ys(i)<=eyesRect(3)+(heights(i)/div)) || (ys(i)+heights(i)>=eyesRect(4)-(heights(i)/div) && ys(i)+heights(i)<=eyesRect(4)+(heights(i)/div)))) % if rectangle i belongs to eyesRect
                % add rectangle i to eyesRect
                % pondération inversement proportionnelle à la distance au modèle de visage
                % avec un paramètre mult que l'ont peut changer
                % pondération en fonction de la taille
                param_dist = true;
                param_dist = param_dist*(mult/distances(i)-1)+1;
                eyesRect(1) = (eyesRect(1)+(xs(i)*param_dist))/(1+param_dist);
                eyesRect(2) = (eyesRect(2)+((xs(i)+widths(i))*param_dist))/(1+param_dist);
                eyesRect(3) = (eyesRect(3)+(ys(i)*param_dist))/(1+param_dist);
                eyesRect(4) = (eyesRect(4)+((ys(i)+heights(i)))*param_dist)/(1+param_dist);
                break;
            end
        end
        xs = eyesRect(1);
        widths = eyesRect(2)-eyesRect(1);
        ys = eyesRect(3);
        heights = eyesRect(4)-eyesRect(3);
    end
    
%     if(isempty(xs)==false)
%         xs(length(xs)+1)=minx; ys(length(ys)+1)=miny; widths(length(widths)+1)=minwidth; heights(length(heights)+1)=minheight;
%     end
            
    disp("LBP is done");
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
    for i=1:(sz(2))-1
       if(array(i)~= array(i+1))
           nbTransitions = nbTransitions + 1;
       end
    end
    bool = (nbTransitions <= 2);
end

%% Chi Square distance calculation
function dist = chisq(S, M, m, windows, firstBar, lastBar, weights)
    dist = 0;    
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

%% Sliding-neighborhood LBP function
function res = LBP(neighborhood)
    new_neighborhood = zeros(3,3);
    for j=1:3  % parcours des voisins du pixel (i,j)
        for z=1:3
            new_neighborhood(j,z) = (neighborhood(j,z)>=neighborhood(2,2)); % tresholding the neighbors' values with the pixel's one
        end
    end
    res = new_neighborhood(1,1)*2^7 + new_neighborhood(1,2)*2^6 + new_neighborhood(1,3)*2^5 + new_neighborhood(2,3)*2^4 + new_neighborhood(3,3)*2^3 + new_neighborhood(3,2)*2^2 + new_neighborhood(3,1)*2 + new_neighborhood(2,1)*1; % combining all the neighbors' values into a binary number
%     res = new_neighborhood(2,1)*2^7 + new_neighborhood(3,1)*2^6 + new_neighborhood(3,2)*2^5 + new_neighborhood(3,3)*2^4 + new_neighborhood(2,3)*2^3 + new_neighborhood(1,3)*2^2 + new_neighborhood(1,2)*2 + new_neighborhood(1,1)*1; % combining all the neighbors' values into a binary number
end