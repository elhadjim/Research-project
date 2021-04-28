% On a un ensemble d'images de la base d'apprentissage.
%
% Description : La fonction transforme les images 2D en un vecteur colonne 1D.
% Ensuite, elle met chaque vecteur colonne 1D dans une colonne pour  
% construire la matrice 2D 'M'.
% La dimension de chaque image de la base d'apprentissage est MxN.
% Soit P le nombre total d'images dans la base d'apprentissage.
% Soit C le nombre de classes. 
%                 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill the following lines every time you change database
basedirectory = 'BaseGroupePIR/';
imagetype = '**/*.gif';
I = 5; % number of individuals in the base
P = 20; % number of pictures of each individual
P_apprentissage = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];    % pictures of each individual used for the training set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagefiles = dir(fullfile(basedirectory, imagetype));
M_total = length(imagefiles); % number of images in the base
M_apprentissage = length(P_apprentissage) * I;  % number of pictures used for the training set
sz = size(imread([imagefiles(1).folder '\' imagefiles(1).name]));   % dimensions of the image
M = [];
k=1;
for i=1:M_total
    bool = false;
    for j=1:length(P_apprentissage)
        if(rem(i-P_apprentissage(j),P)==0)
            bool=true;
            break;
        end
    end
    if bool
        currentimage = imread([imagefiles(i).folder '\' imagefiles(i).name]);
        if(strcmp(imagetype,'**/*.jpg'))
            graycurrentimage = rgb2gray(currentimage);
            temp = reshape(graycurrentimage',sz(1)*sz(2),1);   % Transformation des images 2D en vecteurs colonne 1D
        else
            temp = reshape(currentimage',sz(1)*sz(2),1);   % Transformation des images 2D en vecteurs colonne 1D
        end
        M = [M temp];  % Création de la matrice 
    end
end


% %%%%%%%%%%%%%%%%%%%%%%%% Gestion de fichier (Amir Omidvarnia)
% D = dir(DataPath);
% nbimg = 0; % Initialisation de la variable qui comportera le nombre d'image dans la base d'apprentissage.
% 
% for i = 1:size(D,1)
%     if not(strcmp(D(i).name,'.')|strcmp(D(i).name,'..')|strcmp(D(i).name,'Thumbs.db'))
%         nbimg = nbimg + 1; % Nombre d'images dans la base 
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Construction de la matrice 2D 
% M = [];
% for i = 1 : nbimg
%     str = strcat(DataPath,'/',int2str(i),'.gif'); % chemin vers chaque image de la base d'apprentissage
%     img = imread(str); % lecture de l'image
%     [row ,col] = size(img); % récupération des dimensions 
%     temp = reshape(img',row*col,1);   % Transformation des images 2D en vecteurs colonne 1D
%     M = [M temp];  % Création de la matrice                 
% end

M = double(M);

[m_data ,eigenfaces ,fisherfaces ,ProjectedImg_Fisher] = FisherfaceCore(M);

%% Saving the variables
save('Fisherfaces_training', 'm_data', 'imagefiles', 'P', 'eigenfaces', 'fisherfaces', 'ProjectedImg_Fisher');
