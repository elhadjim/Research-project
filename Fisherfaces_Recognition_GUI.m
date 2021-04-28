function Recognized_img = Fisherfaces_Recognition_GUI(TestImage)
% L'étape de la reconnaissance 
%
% Description: La fonction compare deux visages en projettant les images dans le facespace
% tout en mesurant la distance euclidienne entre elles.
%
% Arguments de la fonction :      TestImage              - image test
%
%                                 m_data                 - (M*Nx1) moyenne de l'image d'apprentissage
%                                         
%
%                                 eigenfaces             - (M*Nx(P-1)) vecteurs propres de la matrice de covariance 
%                                                           de la base d'apprentissage

%                                 fisherfaces            - ((P-1)x(C-1))  (C-1) vecteurs propres les plus grands de la matrice T = Sb/Sw 

%                                 ProjectedImg_Fisher    - ((C-1)xP) images d'apprentissage qui sont projettées dans le sous espace linéaire de Fisher
%                                         
% 
% Returns:                        Recognized_img          - Nom de l'image reconnue dans la base d'apprentissage.
                 
disp("Fisherfaces is launched");
load('Fisherfaces_training');
%%%%%%%%%%%%%%%%%%%%%%%% Récupération et Projection de l'Image test
% Image = imread(TestImage);%Lecture de l'image test
Image = TestImage;
% figure;
% imshow(Image);
[row ,col] = size(Image); % Récupération des dimensions
Reshaped_Image = reshape(Image',row*col,1); % Transformation de l'image en un vecteur colonne
ProjectedTestImage = fisherfaces' * eigenfaces' * (double(Reshaped_Image)-m_data); % Projection de l'image test dans le sous espace linéaire de Fisher

%%%%%%%%%%%%%%%%%%%%%%%% Calcul des distances euclidiennes : 
% On calcule les distances euclidiennes entre l'image test projettée et la projection de
% toutes les images d'apprentissage centrées.
% L'image test doit avoir une distance euclidienne minimale avec l'image qui lui correspond
% dans la base d'apprentissage.

Euclide_dist = [];
imgcount = size(ProjectedImg_Fisher,2);
for i = 1 : imgcount
    col = ProjectedImg_Fisher(:,i);
    temp = ( norm( ProjectedTestImage - col ) )^2; % calcul de la distance euclidienne entre l'image test projettée et l'image de la base d'apprentissage
    Euclide_dist = [Euclide_dist temp];
end

[~ , Recognized_index] = min(Euclide_dist); % On récupère l'indexe de l'image correspondant au minimum
Recognized_img = strsplit(imagefiles(Recognized_index).name, {'-','.gif','1','2','3','3','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','0'}, 'CollapseDelimiters', true);
Recognized_img = Recognized_img(1);
% Recognized_img = strcat(int2str(Recognized_index),'.gif');
disp("Fisherfaces is done");
