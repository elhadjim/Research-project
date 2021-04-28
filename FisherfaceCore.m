function [m_data ,eigenfaces ,fisherfaces ,ProjectedImg_Fisher] = FisherfaceCore(M)
% Utilisation de la méthode PCA et FLD.
% 
%
% Description: La fonction prend en argument la matrice 2D qui contient toutes les images colonnes 1D.
% On suppose Mi une image de la base d'apprentissage.
% Soit P le nombre total d'images de la base d'apprentissage de dimension MxN
% Soit C le nombre de classes.
% On projette d'abord les images dans le sous espace linéaire de dimension (P-C)

% La matrice de transfert est: Xi = eigenfaces' * (Mi - m_data). ( cf.Pentland)

% Ensuite Xi est convertie en Yi en projettant les images dans un autre espace de dimension plus petite (C-1) afin 
% afin d'effectuer une bonne classification des classes.

% La matrice de transfert est : Yi = fisherfaces' * Xi = fisherfaces' * eigenfaces' * (Mi - m_data)
%
% Argument de la fonction :      M                          - La matrice 2D de dimensions MNxP.
% 
% Returns:                       m_data                     - (M*Nx1) moyenne de la base d'apprentissage.
%                                eigenfaces                 - (M*Nx(P-C)) Les vecteurs propres de la matrice de covariance des images 
%                                                              de la base d'apprentissage ( Eigen vectors).
%                                fisherfaces                - ((P-C)x(C-1)) Les (C-1) vecteurs propres les plus grands de la matrice T = Sb/Sw.
%                                ProjectedImg_Fisher        - ((C-1)xP) Les images d'apprentissages qui sont projetées dans le sous espace linéaire de Fisher.
%                 

imgcount = size(M,2);
Class_number = imgcount/10; % Nombre de personnes ou de classes
Class_population = 10; % Nombre d'images dans chaque classe
P = Class_population * Class_number; % Le nombre total des images de la base d'apprentissage

%%%%%%%%%%%%%%%%%%%%%%%% Calcul de l'image moyenne

m_data = mean(M,2); 

%%%%%%%%%%%%%%%%%%%%%%%% Calcul de l'écart de chaque image avec l'image moyenne

A = [];
for i=1 : imgcount
    temp = double(M(:,i)) - m_data;
    A = [A temp];
end

%%%%%%%%%%%%%%%%%%%%%%%% La dimension de la matrice de covariance étant très grande (cf.Pentland)
%                        On utilise donc la matrice L = A'*A. ( C= A*A')
%                        Dim ( C ) = (MxN)²
%			             Dim ( L ) = P²
%			             En général le nombre d'image de la base d'apprentissage P est souvent très petit que le nombre de pixel de l'image 
%		                 ==> P<MxN. 
%			             De plus, on sait que le nombre maximal de valeurs propres non nulles d'une matrice de dimension PxQ est min(P-1,Q-1)
%			             Donc en utilisant la transposée de la matrice C qui est la matrice L, il y'aura seulement P-1 vecteurs propres significatifs plutôt que MxN.
%			             Les vecteurs propres restant seront donc assignés et associés à des valeurs nulles.

%%%%%%%%%%%%%%%%%%%%%%%% Méthode PCA

L = A'*A;       % La transposée de la matrice de covariance
[V D] = eig(L); % L et C ont les mêmes valeurs propres. Ce sont les éléments de la diagonale de la matrice D
		        % La matrice V contient les vecteurs propres de la matrice L.

%%%%%%%%%%%%%%%%%%%%%%%% Elimination de quelques valeurs propres.
L_eig_vec = [];
for i = 1 : P-Class_number %% Réduction de l'espace ( cf. Belhumeur)
			               %C'est cette réduction en utilisant la méthode PCA qui permet d'éviter le problème 
                           % de singularité de la matrice de variance intra-classes  Sw                  
    L_eig_vec = [L_eig_vec V(:,i)];
end

%%%%%%%%%%%%%%%%%%%%%%%% Calcul des vecteurs propres de la matrice de covariance'C'

eigenfaces = A * L_eig_vec; %% (cf. Pentland ) Les vecteurs propres de la matrices de covariance C peuvent être calculés à partir des vecteurs
                            % propres de L.

%%%%%%%%%%%%%%%%%%%%%%%% Projection des images centrées dans le "facespace"
% Xi = eigenfaces' * (Mi-m_database)

ProjectedImg_PCA = [];
for i = 1 : P
    temp = eigenfaces'*A(:,i);
    ProjectedImg_PCA = [ProjectedImg_PCA temp]; 
end

%%%%%%%%%%%%%%%%%%%%%%%% Calcul de la moyenne de chaque classe dans
%%%%%%%%%%%%%%%%%%%%%%%% l'eigenspace

mean_PCA = mean(ProjectedImg_PCA,2);          % Moyenne dans l'eigenspace (facespace)
m = zeros(P-Class_number,Class_number);    % Initialisation de la matrice contenant la moyenne de chaque classe 
Sw = zeros(P-Class_number,P-Class_number); % Initialisation de la matrice de variance intra-classe
Sb = zeros(P-Class_number,P-Class_number); % Initialisation de la matrice de variance inter-classe

for i = 1 : Class_number
    m(:,i) = mean( ( ProjectedImg_PCA(:,((i-1)*Class_population+1):i*Class_population) ), 2 )';    
    
    S  = zeros(P-Class_number,P-Class_number); 
    for j = ( (i-1)*Class_population+1 ) : ( i*Class_population )
        S = S + (ProjectedImg_PCA(:,j)-m(:,i))*(ProjectedImg_PCA(:,j)-m(:,i))'; % (cf.Belhumeur)
    end
    
    Sw = Sw + S; % Matrice de variance intra-classe
    Sb = Sb +(m(:,i)-mean_PCA) * (m(:,i)-mean_PCA)'; % Matrice de variance inter-classe (cf.Belhumeur)
end

%%%%%%%%%%%%%%%%%%%%%%%% Méthode FLD
% On souhaite maximiser la matrice de variance inter-classes et minimiser la matrice de variance intra-classe.
% On introduit la matrice suivante: T = Sb/Sw
% La solution à ce problème est donnée par la décomposition généralisée des
% valeurs propres : Sb*V= Sw*V*D  soit T*V = V*D où D est la matrice diagonale des valeurs
% propres généralisées et V est la matrice contenant les vecteurs propres
% associés à chaque valeur propre.

[T_eig_vec, ~] = eig(Sb,Sw); 
T_eig_vec = fliplr(T_eig_vec); %% Tri décroissant pour ne garder que les plus grands vecteurs propres 

%%%%%%%%%%%%%%%%%%%%%%%% On élimine les valeurs propres nulles en réduisant la dimension ( cf.Belhumeur )

for i = 1 : Class_number-1 
    fisherfaces(:,i) = T_eig_vec(:,i); % On ne garde que (C-1) vecteurs propres les plus grands de la matrice T
end

%%%%%%%%%%%%%%%%%%%%%%%% Projection des images dans le sous espace linéaire de Fisher
% Yi = fisherfaces' * eigenfaces' * (Mi - m_data) 

for i = 1 : P
    ProjectedImg_Fisher(:,i) = fisherfaces' * ProjectedImg_PCA(:,i);
end
