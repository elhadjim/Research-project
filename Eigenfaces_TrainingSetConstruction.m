close all;

treshold_analysis = false; % pour voir les boxplots des distances pour pouvoir déterminer un seuil
human_treshold_analysis = false; % pour voir les boxplots des distances à des visages/non visages pour déterminer un seuil entre visage humain et visage non humain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill the following lines every time you change database
basedirectory = 'BaseGroupePIR/';
imagetype = '**/*.gif';
I = 5; % number of individuals in the base
P = 20; % number of pictures of each individual
treshold = 0.17;
facetreshold = 9700;
Mprime = 10;    % number of eigenfaces that construct the face space
P_apprentissage = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];    % pictures of each individual used for the training set
P_test = [];   % pictures of each individual used for testing
% P_apprentissage = [1 3 5 7 9 11 13 15 17 19];    % pictures of each individual used for the training set
% P_test = [2 4 6 8 10 12 14 16 18 20];   % pictures of each individual used for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construction of the training set (as a list of vectors under the form of a matrix)
% Base that is used to construct the face space

imagefiles = dir(fullfile(basedirectory, imagetype));
M = length(imagefiles); % number of images in the base
M_apprentissage = length(P_apprentissage) * I;  % number of pictures used for the training set
M_test = length(P_test)*I;  % number of pictures used for testing
sz = size(imread([imagefiles(1).folder '\' imagefiles(1).name]));   % dimensions of the image
N = sz(1)*sz(2);    % dimension of an image vector
gamma = zeros(N, M_apprentissage); % matrix that will contain the M_apprentissage column image vectors
k=1;
for i=1:M
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
            gamma(:,k) = graycurrentimage(:);
        else
            gamma(:,k) = currentimage(:);
        end
        k=k+1;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    Construction of the eigenfaces                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_psi = (mean(gamma'))';  % average face (column vector of size N)
% imshow(mat2gray(vec2mat(my_psi, sz(1))));    % to display the average face
A = gamma - (my_psi*ones(1,M_apprentissage));    % matrix A=[omega1 omega2 omega3 ... omegaM]
L = A'*A;  % matrix L=A^T.A
[V,D] = eig(L); % V: matrix of eigenvectors (of dimension Nx1), D: diagonal matrix of associated eigenvalues
[d,ind] = sort(diag(D));    % d: column vector containing the eigenvalues sorted in ascending order
                            % ind: column vector containing the corresponding indices
Vs = V(:,ind);  % Vs: matrix containing the eigenvectors but sorted according to the order of the sorted eigenvalues
m10 = zeros(1,Mprime);
for i=1:Mprime
    m10(1,i) = size(V,2)+1-i;   % indices of the last Mprime columns of Vs
end
Vsb = Vs(:,m10);   % Vsb contains the last Mprime columns of Vs so the Mprime eigenvectors with highest associated eigenvalues
Ds = D(:,ind);
Dsb = Ds(:,m10);    % Dsb contains the Mprime highest eigenvalues
dsb = d(m10',:);     % dsb contains the Mprime highest eigenvalues
eigenfacesmatrix = zeros(N,Mprime); % matrix of eigenfaces
% figure;
for i=1:Mprime
    eigenfacesmatrix(:,i) = A*Vsb(:,i);
    eigenfacesmatrix(:,i) = eigenfacesmatrix(:,i)/norm(eigenfacesmatrix);   % normalization of the eigenfaces
%     subplot(1,10,i);
%     imshow(imrotate(mat2gray(vec2mat(eigenfacesmatrix(:,i),sz(1))),-90)); % to display the eigenfaces
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  Construction of the face classes                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patternvectors = zeros(M_apprentissage,Mprime);   % matrix of M_apprentissage pattern vectors - each pattern vector is a line vector of size Mprime
for i=1:M_apprentissage
    patternvectors(i,:) = (eigenfacesmatrix'*(gamma(:,i)-my_psi))';
end
classvectors = zeros(I,Mprime); % matrix of I class vectors - each class vector is a mean of length(P_apprentissage) pattern vectors
% figure;
for i=1:I
    if(length(P_apprentissage)==1)
        classvectors(i,:) = patternvectors(i,:);
    else
        classvectors(i,:) = mean(patternvectors(length(P_apprentissage)*(i-1)+1:length(P_apprentissage)*i,:));
    end
%     subplot(1,I,i);
%     imshow(imrotate(mat2gray(vec2mat((classvectors(i,:)*(eigenfacesmatrix'))',sz(1))),-90)); % to display the eigenfaces
%     title(i);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analysis of the sensibility and specificity of the recognition test %%%
%%%                      To define a treshold                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each image of the base, calculate its distance to every face class
% -> matrix of size M_test*I
% for the distance to the right face class: (hypothèse vérifiée)
% - if it is the smallest distance and it is under the treshold -> vrai positif
% - otherwise -> faux négatif
% for the distances to all the other face classes: (hypothèse non vérifiée)
% - if it is not the smallest distance or it is over the treshold -> vrai négatif
% - if it is the smallest distance and it is under the treshold -> faux positif

if(treshold_analysis)
    M_test = length(P_test)*I;  % number of pictures used for testing
    gamma_test = zeros(N, M_test); % matrix that will contain the M_test column image vectors
    k=1;
    for i=1:M
        bool = false;
        for j=1:length(P_test)
            if(rem(i-P_test(j),P)==0)
                bool=true;
                break;
            end
        end
        if bool
            currentimage = imread([imagefiles(i).folder '\' imagefiles(i).name]);
            if(strcmp(imagetype,'**/*.jpg'))
                graycurrentimage = rgb2gray(currentimage);
                gamma_test(:,k) = graycurrentimage(:);
            else
                gamma_test(:,k) = currentimage(:);
            end
            k=k+1;
        end
    end

    vp=0;   %vrais positifs
    fp=0;   %faux positifs
    vn=0;   %vrais négatifs
    fn=0;   %faux négatifs
    total=0;
    r=0;    %bonne personne
    w=0;    %mauvaise personne

    distances = zeros(I,M_test);   % matrix of M_test columns vectors, each containing the distances to the I face classes
    rightdistances = zeros(1,M_test);   % for each image, distance to the face class of the person it corresponds to
    wrongdistances = zeros(I,M_test);  % for each image, distances to all the other face classes
    indices = zeros(I,M_test);
    
    patternvectors_test = zeros(M_test,Mprime);   % matrix of M_apprentissage pattern vectors - each pattern vector is a line vector of size Mprime
    for i=1:M_test
        patternvectors_test(i,:) = (eigenfacesmatrix'*(gamma_test(:,i)-my_psi))';
        for j=1:I
    %         distances(j,i) = norm(patternvectors(i,:) - classvectors(j,:));   % euclidian distance
            distances(j,i) = sum((patternvectors_test(i,:)-classvectors(j,:)).^2./dsb'); % mahalanobis distance
        end
        [~,indices(i)] = min(distances(:,i));
    end
    
    for i=1:M_test
        for j=1:I
            total=total+1;
            if (j==fix((i-1)/(length(P_test)))+1) % right person
                rightdistances(i) = distances(j,i);
                wrongdistances(j,i)=Inf;
                r = r+1;
                if((distances(j,i)==(min(distances(:,i)))) && (distances(j,i)<treshold))    % vrai positif
                    vp=vp+1;
                elseif((distances(j,i)~=(min(distances(:,i)))) && (distances(j,i)<treshold))   % faux négatif
                    fp=fp+1;
                else
                    fn=fn+1;
                end
            else   % wrong person
                wrongdistances(j,i) = distances(j,i);
%                 w=w+1;
%                 if((distances(j,i)==(min(distances(:,i)))) && (distances(j,i)<treshold))    % faux positif
%                     fp=fp+1;
%                 else   % vrai négatif
%                     vn=vn+1;
%                 end
            end
        end
    end
    [mineps, ind] = min(distances);
    % disp(mineps);
    figure;
    hold on;
    ax(1) = subplot(1,2,1);
    boxplot(rightdistances);
    title('Distances à la bonne classe');
    ax(2) = subplot(1,2,2);
    boxplot(wrongdistances(:));
    title('Distances aux mauvaises classes');
    linkaxes(ax,'y');
    % set(ax,'ylim',[2 5.5e8]);
%     disp(['vrais positifs: ' num2str(vp)]);
%     disp(['faux négatifs: ' num2str(fn)]);
%     disp(['vrais négatifs: ' num2str(vn)]);
%     disp(['faux positifs: ' num2str(fp)]);
%     disp(['sensibilité: ' num2str(vp/(vp+fn)*100) '%']);
%     disp(['spécificité: ' num2str(vn/(vn+fp)*100) '%']);
    disp(['Pourcentage de reconnaissance : ' num2str(vp/M_test*100) '%']);
    disp(['Pourcentage de non-reconnaissances erronées : ' num2str(fn/M_test*100) '%']);
    disp(['Pourcentage de reconnaissances erronées : ' num2str(fp/M_test*100) '%']);
%     disp(['Poucentage de non-reconnaissances correctes : ' num2str(vn/w*100) '%']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Analysis of distance to the face space to define a treshold     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(human_treshold_analysis)
    basedirectory3 = 'FacesVsNonFaces/';
    imagetype3 = '**/*.gif';
    imagefiles3 = dir(fullfile(basedirectory3, imagetype3));
    M3 = length(imagefiles3); % number of images in the base
    H3 = 0; %number of pictures of humans
    NH3 = 0; %number of pictures of non humans
    gamma3 = zeros(N, M3); % matrix that will contain the M3 column image vectors
    for i=1:M3
        currentimage = imread([imagefiles3(i).folder '\' imagefiles3(i).name]);
        if(strcmp(imagetype3,'**/*.jpg'))
            graycurrentimage = rgb2gray(currentimage);
            gamma3(:,i) = graycurrentimage(:);
        else
            gamma3(:,i) = currentimage(:);
        end
        if(contains(imagefiles3(i).name,'human'))
            H3=H3+1;
        else
            NH3=NH3+1;
        end
    end

    patternvectors3 = zeros(M3,Mprime);   % matrix of M3 line vectors of size Mprime
    epsilon3 = zeros(M3,1); % vector of the distances between each image (mean-adjusted) and its projection onto face space
    humandistances = zeros(1,H3);
    nonhumandistances = zeros(1,NH3);
    j=1;
    k=1;
    for i=1:M3
        patternvectors3(i,:) = (eigenfacesmatrix'*(gamma3(:,i)-my_psi))';
        epsilon3(i) = norm((gamma3(:,i)-my_psi)-(patternvectors3(i,:)*(eigenfacesmatrix'))');
        if(contains(imagefiles3(i).name,'human'))
            humandistances(j) = epsilon3(i);
            j=j+1;
        else
%             figure;
%             imshow(imread([basedirectory3 imagefiles3(i).name]), [0 255]);
%             disp(epsilon3(i));
            nonhumandistances(k) = epsilon3(i);
            k=k+1;
        end
    end
    figure;
    hold on;
    ax3(1) = subplot(1,2,1);
    boxplot(humandistances);
    ax3(2) = subplot(1,2,2);
    boxplot(nonhumandistances);
    linkaxes(ax3,'y');
end

%% Saving the variables
save('Eigenfaces_training', 'eigenfacesmatrix', 'my_psi', 'I', 'P', 'classvectors', 'treshold', 'imagefiles', 'facetreshold', 'dsb');
