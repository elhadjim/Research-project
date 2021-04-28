function [RGB,RGB2,matchedPoints1,matchedPoints2] = SIFT_GUI(guiimage1, guiimage2)
disp("SIFT is launched");
% clear all;
% close all;
%1ère image
image = guiimage1;
% figure;
% imshow(guiimage1);
% nameImage = 'face.jpeg';
% image=imread(nameImage);

[~,~,deepth]=size(image);
if deepth==3 %Si l'image est en couleur, on la transforme en gris
    image=rgb2gray(image);
end
image=im2double(image); %Transformer l'image en matrice de double
original=image;

%initialisation des octaves (matrices vides pour l'instant)
octave1=[];
octave2=[];
octave3=[];
tic;

%% Octave 1

image_temporaire=imresize(image,2,'nearest');         %Double la taille de l'image initiale pour de meilleurs resultats
[hight1,width1,~]=size(image_temporaire);                       %Récupération des dimensions de l'image

k2=0;                                       %1ere octave
image_temporaire(hight1:hight1+4,width1:width1+4)=0;  %Zero padding
clear c;

for k1=0:3                                  %Differentes valeurs de sigma pour 1 octave 
    k=sqrt(2);
    sigma=(k^(k1+(2*k2)))*1.6;

    for x=-2:2                              %Generation du filtre gaussien
        for y=-2:2
             h(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma))); 
        end
    end

    for i=1:hight1                          %Convolution image et filtre
        for j=1:width1
            temp1=image_temporaire(i:i+4,j:j+4)'.*h; % Produit de 2 matrices terme à terme
            conv_1(i,j)=sum(sum(temp1));   % 2 fois car sum fonctionne sur un vecteur
        end
    end

    octave1=[octave1 conv_1];               %Ajout de l'image filtree a l'octave
end
%figure(1);        
%imshow(octave1);

%% Octave 2

%clear I_temp conv_1 h;
image_temporaire2=original;
[hight2,width2,~]=size(image_temporaire2);

k2=1; %deuxième octave
image_temporaire2(hight2:hight2+4,width2:width2+4)=0; %zéro padding
clear c;

for k1=0:3
    k=sqrt(2);
    sigma=(k^(k1+(2*k2)))*1.6;

    for x=-2:2
        for y=-2:2
            h1(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma)));
        end
    end

    for i=1:hight2
        for j=1:width2
            temp2=image_temporaire2(i:i+4,j:j+4)'.*h1;       
            conv_2(i,j)=sum(sum(temp2));
        end
    end

    octave2=[octave2 conv_2];                     
end
%figure(2);
%imshow(octave2);

%% Octave 3
clear I_temp2 conv2 h1;
image_temporaire3=imresize(original,0.5,'nearest'); %taille divisée par 2      
[hight3,width3,~]=size(image_temporaire3);

k2=2; %3ème octave
image_temporaire3(hight3:hight3+4,width3:width3+4)=0;
clear c;

for k1=0:3
    k=sqrt(2);
    sigma=(k^(k1+(2*k2)));

    for x=-2:2
        for y=-2:2
            h2(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma)));
        end
    end


     for i=1:hight3
         for j=1:width3
             temp3=image_temporaire3(i:i+4,j:j+4)'.*h2;   
             conv_3(i,j)=sum(sum(temp3));
         end
     end
     octave3=[octave3 conv_3];
end
%figure (3);
%imshow(octave3);

fprintf('\n Temps pour construction des octaves %.3f s\n',toc) ;
% SORTIE :  3 matrices representants chacune une octave


%% Differences de gaussiennes
DoG1=[];
DoG2=[];
DoG3=[];

for k=1:3
    DoG1=[DoG1 octave1(1:hight1,(k-1)*width1+1:k*width1)-octave1(1:hight1,k*width1+1:(k+1)*width1)];
    DoG2=[DoG2 octave2(1:hight2,(k-1)*width2+1:k*width2)-octave2(1:hight2,k*width2+1:(k+1)*width2)];
    DoG3=[DoG3 octave3(1:hight3,(k-1)*width3+1:k*width3)-octave3(1:hight3,k*width3+1:(k+1)*width3)];
end

%figure(4);
%imshow(DoG1);
%figure(5);
%imshow(DoG2);
%figure(6);
%imshow(DoG3);


% SORTIE :  3 matrices regroupants les 3 resultats des differences de
% gaussiennes pour chaque octave 

%% Localisation des extremas

Localisations=[];

for w=width1+4:2:2*width1-4% supprimer les :2: pour petites images?
    for h=6:2:hight1-4 %on se balade sur l'image, on commence a 6 pour eviter les pb de borne
        HighTab=DoG1(h-3:h+3,w+width1-3:w+width1+3);

        HtmpM=int32(h/2);%On trouve la position correspondante sur l'octave 2 
        WtmpM=int32(w/2)+1;
        MiddleTab=DoG2(HtmpM-2:HtmpM+2,WtmpM-2:WtmpM+2);

        HtmpL=int32(h/4);
        WtmpL=int32(w/4)+1;
        LowTab=DoG3(HtmpL-1:HtmpL+1,WtmpL-width3-1:WtmpL-width3+1);


        maxMT=max(max(MiddleTab));%On cherche les extremums pour chaque octave
        maxLT=max(max(LowTab));
        maxHT=max(max(HighTab));
        minMT=min(min(MiddleTab));
        minLT=min(min(LowTab));
        minHT=min(min(HighTab));


        a=sum(sum(MiddleTab(:,:)==maxMT));%Compte le nombre de fois que maxMT est présent dans la région 5x5 autour du point
        b=sum(sum(MiddleTab(:,:)==minMT));

        maximum = max([maxMT,maxLT,maxHT]);
        minimum = min([minMT,minLT,minHT]);
        %Le pixel courant et soit un max soit un min, et il est le
        %seul de sa région à chaque fois
        if(a==1 && DoG2(HtmpM,WtmpM)==maximum) ||(b==1 && DoG2(HtmpM,WtmpM)==minimum)
            if HtmpM>2 && WtmpM>width2+2 && HtmpM<hight2-2 && WtmpM<2*width2-2
                Localisations = [Localisations HtmpM WtmpM-width2];
            end
        end
    end
end

% SORTIE :  1 vecteur contenant les coordonnées de chaque point d'interet


%% Ecrémation des PI par des seuils en ne selectionnant que les coins
[Gx,Gy]=imgradientxy(DoG1(:,width1+1:2*width1));%(diff_12(:,:));% retour la magnitude des gradient selon x et Y
[Gmag, Gdir]=imgradient(DoG1(:,width1+1:2*width1));

[~, totalPI]=size(Localisations);
newLocalisations=[];
j=1;
for i=2:2:totalPI
    if (abs(Gx(Localisations(1,i-1),Localisations(1,i)))>=0.0015 && abs(Gy(Localisations(1,i-1),Localisations(1,i)))>=0.0015)
        newLocalisations=[newLocalisations j Localisations(1,i-1) Localisations(1,i)];%Adapter le seuil pour qu'il ne me reste plus que X KP
        j=j+1;
    end
end
 % SORTIE : % SORTIE :  1 vecteur contenant les coordonnées de chaque point
 % d'interet triés


%% Orientation des points d'intérêt

[~, totalPI]=size(newLocalisations);
keypoints=[];% stocke mag et dir gradient pour chaque pt
kp=0;
for w=3:3:totalPI % parcours de 3 en 3 car les points d'intérêts occupent 3 places (n°, coordX, coordY)
    kp=kp+1;
    histogram = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    %pour tous les points d'intérêt : permet d'accéder à une zone 5*5
    %autour du point d'intérêt
    %pourMoicoordX=newLocalisations(1,w-1);
    %pourMoicoordY=newLocalisations(1,w);
    for i=newLocalisations(1,w)-2:newLocalisations(1,w)+2%width : coordonnée Y
        for j=newLocalisations(1,w-1)-2:newLocalisations(1,w-1)+2%hight : coordonnée X
            %pourMoiGdir=Gdir(j,i);
            %pourMoiFloor=floor((Gdir(j,i)+180)/20)+1; % pour savoir quelle colonne de l'histogramme remplir : on ajoute 180° pour se ramener à l'intervalle [0;360] puis on divise par 20 pour former des paquets de 20
            %pourMoiGmag=Gmag(j,i)
            histogram(1,floor((Gdir(j,i)+180)/20)+1) = histogram(1,floor((Gdir(j,i)+180)/20)+1) + Gmag(j,i);%On range les differents angles par paquet de 20
        end
    end
    [histMaxMag, indexOfMax] = max(histogram); %on récupère le maximum de l'histogramme (le pic) pour connaitre l'orientation principale du point d'intérêt
    keypoints(1,w) = indexOfMax; %à la place de la coordonnée Y de NewLocalisations, on stocke l'indice correspondant à l'orientation principale
    keypoints(1,w-1) = histMaxMag; %à la place de coordX, on stocke la valeur de cette orientation principale 
    keypoints(1,w-2) = kp;
end
% SORTIE :  1 vecteur contenant pour chaque point d'interet sa magnitude et
% son orientation. L'ordre de rangement est le meme que pour les
% localisation


%% Affichage des points d'intérêt détectés
coordAffichage=[];
for w=3:3:totalPI
    coordAffichage=[coordAffichage; newLocalisations(1,w-1) newLocalisations(1,w)];
    original(newLocalisations(1,w-1),newLocalisations(1,w))=255;
end
RGBA = insertMarker(image,coordAffichage,'color','red');
% figure,imshow(RGBA);
%figure,imshow(original);



%% Creation de la signature des PI
Descripteur1=[];
Descripteur1AvecCoordonnees=[];
KPSignature=[];
dEyestmp=[];
dNosetmp=[];
dBottomtmp=[];
j=0;
id=0;
for w=3:3:totalPI
    id=id+1;
    Wcurrent=newLocalisations(1,w);
    Hcurrent=newLocalisations(1,w-1); % On récupère les localisations des point d'intérêt
    %calcule de l'orientation et de l'amplitude de chaque gradient de chaque pixel dans la zone 16x16 autour du point d'intérêt
    if Wcurrent>7 && Hcurrent>7 && Wcurrent<=width1-8 && Hcurrent<hight1-8%on vérifie qu'il n'y a pas de soucis de bordure
        j=j+1;
        regionMag=Gmag(Hcurrent-7:Hcurrent+8,Wcurrent-7:Wcurrent+8);%on génère les maps des gradients sur la région
        regionDir=Gdir(Hcurrent-7:Hcurrent+8,Wcurrent-7:Wcurrent+8);
    else
        j=j+1;
        continue;
    end

    KPSignature=[];

    %% Division en blocs de 4x4
    for k=1:4
        for l=1:4
            %construction des sous régions (détermination des bordures)
            %une sous région pour l'amplitude, une autre pour l'orientation
            subRegionMag=regionMag(1+(k-1)*4:k*4,1+(l-1)*4:l*4);
            subRegionDir=regionDir(1+(k-1)*4:k*4,1+(l-1)*4:l*4);
            hist=[0 0 0 0 0 0 0 0]; %histogramme de 8 alvéoles
            for m=1:4
                for n=1:4
                    %pour chaque sous région, on construit l'histogramme
                    pourVisuelCoordH=floor((abs(subRegionDir(m,n)-Gdir(Hcurrent,Wcurrent)))/45)+1;
                    hist(1,floor((abs(subRegionDir(m,n)-Gdir(Hcurrent,Wcurrent)))/45)+1) =...
                    hist(1,floor((abs(subRegionDir(m,n)-Gdir(Hcurrent,Wcurrent)))/45)+1) +...
                    subRegionMag(m,n);%On ajoute la magnitude de chaque point dans la partie de l'histogrmme qui correspond a l'orientation relative du point courant avec le point d'interet
                end
            end
            KPSignature=[KPSignature hist]; 
        end
    end
    %save('descriptor.txt','KPSignature','-ascii','-append'); %ok
    %Descripteur1(id,:)=KPSignature; %Descripteurs SIFT pour chaque point clé
    Descripteur1=[Descripteur1;KPSignature];
    Descripteur1AvecCoordonnees=[Descripteur1AvecCoordonnees;double(newLocalisations(1,w-2)) ...
        double(Hcurrent) double(Wcurrent) KPSignature];
end
%newName = input('I want to save the variable under the name:', 's');
%save('Data.mat', '-struct', 'S')  % EDITED
%save('descriptor.m','KPSignature','-ascii','-append'); %ok
%save descriptor.txt Descripteur -ascii -append;
% clearvars -except Descripteur1 Descripteur1AvecCoordonnees image

%2eme image

image2 = guiimage2;
% nameImage2 = 'Noeline-09.jpg';
% image2=imread(nameImage2);
%image2=imtranslate(image,[50,20]);
%image2=imrotate(image,10);
%imshow(image2);


[~,~,deepth]=size(image2);
if deepth==3 %Si l'image est en couleur, on la transforme en gris
    image2=rgb2gray(image2);
end
image2=im2double(image2); %Transformer l'image en matrice de double
original=image2;

%initialisation des octaves (matrices vides pour l'instant)
octave1=[];
octave2=[];
octave3=[];
tic;

%% Octave 1

image_temporaire=imresize(image2,2,'nearest');         %Double la taille de l'image initiale pour de meilleurs resultats
[hight1,width1,~]=size(image_temporaire);                       %Récupération des dimensions de l'image

k2=0;                                       %1ere octave
image_temporaire(hight1:hight1+4,width1:width1+4)=0;  %Zero padding
clear c;

for k1=0:3                                  %Differentes valeurs de sigma pour 1 octave 
    k=sqrt(2);
    sigma=(k^(k1+(2*k2)))*1.6;

    for x=-2:2                              %Generation du filtre gaussien
        for y=-2:2
             h(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma))); 
        end
    end

    for i=1:hight1                          %Convolution image et filtre
        for j=1:width1
            temp1=image_temporaire(i:i+4,j:j+4)'.*h; % Produit de 2 matrices terme à terme
            conv_1(i,j)=sum(sum(temp1));   % 2 fois car sum fonctionne sur un vecteur
        end
    end

    octave1=[octave1 conv_1];               %Ajout de l'image filtree a l'octave
end
%figure(1);        
%imshow(octave1);

%% Octave 2

%clear I_temp conv_1 h;
image_temporaire2=original;
[hight2,width2,~]=size(image_temporaire2);

k2=1; %deuxième octave
image_temporaire2(hight2:hight2+4,width2:width2+4)=0; %zéro padding
clear c;

for k1=0:3
    k=sqrt(2);
    sigma=(k^(k1+(2*k2)))*1.6;

    for x=-2:2
        for y=-2:2
            h1(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma)));
        end
    end

    for i=1:hight2
        for j=1:width2
            temp2=image_temporaire2(i:i+4,j:j+4)'.*h1;       
            conv_2(i,j)=sum(sum(temp2));
        end
    end

    octave2=[octave2 conv_2];                     
end
%figure(2);
%imshow(octave2);

%% Octave 3
clear I_temp2 conv2 h1;
image_temporaire3=imresize(original,0.5,'nearest'); %taille divisée par 2      
[hight3,width3,~]=size(image_temporaire3);

k2=2; %3ème octave
image_temporaire3(hight3:hight3+4,width3:width3+4)=0;
clear c;

for k1=0:3
    k=sqrt(2);
    sigma=(k^(k1+(2*k2)));

    for x=-2:2
        for y=-2:2
            h2(x+3,y+3)=(1/((2*pi)*((k*sigma)*(k*sigma))))*exp(-((x*x)+(y*y))/(2*(k*k)*(sigma*sigma)));
        end
    end


     for i=1:hight3
         for j=1:width3
             temp3=image_temporaire3(i:i+4,j:j+4)'.*h2;   
             conv_3(i,j)=sum(sum(temp3));
         end
     end
     octave3=[octave3 conv_3];
end
%figure (3);
%imshow(octave3);

fprintf('\n Temps pour construction des octaves %.3f s\n',toc) ;
% SORTIE :  3 matrices representants chacune une octave


%% Differences de gaussiennes
DoG1=[];
DoG2=[];
DoG3=[];

for k=1:3
    DoG1=[DoG1 octave1(1:hight1,(k-1)*width1+1:k*width1)-octave1(1:hight1,k*width1+1:(k+1)*width1)];
    DoG2=[DoG2 octave2(1:hight2,(k-1)*width2+1:k*width2)-octave2(1:hight2,k*width2+1:(k+1)*width2)];
    DoG3=[DoG3 octave3(1:hight3,(k-1)*width3+1:k*width3)-octave3(1:hight3,k*width3+1:(k+1)*width3)];
end

%figure(4);
%imshow(DoG1);


% SORTIE :  3 matrices regroupants les 3 resultats des differences de
% gaussiennes pour chaque octave 

%% Localisation des extremas

Localisations=[];

for w=width1+4:2:2*width1-4% supprimer les :2: pour petites images?
    for h=6:2:hight1-4 %on se balade sur l'image, on commence a 6 pour eviter les pb de borne
        HighTab=DoG1(h-3:h+3,w+width1-3:w+width1+3);

        HtmpM=int32(h/2);%On trouve la position correspondante sur l'octave 2 
        WtmpM=int32(w/2)+1;
        MiddleTab=DoG2(HtmpM-2:HtmpM+2,WtmpM-2:WtmpM+2);

        HtmpL=int32(h/4);
        WtmpL=int32(w/4)+1;
        LowTab=DoG3(HtmpL-1:HtmpL+1,WtmpL-width3-1:WtmpL-width3+1);


        maxMT=max(max(MiddleTab));%On cherche les extremums pour chaque octave
        maxLT=max(max(LowTab));
        maxHT=max(max(HighTab));
        minMT=min(min(MiddleTab));
        minLT=min(min(LowTab));
        minHT=min(min(HighTab));


        a=sum(sum(MiddleTab(:,:)==maxMT));%Compte le nombre de fois que maxMT est présent dans la région 5x5 autour du point
        b=sum(sum(MiddleTab(:,:)==minMT));

        maximum = max([maxMT,maxLT,maxHT]);
        minimum = min([minMT,minLT,minHT]);
        %Le pixel courant et soit un max soit un min, et il est le
        %seul de sa région à chaque fois
        if(a==1 && DoG2(HtmpM,WtmpM)==maximum) ||(b==1 && DoG2(HtmpM,WtmpM)==minimum)
            if HtmpM>2 && WtmpM>width2+2 && HtmpM<hight2-2 && WtmpM<2*width2-2
                Localisations = [Localisations HtmpM WtmpM-width2];
            end
        end
    end
end

% SORTIE :  1 vecteur contenant les coordonnées de chaque point d'interet


%% Ecrémation des PI par des seuils en ne selectionnant que les coins
[Gx,Gy]=imgradientxy(DoG1(:,width1+1:2*width1));%(diff_12(:,:));% retour la magnitude des gradient selon x et Y
[Gmag, Gdir]=imgradient(DoG1(:,width1+1:2*width1));

[~, totalPI]=size(Localisations);
newLocalisations=[];
j=1;
for i=2:2:totalPI
    if (abs(Gx(Localisations(1,i-1),Localisations(1,i)))>=0.0015 && abs(Gy(Localisations(1,i-1),Localisations(1,i)))>=0.0015)
        newLocalisations=[newLocalisations j Localisations(1,i-1) Localisations(1,i)];%Adapter le seuil pour qu'il ne me reste plus que X KP
        j=j+1;
    end
end
 % SORTIE : % SORTIE :  1 vecteur contenant les coordonnées de chaque point
 % d'interet triés


%% Orientation des points d'intérêt

[~, totalPI]=size(newLocalisations);
keypoints=[];% stocke mag et dir gradient pour chaque pt
kp=0;
for w=3:3:totalPI % parcours de 3 en 3 car les points d'intérêts occupent 3 places (n°, coordX, coordY)
    kp=kp+1;
    histogram = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    %pour tous les points d'intérêt : permet d'accéder à une zone 5*5
    %autour du point d'intérêt
    %pourMoicoordX=newLocalisations(1,w-1);
    %pourMoicoordY=newLocalisations(1,w);
    for i=newLocalisations(1,w)-2:newLocalisations(1,w)+2%width : coordonnée Y
        for j=newLocalisations(1,w-1)-2:newLocalisations(1,w-1)+2%hight : coordonnée X
            %pourMoiGdir=Gdir(j,i);
            %pourMoiFloor=floor((Gdir(j,i)+180)/20)+1; % pour savoir quelle colonne de l'histogramme remplir : on ajoute 180° pour se ramener à l'intervalle [0;360] puis on divise par 20 pour former des paquets de 20
            %pourMoiGmag=Gmag(j,i)
            histogram(1,floor((Gdir(j,i)+180)/20)+1) = histogram(1,floor((Gdir(j,i)+180)/20)+1) + Gmag(j,i);%On range les differents angles par paquet de 20
        end
    end
    [histMaxMag, indexOfMax] = max(histogram); %on récupère le maximum de l'histogramme (le pic) pour connaitre l'orientation principale du point d'intérêt
    keypoints(1,w) = indexOfMax; %à la place de la coordonnée Y de NewLocalisations, on stocke l'indice correspondant à l'orientation principale
    keypoints(1,w-1) = histMaxMag; %à la place de coordX, on stocke la valeur de cette orientation principale 
    keypoints(1,w-2) = kp;
end
% SORTIE :  1 vecteur contenant pour chaque point d'interet sa magnitude et
% son orientation. L'ordre de rangement est le meme que pour les
% localisation


%% Affichage des points d'intérêt détectés
for w=3:3:totalPI
    original(newLocalisations(1,w-1),newLocalisations(1,w))=255;
    %original(keypoints(1,w-1),keypoints(1,w))=255;
end
% figure,imshow(original);



%% Creation de la signature des PI
Descriptors=[];
Descripteur2=[];
Descripteur2AvecCoordonnees=[];
KPSignature=[];
dEyestmp=[];
dNosetmp=[];
dBottomtmp=[];
j=0;
id=0;
for w=3:3:totalPI
    id=id+1;
    Wcurrent=newLocalisations(1,w);
    Hcurrent=newLocalisations(1,w-1); % On récupère les localisations des point d'intérêt
    %calcule de l'orientation et de l'amplitude de chaque gradient de chaque pixel dans la zone 16x16 autour du point d'intérêt
    if Wcurrent>7 && Hcurrent>7 && Wcurrent<=width1-8 && Hcurrent<hight1-8%on vérifie qu'il n'y a pas de soucis de bordure
        j=j+1;
        regionMag=Gmag(Hcurrent-7:Hcurrent+8,Wcurrent-7:Wcurrent+8);%on génère les maps des gradients sur la région
        regionDir=Gdir(Hcurrent-7:Hcurrent+8,Wcurrent-7:Wcurrent+8);
    else
        j=j+1;
        continue;
    end

    KPSignature=[];

    %% Division en blocs de 4x4
    for k=1:4
        for l=1:4
            %construction des sous régions (détermination des bordures)
            %une sous région pour l'amplitude, une autre pour l'orientation
            subRegionMag=regionMag(1+(k-1)*4:k*4,1+(l-1)*4:l*4);
            subRegionDir=regionDir(1+(k-1)*4:k*4,1+(l-1)*4:l*4);
            hist=[0 0 0 0 0 0 0 0]; %histogramme de 8 alvéoles
            for m=1:4
                for n=1:4
                    %pour chaque sous région, on construit l'histogramme
                    pourVisuelCoordH=floor((abs(subRegionDir(m,n)-Gdir(Hcurrent,Wcurrent)))/45)+1;
                    hist(1,floor((abs(subRegionDir(m,n)-Gdir(Hcurrent,Wcurrent)))/45)+1) =...
                    hist(1,floor((abs(subRegionDir(m,n)-Gdir(Hcurrent,Wcurrent)))/45)+1) +...
                    subRegionMag(m,n);%On ajoute la magnitude de chaque point dans la partie de l'histogrmme qui correspond a l'orientation relative du point courant avec le point d'interet
                end
            end
            KPSignature=[KPSignature hist]; 
        end
    end
    %save('descriptor.txt','KPSignature','-ascii','-append'); %ok
    %Descripteur2(id,:)=KPSignature; %Descripteurs SIFT pour chaque point clé
    Descripteur2=[Descripteur2;KPSignature];
    Descripteur2AvecCoordonnees=[Descripteur2AvecCoordonnees;double(newLocalisations(1,w-2)) ...
        double(Hcurrent) double(Wcurrent) KPSignature];
end

%Descripteur1AvecCoordonnees(120:length(Descripteur1),:)=[];
%Descripteur2AvecCoordonnees(120:length(Descripteur2),:)=[];


Indices=[];
posA=[];
posB=[];
Distance= pdist2(Descripteur1,Descripteur2,'euclidean');
Distance2= pdist2(Descripteur1AvecCoordonnees(:,4:131),Descripteur2AvecCoordonnees(:,4:131),'euclidean');
[ligneD1,colonneD1]=size(Descripteur1AvecCoordonnees);
[ligneD2,colonneD2]=size(Descripteur2AvecCoordonnees);

%%Méthode ratio
%{
tic;
for i=1:ligneD1
    [min coordmin]=mink(Distance2(i,:),2); %on récupère les 2 plus petites valeurs
    ratio=min(1,1)/min(1,2);
    if ratio<1
        Indices=[Indices; i coordmin(1,1)];%récupère les indices des vecteurs
        xa=Descripteur1AvecCoordonnees(i,2);
        ya=Descripteur1AvecCoordonnees(i,3);
        xb=Descripteur2AvecCoordonnees(coordmin(1,1),2);
        yb=Descripteur2AvecCoordonnees(coordmin(1,1),3);
        shift=10;
         if(xb<=xa+shift && xb>=xa-shift && yb<=ya+shift && yb>=ya-shift)
            posA=[posA;xa ya];
            posB=[posB;xb yb];
         end
    end
end
fprintf('\n Temps pour méthode ratio %.3f s\n',toc) ;
RGB = insertMarker(image,posA);
RGB2 = insertMarker(image2,posB);
figure,imshow([RGB,RGB2]);
figure,showMatchedFeatures(RGB,RGB2,posA,posB,'montage')
%}

%%Méthode distance euclidienne

tic;
for i=1:ligneD1
    for j=1:ligneD2
        if Distance2(i,j)<0.8
            Indices=[Indices; i j];%récupère les indices des vecteurs
            %forMe = Descripteur1AvecCoordonnees(i,2);
            %forMe2 = Descripteur1AvecCoordonnees(i,3);
            xa=Descripteur1AvecCoordonnees(i,2);
            ya=Descripteur1AvecCoordonnees(i,3);
            xb=Descripteur2AvecCoordonnees(j,2);
            yb=Descripteur2AvecCoordonnees(j,3);
            
            %tester si dans la zone
            shift=10;
             if(xb<=xa+shift && xb>=xa-shift && yb<=ya+shift && yb>=ya-shift)
                posA=[posA;xa ya];
                posB=[posB;xb yb];
            end
            
            %{    
            %cas ou l'image est translatée
            if(xb-50<=xa+10 && xb-50>=xa-10 && yb-20<=ya+10 && yb-20>=ya-10)
                posA=[posA;xa ya];
                posB=[posB;xb yb];
            end
            %}
        end
    end
end
fprintf('\n Temps pour méthode distance euclidienne %.3f s\n',toc) ;

RGB = insertMarker(image,posA);
RGB2 = insertMarker(image2,posB);
matchedPoints1 = posA;
matchedPoints2 = posB;
% figure,imshow([RGB,RGB2]);
% figure,showMatchedFeatures(RGB,RGB2,posA,posB,'montage')
%}



%newName = input('I want to save the variable under the name:', 's');
%save('Data.mat', '-struct', 'S')  % EDITED
%save('descriptor.m','KPSignature','-ascii','-append'); %ok
%save descriptor.txt Descripteur -ascii -append;


%%Méthode KDTree
%{
tic;
Mdl = KDTreeSearcher(Descripteur1);
[Idx,D]=knnsearch(Mdl,Descripteur2);
%Idx=kmeans(Descripteur1,10);
coordImage1=[];
for i=1:length(Idx)
    coordImage1=[coordImage1;Descripteur1AvecCoordonnees(Idx(i),2) ...
       Descripteur1AvecCoordonnees(Idx(i),3) ];
end

coordImage2=[];
for i=1:ligneD2
    coordImage2=[coordImage2; Descripteur2AvecCoordonnees(i,2) Descripteur2AvecCoordonnees(i,3)];
end

matchedPoints1=[];
matchedPoints2=[];
for i=1:ligneD2
    if(coordImage2(i,1)<=coordImage1(i,1)+10 && coordImage2(i,1)>=coordImage1(i,1)-10 && ...
            coordImage2(i,2)>=coordImage1(i,2)-10 && coordImage2(i,2)<=coordImage1(i,2)+10)
        matchedPoints1=[matchedPoints1; coordImage1(i,1) coordImage1(i,2)];
        matchedPoints2=[matchedPoints2; coordImage2(i,1) coordImage2(i,2)];
    end
end
fprintf('\n Temps pour méthode KDTree %.3f s\n',toc) ;       

RGB = insertMarker(image,matchedPoints1);
RGB2 = insertMarker(image2,matchedPoints2);
% figure,imshow([RGB,RGB2]);
% figure,showMatchedFeatures(RGB,RGB2,matchedPoints1,matchedPoints2,'montage');
%}
disp("SIFT is done");

