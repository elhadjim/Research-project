function result = Eigenfaces_GUI(GUIimage)
    load('Eigenfaces_training');
    disp("Eigenfaces is launched!");
%     figure;
%     imshow(GUIimage);
    imagematrix = GUIimage;
    V = double(imagematrix(:));   % convert NxN matrix into vector of dimension N²    

    patternvector = (eigenfacesmatrix'*(V-my_psi))';   % pattern vector (line vector of size Mprime)

    epsilon = zeros(I,1);   % euclidian distances to each of the I class vectors (column vector of size I)
    for i=1:I
    %     epsilon(i) = norm(patternvector - classvectors(i,:));
        epsilon(i) = sum((patternvector - classvectors(i,:)).^2./dsb'); % mahalanobis distance
    end
    [mineps, ind] = min(epsilon);
%     disp(mineps);

    epsilon2 = norm((V-my_psi)-(patternvector*(eigenfacesmatrix'))');  % euclidian distance to the face space
    % (euclidian distance between the mean-adjusted image (V-my_psi) and its projection onto face space
%     disp(epsilon2);

%     for i=1:I
%         if(epsilon(i)<treshold)
%             disp([num2str(epsilon(i)) ': ' imagefiles(i*P).name]);
%         end
%     end
%         disp(mineps);

    if(mineps < treshold && epsilon2<facetreshold)
%         figure;
%         imshow(imrotate(mat2gray(vec2mat((classvectors(ind,:)*(eigenfacesmatrix'))',256)),-90)); % to display the eigenfaces
        str = strsplit(imagefiles(ind*P).name, {'-','.gif','1','2','0'}, 'CollapseDelimiters', true);
%         disp(['This image belongs to the face class n°' num2str(ind) ' (face class of' chaine ').']);
        result = str{1};
        if(strcmp(str{1}, 'Penelope'))
            result = "Pénélope";
        elseif(strcmp(str{1}, 'Noeline'))
            result = "Noëline";
        elseif(strcmp(str{1}, 'Stephan'))
            result = "Stéphan";
        end
    elseif (mineps >= treshold && epsilon2<facetreshold)
%     elseif (mineps >= treshold && mineps<0.3)
        result = "Personne inconnue";
    else
        result = "Personne inconnue";
    end
    disp('Eigenfaces is done');

end

    