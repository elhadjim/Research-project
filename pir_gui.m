%%
function varargout = pir_gui(varargin)
% PIR_GUI MATLAB code for pir_gui.fig
%      PIR_GUI, by itself, creates a new PIR_GUI or raises the existing
%      singleton*.
%
%      H = PIR_GUI returns the handle to a new PIR_GUI or the handle to
%      the existing singleton*.
%
%      PIR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIR_GUI.M with the given input arguments.
%
%      PIR_GUI('Property','Value',...) creates a new PIR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pir_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pir_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pir_gui

% Last Modified by GUIDE v2.5 16-Jun-2020 21:49:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pir_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pir_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pir_gui is made visible.
function pir_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pir_gui (see VARARGIN)

% Choose default command line output for pir_gui
handles.output = hObject;

axes(handles.camera);
handles.cam = webcam;
imWidth = handles.cam.AvailableResolutions(1);
hImage = image(zeros(480,640,3), 'Parent', handles.camera);
handles.totalWidth = handles.camera.XLim(2) - handles.camera.XLim(1);

% view(axes, [0 90]);
% setappdata(hImage, 'UpdatePreviewWindowFcn', @mypreview_fcn);

camorbit(180,180);
preview(handles.cam, hImage);
% closePreview(handles.cam);
set(hObject, 'HandleVisibility', 'off');
close all;
set(hObject, 'HandleVisibility', 'on');
handles.detectedfaces = {};
handles.detectedfaces_coordinates = [];
handles.lbp_scale = 5;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pir_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pir_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Capture image button
% --- Executes on button press in capturebutton.
function capturebutton_Callback(hObject, eventdata, handles)
% hObject    handle to capturebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject, 'String') == "Prendre la photo")
    set(handles.text2, 'String', 'Photo prise');
    handles.img = flip(snapshot(handles.cam),2);
    imshow(handles.img);
    delete(handles.cam);
    set(hObject, 'String', "Prendre une autre photo");
else
    axes(handles.camera);
    handles.cam = webcam;
    imWidth = handles.cam.AvailableResolutions(1);
    disp(imWidth);
    hImage = image(zeros(480,640,3), 'Parent', handles.camera);
    camorbit(180,180);
    preview(handles.cam, hImage);
    set(hObject, 'String', "Prendre la photo");
    handles.img = [];
    handles.detectedfaces = {};
    handles.detectedfaces_coordinates = [];
    handles.siftind=0;
end
guidata(hObject, handles);


%% Select File Button
% --- Executes on button press in selectfilebutton.
function selectfilebutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectfilebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile({'*.jpg; *.gif; *.png'}, 'Choisir un fichier');
handles.detectedfaces = {};
handles.detectedfaces_coordinates = [];
if(file~=0)
    set(handles.text2, 'String', 'Image sélectionnée');
    handles.img = imread(fullfile(path, file));
    handles.detectedfaces{1} = imresize(imread(fullfile(path, file)),[256 192]);
    handles.detectedfaces_coordinates = [1 1 size(handles.img,2) size(handles.img,1)];
    imshow(handles.img);
    delete(handles.cam);
    set(handles.capturebutton, 'String', "Prendre une photo");
else
    handles.img = [];
end
guidata(hObject, handles);


%% Selection of face detection algorithm
% --- Executes on selection change in facedetectionchoice.
function facedetectionchoice_Callback(hObject, eventdata, handles)
% hObject    handle to facedetectionchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns facedetectionchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from facedetectionchoice

contents = cellstr(get(hObject,'String'));
if(contents{get(hObject,'Value')}=="SURF")
    set(handles.text2, 'String', "Cliquez sur 'Photo n°2' pour prendre une photo de l'objet à détecter dans la photo actuelle");
    set(handles.siftbutton1, 'String', 'Photo n°2');
    set(handles.siftbutton1, 'Visible', 'on');
    set(handles.edit1, 'Visible', 'off');    
elseif(contents{get(hObject,'Value')}=="SIFT")
    set(handles.text2, 'String', "Cliquez sur 'Photo n°2' pour prendre une deuxième photo, qui sera comparée à la première à l'aide du descripteur SIFT");
    set(handles.siftbutton1, 'String', 'Photo n°2');
    set(handles.siftbutton1, 'Visible', 'on');
    set(handles.edit1, 'Visible', 'off');
elseif(contents{get(hObject,'Value')}=="LBP")
    set(handles.edit1, 'Visible', 'on');
    set(handles.text2, 'String', "Entrez une valeur d'échelle (approximation de taille de l'image/taille d'un visage) pour optimiser l'exécution de l'algorithme :");
    set(handles.siftbutton1, 'Visible', 'off');
    set(handles.siftbutton2, 'Visible', 'off');
else
    set(handles.text2, 'String', "");
    set(handles.siftbutton1, 'Visible', 'off');
    set(handles.siftbutton2, 'Visible', 'off');
    set(handles.edit1, 'Visible', 'off');
end


%% Face detection popup menu settings
% --- Executes during object creation, after setting all properties.
function facedetectionchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to facedetectionchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'LBP'; 'SIFT'; 'SURF'; 'Réseau de neurones'});


%% Selection of face recognition algorithm
% --- Executes on selection change in facerecognitionchoice.
function facerecognitionchoice_Callback(hObject, eventdata, handles)
% hObject    handle to facerecognitionchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns facerecognitionchoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from facerecognitionchoice

contents = cellstr(get(hObject,'String'));
if(contents{get(hObject,'Value')}=="Réseau de neurones n°1" || contents{get(hObject,'Value')}=="Réseau de neurones AlexNet")
    set(handles.text2, 'String', "Attention, pour le réseau de neurones, l'image doit être en couleurs.");
else
    set(handles.text2, 'String', "");
end
set(handles.siftbutton1, 'Visible', 'off');
set(handles.siftbutton2, 'Visible', 'off');
set(handles.edit1, 'Visible', 'off');

%% Face recognition popup menu settings
% --- Executes during object creation, after setting all properties.
function facerecognitionchoice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to facerecognitionchoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Eigenfaces'; 'Fisherfaces'; 'Réseau de neurones n°1'; 'Réseau de neurones AlexNet'});


%% Launch detection button
% --- Executes on button press in launchdetection.
function launchdetection_Callback(hObject, eventdata, handles)
% hObject    handle to launchdetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA

try
    if(isempty(handles.img))
        set(handles.text2, 'String', 'Veuillez prendre une photo.');
    else
        set(handles.text2, 'String', 'La détection de visages est en cours...');
        handles.detectedfaces = {};
        handles.detectedfaces_coordinates = [];
        contents = cellstr(get(handles.facedetectionchoice,'String'));
        detectionchoice = contents{get(handles.facedetectionchoice,'Value')};
        imshow(handles.img);
        drawnow;
        
        % Détection avec LBP
        if(detectionchoice=="LBP")
            disp(handles.lbp_scale);
            [xs, ys, widths, heights] = LBP_GUI(handles.lbp_scale, handles.img);
            if(isempty(xs))
                set(handles.text2, 'String', "Aucun visage n'a été détecté");
            else
                k=1;
                set(handles.text2, 'String', "");
                for i=1:length(xs)
                    % rectangle bleu autour du visage détecté
%                     rectangle('Position', [xs(i) ys(i) widths(i) heights(i)], 'EdgeColor', 'blue', 'LineWidth', 1);
                    % Détection des yeux puis recadrage en fonction de leur position
                    [eyex, eyey, eyew, eyeh] = LBP_CadrageYeux(imcrop(handles.img, [xs(i) ys(i) widths(i) heights(i)]));
                    if(isempty(eyex)==false)   
                        % rectangle vert autour des yeux
%                         rectangle('Position', [xs(i)+eyex ys(i)+eyey eyew eyeh], 'EdgeColor', 'green', 'LineWidth', 1);
                        % Calcul des nouvelles coordonnées x et y en fonction de la position des yeux
                        newx = xs(i)+eyex+(eyew/2)-(widths(i)/2); % horizontalement : centrer
                        newy = ys(i)+eyey+(eyeh/2)-(0.4*heights(i)); % verticalement : mettre les yeux à 40% du haut du visage
                        % rectangle rouge autour de l'image recadrée
                        rectangle('Position', [newx newy widths(i)-1 heights(i)-1], 'EdgeColor', 'red', 'LineWidth', 1);
                        handles.detectedfaces{end+1} = imcrop(handles.img, [newx newy widths(i) heights(i)]);
                        handles.detectedfaces_coordinates(k,1) = newx;
                        handles.detectedfaces_coordinates(k,2) = newy;
                        handles.detectedfaces_coordinates(k,3) = widths(i);
                        handles.detectedfaces_coordinates(k,4) = heights(i);
                        k=k+1;
                    end
                end
%                 rectangle('Position', [handles.totalWidth-xs(end)-widths(end) ys(end) widths(end) heights(end)], 'EdgeColor', 'blue', 'LineWidth', 1);
            end
            if(isempty(handles.detectedfaces))
                set(handles.text2, 'String', "Aucun visage n'a été détecté");
            end
            
        % Détection avec SIFT
        elseif(detectionchoice=="SIFT")
            if(isempty(handles.img2))
                set(handles.text2, 'String', "Veuillez prendre une deuxième photo");
            else
                set(handles.text2, 'String', "La détection de visages est en cours...");
                disp(size(handles.img));
                disp(size(handles.img2));
                minimage = min(size(handles.img, 1), size(handles.img2, 1));
                disp(minimage);
                image1 = imresize(handles.img, [minimage size(handles.img,2)*minimage/size(handles.img,1)]);
                image2 = imresize(handles.img2, [minimage size(handles.img2,2)*minimage/size(handles.img2,1)]);
                if(size(image1,3)~=size(image2,3))
                    if(size(image1,3)==3)
                        image1=rgb2gray(image1);
                    elseif(size(image2,3)==3)
                        image2=rgb2gray(image2);
                    end
                end
                disp(size(image1));
                disp(size(image2));
                newimage = [image1 image2];
                imshow(newimage);
                drawnow;
                [RGB,RGB2,matchedPoints1,matchedPoints2] = SIFT_GUI(image1, image2);
                newrgb = [RGB RGB2];
                imshow(newrgb);
                showMatchedFeatures(RGB,RGB2,matchedPoints1,matchedPoints2,'montage');
                set(handles.text2, 'String', "");
            end
            
        % Détection avec SURF
        elseif(detectionchoice=="SURF")
            if(isempty(handles.img2))
                set(handles.text2, 'String', "Veuillez prendre une deuxième photo");
            else
                set(handles.text2, 'String', "La détection de visages est en cours...");
                drawnow;
                minimage = min(size(handles.img, 1), size(handles.img2, 1));
                image1 = imresize(handles.img, [minimage size(handles.img,2)*minimage/size(handles.img,1)]);
                image2 = imresize(handles.img2, [minimage size(handles.img2,2)*minimage/size(handles.img2,1)]);
                if(size(image1,3)~=size(image2,3))
                    if(size(image1,3)==3)
                        image1=rgb2gray(image1);
                    elseif(size(image2,3)==3)
                        image2=rgb2gray(image2);
                    end
                end
                newimage = [image1 image2];
                imshow(newimage);
                drawnow;
                newBoxPolygon = SURF_GUI(image2, image1);
                imshow(image1);
                line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'red', 'LineWidth', 1);
%                 rectangle('Position', [newBoxPolygon(1,1) newBoxPolygon(1,2) newBoxPolygon(3,1)-newBoxPolygon(1,1) newBoxPolygon(3,2)-newBoxPolygon(1,2)]);
                handles.detectedfaces{1} = imcrop(handles.img, [newBoxPolygon(1,1) newBoxPolygon(1,2) newBoxPolygon(3,1)-newBoxPolygon(1,1) newBoxPolygon(3,2)-newBoxPolygon(1,2)]);
                handles.detectedfaces_coordinates(1,1) = newBoxPolygon(1,1);
                handles.detectedfaces_coordinates(1,2) = newBoxPolygon(1,2);
                handles.detectedfaces_coordinates(1,3) = newBoxPolygon(3,1)-newBoxPolygon(1,1);
                handles.detectedfaces_coordinates(1,4) = newBoxPolygon(3,2)-newBoxPolygon(1,2);
                set(handles.text2, 'String', "");
            end

        % Détection avec les Eigenfaces
%         elseif(detectionchoice=="Eigenfaces")
%             disp(handles.totalWidth);
%             [xs, ys, widths, heights] = EigenfacesForDetection(2, handles.img);
%             clear('newface');
%             if(isempty(xs))
%                 set(handles.text2, 'String', "Aucun visage détecté");
%             else
%                 handles.detectedfaces_coordinates(:,1) = xs;
%                 handles.detectedfaces_coordinates(:,2) = ys;
%                 handles.detectedfaces_coordinates(:,3) = widths;
%                 handles.detectedfaces_coordinates(:,4) = heights;
%                 set(handles.text2, 'String', "");
%                 for i=1:length(xs)
%                     rectangle('Position', [xs(i) ys(i) widths(i)-1 heights(i)-1], 'EdgeColor', 'red', 'LineWidth', 1);
%                     handles.detectedfaces{end+1} = imcrop(handles.img, [xs(end) ys(end) widths(end) heights(end)]);
%                 end
%     %             rectangle('Position', [xs(end) ys(end) widths(end) heights(end)], 'EdgeColor', 'blue', 'LineWidth', 1);
%             end
        elseif(detectionchoice=="Réseau de neurones")
            FaceDetector = vision.CascadeObjectDetector();
            Bounding_Box = step(FaceDetector, handles.img);
            n = size(Bounding_Box, 1);
            if(n==0)
                set(handles.text2, 'String', "Aucun visage n'a été détecté");
            else
                set(handles.text2, 'String', "");
                for ii=1:n
                    rectangle('Position', Bounding_Box(ii,:), 'EdgeColor', 'red', 'LineWidth', 1);
                    handles.detectedfaces{end+1} = imcrop(handles.img, Bounding_Box(ii,:));
                    handles.detectedfaces_coordinates(ii,1) = Bounding_Box(ii,1);
                    handles.detectedfaces_coordinates(ii,2) = Bounding_Box(ii,2);
                    handles.detectedfaces_coordinates(ii,3) = Bounding_Box(ii,3);
                    handles.detectedfaces_coordinates(ii,4) = Bounding_Box(ii,4);
                end
            end
        end
    end
catch ME
    if (strcmp(ME.message,"Reference to non-existent field 'img'."))
        set(handles.text2, 'String', 'Veuillez prendre une photo.');
        disp(ME.message);
    elseif (strcmp(ME.message,"Reference to non-existent field 'img2'."))
        set(handles.text2, 'String', 'Veuillez prendre une deuxième photo.');
        disp(ME.message);
    elseif(strcmp(ME.identifier, "vision:points:notEnoughMatchedPts"))
        set(handles.text2, 'String', "Il n'y a pas assez de points de correspondance entre ces deux images.")
        disp(ME.message);
    elseif(strcmp(ME.identifier, "MATLAB:insertMarker:incorrectSize"))
        set(handles.text2, 'String', "Il n'y a pas assez de points de correspondance entre ces deux images.")
        disp(ME.message);
    else
        rethrow(ME);
    end
end

guidata(hObject, handles);


%% Launch recognition button
% --- Executes on button press in launchrecognition.
function launchrecognition_Callback(hObject, eventdata, handles)
% hObject    handle to launchrecognition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if(isempty(handles.img))
        set(handles.text2, 'String', 'Veuillez prendre une photo.');
    else
        if(isempty(handles.detectedfaces))
            set(handles.text2, 'String', 'Veuillez lancer la détection de visages ou importer une image de visage');
        else
            contents = cellstr(get(handles.facerecognitionchoice,'String'));
            recognitionchoice = contents{get(handles.facerecognitionchoice,'Value')};
            set(handles.text2, 'String', 'La reconnaissance du visage est en cours...');
            drawnow;
            % Tracé des rectangles reconnus
            imshow(handles.img);
            for i=1:size(handles.detectedfaces_coordinates, 1)
                rectangle('Position', [handles.detectedfaces_coordinates(i,1) handles.detectedfaces_coordinates(i,2) handles.detectedfaces_coordinates(i,3)-1 handles.detectedfaces_coordinates(i,4)-1], 'EdgeColor', 'red', 'LineWidth', 1);
            end
            % Reconnaissance
            % Eigenfaces
            if(recognitionchoice=="Eigenfaces")      
                results = strings([1,length(handles.detectedfaces)]);
                for i=1:length(handles.detectedfaces)
                    newface = handles.detectedfaces{i};
                    if(size(newface,1)==size(newface,2))
                        newwidth = size(newface,2)/1.1;
                        newx = handles.detectedfaces_coordinates(i,1) + (size(newface,2)-newwidth);
                        newheight = 4/3*newwidth;
                        newy = handles.detectedfaces_coordinates(i,2) - (newheight-size(newface,1));
                        newface = imcrop(handles.img, [newx, newy, newwidth, newheight]);
                    end
                    if(size(newface,3)==3)
                        newface = rgb2gray(newface);
                    end
                    newface = imresize(newface, [256, 192]);
                    results(i) = Eigenfaces_GUI(newface);
                    text(handles.detectedfaces_coordinates(i,1)+handles.detectedfaces_coordinates(i,3)/2, handles.detectedfaces_coordinates(i,2)-30, results(i), 'Color', 'red', 'HorizontalAlignment', 'center');
                end
                set(handles.text2, 'String', results);
            % Fisherfaces
            elseif(recognitionchoice=="Fisherfaces")
                results = strings([1,length(handles.detectedfaces)]);
                for i=1:length(handles.detectedfaces)
                    newface = handles.detectedfaces{i};
                    if(size(newface,1)==size(newface,2))
                        newwidth = size(newface,2)/1.1;
                        newx = handles.detectedfaces_coordinates(i,1) + (size(newface,2)-newwidth);
                        newheight = 4/3*newwidth;
                        newy = handles.detectedfaces_coordinates(i,2) - (newheight-size(newface,1));
                        newface = imcrop(handles.img, [newx, newy, newwidth, newheight]);
                    end
                    if(size(newface,3)==3)
                        newface = rgb2gray(newface);
                    end
                    newface = imresize(newface, [256, 192]);
                    results(i) = Fisherfaces_Recognition_GUI(newface);
                    text(handles.detectedfaces_coordinates(i,1)+handles.detectedfaces_coordinates(i,3)/2, handles.detectedfaces_coordinates(i,2)-30, results(i), 'Color', 'red', 'HorizontalAlignment', 'center');
                end
                set(handles.text2, 'String', results);
            % Réseau Elhadji
            elseif(recognitionchoice=="Réseau de neurones n°1")
                results = strings([1,length(handles.detectedfaces)]);
                for i=1:length(handles.detectedfaces)
                    newface = handles.detectedfaces{i};
                    if(size(newface,1)~=size(newface,2))
                        newx = handles.detectedfaces_coordinates(i,1) - (size(newface,1)-size(newface,2))/2;
                        newface = imcrop(handles.img, [newx, handles.detectedfaces_coordinates(i,2), size(newface,1), size(newface,1)]);
                    end
                    if(size(newface,3)~=3)
                        I = repmat(newface, [1 1 3]);
                    else
                        I = newface;
                    end
                    I=imresize(I,[224 224]);
                    load('net_training');
                    [label,score] = classify(net,I);
                    results(i) = string(label) + " avec "  + 100*max(score) + "%";
                    text(handles.detectedfaces_coordinates(i,1)+handles.detectedfaces_coordinates(i,3)/2, handles.detectedfaces_coordinates(i,2)-30, string(label), 'Color', 'red', 'HorizontalAlignment', 'center');
                end
                set(handles.text2, 'String', results);
            % AlexNet
            elseif(recognitionchoice=="Réseau de neurones AlexNet")
                results = strings([1,length(handles.detectedfaces)]);
                for i=1:length(handles.detectedfaces)
                    newface = handles.detectedfaces{i};
                    if(size(newface,1)~=size(newface,2))
                        newx = handles.detectedfaces_coordinates(i,1) - (size(newface,1)-size(newface,2))/2;
                        newface = imcrop(handles.img, [newx, handles.detectedfaces_coordinates(i,2), size(newface,1), size(newface,1)]);
                    end
                    if(size(newface,3)~=3)
                        I = repmat(newface, [1 1 3]);
                    else
                        I = newface;
                    end
                    I=imresize(I,[227 227]);
                    load('alexnet_training');
                    [label,score] = classify(netTransfer,I);
                    results(i) = string(label) + " avec "  + 100*max(score) + "%";
                    text(handles.detectedfaces_coordinates(i,1)+handles.detectedfaces_coordinates(i,3)/2, handles.detectedfaces_coordinates(i,2)-30, string(label), 'Color', 'red', 'HorizontalAlignment', 'center');
                end
                set(handles.text2, 'String', results);
            end
        end
    end
catch ME
    if (strcmp(ME.identifier,'MATLAB:nonExistentField'))
        set(handles.text2, 'String', 'Veuillez prendre une photo.');
        disp(ME.message);
    else
        rethrow(ME);
    end
end
guidata(hObject, handles);


% --- Executes on button press in siftbutton1.
function siftbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to siftbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(hObject, 'String') == "Photo n°2" || get(hObject, 'String') == "Prendre une autre photo" || get(hObject, 'String') == "Prendre une photo")
    set(handles.siftbutton2, 'Visible', 'on');
    axes(handles.camera);
    handles.cam = webcam;
    imWidth = handles.cam.AvailableResolutions(1);
    disp(imWidth);
    hImage = image(zeros(480,640,3), 'Parent', handles.camera);
    camorbit(180,180);
    preview(handles.cam, hImage);
    set(hObject, 'String', "Prendre la photo");
    handles.img2 = [];
else
    set(handles.text2, 'String', 'Photo prise');
    handles.img2 = flip(snapshot(handles.cam),2);
    imshow(handles.img2);
    delete(handles.cam);
    set(hObject, 'String', "Prendre une autre photo");
end
guidata(hObject, handles);



% --- Executes on button press in siftbutton2.
function siftbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to siftbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile({'*.jpg; *.gif; *.png'}, 'Choisir un fichier');
if(file~=0)
    set(handles.text2, 'String', 'Image sélectionnée');
    handles.img2 = imread(fullfile(path, file));
    imshow(handles.img2);
    delete(handles.cam);
    set(handles.siftbutton1, 'String', "Prendre la photo");
else
    handles.img2 = [];
end
guidata(hObject, handles);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

handles.lbp_scale = str2double(get(hObject, 'String'));
if(handles.lbp_scale==0 || isnan(handles.lbp_scale))
    handles.lbp_scale=5;
end
disp(['LBP scale changed to ' num2str(handles.lbp_scale)]);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
