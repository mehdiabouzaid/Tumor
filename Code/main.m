%% Postulat

% Comme les mask114.png, mask121.png, mask129.png ne contiennent qu'une
% seule tumeur, on suppose qu'on a qu'une seule tumeur a detecter par image

%%
clc
close all
clear
image = [];

%% Ajout du repertoire des donnees

path=which('main.m');
path=[path(1:(end-length('Code/main.m'))) 'Donnees/tumeur'];

if ~exist(path, 'dir')
    disp('Pour des raisons de confidentialite, les images medicales ne sont pas stockees en ligne.')
    error('Creez le repertoire tumeur dans le repertoire Donnees et mettez les images dans le repertoire tumeur')
end

addpath(path)

%% Image Originale

nomImage = 'dcm114.png';

image = im2double(imread(nomImage));
[n,p,z] = size(image);

disp(['Dimensions Image Originale : n=', num2str(n), ' p=', num2str(p), ' z=', num2str(z)])

if (z>1)
    image = rgb2gray(image);
end

i_subplot = 1;

figure('Renderer','painters','Position',[10 10 1600 800])
subplot(4,4,i_subplot)
imagesc(image)
title(['Image Originale ' nomImage])
i_subplot = i_subplot+1;

%% Mask : verite terrain

nomMask = 'mask114.png';
mask = im2double(imread(nomMask)); % verite terrain

subplot(4,4,i_subplot)
imagesc(mask)
title(['Mask ' nomMask ' : verite terrain'])
i_subplot = i_subplot+1;

%% Histogramme 

subplot(4,4,i_subplot)
h=histogram(image);
grid on
title('Histogramme Image Originale')
i_subplot = i_subplot+1;

% On remarque que l'image originale contient enormement de pixels noirs
% (arriere plan) qui peuvent fausser la determination automatique du seuil
% pour le seuillage

%% Pre-traitement : Filtrage Gaussien

subplot(4,4,i_subplot)
imageFiltrageGaussien = imgaussfilt(image);
imagesc(image)
title('Image Originale Apres filtrage gaussien')
i_subplot = i_subplot+1;

%% Determination du seuil avec Otsu

seuilOtsuImageOriginale = graythresh(imageFiltrageGaussien)
% la fonction graythresh de Matlab utilise la methode d'Otsu vue en TP pour
% determiner le seuil de maniere automatique

%% Segmentation (seuillage)

imageSeuillee = imbinarize(imageFiltrageGaussien, seuilOtsuImageOriginale);

subplot(4,4,i_subplot)
imagesc(imageSeuillee)
title(['Apres segmentation avec Otsu, s=' num2str(seuilOtsuImageOriginale)])
i_subplot = i_subplot+1;

% On remarque que la methode d'Otsu n'est pas satisfaisante dans notre cas
% car le seuil est fausse par l'arriere plan noir. La tumeur et l'organe se
% trouve alors dans la meme region et on ne pas les dissocier

%% Image Rognee
%% Rogner l'image originale (retirer les colonnes et les lignes qui ne 
%% contiennent des valeurs entre 0 et 0.1)

% Retirer les lignes et les colonnes qui ne contiennent que des 0 ne
% suffisent pas pour avoir des resultats satisfaisants avec Otsu. C'est
% pourquoi, on retire les lignes et les colonnes de l'image originale qui
% ont des valeurs entre 0 et 0.1. On a choisi 0.1 de maniere empirique

seuilMax = 0.1;
[imageRognee,indicesLignesNoires,indicesColonnesNoires] = rognerImage(image,seuilMax,n,p);

subplot(4,4,i_subplot)
imagesc(imageRognee);
title('Image Rognee')
i_subplot = i_subplot+1;

subplot(4,4,i_subplot)
h=histogram(imageRognee);
grid on
title('Histogramme Image Rognee')
i_subplot = i_subplot+1;

%% Pre-traitement : Filtrage Gaussien

imageRogneeFitlrageGaussien = imgaussfilt(imageRognee);

subplot(4,4,i_subplot)
imagesc(imageRogneeFitlrageGaussien)
title('Image Rognee Apres filtrage gaussien')
i_subplot = i_subplot+1;

%% Determination du seuil avec Otsu

seuilOtsuImageRognee = graythresh(imageRogneeFitlrageGaussien)

% Le seuil trouve avec Otsu a partir de l'image rognee est tres different
% de celui trouve a partir de l'image originale

%% Segmentation (seuillage)

imageRogneeSeuillee = imbinarize(imageRogneeFitlrageGaussien, seuilOtsuImageRognee);

subplot(4,4,i_subplot)
imagesc(imageRogneeSeuillee)
title(['Apres segmentation avec Otsu, s=',num2str(seuilOtsuImageRognee)])
i_subplot = i_subplot+1;

% imageRogneeSeuillee n'a pas les memes dimensions que le mask (verite
% terrain) donc on peut pas utiliser la fonction dice de Matlab directement
% pour calculer la similitude.

% Nous devons alors reconstruire l'imageSeuillee a partir de l'imageRogneeSeuillee

%% Reconstruction

resOtsu = ones(size(image));

for j = 1:length(indicesColonnesNoires)
    resOtsu(:,indicesColonnesNoires(j)) = 0;
end

for i = 1:length(indicesLignesNoires)
    resOtsu(indicesLignesNoires(i),:) = 0;
end

resOtsu(resOtsu==1) = imageRogneeSeuillee;
    
i_subplot = afficherResultats(resOtsu,image,mask,'Otsu',i_subplot);

%% K-Means

K = 3; % nombre de clusters
% 3 clusters : arriere plan, organe, tumeur

resKmeans = reshape(kmeans(imageFiltrageGaussien(:),K),size(imageFiltrageGaussien));

subplot(4,4,i_subplot)
imagesc(resKmeans)
title(['K-Means, K=' num2str(K)])
i_subplot = i_subplot+1;

i_subplot = afficherResultats(resKmeans,image,mask,'K-Means',i_subplot);

%% Interpretation des resultats

% Otsu ne donne pas de bons resultats avec l'image originale car l'arriere
% plan contient trop de pixels noirs

% Une fois qu'on a rogne l'image (retire les pixels entre 0 et 0.1 :
% pixels noirs, fonces), Otsu donne de bons resultats.

% La methode de K-Means donne de bons resultats sur l'image originale.

% Les methodes d'Otsu sur l'image rognee et K-Means sur l'image originale
% donnent les memes valeurs de similitudes pour dcm114.png (0.70175) et dcm121.png (0.8046).
% Pour dcm129.png, :
% - parfois K-Means donne une similitude de 0.81818, un peu plus elevee qu'Otsu (0.80597)
% - parfois K-Means donne la meme similitude qu'Otsu (0.80597)
% Ces similitudes qui changent pour dcm129.png sont dues a l'initialisation
% des centres qui peut changer dans l'algorithme de K-Means 

% Ces similitudes sont proches de 1 donc ce sont de bons resultats.

% Quelle que soit la methode, quelle que soit l'image, le centre de masse
% de la tumeur predite et du mask ne sont pas les memes.

% Limites du projet
% Notre projet fonctionne sur les 3 images qui nous ont ete donnees mais il
% pourrait ne pas fonctionner sur d'autres images


%% Fonctions

%% afficherResultats

function [i_subplot] = afficherResultats(res,image,mask,methode,i_subplot)
    %% Recuperation de la tumeur

    % Combien y a-t-il de clusters ?
    modaliteClusters = unique(res)';

    % Combien y a-t-il de pixels dans chaque cluster ?
    tailleClusters = zeros(length(modaliteClusters),1);
    indice=1;
    for i = modaliteClusters
        tailleClusters(indice) = length(res(res==i));
        indice = indice+1;
    end

    % Quel cluster est le plus petit (ce cluster correspond alors a la tumeur)
    [~, indiceTailleClusterMin] = min(tailleClusters);

    %% Isoler la tumeur

    indicesTumeur = find(res==modaliteClusters(indiceTailleClusterMin));
    imageBinarisee = zeros(size(image));
    imageBinarisee(indicesTumeur) = 1;
    
    subplot(4,4,i_subplot)
    imagesc(imageBinarisee)
    title(['Isoler la tumeur, methode : ' methode])
    i_subplot = i_subplot+1;
    
    subplot(4,4,i_subplot)
    imagesc(image)
    hold on
    title(['Localiser la tumeur, methode : ' methode])
    i_subplot = i_subplot+1;

    %% Localiser la tumeur dans l'image originale

    s_predit = regionprops(bwlabel(imageBinarisee)); % determine entre autres le centre de masse de la tumeur predite

    rectangle('Position',s_predit.BoundingBox,'EdgeColor', 'r','LineWidth', 3)
    hold off

    %% Distance entre la tumeur dectectee et le masque
    
    disp([methode ' :'])
    
    s_mask = regionprops(bwlabel(mask)); % determine entre autres le centre de masse de la tumeur du mask
    distance_s_predit_s_mask = norm(s_predit.Centroid-s_mask.Centroid,2); % norme 2
    disp(['Distance entre le centre de masse de la tumeur predite et le centre de masse de la tumeur dans le mask : ' num2str(distance_s_predit_s_mask)])
    
    similarite = dice(imageBinarisee, mask); % quantifie la similarite entre deux images
    disp(['Similarite entre notre prediction et le mask = ' num2str(similarite) ])
    fprintf('\n')
    
end

%% rognerImage

function [imageRognee,indicesLignesNoires,indicesColonnesNoires] = rognerImage(image,seuilMax,n,p)
    
    indicesLignesNoires = [];
    for i = 1:n
        valeurs = unique(image(i,:));

        if ((valeurs>=0) & (valeurs<seuilMax))
            indicesLignesNoires = [indicesLignesNoires; i];
        end
    end

    indicesColonnesNoires = [];
    for i = 1:p
        valeurs = unique(image(:,i));

        if (valeurs>=0 & valeurs<seuilMax)
            indicesColonnesNoires = [indicesColonnesNoires; i];
        end
    end

    image(indicesLignesNoires,:) = [];
    image(:,indicesColonnesNoires) = [];
    imageRognee = image;
end