Iun = im2double(imread('results\exper8-UN-Frm6.jpg'));
Isvf = im2double(imread('results\exper8-SVF-Frm6.jpg'));
Ionmca = im2double(imread('results\exper8-ONMCA-Frm6.jpg'));
Ioffmca = im2double(imread('results\exper8-OFFMCA-Frm6.jpg'));
Itamca = im2double(imread('results\exper8-TAMCA-Frm6.jpg'));

close all

% figure, imshow(Isvf);
% figure, imshow(Ionmca);
% figure, imshow(Ioffmca);
% figure, imshow(Itamca);

figure, imshow(Iun);

color = [1, 0, 0];

Isvf2 = repmat(Isvf, [1,1,3]);
Isvf2 = drawArrow(Isvf2, [250 203], [253 184], color);
Isvf2 = drawArrow(Isvf2, [150 350], [149 368], color);
Isvf2 = drawArrow(Isvf2, [179 126], [179 110], color);
Isvf2 = drawArrow(Isvf2, [129 303], [114 299], color);
figure, imshow(Isvf2);
imwrite(Isvf2, 'results\figureExampleRealSVF.png');

Ionmca2 = repmat(Ionmca, [1,1,3]);
Ionmca2 = drawArrow(Ionmca2, [69 315], [76 331], color);
Ionmca2 = drawArrow(Ionmca2, [140 278], [164 278], color);
Ionmca2 = drawArrow(Ionmca2, [133 303], [114 299], color);
figure, imshow(Ionmca2);
imwrite(Ionmca2, 'results\figureExampleRealONMCA.png');

Ioffmca2 = repmat(Ioffmca, [1,1,3]);
% Ioffmca2 = drawArrow(Ioffmca2, [69 315], [76 331], color);
Ioffmca2 = drawArrow(Ioffmca2, [246 297], [225 297], color);
Ioffmca2 = drawArrow(Ioffmca2, [140 361], [144 378], color);
figure, imshow(Ioffmca2);
imwrite(Ioffmca2, 'results\figureExampleRealOFFMCA.png');

Itamca2 = repmat(Itamca, [1,1,3]);
% Itamca2 = drawArrow(Itamca2, [69 315], [76 331], color);
Itamca2 = drawArrow(Itamca2, [246 297], [225 297], color);
Itamca2 = drawArrow(Itamca2, [140 361], [144 378], color);
figure, imshow(Itamca2);
imwrite(Itamca2, 'results\figureExampleRealTAMCA.png');
