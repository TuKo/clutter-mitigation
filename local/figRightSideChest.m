% Figure: create an image of a scan from the right side of the chest.

addpath('ompbox\');
addpath('..\helpers\');
addpath('..\');


[I, nMLAs] = readDataset('..\datasets_new\p1s1\');
frame = 5;

writefolder = '.\results\';


im = scanConvertMovie(tissueProcessing(stb(I(:,:,frame),nMLAs, nMLAs/2),30,-50),400,400);

imwrite(fliplr(im), [writefolder 'figureClutterExample.jpg']);

rmpath('..\ompbox\');
rmpath('..\helpers\');
rmpath('..\');
