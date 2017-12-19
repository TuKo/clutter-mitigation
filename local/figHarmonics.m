% Figure: create an image of the 1st and 2nd harmonics of the same
% acquisition.

addpath('ompbox\');
addpath('..\helpers\');
addpath('..\');


[Ifst, Isnd, nMLAs] = readDatasetDualHarm('..\datasets\p9s2\');
frame = 15;

writefolder = '.\results\';


im1st = scanConvertMovie(tissueProcessing(stb(Ifst(:,:,frame),nMLAs, nMLAs/2),40,-30),400,500);
im2nd = scanConvertMovie(tissueProcessing(stb(Isnd(:,:,frame),nMLAs, nMLAs/2),30,-20),400,500);

imwrite(fliplr(im1st), [writefolder 'figureFstExample.jpg']);
imwrite(fliplr(im2nd), [writefolder 'figureSndExample.jpg']);

rmpath('..\ompbox\');
rmpath('..\helpers\');
rmpath('..\');
