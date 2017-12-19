function [out1,out2,D1C,D2C,D2,tissue1,tissue2] = experiment1(I1, I2, I1clutter,I2clutter)
% addpath('ompbox/');

% setup for the problem and algorithms
N = 15; % number of frames
M = 15; % number of axial (vertical) elements

% normalize I1
I1 = (max(abs(I2(:)))/ max(abs(I1(:))) )*I1;

% Learn the clutter dictionary on 1st Harmonic
% [I1clutter, nMLAs] = readDataset('data\clutter1-1\');

params.M = M-1;
params.N = N-1;
params.OL_M = params.M;
params.partition = 1.0;
params.parallel = 1;

params.OMP_ERROR = 0;
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.DICT_SIZE = round((params.M+1)*(params.N+1)*1.1);
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));

D1C = LearnDictionary(I1clutter, params);

% Learn the clutter dictionary on 2nd Harmonic
% [I2clutter, nMLAs] = readDataset('data\clutter2-2\');
D2C = LearnDictionary(I2clutter, params);


% Learn the tissue dictionary on clean 2nd Harmonic
% This dictionary is bigger to clean up the clutter part if there is such
params.DICT_SIZE = round((params.M+1)*(params.N+1)*2);
% [I2tissue, nMLAs] = readDataset('data\tissue1-2\');
% D2 = LearnDictionary(I2tissue, params);
% Clean the clutter part from D2 (TA-MCA)
% A2 = comp(D2C'*D2, D2C'*D2C, params.M+1);

% Clean and fuse the harmonics
% [I1, I2, nMLAs] = readDatasetDualHarm('data\dual1\');

D2 = LearnDictionary(I2, params);
% Clean the clutter part from D2 (TA-MCA)
A2 = comp(D2C'*D2, D2C'*D2C, params.M+1);

[height,width,nFrames] = size(I1);    
Mask = true(height,width);
startFrame = (N-1)/2+1;
stopFrame = nFrames - (N-1)/2;

% params.CUTTOFFs = 0.35;
params.CUTTOFFs = 0.2;
params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
[out1,out2,tissue1,tissue2] = FusionTA(D1C, D2C, D2, A2, I1, I2, params, startFrame, stopFrame, Mask);

% GAIN = -30;
% H = 400;
% W = 400;
% I1det = tissueProcessing(I1(:,:,startFrame:stopFrame),30,GAIN);
% I2det = tissueProcessing(I2(:,:,startFrame:stopFrame),30,GAIN);
% O1det = tissueProcessing(out1,30,GAIN);
% O2det = tissueProcessing(out2,30,GAIN);

%writeVideo(['results\output.avi'], convertMovie([scanConvertMovie(I1det,H,W) scanConvertMovie(I2det,H,W); scanConvertMovie(O1det,H,W) scanConvertMovie(O2det,H,W)]), 20);
%writeVideo(['results\output2-mix.avi'], convertMovie([scanConvertMovie(0.5*stb(I1det,nMLAs,2)+0.5*stb(I2det,nMLAs,2),H,W), scanConvertMovie(0.5*stb(r{1},nMLAs,2)+0.5*stb(r{2},nMLAs,2),H,W)]), 20);
