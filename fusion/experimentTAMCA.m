function [out1,out2,D1C,D2C,D1T,D2T,tissue1,tissue2] = experimentTAMCA(I1, I2, I1clutter, I2clutter, D1C, D2C, D1T, D2T)
% Learn a clutter dictionary for 1st and one for 2nd harmonics.
% Learn an adaptive tissue dictionary for each harmonic separately (as
% TAMCA - Eilat paper) (after normalization of the first harmonic).

out1=[];out2=[];D1C=[];D2C=[];D1T=[];D2T=[];tissue1=[];tissue2=[];

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
params.frames = false;
if isempty(D1C)
    D1C = LearnDictionary(I1clutter, params);
end

% Learn the clutter dictionary on 2nd Harmonic
% [I2clutter, nMLAs] = readDataset('data\clutter2-2\');
if isempty(D2C)
    D2C = LearnDictionary(I2clutter, params);
end

% Learn the tissue dictionary on clean 2nd Harmonic
% This dictionary is bigger to clean up the clutter part if there is such
params.DICT_SIZE = round((params.M+1)*(params.N+1)*2);
% [I2tissue, nMLAs] = readDataset('data\tissue1-2\');
% D2 = LearnDictionary(I2tissue, params);
% Clean the clutter part from D2 (TA-MCA)

% Clean and fuse the harmonics
if isempty(D1T)
    params.frames = true; % learn 6 frames and not 3 as before
    params.iternum = 10; % more iterations for KSVD in TISSUE
    D1T = LearnDictionary(I1, params);
end

if isempty(D2T)
    params.frames = true; % learn 6 frames and not 3 as before
    params.iternum = 10; % more iterations for KSVD in TISSUE
    D2T = LearnDictionary(I2, params);
end

% Clean the clutter part from D2 (TA-MCA)
A1 = comp(D1C'*D1T, D1C'*D1C, params.M+1);
A2 = comp(D2C'*D2T, D2C'*D2C, params.M+1);

[height,width,nFrames] = size(I1);    
Mask = true(height,width);
startFrame = (N-1)/2+1;
stopFrame = nFrames - (N-1)/2;

% params.CUTTOFFs = 0.35;
params.CUTTOFFs = 0.2;
params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
[out1,tissue1] = TAMCA(D1C, D1T, A1, I1, params, startFrame, stopFrame, Mask);
[out2,tissue2] = TAMCA(D2C, D2T, A2, I2, params, startFrame, stopFrame, Mask);
