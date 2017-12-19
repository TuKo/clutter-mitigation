%% send dataset1

c = parcluster('GIPCluster');
test = 'dual1';

N=15;

[I1clutter, nMLAs] = readDataset('data\clutter1-1\');
I1clutter = stb(I1clutter, nMLAs, nMLAs/2);
[I2clutter, nMLAs] = readDataset('data\clutter2-2\');
I2clutter = stb(I2clutter, nMLAs, nMLAs/2);
[I1, I2, nMLAs] = readDatasetDualHarm(['data\' test '\']);
I1 = stb(I1, nMLAs, nMLAs/2);
I2 = stb(I2, nMLAs, nMLAs/2);


D1C = [];
D2C = [];
D2T = [];
D1T = [];

j1 = batch(c,@experimentTAMCA,nargout(@experimentTAMCA),{I1, I2, I1clutter, I2clutter, D1C, D2C, D1T, D2T},'Pool',7,'AttachedFiles',{'ompbox\'});


%% send dataset2

test = 'dual2';

[I1, I2, nMLAs] = readDatasetDualHarm(['data\' test '\']);
I1 = stb(I1, nMLAs, nMLAs/2);
I2 = stb(I2, nMLAs, nMLAs/2);

D1C = [];
D2C = [];
D2T = [];
D1T = [];

j2 = batch(c,@experimentTAMCA,nargout(@experimentTAMCA),{I1, I2, I1clutter, I2clutter, D1C, D2C, D1T, D2T},'Pool',7,'AttachedFiles',{'ompbox\'});

%% wait them to finish

wait(j1);
wait(j2);

%% read dataset1
% OUTPUTS: out1,out2,D1C,D2C,D1T,D2T,tissue1,tissue2

test = 'dual1';

% Convert input sizes
[I1, I2, nMLAs] = readDatasetDualHarm(['data\' test '\']);
I1 = stb(I1, nMLAs, nMLAs/2);
I2 = stb(I2, nMLAs, nMLAs/2);
I1N = (max(abs(I2(:)))/ max(abs(I1(:))) )*I1;
GAIN1 = -40;
GAIN2 = GAIN1+20;
H = 400;
W = 400;
[height,width,nFrames] = size(I1);
startFrame = (N-1)/2+1;
stopFrame = nFrames - (N-1)/2;
I1N = I1N(:,:,startFrame:stopFrame);
I2N = I2(:,:,startFrame:stopFrame);

r = fetchOutputs(j1);
out1 = r{1};
out2 = r{2};
D1C = r{3};
D2C = r{4};
D1T = r{5};
D2T = r{6};
clear r;
clear mixBefore mixAfter originals results allmix allresult

save(['results/mat-nosync/TAMCA-' test '.mat']);

mixBefore = scanConvertMovie(tissueProcessing(0.5*I1N+0.5*I2N,30,GAIN2),H,W);
mixAfter  = scanConvertMovie(tissueProcessing(0.5*out1+0.5*out2,30,GAIN2),H,W);
originals = [scanConvertMovie(tissueProcessing(I1N,30,GAIN2),H,W), scanConvertMovie(tissueProcessing(I2N,30,GAIN2),H,W)];
results   = [scanConvertMovie(tissueProcessing(out1,30,GAIN2),H,W), scanConvertMovie(tissueProcessing(out2,30,GAIN2),H,W)];
allmix    = [originals, mixBefore];
allresult = [results,   mixAfter];

writeVideo(['results\' test '-TAMCA-mix-nostb.avi'], convertMovie([mixBefore, mixAfter]), 20);
writeVideo(['results\' test '-TAMCA-results-nostb.avi'], convertMovie([originals; results]), 20);
writeVideo(['results\' test '-TAMCA-all-nostb.avi'], convertMovie([allmix; allresult]), 20);


%% read dataset2
% OUTPUTS: out1,out2,D1C,D2C,D1T,D2T,tissue1,tissue2

test = 'dual2';

% Convert input sizes
[I1, I2, nMLAs] = readDatasetDualHarm(['data\' test '\']);
I1 = stb(I1, nMLAs, nMLAs/2);
I2 = stb(I2, nMLAs, nMLAs/2);
I1N = (max(abs(I2(:)))/ max(abs(I1(:))) )*I1;
GAIN1 = -40;
GAIN2 = GAIN1+20;
H = 400;
W = 400;
[height,width,nFrames] = size(I1);
startFrame = (N-1)/2+1;
stopFrame = nFrames - (N-1)/2;
I1N = I1N(:,:,startFrame:stopFrame);
I2N = I2(:,:,startFrame:stopFrame);

r = fetchOutputs(j2);
out1 = r{1};
out2 = r{2};
D1C = r{3};
D2C = r{4};
D1T = r{5};
D2T = r{6};
clear r;
clear mixBefore mixAfter originals results allmix allresult

save(['results/mat-nosync/TAMCA-' test '.mat']);

mixBefore = scanConvertMovie(tissueProcessing(0.5*I1N+0.5*I2N,30,GAIN2),H,W);
mixAfter  = scanConvertMovie(tissueProcessing(0.5*out1+0.5*out2,30,GAIN2),H,W);
originals = [scanConvertMovie(tissueProcessing(I1N,30,GAIN2),H,W), scanConvertMovie(tissueProcessing(I2N,30,GAIN2),H,W)];
results   = [scanConvertMovie(tissueProcessing(out1,30,GAIN2),H,W), scanConvertMovie(tissueProcessing(out2,30,GAIN2),H,W)];
allmix    = [originals, mixBefore];
allresult = [results,   mixAfter];

writeVideo(['results\' test '-TAMCA-mix-nostb.avi'], convertMovie([mixBefore, mixAfter]), 20);
writeVideo(['results\' test '-TAMCA-results-nostb.avi'], convertMovie([originals; results]), 20);
writeVideo(['results\' test '-TAMCA-all-nostb.avi'], convertMovie([allmix; allresult]), 20);

