addpath('ompbox/');
addpath('../helpers/');
addpath('../');

addpath('ompbox\');
addpath('..\helpers\');
addpath('..\');

datasets = {'..\datasets\set1\','..\datasets\set5\','..\datasets\set6\',...
    '..\datasets\set7\','..\datasets\set8\'};

% datasets = {'..\datasets_new\p3s25\','..\datasets_new\p3s27\','..\datasets_new\p3s6\','..\datasets_new\p3s8\'};
% datasets = {'..\datasets\set8\'};
% datasets = {'..\datasets_new\p3s6\'};
% clutterset = '..\datasets\p5s14\';
% tissueset = '..\datasets\p4s1\';

% datasets = {'..\datasets_new\p3s25\'};
clutterset = '..\datasets_new\p1s1\';
tissueset = '..\datasets_new\p2s1\';

Nvalue = 15;
Mvalue = 15;
BETAs = 0.05:0.025:0.65;
CUTOFFs = 0.05:0.025:0.95;
TAUs = 0.05:0.025:0.65;
% BETAs = [0.3 0.35];
% TAUs = BETAs;
% CUTOFFs = [0.2 0.4];

%%

tic
[CNR_ONMCA, CNR_OFFMCA, CNR_TAMCA, CNR_SVF, CNR_UN] = ...
    figRealData(datasets, clutterset, tissueset, Mvalue, Nvalue, TAUs, BETAs, CUTOFFs);
toc

save('results\realdata.mat');

%% merge results of 2 runs
load('results\realdata.mat');
PSNR_OFFMCA2 = PSNR_OFFMCA;
PSNR_ONMCA2 = PSNR_ONMCA;
PSNR_TAMCA2 = PSNR_TAMCA;
PSNR_SVF2 = PSNR_SVF;
CNR_OFFMCA2 = CNR_OFFMCA;
CNR_ONMCA2 = CNR_ONMCA;
CNR_TAMCA2 = CNR_TAMCA;
CNR_SVF2 = CNR_SVF;
CNR_UN2 = CNR_UN;
datasets2 = datasets;

load('results\realdata2.mat');
PSNR_OFFMCA = [PSNR_OFFMCA; PSNR_OFFMCA2];
PSNR_ONMCA  = [PSNR_ONMCA;  PSNR_ONMCA2];
PSNR_TAMCA  = [PSNR_TAMCA; PSNR_TAMCA2];
PSNR_SVF    = [PSNR_SVF; PSNR_SVF2];
CNR_OFFMCA = [CNR_OFFMCA; CNR_OFFMCA2];
CNR_ONMCA  = [CNR_ONMCA; CNR_ONMCA2];
CNR_TAMCA  = [CNR_TAMCA; CNR_TAMCA2];
CNR_SVF    = [CNR_SVF; CNR_SVF2];
CNR_UN     = [CNR_UN; CNR_UN2];
datasets = [datasets, datasets2];

save('results\realdatamixed.mat');

%% Analyze results
% ONMCA -> 17
% SVF   -> 16
% TAMCA -> 13
mean(PSNR_TAMCA - repmat(PSNR_OFFMCA,[1,numel(CUTOFFs)]))
mean(PSNR_ONMCA - repmat(PSNR_OFFMCA,[1,numel(BETAs)]))
mean(PSNR_SVF - repmat(PSNR_OFFMCA,[1,numel(TAUs)]))

fprintf('MEAN error with similar PSNR\n');
mean(CNR_OFFMCA      - CNR_UN)
mean(CNR_TAMCA(:,13) - CNR_UN)
mean(CNR_ONMCA(:,17) - CNR_UN)
mean(CNR_SVF(:,16)   - CNR_UN)

%% 
figure;

fontsz = 20;

plot(mean(PSNR_ONMCA), mean(CNR_ONMCA - repmat(CNR_UN,[1,size(CNR_ONMCA,2)])),'-r','linewidth',2,'markersize',15)
hold on
plot(mean(PSNR_TAMCA), mean(CNR_TAMCA - repmat(CNR_UN,[1,size(CNR_TAMCA,2)])),'-B','linewidth',2,'markersize',15)
plot(mean(PSNR_SVF), mean(CNR_SVF - repmat(CNR_UN,[1,size(CNR_SVF,2)])),'-g','linewidth',2,'markersize',15)
plot(mean(PSNR_OFFMCA),mean(CNR_OFFMCA-CNR_UN),'+-k','linewidth',2,'markersize',15);
xlabel('PSNR','FontSize',fontsz)
ylabel('CNR','FontSize',fontsz)
legend('ON-MCA','TA-MCA','SVF','OFF-MCA','FontSize',fontsz);
set(gca,'FontSize',fontsz);
print(['results\PSNR-vs-CNR.png'],'-dpng');

%% Create the output for this values

algo =[];
k = 0;

k = k+1;
algo{k}.method = 'UN';

k = k+1;
algo{k}.method = 'ONMCA';
algo{k}.best_BETA = BETAs(17);

k = k+1;
algo{k}.method = 'SVF';
algo{k}.best_TAU = TAUs(16);

k = k+1;
algo{k}.method = 'OFFMCA';

k = k+1;
algo{k}.method = 'TAMCA';
algo{k}.best_CUTOFF = CUTOFFs(13);


Nvalue = 15;
Mvalue = 15;

clutterset = '..\datasets_new\p1s1\';
tissueset = '..\datasets_new\p2s1\';

datasets = {'..\datasets\set1\','..\datasets\set5\','..\datasets\set6\',...
    '..\datasets\set7\','..\datasets\set8\','..\datasets_new\p3s25\','..\datasets_new\p3s27\','..\datasets_new\p3s6\','..\datasets_new\p3s8\'};

output.GAIN = [-55, -20, -55, -55, -55, -30, -30, -30, -30];
output.visualFrames = [];
%output.GAIN = [-55, -20, -55, -55, -55, -30, -30, -30, -30];
% output.visualFrames = 1:10;
output.folder = 'results\';

tic
figExamplesRealData(datasets, clutterset, tissueset, output, Nvalue, Mvalue, algo);
toc

%% write example of UN with IN-OUT boxes
addpath('ompbox\');
addpath('..\helpers\');
addpath('..\');

[seq,nMLAs] = readDataset('..\datasets_new\p3s6\');
Nvalue = 15; Mvalue = 15;
[height,width,nFrames] = size(seq);
startFrame = max(Nvalue-1)/2+1;
stopFrame = nFrames - max(Nvalue-1)/2;
seq_out = stb(seq(:,:,startFrame:stopFrame),nMLAs,2);
frame = 6;

[InX, InY, OutX, OutY] = readMeasurementWindow('..\datasets_new\p3s6');
InL = min(InX(:))/2;
InR = max(InX(:))-min(InX(:))/2;
InT = min(InY(:));
InB = max(InY(:));
OutL = min(OutX(:))/2;
OutR = max(OutX(:))-min(OutX(:))/2;
OutT = min(OutY(:));
OutB = max(OutY(:));

sampleImage = tissueProcessing(seq_out(:,:,frame),30,-30);
sampleImage(InT:InB,InL) = 1;
sampleImage(InT:InB,InR) = 1;
sampleImage(InT,InL:InR) = 1;
sampleImage(InB,InL:InR) = 1;
sampleImage(OutT:OutB,OutL) = 0;
sampleImage(OutT:OutB,OutR) = 0;
sampleImage(OutT,OutL:OutR) = 0;
sampleImage(OutB,OutL:OutR) = 0;
sampleImage(InT-1:InB+1,InL-1) = 1;
sampleImage(InT-1:InB+1,InR+1) = 1;
sampleImage(InT-1,InL-1:InR+1) = 1;
sampleImage(InB+1,InL-1:InR+1) = 1;
sampleImage(OutT-1:OutB+1,OutL-1) = 0;
sampleImage(OutT-1:OutB+1,OutR+1) = 0;
sampleImage(OutT-1,OutL-1:OutR+1) = 0;
sampleImage(OutB+1,OutL-1:OutR+1) = 0;
imshow(sampleImage)
imwrite(scanConvert(sampleImage,400,400),[ 'results\exper8-Frm6-UN2.jpg']);

%%
figure
plot(CNR_OFFMCA-CNR_UN,'-b');
hold on
plot(max(CNR_SVF,[],2)-CNR_UN,'-c');
plot(max(CNR_ONMCA,[],2)-CNR_UN,'-r');
plot(max(CNR_TAMCA,[],2)-CNR_UN,'-g');
% plot(CNR_UN,'--k');

mean(CNR_OFFMCA-CNR_UN)
mean(max(CNR_SVF,[],2)-CNR_UN)
mean(max(CNR_ONMCA,[],2)-CNR_UN)
mean(max(CNR_TAMCA,[],2)-CNR_UN)

[~,idx] = max(CNR_SVF,[],2)
TAUs(12)

[~,idx] = max(CNR_ONMCA,[],2)
BETAs(12)

[~,idx] = max(CNR_TAMCA,[],2)
CUTOFFs(24)

figure
plot(CNR_OFFMCA-CNR_UN,'-b');
hold on
plot(CNR_SVF(:,12)-CNR_UN,'-c');
plot(CNR_ONMCA(:,12)-CNR_UN,'-r');
plot(CNR_TAMCA(:,24)-CNR_UN,'-g');



