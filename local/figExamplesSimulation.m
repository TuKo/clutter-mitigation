%% This file creates an example image of each method using the simulations
clear all;

% load data
load ..\datasets\Simulated\simulatedData40.mat
artIQ = artifactIQdata;
IQdata = IQdata(:,:,1:20);
clear artifactIQdata;

%% Simulation setup
addpath('ompbox\');
addpath('..\helpers\');

shift = 80; % tissue displacement = 1
art_shift = 10; % clutter displacement = 1/8
downsampleStep = 10; % downsampling for sequence creation
rho = 0.98;
SNR = 30;

M = 33-1;
N = 9-1;

% MCA methods setup
CLUTTER_METHODs = 1; % use svd for selecting the atoms in  fully-online MCA
OMP_ERROR = 2.3;
ATOM_THR = 0.775;
CUTOFF = 0.35;
DICT_SIZE = 2;

% SVF methods setup
tau = 0.75;
alpha = 30;

% Compute the Input and Output windows for all the methods
[zone1Y,zone1X] = meshgrid((565:765)-(shift./downsampleStep),6:15);
[zone2Y,zone2X] = meshgrid((565:765)-(shift./downsampleStep),65:74);
[zoneInY,zoneInX] = meshgrid((565:765)-(shift./downsampleStep),31:50);
zoneInX = zoneInX(:);
zoneInY = zoneInY(:);    
zoneOutX = [zone1X(:); zone2X(:)];
zoneOutY = [zone1Y(:); zone2Y(:)];

[PSNRzoneInY,PSNRzoneInX] = meshgrid((400:900)-(shift./downsampleStep),6:74);
PNSRzoneX = PSNRzoneInX(:);
PNSRzoneY = PSNRzoneInY(:);

LinesSaved = 250:(1090-50);

% Data Sample for learning the clutter dictionary
[~, ~, RFDataSample, RFData, ~] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
ClutterDataSample = RFData - RFDataSample;
clear RFData

% Original Data
[FinalData, Noise, RFdataClean, RFdata, NSTD] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
[~, ~, T] = size(FinalData);

RFCleanImg = RFdataClean(:,:,(T+1)/2);
NoiseImg = Noise(:,:,(T+1)/2);

% Unfiltered image
UNImage = abs(FinalData(:,:,(T+1)/2));
[~,~,CNR_UN,~] = computeSNRandCNR(UNImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_UN, ~] = computeMSE(UNImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);

% Expected image
PFImage = abs(RFCleanImg+NoiseImg);
[~,~,CNR_PF,~] = computeSNRandCNR(PFImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_PF, ~] = computeMSE(PFImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);

% Clean image
CImage = abs(RFCleanImg);

%% Method 1: OFF-MCA
params.P = 0;
params.N = N;
params.M = M;
params.OL_M = params.M;
params.partition = 1.0;
params.parallel = 1;
params.DICT_SIZE = round((params.M+1)*(params.N+1)*1.1);
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
params.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*2)*NSTD*OMP_ERROR;
% learn the clutter dictionary for this size of patches
Dc = KSVDCDoubleField2(ClutterDataSample, params);
% learn the tissue dictionary for this size of patches
Dt = KSVDCDoubleField2(RFDataSample, params);
OFFMCAImage = KSVDCField2SupressClutterGlobalDict(Dt, Dc, FinalData, params);
% Measure performance
[~,~,CNR_OFFMCA,~] = computeSNRandCNR(OFFMCAImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_OFFMCA, ~] = computeMSE(OFFMCAImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);

%% Method 2: TA-MCA
params.P = 0;
params.N = N;
params.M = M;
params.OL_M = params.M;
params.partition = 1.0;
params.parallel = 1;
params.DICT_SIZE = round((params.M+1)*(params.N+1)*2);
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
params.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*2)*NSTD*OMP_ERROR; 
params.cutoff = CUTOFF;
% The clutter dictionary is learned from Method 1.
DtA = KSVDCDoubleField2(FinalData, params);
At = comp(Dc'*DtA, Dc'*Dc, params.M+1);
TAMCAImage = KSVDCField2SupressClutterAdaptiveDict(DtA, Dc, At, FinalData, params);

% Measure performance
[~,~,CNR_TAMCA,~] = computeSNRandCNR(TAMCAImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_TAMCA, ~] = computeMSE(TAMCAImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);


%% Method 3: ON-MCA
params.P = 0;
params.N = N;
params.M = M;
params.OL_M = params.M;
params.partition = 0.8; % Perform Two KSVDs
params.parallel = 1;
params.DICT_SIZE = round((params.M+1)*(params.N+1)*DICT_SIZE);
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
params.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*2)*NSTD*OMP_ERROR; 
params.CLUTTER_METHOD = CLUTTER_METHODs;
params.ATOM_THRs = ATOM_THR;
DtFO = KSVDCDoubleField2(FinalData, params); 
DtDFO = DtFO'*DtFO;
ONMCAImage = KSVDCField2SupressClutter_2(DtFO, DtDFO, FinalData, params);

% Measure performance
[~,~,CNR_ONMCA,~] = computeSNRandCNR(ONMCAImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_ONMCA, ~] = computeMSE(ONMCAImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);


%% Method 4: SVF
midFrame = (size(FinalData,3)-1)/2 +1;
SVFImage = SVF_oneFrame(FinalData, M, N, alpha, tau);
SVFImage = abs(SVFImage(:,:,midFrame));

% Measure performance
[~,~,CNR_SVF,~] = computeSNRandCNR(SVFImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_SVF, ~] = computeMSE(SVFImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);

%% FIR image
FIRImage = abs(FinalData(:,:,midFrame)-(FinalData(:,:,midFrame-1)));
[~,~,CNR_FIR,~] = computeSNRandCNR(FIRImage, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_FIR, ~] = computeMSE(FIRImage, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);


%% save  results for prosperity
save('temp/phd_seminar.mat');

%% Save sample image with CNR areas

UNi = tissueProcessing(UNImage(LinesSaved,:)*1e25,30,0);
lines = min(LinesSaved);
color = 1.0;
min1Y = min(zone1Y(:))-lines;
min1X = min(zone1X(:));
max1Y = max(zone1Y(:))-lines;
max1X = max(zone1X(:));
min2Y = min(zone2Y(:))-lines;
min2X = min(zone2X(:));
max2Y = max(zone2Y(:))-lines;
max2X = max(zone2X(:));
minInY = min(zoneInY(:))-lines;
minInX = min(zoneInX(:));
maxInY = max(zoneInY(:))-lines;
maxInX = max(zoneInX(:));

UNi = drawLine(UNi, [min1X min1Y], [max1X min1Y], color);
UNi = drawLine(UNi, [min1X min1Y-1], [max1X min1Y-1], color);
UNi = drawLine(UNi, [min1X min1Y+1], [max1X min1Y+1], color);
UNi = drawLine(UNi, [max1X min1Y], [max1X max1Y], color);
UNi = drawLine(UNi, [min1X max1Y], [max1X max1Y], color);
UNi = drawLine(UNi, [min1X max1Y+1], [max1X max1Y+1], color);
UNi = drawLine(UNi, [min1X max1Y-1], [max1X max1Y-1], color);
UNi = drawLine(UNi, [min1X min1Y], [min1X max1Y], color);

UNi = drawLine(UNi, [min2X min2Y], [max2X min2Y], color);
UNi = drawLine(UNi, [min2X min2Y-1], [max2X min2Y-1], color);
UNi = drawLine(UNi, [min2X min2Y+1], [max2X min2Y+1], color);
UNi = drawLine(UNi, [max2X min2Y], [max2X max2Y], color);
UNi = drawLine(UNi, [min2X max2Y], [max2X max2Y], color);
UNi = drawLine(UNi, [min2X max2Y+1], [max2X max2Y+1], color);
UNi = drawLine(UNi, [min2X max2Y-1], [max2X max2Y-1], color);
UNi = drawLine(UNi, [min2X min2Y], [min2X max2Y], color);

color = 0.0;
UNi = drawLine(UNi, [minInX minInY], [maxInX minInY], color);
UNi = drawLine(UNi, [minInX minInY-1], [maxInX minInY-1], color);
UNi = drawLine(UNi, [minInX minInY+1], [maxInX minInY+1], color);
UNi = drawLine(UNi, [maxInX minInY], [maxInX maxInY], color);
UNi = drawLine(UNi, [minInX maxInY], [maxInX maxInY], color);
UNi = drawLine(UNi, [minInX maxInY+1], [maxInX maxInY+1], color);
UNi = drawLine(UNi, [minInX maxInY-1], [maxInX maxInY-1], color);
UNi = drawLine(UNi, [minInX minInY], [minInX maxInY], color);

showSimulationImage(UNi);
print([folder '\figureExampleUNboxes.jpg'],'-djpeg');

%% Save images and print performance results
close all;

folder = 'results';

fprintf('PSNR Results:\n');
fprintf('OFF-MCA = %g, \tTA-MCA = %g, \tON-MCA = %g, \tSVF = %g, \tUN = %g, \tPF = %g\n',PSNR_OFFMCA,PSNR_TAMCA,PSNR_ONMCA,PSNR_SVF,PSNR_UN,PSNR_PF);
fprintf('CNR Results:\n');
fprintf('OFF-MCA = %g, \tTA-MCA = %g, \tON-MCA = %g, \tSVF = %g, \tUN = %g, \tPF = %g\n',CNR_OFFMCA,CNR_TAMCA,CNR_ONMCA,CNR_SVF,CNR_UN,CNR_PF);

showSimulationImage(tissueProcessing(CImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleClean.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(UNImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleUN.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(PFImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExamplePF.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(OFFMCAImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleOFFMCA.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(TAMCAImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleTAMCA.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(ONMCAImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleONMCA.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(SVFImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleSVF.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(FIRImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleFIR.jpg'],'-djpeg');

showSimulationImage(tissueProcessing(OFFMCAImage(LinesSaved,:)*1e25-PFImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleDiffOFFMCA.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(TAMCAImage(LinesSaved,:)*1e25-PFImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleDiffTAMCA.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(ONMCAImage(LinesSaved,:)*1e25-PFImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleDiffONMCA.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(SVFImage(LinesSaved,:)*1e25-PFImage(LinesSaved,:)*1e25,30,0));
print([folder '\figureExampleDiffSVF.jpg'],'-djpeg');

%% Save the dictionaries

%MCA Methods
dN = ceil(sqrt(size(Dc,2)*3))-2;
dM = ceil(size(Dc,2)/dN);
imDC = showdict(tissueProcessing(Dc,3.0,0),[M+1 N+1],dM, dN,'highcontrast','whitelines');
imwrite(imDC,'results/Dictionary-OFFMCA-DC.png');


dN = ceil(sqrt(size(Dt,2)*3))-2;
dM = ceil(size(Dt,2)/dN);
imDT = showdict(tissueProcessing(Dt,3.0,0),[M+1 N+1],dM, dN,'highcontrast','whitelines');%
imwrite(imDT, 'results/Dictionary-OFFMCA-DT.png');


dN = ceil(sqrt(size(DtA,2)*3))-2;
dM = ceil(size(DtA,2)/dN);
imDTA = showdict(tissueProcessing(DtA,3.0,0),[M+1 N+1],dM, dN-1,'highcontrast','whitelines');%
imwrite(imDTA, 'results/Dictionary-TAMCA-DT-Full.png');

DtAA = DtA(:,~(sum(abs(DtA - Dc*At).^2)<CUTOFF));
dN = ceil(sqrt(size(DtAA,2)*3))-2;
dM = ceil(size(DtAA,2)/dN);
imDTA = showdict(tissueProcessing(DtAA,3.0,0),[M+1 N+1],dM, dN,'highcontrast','whitelines');%
imwrite(imDTA, 'results/Dictionary-TAMCA-DT-Clean.png');

dN = ceil(sqrt(size(DtFO,2)*3))-2;
dM = ceil(size(DtFO,2)/dN);
imDTON = showdict(tissueProcessing(DtFO,3.0,0),[M+1 N+1],dM, dN-1,'highcontrast','whitelines');
imwrite(imDTON,'results/Dictionary-ONMCA-Full.png');

staticAtom = findClutterAtoms(DtFO, ATOM_THR, CLUTTER_METHODs, M+1, N+1, 1);
DtON = DtFO(:,~staticAtom);
dN = ceil(sqrt(size(DtON,2)*3))-2;
dM = ceil(size(DtON,2)/dN);
imDTON = showdict(tissueProcessing(DtON,3.0,0),[M+1 N+1],dM, dN-1,'highcontrast','whitelines');
imwrite(imDTON,'results/Dictionary-ONMCA-DT.png');

DcON = DtFO(:,staticAtom);
dN = ceil(sqrt(size(DcON,2)*3))-2;
dM = ceil(size(DcON,2)/dN);
imDCON = showdict(tissueProcessing(DcON,3.0,0),[M+1 N+1],dM, dN,'highcontrast','whitelines');
imwrite(imDCON,'results/Dictionary-ONMCA-DC.png');

% SVF
k = 650; % 565:765
j = 40;
l = (T-1)/2+1;
patch = reshape(1e25*FinalData(k+(-M/2:M/2),j,l+ ((-N/2):(N/2))), [M+1, N+1]); %#ok<IJCL>
[u,s,v] = svd(patch,0);
s1 = s;
DCSVF(:,1) = patch(:);
for i = 1:size(u,2)
    ccc = u(:,i)*v(i,:);
    DCSVF(:,i+1) = ccc(:);
end
imSVFC = showdict(tissueProcessing(DCSVF,3.0,0),[M+1 N+1],1, size(u,2),'highcontrast','whitelines');
imwrite(imSVFC,'results/Dictionary-SVF-Clutter.png');

k = 650; % 565:765
j = 10;
patch = reshape(1e25*FinalData(k+(-M/2:M/2),j,l+ ((-N/2):(N/2))), [M+1, N+1]); %#ok<IJCL>
[u,s,v] = svd(patch,0);
s2 = s;
DTSVF(:,1) = patch(:);
for i = 1:size(u,2)
    ccc = u(:,i)*v(i,:);
    DTSVF(:,i+1) = ccc(:);
end
imSVFT = showdict(tissueProcessing(DTSVF,3.0,0),[M+1 N+1],1, size(u,2),'highcontrast','whitelines');
imwrite(imSVFT,'results/Dictionary-SVF-Tissue.png');

fontsz = 20;
figure
plot(diag(s1)./trace(abs((s1))),'-b','linewidth',2,'markersize',15);
hold on;
plot(diag(s2)./trace(abs((s2))),'-r','linewidth',2,'markersize',15);
plot(tau*ones(size(u,2),1),'--g','linewidth',2,'markersize',15);
xlabel('i','FontSize',fontsz);
ylabel('Normalized Singular Value \sigma_i','FontSize',fontsz);
axis([1,numel(diag(s)) 0,1]);
legend({'Clutter','Tissue'},'FontSize',fontsz,'Location','NorthEast');
set(gca,'FontSize',fontsz);
print(['results\figureSVDs.png'],'-dpng');

%% Different CUTOFFs TA-MCA results
folder = 'results';
params.P = 0;
params.N = N;
params.M = M;
params.OL_M = params.M;
params.partition = 1.0;
params.parallel = 1;
params.DICT_SIZE = round((params.M+1)*(params.N+1)*2);
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
params.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*2)*NSTD*OMP_ERROR; 



TAMCAImageCUT1 = TAMCAImage;
params.cutoff = 0.4;
TAMCAImageCUT2 = KSVDCField2SupressClutterAdaptiveDict(DtA, Dc, At, FinalData, params);
params.cutoff = 0.6;
TAMCAImageCUT3 = KSVDCField2SupressClutterAdaptiveDict(DtA, Dc, At, FinalData, params);
params.cutoff = 0.9;
TAMCAImageCUT4 = KSVDCField2SupressClutterAdaptiveDict(DtA, Dc, At, FinalData, params);

% Measure performance
showSimulationImage(tissueProcessing(TAMCAImageCUT1(LinesSaved,:)*1e25,30,0));
print([folder '\figureTAMCAtau1.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(TAMCAImageCUT2(LinesSaved,:)*1e25,30,0));
print([folder '\figureTAMCAtau2.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(TAMCAImageCUT3(LinesSaved,:)*1e25,30,0));
print([folder '\figureTAMCAtau3.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(TAMCAImageCUT4(LinesSaved,:)*1e25,30,0));
print([folder '\figureTAMCAtau4.jpg'],'-djpeg');

%% Different CUTOFFs ON-MCA results
params.P = 0;
params.N = N;
params.M = M;
params.OL_M = params.M;
params.partition = 0.8; % Perform Two KSVDs
params.parallel = 1;
params.DICT_SIZE = round((params.M+1)*(params.N+1)*DICT_SIZE);
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
params.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*2)*NSTD*OMP_ERROR; 
params.CLUTTER_METHOD = CLUTTER_METHODs;
params.ATOM_THRs = ATOM_THR;


params.ATOM_THRs = 0.6;
ONMCAImageCUT1 = KSVDCField2SupressClutter_2(DtFO, DtDFO, FinalData, params);
ONMCAImageCUT2 = ONMCAImage;
params.ATOM_THRs = 0.85;
ONMCAImageCUT3 = KSVDCField2SupressClutter_2(DtFO, DtDFO, FinalData, params);

% Measure performance
showSimulationImage(tissueProcessing(ONMCAImageCUT1(LinesSaved,:)*1e25,30,0));
print([folder '\figureONMCAbeta1.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(ONMCAImageCUT2(LinesSaved,:)*1e25,30,0));
print([folder '\figureONMCAbeta2.jpg'],'-djpeg');
showSimulationImage(tissueProcessing(ONMCAImageCUT3(LinesSaved,:)*1e25,30,0));
print([folder '\figureONMCAbeta3.jpg'],'-djpeg');

