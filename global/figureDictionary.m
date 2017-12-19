%% load data
load datasets\Simulated\simulatedData40
artIQ = artifactIQdata;
IQdata = IQdata(:,:,1:20);
clear artifactIQdata;

%%
Mvalue = 33;
Nvalue = 9;
OMP_ERROR = 2.3;

shift = 80; % tissue displacement = 1
art_shift = 10; % clutter displacement = 1/8
downsampleStep = 10; % downsampling for sequence creation
rho = 0.98;
SNR = 30;

DICT_SIZEs = 2; % dictionary redundancy
ATOMS_KSVDs = 30;
OMP_MAX_ATOMSs = 40; % max atoms for OMP with error threshold
CLUTTER_METHOD = 1; % use svd for selecting the atoms


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

params.partition =0.8;
params.parallel = 1;
params.P = 0; % radial size of patches (usually 0, only axial-time)
params.OL_P = params.P/2;
params.ATOMS_KSVD = ATOMS_KSVDs;
params.N = Nvalue-1;
params.M = Mvalue-1;
params.OL_M = params.M;
params.ATOMS_KSVD = round(0.10*(params.M+1)*(params.N+1)*(params.P+1));
params.DICT_SIZE = round((params.M+1)*(params.N+1)*(params.P+1)*DICT_SIZEs);

[FinalData, Noise, RFdataClean, RFdata, NSTD] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
[~, ~, T] = size(FinalData);


IQDictionary = KSVDCDoubleField2(FinalData, params); % Two KSVDs

%% clean the sequence with MCA
params.ATOM_THRs = 0.75;
params.OMP_MAX_ATOMSs = OMP_MAX_ATOMSs;
params.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*(params.P+1)*2)*NSTD*OMP_ERROR; %*std(FinalData(:));

DtD = IQDictionary'*IQDictionary;
imgMCA = KSVDCField2SupressClutter_2(IQDictionary, DtD, FinalData, params);
[~,~,CNR_MCA,~] = computeSNRandCNR(imgMCA, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_MCA, ~] = computeMSE(imgMCA, abs(RFdataClean(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);

%% clean the sequence with SVF
MvalSVF = 33-1;
NvalSVF = 7-1;
alpha = 25;
tau = 0.75;
midFrame = 10;
imgSVF = abs(SVF_oneFrame(FinalData, MvalSVF, NvalSVF, alpha, tau));
imgSVF = imgSVF(:,:,midFrame);
[~,~,CNR_SVF,~] = computeSNRandCNR(imgSVF, zoneInX, zoneInY, zoneOutX, zoneOutY);
[PSNR_SVF, ~] = computeMSE(imgSVF, abs(RFdataClean(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);

%% save images and diff images


%% save the dictionary image
ATOM_THR = 0.75;
stationaryAtom = findClutterAtoms(IQDictionary, ATOM_THR, CLUTTER_METHOD, params.M+1, params.N+1, params.P+1);

ImageD = showDictionary(IQDictionary, params.M+1, params.N+1, params.P+1,[13, 46]);
ImageDc = showDictionary(IQDictionary(:,stationaryAtom), params.M+1, params.N+1, params.P+1,[2, 6]);
ImageDt = showDictionary(IQDictionary(:,~stationaryAtom), params.M+1, params.N+1, params.P+1,[13, 46]);

imwrite(ImageD,'results\figureD.png');
imwrite(ImageDt,'results\figureDt.png');
imwrite(ImageDc,'results\figureDc.png');

%% save atom images
Dc = IQDictionary(:,stationaryAtom);
atom = reshape(tissueProcessing((Dc(:,3)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomClutter1.png');
atom = reshape(tissueProcessing((Dc(:,8)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomClutter2.png');
atom = reshape(tissueProcessing((Dc(:,9)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomClutter3.png');
atom = reshape(tissueProcessing((Dc(:,2)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomClutter4.png');

Dt = IQDictionary(:,~stationaryAtom);
atom = reshape(tissueProcessing((Dt(:,57)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomTissue1.png');
atom = reshape(tissueProcessing((Dt(:,end-1)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomTissue2.png');
atom = reshape(tissueProcessing((Dt(:,46*3+1)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomTissue3.png');
atom = reshape(tissueProcessing((Dt(:,46*2+2)),1,0),params.M+1, params.N+1);
imwrite(atom,'results\figureAtomTissue4.png');
%% Intentar mostrar el diccionario en una misma imagen
close all
ImageD = showDictionary([IQDictionary(:,~stationaryAtom),IQDictionary(:,stationaryAtom)] , params.M+1, params.N+1, params.P+1,[13, 46],sum(~stationaryAtom));
imwrite(ImageD,'results\figureD.png');