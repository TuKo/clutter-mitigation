function [FinalData, Noise, RFdataClean, RFdata, NSTD] = createSequence(shift, art_shift, downsampleStep, varargin)

% addpath('.\FieldII');
% addpath('.\ompbox');

nLines = 80;
fs = 400e6;
freq = 5e6;
BW = 0.5;
taps = 50;

% load simulatedData
% load D:\simulatedData40
% artIQ = artifactIQdata;
% clear artifactIQdata;

% default values for the parameters

if nargin <4
    rho = 0.98;
    SNR = 30;
elseif nargin <5
    rho = varargin{1};
    SNR = 30;
else
    rho = varargin{1};
    SNR = varargin{2};
end

if nargin > 6
    IQdata = varargin{3};
    artIQ = varargin{4};
else
    error('data missing');
%     load D:\simulatedData40
%     artIQ = artifactIQdata;
%     clear artifactIQdata;    
end

% caculate the noise STD
kern = IQdata(4000:10000,1:20,1);
v = kern(:);
NSTD = (std(real(v) * 1e24)^20 / (10^SNR))^(1/20); % Noise standard deviation
NSTD = NSTD * 1e-24; 


decorrData = decorrelateData(IQdata, rho);
[lineLength, nElements, mFrames] = size(decorrData);

decorrData = decorrData(:,:,1:19);

% Generate the displaced ensembles with the clutter artifact
lineSamples = 1090;
[RFdata, RFdataClean] = displacedData_uponly(decorrData, artIQ, shift, art_shift, downsampleStep, lineSamples);
RFdata = RFdata(1:lineSamples,:,:);
RFdataClean = RFdataClean(1:lineSamples,:,:);

% Add white gaussian noise to the signal
Noise = NSTD * randn(size(RFdata)) + 1i * NSTD * randn(size(RFdata));
FinalData = RFdata + Noise;
% [R, S, T] = size(FinalData);

% figure(1); 
% imagesc(abs(FinalData(:,:,1)));
% colormap gray
% % axis image;
% drawnow();
% for i = 1:T
%     imagesc(abs(FinalData(:,:,i)));
%     colormap gray
% %     axis image;
%     drawnow();
% end

end