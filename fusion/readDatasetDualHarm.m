function [Ifst, Isnd, nMLAs] = readDatasetDualHarm(folder)
% [FST, SND, MLAS] = READDATASETDUALHARM(FOLDER)
%   FOLDER - a string with a relative or absolute path with the binary
%   files of a IQ sequence.
%
%   Output:
%       FST  - The sequence for the fundamental harmonic
%       SND  - The sequence for the second harmonic
%       MLAS - Number of lines per acquired radius


% Copyright Javier Turek (c)
% Technion, Israel 3200003
% All rights reserved

if ~exist(folder,'file')
    error('Unknown folder');
end

%read info file
[results]=inifile([folder 'info.txt'],'read',{'dump','','numcycles';'dump','','numvectors';'dump','','numsamples';'dump','','nummlas'});
nFiles = str2num(results{1});
nVectors = str2num(results{2});
nSamples = str2num(results{3});
nMLAs = str2num(results{4});

%TODO: need to read from the ini file
SR = 4.0404e6;
f1st = 1.7e6;
f2nd = 3.4e6;
%fMid = 2.55e6; %or? 
fMid = 2.7e6;
BW1st = 0.7e6;
BW2nd = 1e6;

if isempty(nMLAs)
    nMLAs = 1;
end
    
Ifst = zeros(nSamples,nVectors*nMLAs,nFiles);
Isnd = zeros(nSamples,nVectors*nMLAs,nFiles);
for j = 1:nFiles
    file = [folder num2str(j-1) '.bin'];
    fid = fopen(file);
    Comp = fread(fid, [2*nSamples, nVectors*nMLAs], 'float');
    [frame1st, frame2nd] = getHarmonics(complex(Comp(1:2:end,:),Comp(2:2:end,:)), SR, BW1st, BW2nd, f1st, f2nd, fMid);
    Ifst(:,:,j) = frame1st;
    Isnd(:,:,j) = frame2nd;
    fclose(fid);
end

end

function [fst, snd] = getHarmonics(data, samplingRate, bandWidth1st, bandWidth2nd, freq1st, freq2nd, freqMid)

fst = zeros(size(data));
snd = zeros(size(fst));

%SR=4.0404e6;
%freqMid = 2.55e6;
%bandWidth = 0.7e6
%freq1st = 1e6
%freq2nd = 0.7e6

t = (0:length(data(:,1))-1)/samplingRate;
b1 = MyLPF(160, bandWidth1st/(samplingRate/2));
b2 = MyLPF(160, bandWidth2nd/(samplingRate/2));

for i = 1:size(data,2)
    %shift freq to -(mid-first) MHz
    IQ = data(:,i);
    % close all
    % DisplaySpectrum(IQ, samplingRate, 0, 'fist vector', '');
    % DisplaySpectrum(IQ, samplingRate, freqMid, 'fist vector', '');
    IQ1st = IQ.*exp(complex(0,(freqMid-freq1st)*2*pi*t'));
    IQ1stfilt = filter(b1,1,IQ1st);
    IQ2nd = IQ.*exp(complex(0,-(freq2nd-freqMid)*2*pi*t'));
    IQ2ndfilt = filter(b2,1,IQ2nd);
    
    fst(:,i) = IQ1stfilt;
    snd(:,i) = IQ2ndfilt;

%     if 0 
%         DisplaySpectrum(IQ1st, samplingRate, 0, 'IQ1st', '');
%         DisplaySpectrum(IQ1stfilt, samplingRate, 0, 'IQ1srfilt', '');
%         DisplaySpectrum(IQ2nd, samplingRate, 0, 'IQ2nd', '');
%         DisplaySpectrum(IQ2ndfilt, samplingRate, 0, 'IQ2ndfilt', '');
%     end
end

end