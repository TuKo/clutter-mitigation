function [I, nMLAs] = readDataset(folder)

if ~exist(folder,'file')
    error('Unknown folder');
end

%read info file
[results]=inifile([folder 'info.txt'],'read',{'dump','','numcycles';'dump','','numvectors';'dump','','numsamples';'dump','','nummlas'});
nFiles = str2num(results{1});
nVectors = str2num(results{2});
nSamples = str2num(results{3});
nMLAs = str2num(results{4});

if isempty(nMLAs)
    nMLAs = 1;
end
    
% I = zeros(floor(nSamples/nMLAs),nVectors*nMLAs,nFiles);
I = zeros(floor(nSamples),nVectors*nMLAs,nFiles);
for j = 1:nFiles
    file = [folder num2str(j-1) '.bin'];
    fid = fopen(file);
%     Comp = fread(fid, [2*nSamples/nMLAs, nVectors*nMLAs], 'float');
    Comp = fread(fid, [2*nSamples, nVectors*nMLAs], 'float');
    I(:,:,j) = Comp(1:2:end,:) + 1i*Comp(2:2:end,:);
    fclose(fid);
end


end