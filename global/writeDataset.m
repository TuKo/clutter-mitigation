function writeDataset(folder,I,nMLAs)

if ~exist(folder,'file')
    [success, ~,~ ] = mkdir(folder);
    if ~success
        error('Cannot create folder');
    end
end

[nSamples,nVectors,nFiles] = size(I);

%write info file
inifile([folder 'info.txt'],'write',{'dump','','numcycles', nFiles;'dump','','numvectors',nVectors/nMLAs;'dump','','numsamples',nSamples*nMLAs;'dump','','nummlas',nMLAs});

Comp = zeros(2*nSamples*nMLAs,nVectors/nMLAs);
for j = 1:nFiles
    file = [folder num2str(j-1) '.bin'];
    fid = fopen(file,'w');
    Comp(1:2:end,:) = reshape(real(I(:,:,j)),[nSamples*nMLAs, nVectors/nMLAs]);
    Comp(2:2:end,:) = reshape(imag(I(:,:,j)),[nSamples*nMLAs, nVectors/nMLAs]);
    fwrite(fid, Comp, 'float');
    fclose(fid);
end


end