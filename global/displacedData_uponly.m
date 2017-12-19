function [RFdata, RFdataClean] = displacedData_uponly(dataIQ, artifactIQ, shift, art_shift, downsampleStep,lineSamples)

[lineLength, nElements, mFrames] = size(dataIQ);

% lineSamples = floor(lineLength/downsampleStep);%-max(downsampleStep*mFrames,downsampleStep*mFrames);
RFdata = zeros(lineSamples,nElements,mFrames);
RFdataClean = zeros(lineSamples,nElements,mFrames);


T = (mFrames-1)/2;
extraShift = T*shift/downsampleStep;
for k = 1:mFrames
    tempData = dataIQ((1+shift*(k-1)):downsampleStep:end,:,k);
    auxData = zeros(lineSamples,nElements);
    auxData((1+extraShift):lineSamples,:) = tempData(1:(lineSamples-extraShift),:);
%     auxData = dataIQ((1+shift*(k-1)):downsampleStep:end,:,k);
    auxArtifact = artifactIQ((1+art_shift*(k-1)):downsampleStep:end,:);
    RFdata(1:lineSamples,:,k) = auxData(1:lineSamples,:) + auxArtifact(1:lineSamples,:);
    RFdataClean(1:lineSamples,:,k) = auxData(1:lineSamples,:);
end

end