function [PSNR] = computePSNRsequence(detData, detDataOrig, zoneInX, zoneInY)

MAX = 10;

nFrames = size(detData,3);
PSNR = zeros(nFrames,1);

for frame = 1:nFrames
    imgDest = detData(sub2ind(size(detData), zoneInY, zoneInX, frame*ones(numel(zoneInX),1)));
    imgOrig = detDataOrig(sub2ind(size(detData), zoneInY, zoneInX, frame*ones(numel(zoneInX),1)));
    MSE = sum((imgOrig(:)-imgDest(:)).^2)/numel(imgOrig);
    PSNR = 10*log10(MAX^2/MSE);
end

PSNR = mean(PSNR);

end