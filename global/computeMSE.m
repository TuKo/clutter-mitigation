function [PSNR, MSE] = computeMSE(detData, detDataOrig, zoneInX, zoneInY)

imgDest = detData(sub2ind(size(detData), zoneInY, zoneInX))*1e24;

imgOrig = detDataOrig(sub2ind(size(detData), zoneInY, zoneInX))*1e24;

MSE = sum((imgOrig(:)-imgDest(:)).^2)/numel(imgOrig);

MAX = 10;
PSNR = 10*log10(MAX^2/MSE);

MSE = sqrt(MSE*(1e-24)^2);

end