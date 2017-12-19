function [SNRIn, SNROut, CNR, CNR2, imgIn, imgOut] = computeSNRandCNR(detData, zoneInX, zoneInY, zoneOutX, zoneOutY)

imgIn = detData(sub2ind(size(detData), zoneInY, zoneInX));
imgOut = detData(sub2ind(size(detData), zoneOutY, zoneOutX));

SNRIn = 20*log10(mean(imgIn(:))./std(imgIn(:)));
SNROut = 20*log10(mean(imgOut(:))./std(imgOut(:)));

CNR = 20*log10(abs(mean(imgOut(:))-mean(imgIn(:)))./std(imgOut(:)));

CNR2 = 20*log10(mean(imgOut(:))/mean(imgIn(:)));
% CNR = log10(mean(imgOut(:))/mean(imgIn(:)))/log10(std(imgOut(:)));

% width = lesion_width(IQdata(:,40),fs);

end