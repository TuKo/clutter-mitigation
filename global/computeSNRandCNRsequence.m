function [CNR, CNR2, varargout] = computeSNRandCNRsequence(detData, zoneInX, zoneInY, zoneOutX, zoneOutY)
nFrames = size(detData,3);

CNR = zeros(nFrames,1);
CNR2 = zeros(nFrames,1);
SNRout = zeros(nFrames, numel(zoneOutX(:)));

for frame = 1:nFrames
    imgIn = detData(sub2ind(size(detData), zoneInY, zoneInX, frame*ones(numel(zoneInX),1)));
    imgOut = detData(sub2ind(size(detData), zoneOutY, zoneOutX, frame*ones(numel(zoneOutX),1)));

%     SNRIn(frame) = 20*log10(mean(imgIn(:))./std(imgIn(:)));
%     SNROut(frame) = 20*log10(mean(imgOut(:))./std(imgOut(:)));
    SNRout(frame,:) = imgOut(:);
    CNR(frame) = 20*log10(abs(mean(imgOut(:))-mean(imgIn(:)))./std(imgOut(:)));
    CNR2(frame) = 20*log10(mean(imgOut(:))/mean(imgIn(:)));
%     CNR2(frame) = 20*log10( abs(mean(imgOut(:))-mean(imgIn(:)))./sqrt(std(imgOut(:))^2+std(imgIn(:))^2) );

% CNR = log10(mean(imgOut(:))/mean(imgIn(:)))/log10(std(imgOut(:)));

% width = lesion_width(IQdata(:,40),fs);
end

CNR = mean(CNR);
CNR2 = mean(CNR2);
SNR = mean(SNRout(:))/std(SNRout(:));
if max(nargout,1)-1 >=1
    varargout{1} = SNR;
end

end