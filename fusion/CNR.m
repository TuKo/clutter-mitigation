function [CNR, stdCNR] = CNR(detData, zoneInX, zoneInY, zoneOutX, zoneOutY)

nn = numel(zoneInX);
mm = numel(zoneOutX);
CNR = zeros(size(detData,3),1);
for i = 1:size(detData,3)
    imgIn = detData(sub2ind(size(detData), zoneInY(:), zoneInX(:), repmat(i,[nn,1])));
    imgOut = detData(sub2ind(size(detData), zoneOutY(:), zoneOutX(:), repmat(i,[mm,1])));
    CNR(i) = 20*log10(abs(mean(imgOut(:))-mean(imgIn(:)))./std(imgOut(:)));
end

stdCNR = std(CNR);
CNR = mean(CNR);

end