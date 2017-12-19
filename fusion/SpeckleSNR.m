function [SNR, stdSNR, stdspace,stdtime] = SpeckleSNR(seq, RegX, RegY)

[RY, RX] = meshgrid(RegY(1):RegY(2), RegX(1):RegX(2));
RMN = size(RX(:),1);
nFrames = size(seq,3);

% seqreg = seq(RY,RX,:);
if nFrames == 1
    seqreg = seq(sub2ind(size(seq), RY(:),RX(:)));
else
    v = zeros(RMN,1);
    v(1) = 1;
    v = cumsum(repmat(v, [nFrames, 1]));
    seqreg = seq(sub2ind(size(seq), repmat(RY(:),[nFrames 1]), repmat(RX(:), [nFrames 1]), v));
    seqreg = reshape(seqreg,[RMN,nFrames]);
end

SNR = mean(mean(seqreg)./std(seqreg));
stdSNR = std(mean(seqreg)./std(seqreg));

stdspace = mean(std(seqreg));
stdtime = mean(std(seqreg,0,2));
% sss = mean(seqreg(:));
end