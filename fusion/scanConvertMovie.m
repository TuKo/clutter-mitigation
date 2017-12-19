function finalMov = scanConvertMovie(mov, imgH, imgW)

% finalMov = zeros(imgH,imgW,size(mov,3));
% parfor k=1:size(mov,3)
%     finalMov(:,:,k) = scanconvert(mov(:,:,k),imgH,imgW);
% end

finalMov = scanConvert(mov,imgH,imgW);

end