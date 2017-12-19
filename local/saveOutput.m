function [] = saveOutput(seq, filename, visualFrames)

height = 400;
width = 400;

% for reinput
% writeDataset([filename '\'],seq,1);

% video
writeVideo([filename '.avi'], convertMovie(scanConvertMovie(seq,height,width)), 20);

% examples
for j = 1:numel(visualFrames)
    imwrite(scanConvert(seq(:,:,visualFrames(j)),height,width),[filename '-Frm' num2str(j) '.jpg']);
end

end