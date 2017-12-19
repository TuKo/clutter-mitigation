function frames = convertMovie(images)

[videoH, videoW, nFrames] = size(images);

frames = struct('cdata', zeros(videoH, videoW, 3, 'uint8'), 'colormap', []);

for k = 1 : nFrames
    frames(k) = im2frame(repmat(im2uint8(images(:,:,k)),[1,1,3]));
end

end