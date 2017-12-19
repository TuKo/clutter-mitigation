function writeVideo(filename, movie, frameRate)

obj = VideoWriter(filename,'Motion JPEG AVI');
% obj = VideoWriter(filename,'Uncompressed AVI');

obj.FrameRate = frameRate;

open(obj);


for k = 1:size(movie,2)
    writeVideo(obj,movie(k));
end

close(obj);

end
