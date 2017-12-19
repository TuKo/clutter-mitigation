addpath('helpers/');
addpath('ompbox/');

j=1;
datasets = {'datasets_new/p2s3/'};

[seq, nMLAs] = readDataset(datasets{j});

DR = 30;
GAIN = -40;

seq = tissueProcessing(stb(seq,nMLAs,nMLAs/2),DR,GAIN);
seq_out = seq;

% 31, 287 --> 51, 320  (size = 10 by 34)
% x0 = 31; y0 = 287;
% dX0 = 51; dY0 = 320;
% 
% x0 = 25; y0 = 257;
% dX0 = 55; dY0 = 320;
% wx = dX0 - x0;
% wy = dY0 - y0;

wx = 30;
wy = 30;

half_wx = floor(wx/2);
half_wy = floor(wy/2);

dWindowX = 20; % +/- pixels in X
dWindowY = 20; % +/- pixels in Y

dirX = zeros(size(seq_out));
dirY = zeros(size(seq_out));

[maxY,maxX]  = size(seq_out(:,:,1));
maxX = maxX - wx - 1; 
maxY = maxY/2 - wy - 1;

for i = 1:size(seq,3)-1
    fprintf('Frame %g\n',i);
    
    % x0 = 41;
    % y0 = (320+287)/2;
    % dX0 = 10;
    % dY0 = 17;
for x0 = 1:wx:maxX
    for y0 = 1:wy:maxY
    [x1,y1] = correlateWindow(seq(:,:,i), seq(:,:,i+1), x0, y0, x0+wx, y0+wy, dWindowX, dWindowY);
    
    dirX(y0+half_wy,x0+half_wx,i) = x1-x0;
    dirY(y0+half_wy,x0+half_wx,i) = y1-y0;
    end
end

%     x0 = x1;
%     y0 = y1;
%     dX0 = x1 + wx;
%     dY0 = y1 + wy;
%     seq_out(y0,x0:dX0,i+1) = 1;
%     seq_out(dY0,x0:dX0,i+1) = 1;
%     seq_out(y0:dY0,x0,i+1) = 1;
%     seq_out(y0:dY0,dX0,i+1) = 1;
    
    
%     x1
%     y1
end


% implay(seq)
% implay(seq_out)

return
%%
for i= 1:size(seq,3)-1
    figure(1);
    quiver(dirX(:,:,i),dirY(:,:,i));
    drawnow(); 
    pause();
end

