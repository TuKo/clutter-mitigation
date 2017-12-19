function [x1, y1] = correlateWindow(reference, target, x0, y0, dx0, dy0, dWindowX, dWindowY)

% rho = zeros(dWindowY+1,dWindowX+1);
% 
% orig = reference(y0:dy0, x0:dx0);
% orig = orig - mean(orig(:));
% origDivisor = sum(sum(orig.^2));
% 
% dwx = floor((x0+dx0)/2- dWindowX/2);
% dwy = floor((y0+dy0)/2- dWindowY/2);
% 
% wx0 = dx0-x0+1;
% wy0 = dy0-y0+1;
% 
% for alpha = 1:dWindowY+1
%     for beta = 1:dWindowX+1
%         tgt = target(alpha + dwy + (1:wy0), beta + dwx + (1:wx0));
%         tgt = tgt - mean(tgt(:));
%         tgtDivisor = sum(sum(tgt.^2));
%         rho(alpha,beta) = sum(sum(orig.*tgt)) / sqrt(origDivisor * tgtDivisor);
%     end
% end
% 
% [y1,x1] = find(rho == max(rho(:)),1);
% y1 = y1 + dwy;
% x1 = x1 + dwx;

[maxY, maxX] = size(target);
rho = -ones(size(target));

orig = reference(y0:dy0, x0:dx0);
orig = orig - mean(orig(:));
origDivisor = sum(sum(orig.^2));

wx0 = dx0-x0+1;
wy0 = dy0-y0+1;

maxY = maxY - wy0 - dWindowY-1;
maxX = maxX - wx0 -dWindowX-1;

u0 = max(1,y0-dWindowY);
v0 = max(1,x0-dWindowX);

for alpha = u0:min(dy0+dWindowY, maxY)
    for beta = v0:min(dx0+dWindowX, maxX)
        tgt = target(alpha + (0:(wy0-1)), beta + (0:(wx0-1)));
        tgt = tgt - mean(tgt(:));
        tgtDivisor = sum(sum(tgt.^2));
        rho(alpha,beta) = sum(sum(orig.*tgt)) / sqrt(origDivisor * tgtDivisor);
    end
end

[y1,x1] = find(rho == max(rho(:)),1);