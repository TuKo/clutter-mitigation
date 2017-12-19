function [mov, cnt] = joinPatches(patches, H,W, N,P,M, OL_P,OL_M)

mov = zeros(H,W,N+1);
cnt = zeros(H,W,N+1);

cc = 1;
x = zeros(M+1,1,N+1);

for k = (1+M/2):M+1-OL_M:(H-M/2)
    for j = 1:W
        x(:) = patches(:,cc);
        mov(k+(-M/2:M/2),j,1:N+1) = mov(k+(-M/2:M/2),j,1:N+1) + x;
        cnt(k+(-M/2:M/2),j,1:N+1) = cnt(k+(-M/2:M/2),j,1:N+1) + 1;
        cc = cc + 1;
    end
end


end