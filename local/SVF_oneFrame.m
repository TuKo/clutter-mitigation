function mov = SVF_oneFrame(seq, M, N, alpha, tau)

[H,W,Frames] = size(seq);

mov = zeros(size(seq));

% parfor l = (1+N/2):(Frames-N/2)
l = (Frames-1)/2+1;
zz = zeros(H,W);
%     seq_aux = seq(:,:,l+ ((-N/2):(N/2)));
for k =(1+M/2):(H-M/2)
    for j =1:W
%             x = reshape(seq_aux(k+(-M/2:M/2),j,:), [M+1, N+1]);
        x = reshape(seq(k+(-M/2:M/2),j,l+ ((-N/2):(N/2))), [M+1, N+1]);
        [u,s,v] = svd(x,0);
        w = 1 - 1./(1+exp(-alpha*(diag(s)./sum(diag(s))-tau)));
        x = u*diag(w)*s*v';
        zz(k,j) = x(M/2+1,N/2+1);
    end
end
mov(:,:,l) = zz;

end