function patches = extractPatches(seq, offsetFrame,  N,P,M, OL_P,OL_M)

[H,W,~] = size(seq);

if (mod(N,2) == 1 || mod(P,2) == 1 || mod(M,2) == 1 || N < 0 || M < 0 || P < 0)
    error ('Values of M, N, and P should be non-negative and even.');
end

patches = zeros((M+1)*(N+1)*(P+1),round((W-P)*(H-M)/(M+1-OL_M)/(P+1-OL_P)));

cc = 1;
for k = (1+M/2):M+1-OL_M:(H-M/2)
    for j = (1+P/2):P+1-OL_P:(W-P/2)
%         if k-M/2 <= 419 && k+M/2 >=419  && j == 31
%             fprintf('x=%d, y=%d, cc=%d\n',k,j,cc);
%         end
%         x = reshape(seq(k+(-M/2:M/2),j+(-P/2:P/2),offsetFrame+(-N/2:N/2)), [(M+1)*(P+1), (N+1)]);
        x = seq(k+(-M/2:M/2),j+(-P/2:P/2),offsetFrame+(-N/2:N/2));
        patches(:,cc) = x(:);
        cc = cc + 1;
    end
end

%old code
%     patches = zeros((M+1)*(N+1),round(W*(H-M)/(M+1-OL_AX)));
%     patches = zeros((M+1)*(N+1),round(W*(H-M)));
%     cc = 1;
%     for j = 1:W
%         for k = (1):M+1-OL_AX:(H-M)
% %         for k = (1+M/2):(H-M/2)
%             x = reshape(seq(k+(0:M),j,startFrame+(-N/2:N/2)), [M+1, N+1]);
% %             x = reshape(seq(k+(-M/2:M/2),j,startFrame+(-N/2:N/2)), [M+1, N+1]);
%             patches(:,cc) = x(:);
%             cc = cc + 1;
%         end
%     end
%     patches = patches(:,1:cc-1);

end