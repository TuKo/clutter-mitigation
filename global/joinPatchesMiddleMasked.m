function [mov, cnt] = joinPatchesMiddleMasked(patches, H,W, N,P,M, OL_P,OL_M,mid,Mask)


% mov = zeros(H,W,N+1);
% cnt = zeros(H,W,N+1);
% 
% cc = 1;
% for k = (1+M/2):M+1-OL_M:(H-M/2)
%     for j = (1+P/2):P+1-OL_P:(W-P/2)
%         x = zeros(M+1,P+1,N+1);
%         x(:) = patches(:,cc);
%         mov(k+(-M/2:M/2),j+(-P/2:P/2),1:N+1) = mov(k+(-M/2:M/2),j+(-P/2:P/2),1:N+1) + x;
%         cnt(k+(-M/2:M/2),j+(-P/2:P/2),1:N+1) = cnt(k+(-M/2:M/2),j+(-P/2:P/2),1:N+1) + 1;
%         cc = cc + 1;
%     end
% end

if (P > 0)
%     mov = zeros(H,W);
%     cnt = zeros(H,W);
% 
%     x = zeros(M+1,P+1,N+1);
%     cc = 1;
%     for k = (1+M/2):M+1-OL_M:(H-M/2)
%         for j = (1+P/2):P+1-OL_P:(W-P/2)
%             x(:) = patches(:,cc);
%             mov(k+(-M/2:M/2),j+(-P/2:P/2)) = mov(k+(-M/2:M/2),j+(-P/2:P/2)) + x(:,:,mid);
%             cnt(k+(-M/2:M/2),j+(-P/2:P/2)) = cnt(k+(-M/2:M/2),j+(-P/2:P/2)) + 1;
%             cc = cc + 1;
%         end
%     end
    error('not implemented');
else
    mov = zeros(H,W);
    cnt = zeros(H,W);

    x = zeros(M+1,N+1);
    cc = 1;
    for k = (1+M/2):M+1-OL_M:(H-M/2)
        for j = 1:W
            if sum(sum(Mask(k+(-M/2:M/2),j+(-P/2:P/2))))<1
                continue;
            end
            x(:) = patches(:,cc);
            mov(k+(-M/2:M/2),j) = mov(k+(-M/2:M/2),j) + x(:,mid);
            cnt(k+(-M/2:M/2),j) = cnt(k+(-M/2:M/2),j) + 1;
            cc = cc + 1;
        end
    end
end

% cc = 1;
% for k = (1+M/2):M+1-OL_M:(H-M/2)
%     for j = (1+P/2):P+1-OL_P:(W-P/2)
%         x = zeros(M+1,P+1,N+1);
%         x(:) = patches(:,cc);
%         mov(k,j) = mov(k,j) + x(M/2+1,1,mid);
%         cnt(k,j) = cnt(k,j) + 1;
%         cc = cc + 1;
%     end
% end

end