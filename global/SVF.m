function mov = SVF(seq, M, N, alphas, tau, startFrame, stopFrame, Mask)

[H,W,~] = size(seq);

if numel(alphas) >1
    mov = cell(numel(alphas),1);

    parfor l =1:numel(alphas)
        mov{l} = zeros(H,W,stopFrame-startFrame+1);
    end

    rows = sum(Mask,2)>0;
    for l = startFrame:stopFrame
    %     zz = zeros(H,W);
    %     seq_aux = seq(:,:,l+ ((-N/2):(N/2)));
        for k =(1+M/2):(H-M/2)
            if ~rows(k)
                continue;
            end
            for j =1:W
                if ~Mask(k,j)
                    continue;
                end
    %             x = reshape(seq_aux(k+(-M/2:M/2),j,:), [M+1, N+1]);
                x = reshape(seq(k+(-M/2:M/2),j,l+ ((-N/2):(N/2))), [M+1, N+1]);
                [u,s,v] = svd(x,0);
                for m = 1:numel(alphas)
                    w = 1 - 1./(1+exp(-alphas(m)*(diag(s)./sum(diag(s))-tau)));
                    x = u*diag(w)*s*v';
                    mov{m}(k,j,l-startFrame+1) = x(M/2+1,N/2+1);
                end
    %             zz(k,j) = x(M/2+1,N/2+1);
            end
        end
    %     mov(:,:,l-startFrame+1) = zz;
    end
else
    mov = zeros(H,W,stopFrame-startFrame+1);
    alpha = alphas(1);
    
    rows = sum(Mask,2)>0;
    for l = startFrame:stopFrame
        zz = zeros(H,W);
    %     seq_aux = seq(:,:,l+ ((-N/2):(N/2)));
        parfor k =(1+M/2):(H-M/2)
            if ~rows(k)
                continue;
            end
            for j =1:W
                if ~Mask(k,j)
                    continue;
                end
    %             x = reshape(seq_aux(k+(-M/2:M/2),j,:), [M+1, N+1]);
                x = reshape(seq(k+(-M/2:M/2),j,l+ ((-N/2):(N/2))), [M+1, N+1]);
                [u,s,v] = svd(x,0);
                w = 1 - 1./(1+exp(-alpha*(diag(s)./sum(diag(s))-tau)));
                x = u*diag(w)*s*v';
%                 mov(k,j,l-startFrame+1) = x(M/2+1,N/2+1);
                zz(k,j) = x(M/2+1,N/2+1);
            end
        end
        mov(:,:,l-startFrame+1) = zz;
    end
    
end

end