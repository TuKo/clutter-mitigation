function [mov1,mov2,tissue1, tissue2] = FusionTA(D1C, D2C, Dt, At, seq1, seq2, params, startFrame, stopFrame, Mask)


if (isfield(params,'M'))
    M = params.M;
else
    M = 24; % default axial
end
if (isfield(params,'N'))
    N = params.N;
else
    N = min(8,size(seq,3)); %default A-lines
end
if (isfield(params,'P'))
    P = params.P;
else
    P = 0; % default radial
end
if (isfield(params,'OMP_ERROR'))
    OMP_ERROR = params.OMP_ERROR;
    ERR_THR_CRITERION = 1;    
else
    OMP_ERROR = 1;
    ERR_THR_CRITERION = 0;
end
if (isfield(params,'CUTTOFFs'))
    THRs = params.CUTTOFFs;
else
    error('No threshold parameter is defined for choosing the clutter atoms (CUTTOFFs)');
end

if (isfield(params,'OL_M'))
    OL_AX = params.OL_M;
else
    OL_AX = round(M/2); % Maximum axial overlap by default (from 0 to M)
end
if (isfield(params,'OL_P'))
    OL_R = params.OL_P;
else
    OL_R = 0; % overlap in angle axis (from 0 to P/2)
end
if (isfield(params,'OMP_MAX_ATOMS'))
    OMP_MAX_ATOMS = params.OMP_MAX_ATOMS;
else
    OMP_MAX_ATOMS = 40; % Maximum number of atoms when using OMPerr.
end

[H,W,totFrames] = size(seq1);

nFrames = stopFrame-startFrame+1;
totFrames = nFrames;

% mov1 = cell(numel(THRs),1);
% mov2 = cell(numel(THRs),1);
% parfor t=1:numel(THRs)
%     mov1{t} = zeros(H,W,nFrames);
%     mov2{t} = zeros(H,W,nFrames);
% end

mov1 = zeros(H,W,totFrames);
mov2 = zeros(H,W,totFrames);
cnt1 = zeros(H,W,totFrames);
cnt2 = zeros(H,W,totFrames);
tissue1 = zeros(H,W,totFrames);
tissue2 = zeros(H,W,totFrames);
[n1, m1] = size(D1C);
[n2, m2] = size(D2C);

t=1;
% parfor t = 1:numel(THRs)
    
    % Then clean the clutter atoms from the adaptive dictionary to obtain the final Dt
    Dt2 = Dt(:,~(sum(abs(Dt - D2C*At).^2)<THRs(t)));
    
    % Then merge Dt2 and D1C and D2C
    D = [Dt2./sqrt(2), D1C, zeros(n1,m2); Dt2./sqrt(2), zeros(n2,m1), D2C];
    DtD = D'*D;
    stationaryAtom = false(size(D,2),1);
    stationaryAtom((size(Dt2,2)+1):end) = true;    
%     stationaryAtom1 = false(size(D,2),1);
%     stationaryAtom1(size(Dt2,2)+(1:m1)) = true;
%     stationaryAtom2 = false(size(D,2),1);
%     stationaryAtom2(size(Dt2,2)+m1+(1:m2)) = true;
    
    parfor frame = 1:nFrames
        midFrame = startFrame + frame -1;
        tic
        fprintf('Frame = %d ', frame);

        posMidFrame = 1+N/2;
        patches1 = extractPatchesMasked(seq1, midFrame,  N, P, M, OL_R, OL_AX, Mask);
        patches2 = extractPatchesMasked(seq2, midFrame,  N, P, M, OL_R, OL_AX, Mask);
        patches = [patches1; patches2];
        
        if ERR_THR_CRITERION
            G2 = comp2(D'*patches, sum(conj(patches).*patches), DtD, OMP_ERROR, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)
        else
            G2 = comp(D'*patches, DtD, OMP_MAX_ATOMS); % maximum of atoms (last param)
        end
        fprintf('Avg atoms used: %g', mean(full(sum(G2~=0))));
        
%         G = G2;
%         G(stationaryAtom,:) = 0;
%         patchesClean = D*G;
%         [tN1, cntN1] =  joinPatchesMiddleMasked(patchesClean(1:size(patches1,1),:), H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
%         tN1(cntN1 ~=0) = tN1(cntN1 ~=0) ./ cntN1(cntN1 ~=0);
%         tissue1(:,:,frame) = tN1;
%         clear tN1 cntN1

%         [tN2, cntN2] =  joinPatchesMiddleMasked(patchesClean((1+size(patches1,1)):end,:), H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
%         tN2(cntN2 ~=0) = tN2(cntN2 ~=0) ./ cntN2(cntN2 ~=0);
%         tissue2(:,:,frame) = tN2;
%         clear tN2 cntN2
        
        % Extract clutter only with every value of Threshold
        G2(~stationaryAtom,:) = 0;
        patchesClean = patches - D*G2;
        % Average the samples in all patches
        %         [movN, cntN] = joinPatches(patches2, H,W, N,P,M, OL_R,OL_AX);
        %         mov(:,:,t) = mov(:,:,t) + movN(:,:,midFrame-(startFrame-N/2)+1);
        %         cnt(:,:,t) = cnt(:,:,t) + cntN(:,:,midFrame-(startFrame-N/2)+1);
        
%         [movN1, cntN1] =  joinPatches(patchesClean(1:size(patches1,1),:), H,W, N,P,M, OL_R,OL_AX);
%         mov1(:,:,midFrame+(-N/2:N/2)) = mov1(:,:,midFrame+(-N/2:N/2)) + movN1;
%         cnt1(:,:,midFrame+(-N/2:N/2)) = cnt1(:,:,midFrame+(-N/2:N/2)) + cntN1;
%         
%         [movN2, cntN2] =  joinPatches(patchesClean((1+size(patches1,1)):end,:), H,W, N,P,M, OL_R,OL_AX);
%         mov2(:,:,midFrame+(-N/2:N/2)) = mov2(:,:,midFrame+(-N/2:N/2)) + movN2;
%         cnt2(:,:,midFrame+(-N/2:N/2)) = cnt2(:,:,midFrame+(-N/2:N/2)) + cntN2;
        
        [movN1, cntN1] =  joinPatchesMiddleMasked(patchesClean(1:size(patches1,1),:), H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
        movN1(cntN1 ~=0) = movN1(cntN1 ~=0) ./ cntN1(cntN1 ~=0);
        mov1(:,:,frame) = movN1;
%         clear movN1 cntN1
        
        [movN2, cntN2] =  joinPatchesMiddleMasked(patchesClean((1+size(patches1,1)):end,:), H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
        movN2(cntN2 ~=0) = movN2(cntN2 ~=0) ./ cntN2(cntN2 ~=0);
        mov2(:,:,frame) = movN2;
%         clear movN2 cntN2
        
        fprintf(', time = %g\n', toc);
    end

% end

% mov1(cnt1 ~=0) = mov1(cnt1 ~=0) ./ cnt1(cnt1 ~=0);
% mov2(cnt2 ~=0) = mov2(cnt2 ~=0) ./ cnt2(cnt2 ~=0);

end