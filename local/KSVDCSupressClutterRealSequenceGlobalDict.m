function mov = KSVDCSupressClutterRealSequenceGlobalDict(Dt, Dc, seq, params, startFrame, stopFrame, Mask)


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
    ERR_THR_CRITERION = true;    
else
    ERR_THR_CRITERION = false;
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

[H,W,~] = size(seq);

D = [Dc , Dt];
DtD = D'*D;


nFrames = stopFrame-startFrame+1;

stationaryAtom = false(size(D,2),1);
stationaryAtom(1:size(Dc,2)) = true;

mov = zeros(H,W,nFrames);

for frame = 1:nFrames
    % lastFrame = startFrame + N;
    midFrame = startFrame + frame -1;
    
    tic
    
    fprintf('Frame = %d ', frame);

    posMidFrame = 1+N/2;
    patches = extractPatchesMasked(seq, midFrame,  N, P, M, OL_R, OL_AX, Mask);

    if ERR_THR_CRITERION    
        G = comp2(D'*patches, sum(conj(patches).*patches), DtD, OMP_ERROR, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)
    else
        G = comp(D'*patches, DtD, OMP_MAX_ATOMS); % maximum of atoms (last param)
    end
    fprintf('Avg atoms used: %g', mean(full(sum(G~=0))));
    %     G = comp2(D, patches, OMP_ERROR, OMP_MAX_ATOMS); % maximum of atoms (last param)
    %     G = comp2(D, patches*1e25, DtD, OMP_ERROR*1e25, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)

    % Extract clutter only with every value of Threshold
    %     mov = zeros(H,W,numel(THRs));
    %     cnt = zeros(H,W,numel(THRs));
    %     wrong_patches = 0;
    %     patches2 = zeros(size(patches));
    G2 = G;
    G2(~stationaryAtom,:) = 0;
    patches2 = patches - D*G2;

    % Average the samples in all patches
    %         [movN, cntN] = joinPatches(patches2, H,W, N,P,M, OL_R,OL_AX);
    %         mov(:,:,t) = mov(:,:,t) + movN(:,:,midFrame-(startFrame-N/2)+1);
    %         cnt(:,:,t) = cnt(:,:,t) + cntN(:,:,midFrame-(startFrame-N/2)+1);
    [movN, cntN] =  joinPatchesMiddleMasked(patches2, H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
    movN(cntN ~=0) = movN(cntN ~=0) ./ cntN(cntN ~=0);
    mov(:,:,frame) = movN;
    %     fprintf(', wrong energy = %d, time = %g\n',  wrong_patches, toc);
    fprintf(', time = %g\n', toc);
end


end