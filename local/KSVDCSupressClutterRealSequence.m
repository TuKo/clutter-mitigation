function mov = KSVDCSupressClutterRealSequence(D, DtD, seq, params, startFrame, stopFrame, Mask)


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
if (isfield(params,'ATOM_THRs'))
    THRs = params.ATOM_THRs;
else
    error('No threshold parameter is defined for choosing the clutter atoms');
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
if (isfield(params,'CLUTTER_METHOD'))
    CLUTTER_METHOD = params.CLUTTER_METHOD;
else
    CLUTTER_METHOD = 1; % by Default PCA method
end

[H,W,~] = size(seq);

nFrames = stopFrame-startFrame+1;

stationaryAtom = false(size(D,2),numel(THRs));
for t = 1:numel(THRs)
    stationaryAtom(:,t) = findClutterAtoms(D, THRs(t), CLUTTER_METHOD, M+1, N+1, P+1);
end

mov = cell(numel(THRs),1);
for t=1:numel(THRs)
    mov{t} = zeros(H,W,nFrames);
end

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
%         G = comp2(D, patches, OMP_ERROR, OMP_MAX_ATOMS); % maximum of atoms (last param)
    end
    fprintf('Avg atoms used: %g', mean(full(sum(G~=0))));

    % Extract clutter only with every value of Threshold
%     mov = zeros(H,W,numel(THRs));
%     cnt = zeros(H,W,numel(THRs));
%     wrong_patches = 0;
%     patches2 = zeros(size(patches));
    parfor t = 1:numel(THRs)
        G2 = G;
        G2(~stationaryAtom(:,t),:) = 0;
        patches2 = patches - D*G2;
        
        % correct for patches with wrong energy after clutter removal:
        % make sure that Energy(patch)>= Energy(Dt*At)
%         energy = sum(abs(patches2).^2) > sum(abs(patches).^2);
%         patches2(:,energy) = patches(:,energy);
%         wrong_patches = wrong_patches + sum(energy);
        % Average the samples in all patches
        %         [movN, cntN] = joinPatches(patches2, H,W, N,P,M, OL_R,OL_AX);
        %         mov(:,:,t) = mov(:,:,t) + movN(:,:,midFrame-(startFrame-N/2)+1);
        %         cnt(:,:,t) = cnt(:,:,t) + cntN(:,:,midFrame-(startFrame-N/2)+1);
        
        [movN, cntN] =  joinPatchesMiddleMasked(patches2, H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
        movN(cntN ~=0) = movN(cntN ~=0) ./ cntN(cntN ~=0);
        mov{t}(:,:,frame) = movN;
    end
%     fprintf(', wrong energy = %d, time = %g\n',  wrong_patches, toc);
    fprintf(', time = %g\n', toc);
end


end