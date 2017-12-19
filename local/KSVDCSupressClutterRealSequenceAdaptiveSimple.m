function mov_out = KSVDCSupressClutterRealSequenceAdaptiveSimple(Dt, Dc, At, seq, params, startFrame, stopFrame, Mask)


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

[H,W,~] = size(seq);

nFrames = stopFrame-startFrame+1;

% stationaryAtom = false(size(D,2),numel(THRs));
% for t = 1:numel(THRs)
%     stationaryAtom(:,t) = findClutterAtoms(D, THRs(t), CLUTTER_METHOD, M+1, N+1, P+1);
% end

if numel(THRs) > 1
    error('This works with only one THR');
end

mov = zeros(H,W,nFrames);

% Then clean the clutter atoms from the adaptive dictionary to obtain the final Dt
Dt2 = Dt(:,~(sum(abs(Dt - Dc*At).^2)<THRs));
    
% Then merge Dt and Dc
D = [Dc , Dt2];
DtD = D'*D;
stationaryAtom = false(size(D,2),1);
stationaryAtom(1:size(Dc,2)) = true;

parfor frame = 1:nFrames
    midFrame = startFrame + frame -1;
    tic
    fprintf('Frame = %d ', frame);
    
    posMidFrame = 1+N/2;
    patches = extractPatchesMasked(seq, midFrame,  N, P, M, OL_R, OL_AX, Mask);
    
    if ERR_THR_CRITERION
        G2 = comp2(D'*patches, sum(conj(patches).*patches), DtD, OMP_ERROR, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)
        %             G2 = comp2(D, patches, OMP_ERROR, OMP_MAX_ATOMS); % maximum of atoms (last param)
    else
        G2 = comp(D'*patches, DtD, OMP_MAX_ATOMS); % maximum of atoms (last param)
    end
    fprintf('Avg atoms used: %g', mean(full(sum(G2~=0))));
    
    % Extract clutter only with every value of Threshold
    %     mov = zeros(H,W,numel(THRs));
    %     cnt = zeros(H,W,numel(THRs));
    
    %         G2 = G;
    G2(~stationaryAtom,:) = 0;
    patches2 = patches - D*G2;
    patches = [];
    % Average the samples in all patches
    %         [movN, cntN] = joinPatches(patches2, H,W, N,P,M, OL_R,OL_AX);
    %         mov(:,:,t) = mov(:,:,t) + movN(:,:,midFrame-(startFrame-N/2)+1);
    %         cnt(:,:,t) = cnt(:,:,t) + cntN(:,:,midFrame-(startFrame-N/2)+1);
    
    [movN, cntN] =  joinPatchesMiddleMasked(patches2, H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
    movN(cntN ~=0) = movN(cntN ~=0) ./ cntN(cntN ~=0);
    patches2 = [];
    fprintf(', joining, ');
    mov(:,:,frame) = movN;
    %     end
    %     fprintf(', wrong energy = %d, time = %g\n',  wrong_patches, toc);
    fprintf(', time = %g\n', toc);
end

mov_out = cell(1);
mov_out{1} = mov;

end