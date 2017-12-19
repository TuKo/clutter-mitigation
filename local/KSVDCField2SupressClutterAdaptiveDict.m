function mov = KSVDCField2SupressClutterAdaptiveDict(Dt, Dc, At, seq, params)

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
    ERR_THR_CRITERION = 0;
end
% if (isfield(params,'ATOM_THRs'))
%     THRs = params.ATOM_THRs;
% else
%     error('No threshold parameter is defined for choosing the clutter atoms');
% end

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
if (isfield(params,'cutoff'))
    CUTOFF = params.cutoff;
else
    CUTOFF = 0.4; 
end

[H,W,nFrames] = size(seq);

% Then clean the clutter atoms from the adaptive dictionary to obtain the final Dt
Dt = Dt(:,~(sum(abs(Dt - Dc*At).^2)<CUTOFF));

% Then merge Dt and Dc
D = [Dc , Dt];
DtD = D'*D;
stationaryAtom = false(size(D,2),1);
stationaryAtom(1:size(Dc,2)) = true;

mov = zeros(H,W);
cnt = zeros(H,W);

midFrame = (nFrames-1)/2+1;
firstFrame = max(midFrame-N/2 , N/2+1);
lastFrame = min(midFrame+N/2, nFrames-N/2);

% for startFrame = (1+N/2):(nFrames-N/2)
% for startFrame = firstFrame:lastFrame
for startFrame = midFrame
%     fprintf('Frame = %d', startFrame);
%     tic
    posMidFrame = midFrame-startFrame+1+N/2;
    patches = extractPatches(seq, startFrame,  N, P, M, OL_R, OL_AX);

    cc = size(patches,2);
%     G = comp2(D, patches, OMP_ERROR, OMP_MAX_ATOMS); % maximum of atoms (last param)
%     G = comp2(D, patches*1e25, DtD, OMP_ERROR*1e25, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)
    ff = patches*1e25;
    
    if ERR_THR_CRITERION
        G = comp2(D'*ff, sum(conj(ff).*ff), DtD, OMP_ERROR*1e25, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)
    else
        G = comp(D'*ff, DtD, OMP_MAX_ATOMS);
    end
    
    % Extract clutter only with every value of Threshold
%     GG = G * sparse(1:cc,1:cc,1./sqrt(sum(G.^2)));
%         GG(stationaryAtom,:) = 1;
        G2 = G;
        G2(~stationaryAtom,:) = 0;
        patches2 = patches - D*G2*1e-25;
        % Average the samples in all patches
%         [movN, cntN] = joinPatches(patches2, H,W, N,P,M, OL_R,OL_AX);
%         mov(:,:) = mov(:,:) + movN(:,:,midFrame-(startFrame-N/2)+1);
%         cnt(:,:) = cnt(:,:) + cntN(:,:,midFrame-(startFrame-N/2)+1);
        [movN, cntN] =  joinPatchesMiddle(patches2, H,W, N,P,M, OL_R,OL_AX, posMidFrame);
        mov(:,:) = mov(:,:) + movN;
        cnt(:,:) = cnt(:,:) + cntN;
%     fprintf(', time = %g\n', toc);
end

mov(cnt ~=0) = mov(cnt ~=0) ./ cnt(cnt ~=0);
mov = abs(mov);

end