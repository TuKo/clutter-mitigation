function [mov,clutter] = KSVDCField2SupressClutterGlobalSeq(Dc, Dt, seq, params, startFrame, stopFrame, Mask)


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
else
    error('No error threshold is defined for OMPerr');
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

D = [Dc,Dt];
DtD = D'*D;

stationaryAtom = false(size(D,2),1);
stationaryAtom(1:size(Dc,2))= true;

mov = zeros(H,W,nFrames);
clutter = zeros(H,W,nFrames);

for frame = 1:nFrames
    % lastFrame = startFrame + N;
    midFrame = startFrame + frame -1;
    
    tic
    
    fprintf('Frame = %d ', frame);

    posMidFrame = 1+N/2;
    patches = extractPatchesMasked(seq, midFrame,  N, P, M, OL_R, OL_AX, Mask);

%     G = comp2(D'*patches, sum(conj(patches).*patches), DtD, OMP_ERROR, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)

    if params.parallel ~= 1
        G = comp(D'*patches, DtD, OMP_MAX_ATOMS); % maximum of atoms (last param)
    else
        numdata = size(patches,2);
        datachunk = ceil(numdata/4);
        localGamma = cell(1,4);
        localdiv = ones(4,1)*datachunk;
        localdiv(end) = localdiv(end) - (sum(localdiv)-numdata);
        localdata = mat2cell(patches,size(patches,1),localdiv);
        parfor i = 1:4
            localGamma{i} = comp(D'*localdata{i},DtD,OMP_MAX_ATOMS);
        end
        G = cell2mat(localGamma);    
    end
    
    fprintf('Avg atoms used: %g', mean(full(sum(G~=0))));
%     G = comp2(D, patches, OMP_ERROR, OMP_MAX_ATOMS); % maximum of atoms (last param)
%     G = comp2(D, patches*1e25, DtD, OMP_ERROR*1e25, 'maxatoms', OMP_MAX_ATOMS); % maximum of atoms (last param)

    % Extract clutter only with every value of Threshold
    
    G2 = G;
    G2(~stationaryAtom,:) = 0; % G2 marks the atoms with clutter
    clutter(:,:,frame) = joinPixels(full(max(abs(G2))),H,W, N,P,M,OL_R,OL_AX,posMidFrame); 
    patches2 = patches - D*G2;
    % Average the samples in all patches
    %         [movN, cntN] = joinPatches(patches2, H,W, N,P,M, OL_R,OL_AX);
    %         mov(:,:,t) = mov(:,:,t) + movN(:,:,midFrame-(startFrame-N/2)+1);
    %         cnt(:,:,t) = cnt(:,:,t) + cntN(:,:,midFrame-(startFrame-N/2)+1);
    [movN, cntN] =  joinPatchesMiddleMasked(patches2, H,W, N,P,M, OL_R,OL_AX, posMidFrame, Mask);
    movN(cntN ~=0) = movN(cntN ~=0) ./ cntN(cntN ~=0);
    mov(:,:,frame) = movN;

    fprintf(', time = %g\n', toc);
end


end