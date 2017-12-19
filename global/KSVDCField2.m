function D = KSVDCField2(seq, params)
% Cleans clutter using the sparse representations in IQ data.
% SEQ - sequence as a 3D matrix with complex data.
% M - number of axial samples per patch
% N - number of frame samples per patch
% P - number of radial samples per patch
% DICT_SIZE  - Dictionary size
% ATOMS_KSVD - Sparsity for KSVD, the dictionary computation
% OMP_ERROR  - Error Threshold for OMPerr
% ATOM_THR   - Threshold for choosing an atom as clutter/tissue
% [OMP_MAX_ATOMS] - Maximum number of atoms when using OMPerr (optional)
% 
% 
% The algorithms learns a dictionary D using ksvd on the first N frames and
% for every (overlapped) patch.
% Next, the stationary atoms are found by looking at the first singular
% value of each atom as a 2D patch, or by an horizontal (time) Total 
% Variation measure.
% Once the stationary atoms are known, the whole sequence is cleaned frame
% by frame by extracting the overlapping patches, finding the patches'
% sparse representations, and then zeroing the coefficients of the
% stationary atoms. The cleaned patches are merged back to a cleaned
% sequence.


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

if (isfield(params,'DICT_SIZE'))
    DICT_SIZE = params.DICT_SIZE;
else
    DICT_SIZE = (M+1)*(N+1)*(P+1)*2; %redundancy 2 by default
end
if DICT_SIZE < round((N+1)*(M+1)*(P+1)) 
    error('Wrong number of atoms. Please choose a bigger dictionary size.');
end

if (isfield(params,'ATOMS_KSVD'))
    ATOMS_KSVD = params.ATOMS_KSVD;
else
    ATOMS_KSVD = 5; % default # of atoms to use
end
if (isfield(params,'OL_M'))
    OL_AX_DICT = params.OL_M;
else
    OL_AX_DICT = round(M/2); % Maximum axial overlap by default (from 0 to M/2)
end
if (isfield(params,'OL_P'))
    OL_R_DICT = params.OL_P;
else
    OL_R_DICT = 0; % overlap in angle axis (from 0 to P/2)
end

% OL_AX_DICT = M/2; % maximum overlap in axial for dictionary (from 0 to M/2)
% OL_AX_DICT = M; % maximum overlap in axial for dictionary (from 0 to M/2)
% OL_R_DICT = P/2; % overlap in angle axis for dictionary (from 0 to P/2)
offsetFrame = N/2+1;

% Dictionay training with KSVD:
% Took some 2D (spatial-time) patches of N+1 (time) by M+1 (axial) by P+1
% (angular). Then, run KSVD on these patches.
fprintf('Building dictionary\n');
offsetFrame = max(offsetFrame,N/2+1);


offsetFrame = (size(seq,3)+1)/2;
% patches = extractPatches(seq, offsetFrame, N, P, M, OL_R_DICT, OL_AX_DICT);
params.data = extractPatches(seq, offsetFrame, N, P, M, OL_R_DICT, OL_AX_DICT) * 1e25;

% params.exact = 1;
params.exact = 0;
% params.data = patches * 1e25;
params.Tdata = ATOMS_KSVD;
% params.Edata = params.Edata *1e25;
% params.maxatoms = params.OMP_MAX_ATOMS;
params.dictsize = DICT_SIZE;
% params.initdict = ;
params.memusage = 'high';
params.iternum = 5;
% params.iternum = 10;
params.parallel = 1;
% params.muthresh = 0.95;
[D,~] = cksvd(params,'itr');
% DICT_SIZE = size(D,2);

end