function [CNR_ONMCA, CNR_OFFMCA, CNR_TAMCA, CNR_SVF, CNR_UN] = ...
    figRealData(datasets, clutterset, tissueset, Mvalue, Nvalue, TAUs, BETAs, CUTOFFs)

exper = numel(datasets);

alpha = 30;

CNR_OFFMCA = zeros(exper,1);
CNR_TAMCA  = zeros(exper,numel(CUTOFFs));
CNR_ONMCA  = zeros(exper,numel(BETAs));
CNR_SVF    = zeros(exper,numel(TAUs));
CNR_UN     = zeros(exper,1);


PSNR_OFFMCA = zeros(exper,1);
PSNR_TAMCA  = zeros(exper,numel(CUTOFFs));
PSNR_ONMCA   = zeros(exper,numel(BETAs));
PSNR_SVF    = zeros(exper,numel(TAUs));


% learn clutter dictionary for TA-MCA and OFF-MCA
tic
params = [];
params.N = Nvalue-1;
params.M = Mvalue-1;
params.OL_M = params.M;
params.parallel = true;
params.partition = 1.0;
params.exact = 0;
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.DICT_SIZE = 1.1*(params.N+1)*(params.M+1);
[clutterseq, ~] = readDataset(clutterset);
Dc = KSVDCDoubleFieldRealSequence(clutterseq, params);

% learn tissue dictionary for OFF-MCA
[tissueseq, ~] = readDataset(tissueset);
Dt = KSVDCDoubleFieldRealSequence(tissueseq, params);
toc

clear clutterseq tissueseq

for j = 1:numel(datasets)
    fprintf('Dataset %d - %s\n',j,datasets{j});
    % read the dataset
    [seq, nMLAs] = readDataset(datasets{j});

    % read the measurement window data
    [InX, InY, OutX, OutY] = readMeasurementWindow(datasets{j});
    
    [height,width,nFrames] = size(seq);    
    
    Mask = false(height,width);
    Mask(InY, InX) = true;
    Mask(OutY, OutX) = true;
    
    startFrame = max(Nvalue-1)/2+1;
    stopFrame = nFrames - max(Nvalue-1)/2;
    
    % Measure Unfiltered
    fprintf('Unfiltered\n');
    seq_un = seq(:,:,startFrame:stopFrame);
    [CNR_UN(j),~] = computeSNRandCNRsequence(abs(seq_un), InX, InY, OutX, OutY);
    
    % Run SVF
    fprintf('SVF\n');
    tic
    CNR_tmp = zeros(numel(TAUs),1);
    PSNR_tmp = zeros(numel(TAUs),1);
    for l = 1:numel(TAUs)
        seq_svf = SVF(seq, Mvalue-1, Nvalue-1, alpha, TAUs(l), startFrame, stopFrame, Mask);
        [CNR_tmp(l), ~] = computeSNRandCNRsequence(abs(seq_svf), InX, InY, OutX, OutY);
        PSNR_tmp(l) = computePSNRsequence(abs(seq_svf), abs(seq_un), OutX, OutY);
    end
    CNR_SVF(j,:) = CNR_tmp;
    PSNR_SVF(j,:) = PSNR_tmp;
    toc
    
    % Run ON-MCA
    fprintf('ON-MCA\n');
    % use default values for: P, OL_P, DICT_SIZE (2:1),
    params = [];
    params.N = Nvalue-1;
    params.M = Mvalue-1;
    params.OL_M = params.M;
    params.parallel = true;
    params.partition = 0.8;
    params.exact = 0;
    params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
    params.DICT_SIZE = 2*(params.N+1)*(params.M+1);      
    tic
    D = KSVDCDoubleFieldRealSequence(seq, params);
    DtD = D'*D;
    toc
    
    params.ATOM_THRs = BETAs;
    params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
    seq_onmca = KSVDCSupressClutterRealSequence(D, DtD, seq, params, startFrame, stopFrame, Mask);
    CNR_tmp =zeros(numel(BETAs),1);
    PSNR_tmp = zeros(numel(BETAs),1);
    parfor l = 1:numel(BETAs)
        [CNR_tmp(l),~] = computeSNRandCNRsequence(abs(seq_onmca{l}), InX, InY, OutX, OutY);
        PSNR_tmp(l) = computePSNRsequence(abs(seq_onmca{l}), abs(seq_un), OutX, OutY);
    end
    CNR_ONMCA(j,:) = CNR_tmp;
    PSNR_ONMCA(j,:) = PSNR_tmp;
    
    % Run OFF-MCA
    fprintf('OFF-MCA\n');
    params = [];
    params.N = Nvalue-1;
    params.M = Mvalue-1;
    params.OL_M = params.M;
    params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
    seq_offmca = KSVDCSupressClutterRealSequenceGlobalDict(Dt, Dc, seq, params, startFrame, stopFrame, Mask);
    [CNR_OFFMCA(j),~] = computeSNRandCNRsequence(abs(seq_offmca), InX, InY, OutX, OutY);
    PSNR_OFFMCA(j) = computePSNRsequence(abs(seq_offmca), abs(seq_un), OutX, OutY);
    
    % Run TA-MCA % after OFF-MCA to reuse the clutter dictionary
    fprintf('TA-MCA\n');
    % learn adaptive tissue dictionary
    params = [];
    params.N = Nvalue-1;
    params.M = Mvalue-1;
    params.OL_M = params.M;
    params.parallel = true;
    params.partition = 1.0;
    params.exact = 0;
    params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
    params.DICT_SIZE = 2*(params.N+1)*(params.M+1);
    tic
    DtA = KSVDCDoubleFieldRealSequence(seq, params);
    At = comp(Dc'*DtA, Dc'*Dc, params.M+1);
    toc

    params.CUTTOFFs = CUTOFFs;
    params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
    seq_tamca = KSVDCSupressClutterRealSequenceAdaptive(DtA, Dc, At, seq, params, startFrame, stopFrame, Mask);
    CNR_tmp = zeros(numel(CUTOFFs),1);
    PSNR_tmp = zeros(numel(CUTOFFs),1);
    parfor l = 1:numel(CUTOFFs)
        [CNR_tmp(l),~] = computeSNRandCNRsequence(abs(seq_tamca{l}), InX, InY, OutX, OutY);
        PSNR_tmp(l) = computePSNRsequence(abs(seq_tamca{l}), abs(seq_un), OutX, OutY);
    end
    CNR_TAMCA(j,:) = CNR_tmp;
    PSNR_TAMCA(j,:) = PSNR_tmp;
    
    save('results\tempdata.mat');
end

save('results\psnr.mat','PSNR_ONMCA','PSNR_TAMCA','PSNR_OFFMCA','PSNR_SVF');

end