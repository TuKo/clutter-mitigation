function [seq_UNF, seq_FIR, seq_SVF, seq_MCA, nMLAs] = createSequencesRealdata(dataset,Msvf,Nsvf,alpha,tau, Mmca, Nmca, AtomTHR)

seq_UNF = 0;
seq_FIR = 0;
seq_SVF = 0;
seq_MCA = 0;

[seq, nMLAs] = readDataset(dataset);

% CNR areas only
% read the measurement window data
% [InX, InY, OutX, OutY] = readMeasurementWindow(dataset);
% [height,width,nFrames] = size(seq);
% Mask = false(height,width);
% Mask(InY, InX) = true;
% Mask(OutY, OutX) = true;

% full image
[height,width,nFrames] = size(seq);
Mask = true(height,width);

Nvalues = 19;
startFrame = max(Nvalues-1)/2+1;
stopFrame = nFrames - max(Nvalues-1)/2;
    
% Measure Unfiltered
seq_UNF = seq(:,:,startFrame:stopFrame);

% Run FIR
seq_FIR = seq(:,:,startFrame:stopFrame)-seq(:,:,[startFrame:stopFrame]-1);

% Run SVF
fprintf('SVF\n');
tic
seq_SVF = SVF(seq, Msvf-1, Nsvf-1, alpha, tau, startFrame, stopFrame, Mask);
toc
% seq_SVF = seq_SVF{1};

seq_MCA = zeros(size(seq_SVF));

% Run MCA
fprintf('MCA\n');
params = [];
params.N = Nmca-1;
params.M = Mmca-1;
params.OL_M = params.M;
params.parallel = true;
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.partition = 1.0;
params.DICT_SIZE = (params.M+1)*(params.N+1)*4;

% Learn from data
D = KSVDCDoubleFieldRealSequence(seq, params);
DtD = D'*D;

% Remove clutter
params.ATOM_THRs = AtomTHR;
% params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
% params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
% params.OMP_MAX_ATOMS = (params.N+1)*(params.M+1)/2;
params.OMP_MAX_ATOMS = 45;
params.OL_M = params.M;
params.OMP_ERROR = 1; %sqrt(2*(locparams.M+1)*(locparams.N+1))*NSTD*OMP_THRs(l);
seq_MCA = KSVDCField2SupressClutterSeq(D, DtD, seq, params, startFrame, stopFrame, Mask);
seq_MCA = seq_MCA{1};

end
