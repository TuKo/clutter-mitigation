function [first_exper, last_exper, CNR_GMCA, CNR_PF, CNR_BMODE, PSNR_GMCA, PSNR_PF, PSNR_BMODE, CNR_AMCA, PSNR_AMCA, CNR_FOMCA, PSNR_FOMCA, CNR_SVF, PSNR_SVF, CNR_FIR, PSNR_FIR, CNR_GMCA2, PSNR_GMCA2] = ...
    figMotion(IQdata, artIQ, first_exper, last_exper, shifts, Mvalue, Nvalue, OMP_ERROR, cutoffs, taus, ATOM_THRs) 

% Before sending to the cluster

% load simulatedData40.mat
% load D:\javier_temp\simulatedData40.mat
% artIQ = artifactIQdata;
% IQdata = IQdata(:,:,1:20);
% clear artifactIQdata;

% function
max_exper = last_exper-first_exper+1; % total  of experiments in this job

shift = 80; % tissue displacement = 1
art_shift = 10; % clutter displacement = 1/8
downsampleStep = 10; % downsampling for sequence creation
rho = 0.98;
SNR = 30;


CLUTTER_METHODs = 1; % use svd for selecting the atoms in  fully-online MCA
DICT_SIZEs = 2; % dictionary redundancy for FOMCA

% Result matrices
CNR_GMCA   = zeros(max_exper,numel(shifts));
PSNR_GMCA  = zeros(max_exper,numel(shifts));
CNR_GMCA2  = zeros(max_exper,numel(shifts));
PSNR_GMCA2 = zeros(max_exper,numel(shifts));
CNR_AMCA   = zeros(max_exper,numel(shifts),numel(cutoffs));
PSNR_AMCA  = zeros(max_exper,numel(shifts),numel(cutoffs));
CNR_FOMCA  = zeros(max_exper,numel(shifts),numel(ATOM_THRs));
PSNR_FOMCA = zeros(max_exper,numel(shifts),numel(ATOM_THRs));
CNR_SVF    = zeros(max_exper,numel(shifts),numel(taus));
PSNR_SVF   = zeros(max_exper,numel(shifts),numel(taus));
CNR_PF     = zeros(max_exper,numel(shifts));
PSNR_PF    = zeros(max_exper,numel(shifts));
CNR_BMODE  = zeros(max_exper,numel(shifts));
PSNR_BMODE = zeros(max_exper,numel(shifts));
CNR_FIR    = zeros(max_exper,numel(shifts));
PSNR_FIR   = zeros(max_exper,numel(shifts));

[~, ~, RFDataSample, RFData, ~] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
ClutterDataSample = RFData - RFDataSample;
clear RFData

% Run the algorithms
params.P = 0; % radial size of patches (usually 0, only axial-time)
params.OL_P = params.P/2;
params.N = Nvalue-1;
params.M = Mvalue-1;
params.OL_M = params.M;
params.partition = 1.0;
params.parallel = 1;
params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
params.DICT_SIZE = round((params.M+1)*(params.N+1)*1.1);

% Offline dictionary learning:
% learn the clutter dictionary for this size of patches
Dc = KSVDCDoubleField2(ClutterDataSample, params);
% learn the tissue dictionary for this size of patches
Dt = KSVDCDoubleField2(RFDataSample, params);

params.parallel = 0;
params.CLUTTER_METHOD = CLUTTER_METHODs;

for i = first_exper:last_exper
    fprintf('Experiment %d\n',i);

    parfor q = 1:numel(shifts)
        disp(shifts(q));
        sh = shifts(q);
        [FinalData, Noise, RFdataClean, ~, NSTD] = createSequence(sh, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
        [~, ~, T] = size(FinalData);

        RFCleanImg = RFdataClean(:,:,(T+1)/2);
        NoiseImg = Noise(:,:,(T+1)/2);
        
        % Compute the Input and Output windows for all the methods
        [zone1Y,zone1X] = meshgrid((565:765)-(sh./downsampleStep),6:15);
        [zone2Y,zone2X] = meshgrid((565:765)-(sh./downsampleStep),65:74);
        [zoneInY,zoneInX] = meshgrid((565:765)-(sh./downsampleStep),31:50);
        zoneInX = zoneInX(:);
        zoneInY = zoneInY(:);
        zoneOutX = [zone1X(:); zone2X(:)];
        zoneOutY = [zone1Y(:); zone2Y(:)];
        
        [PSNRzoneInY,PSNRzoneInX] = meshgrid((400:900)-(sh./downsampleStep),6:74);
        PNSRzoneX = PSNRzoneInX(:);
        PNSRzoneY = PSNRzoneInY(:);

        
        fprintf('GMCA DL\n');
        locparams = params;
        % the clutter dictionary is learned outside the loops
        % the global tissue dictionary is also learned outside the loops
        
        % learn the tissue dictionary for a sequence with correct
        % tissue velocity
        locparams.DICT_SIZE = round((locparams.M+1)*(locparams.N+1)*1.1);
        [~, ~, RFDataSampleForTissue, ~, ~] = createSequence(sh, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
        Dt2 = KSVDCDoubleField2(RFDataSampleForTissue, locparams);
        
        fprintf('GMCA\n');
        locparams.OMP_ERROR = sqrt((locparams.M+1)*(locparams.N+1)*2)*NSTD*OMP_ERROR(1);
        KSVDCimg = KSVDCField2SupressClutterGlobalDict(Dt, Dc, FinalData, locparams);
        [~,~,cnr_partial,~] = computeSNRandCNR(KSVDCimg, zoneInX, zoneInY, zoneOutX, zoneOutY);
        [psnr_partial, ~] = computeMSE(KSVDCimg, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        CNR_GMCA(i,q) = cnr_partial;
        PSNR_GMCA(i,q) = psnr_partial;

        KSVDCimg2 = KSVDCField2SupressClutterGlobalDict(Dt2, Dc, FinalData, locparams);
        [~,~,cnr_partial,~] = computeSNRandCNR(KSVDCimg2, zoneInX, zoneInY, zoneOutX, zoneOutY);
        [psnr_partial, ~] = computeMSE(KSVDCimg2, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        CNR_GMCA2(i,q) = cnr_partial;
        PSNR_GMCA2(i,q) = psnr_partial;
        
        
        fprintf('AMCA DL\n');
        locparams = params;
        locparams.DICT_SIZE = round((locparams.M+1)*(locparams.N+1)*2);
        % Learn an adaptive Dictionary from the data, for the adaptive tissue
        DtA = KSVDCDoubleField2(FinalData, locparams);
        At = comp(Dc'*DtA, Dc'*Dc, locparams.M+1);
         
        fprintf('AMCA\n');
        locparams.OMP_ERROR = sqrt((locparams.M+1)*(locparams.N+1)*2)*NSTD*OMP_ERROR(2);
        cnr_partial = zeros(numel(cutoffs),1);
        psnr_partial = zeros(numel(cutoffs),1);
        for r = 1:numel(cutoffs)
            locparams.cutoff = cutoffs(r);
            KSVDCimgAdaptive = KSVDCField2SupressClutterAdaptiveDict(DtA, Dc, At, FinalData, locparams);
            [~,~,cnr_partial(r),~] = computeSNRandCNR(KSVDCimgAdaptive, zoneInX, zoneInY, zoneOutX, zoneOutY);
            [psnr_partial(r),~] = computeMSE(KSVDCimgAdaptive, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        end
        CNR_AMCA(i,q,:) = cnr_partial;
        PSNR_AMCA(i,q,:) = psnr_partial;
        
        
        fprintf('FOMCA DL\n');
        locparams = params;
        locparams.DICT_SIZE = round((locparams.M+1)*(locparams.N+1)*DICT_SIZEs);
        locparams.partition = 0.8;
        DtFO = KSVDCDoubleField2(FinalData, locparams); % Two KSVDs
        DtDFO = DtFO'*DtFO;
        
        fprintf('FOMCA\n');
        locparams.OMP_ERROR = sqrt((locparams.M+1)*(locparams.N+1)*2)*NSTD*OMP_ERROR(3);
        locparams.ATOM_THRs = ATOM_THRs;        
        cnr_partial = zeros(numel(ATOM_THRs),1);
        psnr_partial = zeros(numel(ATOM_THRs),1);
        KSVDCimg = KSVDCField2SupressClutter_2(DtFO, DtDFO, FinalData, locparams);
        for t = 1:numel(ATOM_THRs)
            % Compute KSVD Complex
            [~,~,cnr_partial(t),~] = computeSNRandCNR(KSVDCimg(:,:,t), zoneInX, zoneInY, zoneOutX, zoneOutY);
            [psnr_partial(t), ~] = computeMSE(KSVDCimg(:,:,t), abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        end
        CNR_FOMCA(i,q,:) = cnr_partial;
        PSNR_FOMCA(i,q,:) = psnr_partial;
        
        
        fprintf('SVF\n');
        cnr_partial  = zeros(numel(taus),1);
        psnr_partial = zeros(numel(taus),1);
        MvalSVF = params.M;
        NvalSVF = params.N;
        midFrame = (size(FinalData,3)-1)/2 +1;
        tic
        for t = 1:numel(taus)
            tau = taus(t);
            % Compute SVF results
            SVFimg = SVF_oneFrame(FinalData, MvalSVF, NvalSVF, 30, tau);
            SVFimg = abs(SVFimg(:,:,midFrame));
            [~,~,cnr_partial(t),~] = computeSNRandCNR(SVFimg, zoneInX, zoneInY, zoneOutX, zoneOutY);
            [psnr_partial(t),~] = computeMSE(SVFimg, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        end
        toc
        CNR_SVF(i,q,:) = cnr_partial;
        PSNR_SVF(i,q,:) = psnr_partial;
        
        % Compute B-mode (unfiltered image)
        fprintf('Unfiltered\n');
        Bmode = abs(FinalData(:,:,(T+1)/2));
        [~,~,CNR_BMODE(i,q),~] = computeSNRandCNR(Bmode, zoneInX, zoneInY, zoneOutX, zoneOutY);
        [PSNR_BMODE(i,q), ~] = computeMSE(Bmode, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        
        
        % Compute for Perfect filter
        fprintf('Perfect Filtering\n');
        perfectFil = abs(RFCleanImg+NoiseImg);
        [~,~,CNR_PF(i,q), ~] = computeSNRandCNR(perfectFil, zoneInX, zoneInY, zoneOutX, zoneOutY);
        [PSNR_PF(i,q), ~] = computeMSE(perfectFil, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
        
        % Compute FIR method
        fprintf('FIR\n');
        FIR = abs(FinalData(:,:,(T+1)/2+1)-FinalData(:,:,(T+1)/2));
        [~,~,CNR_FIR(i,q),~] = computeSNRandCNR(FIR, zoneInX, zoneInY, zoneOutX, zoneOutY);
        [PSNR_FIR(i,q), ~] = computeMSE(FIR, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
    end
    
end

% save D:\javier_temp\resultsMotion.mat

end