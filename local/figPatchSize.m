function [first_exper, last_exper, CNR_GMCA, CNR_PF, CNR_BMODE, PSNR_GMCA, PSNR_PF, PSNR_BMODE, CNR_AMCA, PSNR_AMCA, CNR_FOMCA, PSNR_FOMCA, CNR_SVF, PSNR_SVF, CNR_FIR, PSNR_FIR] = ...
    figPatchSize(IQdata, artIQ, first_exper, last_exper, Mvalues, Nvalues, OMP_ERRORs, cutoffs, ATOM_THRs, taus)

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

% Default algorithm parameters
CLUTTER_METHODs = 1; % use svd for selecting the atoms in  fully-online MCA
Ps = 0; % use only one lateral line per patch
DICT_SIZEs = 2; % dictionary redundancy for FOMCA
ATOMS_KSVDs = 30; % 
OMP_MAX_ATOMSs = 40; % max atoms for OMP with error threshold


% Compute the Input and Output windows for all the methods
[zone1Y,zone1X] = meshgrid((565:765)-(shift./downsampleStep),6:15);
[zone2Y,zone2X] = meshgrid((565:765)-(shift./downsampleStep),65:74);
[zoneInY,zoneInX] = meshgrid((565:765)-(shift./downsampleStep),31:50);
zoneInX = zoneInX(:);
zoneInY = zoneInY(:);    
zoneOutX = [zone1X(:); zone2X(:)];
zoneOutY = [zone1Y(:); zone2Y(:)];

[PSNRzoneInY,PSNRzoneInX] = meshgrid((400:900)-(shift./downsampleStep),6:74);
PNSRzoneX = PSNRzoneInX(:);
PNSRzoneY = PSNRzoneInY(:);

% Result matrices
CNR_GMCA   = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs));
PSNR_GMCA  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs));
CNR_AMCA   = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(cutoffs));
PSNR_AMCA  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(cutoffs));
CNR_FOMCA  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
PSNR_FOMCA = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
CNR_SVF    = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
PSNR_SVF   = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
CNR_PF     = zeros(max_exper,1);
PSNR_PF    = zeros(max_exper,1);
CNR_BMODE  = zeros(max_exper,1);
PSNR_BMODE = zeros(max_exper,1);
CNR_FIR    = zeros(max_exper,1);
PSNR_FIR   = zeros(max_exper,1);

[~, ~, RFDataSample, RFData, ~] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
ClutterDataSample = RFData - RFDataSample;
clear RFData

for i = first_exper:last_exper
    fprintf('Experiment %d\n',i);
    [FinalData, Noise, RFdataClean, RFdata, NSTD] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
    [~, ~, T] = size(FinalData);
    
    RFCleanImg = RFdataClean(:,:,(T+1)/2);
    NoiseImg = Noise(:,:,(T+1)/2);
    
    params.P = 0;
    for q = 1:numel(Nvalues)
        params.N = Nvalues(q)-1;
        disp(Nvalues(q));
        for p = 1:numel(Mvalues)
            disp(Mvalues(p));
            params.M = Mvalues(p)-1;
            params.OL_M = params.M;
            params.partition = 1.0;
            params.parallel = 1;
            params.DICT_SIZE = round((params.M+1)*(params.N+1)*1.1);
            params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
            
            fprintf('GMCA DL\n');
            % learn the clutter dictionary for this size of patches
            Dc = KSVDCDoubleField2(ClutterDataSample, params);
            % learn the tissue dictionary for this size of patches
            Dt = KSVDCDoubleField2(RFDataSample, params);

            fprintf('AMCA DL\n');
            % Learn an adaptive Dictionary from the data, for the adaptive tissue
            params.DICT_SIZE = round((params.M+1)*(params.N+1)*2);
            params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
            DtA = KSVDCDoubleField2(FinalData, params);
            At = comp(Dc'*DtA, Dc'*Dc, params.M+1);

            fprintf('FOMCA DL\n');
            params.DICT_SIZE = round((params.M+1)*(params.N+1)*(params.P+1)*DICT_SIZEs);
            params.partition = 0.8;
            params.CLUTTER_METHOD = CLUTTER_METHODs;
            DtFO = KSVDCDoubleField2(FinalData, params); % Two KSVDs
            DtDFO = DtFO'*DtFO;
            
            
            % clean the noisy image
            parfor s = 1:numel(OMP_ERRORs)
                fprintf('GMCA ERR=%f\n',OMP_ERRORs(s));
                locparams = params;
                locparams.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
                locparams.parallel = 0;
                locparams.OMP_ERROR = sqrt((params.M+1)*(params.N+1)*(params.P+1)*2)*NSTD*OMP_ERRORs(s); %*std(FinalData(:));
                KSVDCimg = KSVDCField2SupressClutterGlobalDict(Dt, Dc, FinalData, locparams);
                [~,~,partial_res,~] = computeSNRandCNR(KSVDCimg, zoneInX, zoneInY, zoneOutX, zoneOutY);
                [psnr_tmpk, ~] = computeMSE(KSVDCimg, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
                CNR_GMCA(i,q,p,s) = partial_res;
                PSNR_GMCA(i,q,p,s) = psnr_tmpk;
                
                fprintf('AMCA ERR=%f\n',OMP_ERRORs(s));
                cnr_partial = zeros(numel(cutoffs),1);
                psnr_partial = zeros(numel(cutoffs),1);
                for r = 1:numel(cutoffs)
                    locparams.cutoff = cutoffs(r);
                    KSVDCimgAdaptive = KSVDCField2SupressClutterAdaptiveDict(DtA, Dc, At, FinalData, locparams);
                    [~,~,cnr_partial(r),~] = computeSNRandCNR(KSVDCimgAdaptive, zoneInX, zoneInY, zoneOutX, zoneOutY);
                    [psnr_partial(r),~] = computeMSE(KSVDCimgAdaptive, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
                end
                CNR_AMCA(i,q,p,s,:) = cnr_partial;
                PSNR_AMCA(i,q,p,s,:) = psnr_partial; 

                fprintf('FOMCA ERR=%f\n',OMP_ERRORs(s));
                cnr_partial = zeros(numel(ATOM_THRs),1);
                psnr_partial = zeros(numel(ATOM_THRs),1);                
                locparams.ATOM_THRs = ATOM_THRs;
                KSVDCimg = KSVDCField2SupressClutter_2(DtFO, DtDFO, FinalData, locparams);
                for t = 1:numel(ATOM_THRs)
                    % Compute KSVD Complex
                    [~,~,cnr_partial(t),~] = computeSNRandCNR(KSVDCimg(:,:,t), zoneInX, zoneInY, zoneOutX, zoneOutY);
                    [psnr_partial(t), ~] = computeMSE(KSVDCimg(:,:,t), abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
                end
                CNR_FOMCA(i,q,p,s,:) = cnr_partial;
                PSNR_FOMCA(i,q,p,s,:) = psnr_partial;
            end
            
            fprintf('SVF\n');
            cnr_partial = zeros(numel(taus),1);
            psnr_partial = zeros(numel(taus),1);
            MvalSVF = params.M;
            NvalSVF = params.N;
            midFrame = (size(FinalData,3)-1)/2 +1;
            parfor t = 1:numel(taus)
                tau = taus(t);
                % Compute SVF results
                SVFimg = SVF_oneFrame(FinalData, MvalSVF, NvalSVF, 30, tau);
                SVFimg = abs(SVFimg(:,:,midFrame));
                [~,~,cnr_partial(t),~] = computeSNRandCNR(SVFimg, zoneInX, zoneInY, zoneOutX, zoneOutY);
                [psnr_partial(t),~] = computeMSE(SVFimg, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
            end
            CNR_SVF(i,q,p,:) = cnr_partial;
            PSNR_SVF(i,q,p,:) = psnr_partial;        
        end
    end

    % Compute B-mode (unfiltered image)
    fprintf('Unfiltered\n');
    Bmode = abs(FinalData(:,:,(T+1)/2));
    [~,~,CNR_BMODE(i),~] = computeSNRandCNR(Bmode, zoneInX, zoneInY, zoneOutX, zoneOutY);
    [PSNR_BMODE(i), ~] = computeMSE(Bmode, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
    
    
    % Compute for Perfect filter
    fprintf('Perfect Filtering\n');
    perfectFil = abs(RFCleanImg+NoiseImg);
    [~,~,CNR_PF(i), ~] = computeSNRandCNR(perfectFil, zoneInX, zoneInY, zoneOutX, zoneOutY);
    [PSNR_PF(i), ~] = computeMSE(perfectFil, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);

    % Compute FIR method
    fprintf('FIR\n');    
    FIR = abs(FinalData(:,:,(T+1)/2+1)-FinalData(:,:,(T+1)/2));
    [~,~,CNR_FIR(i),~] = computeSNRandCNR(FIR, zoneInX, zoneInY, zoneOutX, zoneOutY);
    [PSNR_FIR(i), ~] = computeMSE(FIR, abs(RFCleanImg), PNSRzoneX, PNSRzoneY);
end

end