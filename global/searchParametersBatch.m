function [first_exper, last_exper, ...
    CNR_KSVDC, CNR2_KSVDC, PSNR_KSVDC, MSE_KSVDC, PSNR2_KSVDC, MSE2_KSVDC, ...
    CNR_SVF, CNR2_SVF, PSNR_SVF, MSE_SVF, PSNR2_SVF, MSE2_SVF, ...
    CNR_FIR, CNR2_FIR, PSNR_FIR, MSE_FIR, PSNR2_FIR, MSE2_FIR, ...
    CNR_BMODE, CNR2_BMODE, PSNR_BMODE, MSE_BMODE, PSNR2_BMODE, MSE2_BMODE, ...
    CNR_PF, CNR2_PF, PSNR_PF, MSE_PF ...
    ] = searchParametersBatch(IQdata, artIQ, first_exper, last_exper, Mvalues, Nvalues, taus, OMP_ERRORs, ATOM_THRs)

% addpath('ompbox/');
% General Parameters
max_exper = last_exper-first_exper+1; % total  of experiments in this job

shift = 80; % tissue displacement = 1
art_shift = 10; % clutter displacement = 1/8
downsampleStep = 10; % downsampling for sequence creation
rho = 0.98;
SNR = 30;

% SVF parameters
% taus = 0.35:0.025:0.95;

% KSVD parameters
Ps = 0;
DICT_SIZEs = 2; % dictionary redundancy
ATOMS_KSVDs = 30;
OMP_MAX_ATOMSs = 40; % max atoms for OMP with error threshold
CLUTTER_METHODs = 1; % use svd for selecting the atoms
% OMP_ERRORs = 0.8:0.1:2.3;
% ATOM_THRs = 0.35:0.025:0.95;

% Results tables
CNR_KSVDC  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
CNR2_KSVDC = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
PSNR_KSVDC = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
MSE_KSVDC  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
PSNR2_KSVDC = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));
MSE2_KSVDC  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(OMP_ERRORs),numel(ATOM_THRs));

CNR_SVF  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
CNR2_SVF = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
PSNR_SVF = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
MSE_SVF  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
PSNR2_SVF = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));
MSE2_SVF  = zeros(max_exper,numel(Nvalues),numel(Mvalues),numel(taus));

CNR_FIR  = zeros(max_exper,1);
CNR2_FIR = zeros(max_exper,1);
PSNR_FIR = zeros(max_exper,1);
MSE_FIR  = zeros(max_exper,1);
PSNR2_FIR = zeros(max_exper,1);
MSE2_FIR  = zeros(max_exper,1);

CNR_BMODE  = zeros(max_exper,1);
CNR2_BMODE = zeros(max_exper,1);
PSNR_BMODE = zeros(max_exper,1);
MSE_BMODE  = zeros(max_exper,1);
PSNR2_BMODE = zeros(max_exper,1);
MSE2_BMODE  = zeros(max_exper,1);

CNR_PF  = zeros(max_exper,1);
CNR2_PF = zeros(max_exper,1);
PSNR_PF = zeros(max_exper,1);
MSE_PF  = zeros(max_exper,1);

% Compute the Input and Output windows for all the methods
[zone1Y,zone1X] = meshgrid((565:765)-(shift./downsampleStep),6:15);
[zone2Y,zone2X] = meshgrid((565:765)-(shift./downsampleStep),65:74);
[zoneInY,zoneInX] = meshgrid((565:765)-(shift./downsampleStep),31:50);
zoneInX = zoneInX(:);
zoneInY = zoneInY(:);    
zoneOutX = [zone1X(:); zone2X(:)];
zoneOutY = [zone1Y(:); zone2Y(:)];

% Windows for measuring PSNR
% dpthB = ((1:size(Bmode,1))./fs).*(1540/2);
% d_x = x_size/(2*nLines); %=.05/2/180
% lat = (-39:40).*d_x;

[PSNRzoneInY,PSNRzoneInX] = meshgrid((400:900)-(shift./downsampleStep),6:74);
PNSRzoneX = PSNRzoneInX(:);
PNSRzoneY = PSNRzoneInY(:);

% Run the algorithms
params.P = 0; % radial size of patches (usually 0, only axial-time)
params.OL_P = params.P/2;
params.OMP_MAX_ATOMS = OMP_MAX_ATOMSs;
params.CLUTTER_METHOD = CLUTTER_METHODs;
params.ATOMS_KSVD = ATOMS_KSVDs;

for i = first_exper:last_exper
    % Create the sequence
    [FinalData, Noise, RFdataClean, RFdata, NSTD] = createSequence(shift, art_shift, downsampleStep, rho, SNR, IQdata, artIQ);
    [~, ~, T] = size(FinalData);
    
    % Run KSVD and SVF
    for q = 1:numel(Nvalues)
        params.N = Nvalues(q)-1;
        disp(Nvalues(q));
        for p = 1:numel(Mvalues)
            disp(Mvalues(p));
            params.M = Mvalues(p)-1;
            params.OL_M = params.M;
            params.ATOMS_KSVD = round(0.10*(params.M+1)*(params.N+1)*(params.P+1));
            params.DICT_SIZE = round((params.M+1)*(params.N+1)*(params.P+1)*DICT_SIZEs);
    %         IQDictionary = KSVDCField2(FinalData, params); % One KSVD
            IQDictionary = KSVDCDoubleField2(FinalData, params); % Two KSVDs
            DtDIQ = IQDictionary'*IQDictionary;

            fprintf('KSVDC\n');
            tic
            parfor s = 1:numel(OMP_ERRORs)
                partial_res = zeros(numel(ATOM_THRs),1);
                partial2 = zeros(numel(ATOM_THRs),1);
                mse_tmpk = zeros(numel(ATOM_THRs),1);
                mse_tmpk2 = zeros(numel(ATOM_THRs),1);
                psnr_tmpk2 = zeros(numel(ATOM_THRs),1);
                psnr_tmpk = zeros(numel(ATOM_THRs),1);
                locparams = params;
                locparams.OMP_ERROR = sqrt((locparams.M+1)*(locparams.N+1)*(locparams.P+1)*2)*NSTD*OMP_ERRORs(s); %*std(FinalData(:));
                locparams.ATOM_THRs = ATOM_THRs;
                KSVDCimg = KSVDCField2SupressClutter_2(IQDictionary, DtDIQ, FinalData, locparams);
                for t = 1:numel(ATOM_THRs)
                    % Compute KSVD Complex
                    [~,~,partial_res(t),partial2(t)] = computeSNRandCNR(KSVDCimg(:,:,t), zoneInX, zoneInY, zoneOutX, zoneOutY);
                    [psnr_tmpk(t), mse_tmpk(t)] = computeMSE(KSVDCimg(:,:,t), abs(RFdataClean(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
                    [psnr_tmpk2(t), mse_tmpk2(t)] = computeMSE(KSVDCimg(:,:,t), abs(RFdata(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
                end
                CNR_KSVDC(i,q,p,s,:) = partial_res;
                CNR2_KSVDC(i,q,p,s,:) = partial2;
                PSNR_KSVDC(i,q,p,s,:) = psnr_tmpk;
                PSNR2_KSVDC(i,q,p,s,:) = psnr_tmpk2;
                MSE_KSVDC(i,q,p,s,:) = mse_tmpk;
                MSE2_KSVDC(i,q,p,s,:) = mse_tmpk2;
            end
            toc

            fprintf('SVF\n');
            partial_res = zeros(numel(taus),1);
            partial2 = zeros(numel(taus),1);
            mse_tmps = zeros(numel(taus),1);
            mse_tmps2 = zeros(numel(taus),1);
            psnr_tmps2 = zeros(numel(taus),1);
            psnr_tmps = zeros(numel(taus),1);
            MvalSVF = params.M;
            NvalSVF = params.N;
            midFrame = (size(FinalData,3)-1)/2 +1;
            parfor t = 1:numel(taus)
                tau = taus(t);
                % Compute SVF results
                SVFimg = SVF_oneFrame(FinalData, MvalSVF, NvalSVF, 30, tau);
                SVFimg = abs(SVFimg(:,:,midFrame));
                [~,~,partial_res(t),partial2(t)] = computeSNRandCNR(SVFimg, zoneInX, zoneInY, zoneOutX, zoneOutY);
                [psnr_tmps(t), mse_tmps(t)] = computeMSE(SVFimg, abs(RFdataClean(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
                [psnr_tmps2(t), mse_tmps2(t)] = computeMSE(SVFimg, abs(RFdata(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
            end
            CNR_SVF(i,q,p,:) = partial_res;
            CNR2_SVF(i,q,p,:) = partial2;
            PSNR_SVF(i,q,p,:) = psnr_tmps;
            PSNR2_SVF(i,q,p,:) = psnr_tmps2;
            MSE_SVF(i,q,p,:) = mse_tmps;
            MSE2_SVF(i,q,p,:) = mse_tmps2;
        end
    end

    aux = squeeze(CNR_KSVDC(i,:,:,:,:));
    fprintf('KSVD: %f\n',max(aux(:)));
    aux = squeeze(CNR2_KSVDC(i,:,:,:,:));
    fprintf('KSVD2: %f\n',max(aux(:)));
    aux = squeeze(CNR_SVF(i,:,:,:));
    fprintf('SVF: %f\n',max(aux(:)));
    aux = squeeze(CNR2_SVF(i,:,:,:));
    fprintf('SVF2: %f\n',max(aux(:)));

    % Compute FIR method
    FIR = abs(RFdata(:,:,(T+1)/2+1)-RFdata(:,:,(T+1)/2));
    [~,~,CNR_FIR(i),CNR2_FIR(i)] = computeSNRandCNR(FIR, zoneInX, zoneInY, zoneOutX, zoneOutY);
    [PSNR_FIR(i), MSE_FIR(i)] = computeMSE(FIR, abs(RFdataClean(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
    [PSNR2_FIR(i), MSE2_FIR(i)] = computeMSE(FIR, abs(RFdata(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
    fprintf('FIR: %f\n',CNR_FIR(i));
    fprintf('FIR2: %f\n',CNR2_FIR(i));

    % Compute B-mode (unfiltered image)
    Bmode = abs(RFdata(:,:,(T+1)/2));
    [~,~,CNR_BMODE(i),CNR2_BMODE(i)] = computeSNRandCNR(Bmode, zoneInX, zoneInY, zoneOutX, zoneOutY);
    [PSNR_BMODE(i), MSE_BMODE(i)] = computeMSE(Bmode, abs(RFdataClean(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
    %[PSNR2_BMODE(i), MSE2_BMODE(i)] = computeMSE(Bmode, abs(RFdata(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
    fprintf('Unfiltered: %f\n',CNR_BMODE(i));
    fprintf('Unfiltered2: %f\n',CNR2_BMODE(i));

    % Compute for Perfect filter
    perfectFil = abs(RFdataClean(:,:,(T+1)/2));
    [~,~,CNR_PF(i),CNR2_PF(i)] = computeSNRandCNR(perfectFil, zoneInX, zoneInY, zoneOutX, zoneOutY);
    [PSNR_PF(i), MSE_PF(i)] = computeMSE(perfectFil, abs(RFdata(:,:,(T+1)/2)), PNSRzoneX, PNSRzoneY);
    fprintf('PerfectFilt: %f\n',CNR_PF(i));
    fprintf('PerfectFilt2: %f\n',CNR2_PF(i));
end

end