function [CNR_MCA, CNR2_MCA, CNR_SVF, CNR2_SVF, ...
          CNR_FIR, CNR2_FIR, CNR_UNF, CNR2_UNF, ...
          PSNR_MCA, PSNR_SVF, PSNR_FIR, ...
          SNR_UNF, SNR_MCA, SNR_SVF, SNR_FIR ...
      ] = searchParametersRealdataBatch(dataset, Mvalue, Nvalues,  taus, alphas, ATOM_THRs,OMP_THRs)

% search for parameters for all the real datasets

warning('full overlap in M or half is enough?');
% if isempty(gcp)
%     parpool('local');
% end

CNR_MCA = zeros(numel(dataset),numel(Nvalues),numel(ATOM_THRs));
CNR2_MCA = zeros(numel(dataset),numel(Nvalues),numel(ATOM_THRs));
PSNR_MCA = zeros(numel(dataset),numel(Nvalues),numel(ATOM_THRs));
SNR_MCA = zeros(numel(dataset),numel(Nvalues),numel(ATOM_THRs));
CNR_SVF = zeros(numel(dataset),numel(Nvalues),numel(taus),numel(alphas));
CNR2_SVF = zeros(numel(dataset),numel(Nvalues),numel(taus),numel(alphas));
PSNR_SVF = zeros(numel(dataset),numel(Nvalues),numel(taus),numel(alphas));
SNR_SVF = zeros(numel(dataset),numel(Nvalues),numel(taus),numel(alphas));
CNR_FIR = zeros(numel(dataset),1);
CNR2_FIR = zeros(numel(dataset),1);
PSNR_FIR = zeros(numel(dataset),1);
SNR_FIR = zeros(numel(dataset),1);
CNR_UNF = zeros(numel(dataset),1);
CNR2_UNF = zeros(numel(dataset),1);
SNR_UNF = zeros(numel(dataset),1);

for j = 1:numel(dataset)
    fprintf('******* Set %g ******\n',j);
    % read the dataset
    [seq, nMLAs] = readDataset(dataset{j});
    seq = stb(seq, nMLAs, nMLAs/2);
    
    % read the measurement window data
    [InX, InY, OutX, OutY] = readMeasurementWindow(dataset{j});
       
    
    [height,width,nFrames] = size(seq);    
    
%    warning('Extended mask! remove these lines for final results');
%     [OutY, OutX] = meshgrid(min(OutY):height,min(OutX):max(OutX));
%     OutY = OutY';
%     OutY = OutY(:);
%     OutX = OutX';
%     OutX = OutX(:);
    
    Mask = false(height,width);
    Mask(InY, InX) = true;
    Mask(OutY, OutX) = true;

%     seq = seq(:,:,1:18);

    startFrame = max(Nvalues-1)/2+1;
    stopFrame = nFrames - max(Nvalues-1)/2;
    
    % Measure Unfiltered
    seq_un = seq(:,:,startFrame:stopFrame);
    [CNR_UNF(j),CNR2_UNF(j),SNR_UNF(j)] = computeSNRandCNRsequence(abs(seq_un), InX, InY, OutX, OutY);
    fprintf('UNF CNR=%g\n',CNR_UNF(j));

    % Run FIR
    seq_fir = abs(seq(:,:,startFrame:stopFrame)-seq(:,:,[startFrame:stopFrame]-1));
    [CNR_FIR(j),CNR2_FIR(j),SNR_FIR(j)] = computeSNRandCNRsequence(seq_fir, InX, InY, OutX, OutY);
    [PSNR_FIR(j)] = computePSNRsequence(seq_fir,  abs(seq_un),  OutX, OutY);
    fprintf('FIR CNR=%g\n',CNR_FIR(j));  
    
    for k = 1:numel(Nvalues)
        N = Nvalues(k)-1;
        fprintf('------------------\nSize N=%d\n',N+1);
        
        % Run SVF
        CNR = zeros(numel(taus),numel(alphas));
        CNR2 = zeros(numel(taus),numel(alphas));
        PSNR3 = zeros(numel(taus),numel(alphas));
        SNR4 = zeros(numel(taus),numel(alphas));
        for l = 1:numel(taus)
%         parfor l = 1:numel(taus)
            tempCNR =zeros(numel(alphas),1);
            tempCNR2 =zeros(numel(alphas),1);
            tempPSNR3 =zeros(numel(alphas),1);
            tempSNR =zeros(numel(alphas),1);
            seq_svf = SVF(seq, Mvalue-1, N, alphas, taus(l), startFrame, stopFrame, Mask);
            for p = 1:numel(alphas)
                [tempCNR(p), tempCNR2(p),tempSNR(p)] = computeSNRandCNRsequence(abs(seq_svf{p}), InX, InY, OutX, OutY);
                tempPSNR3(p) = computePSNRsequence(abs(seq_svf{p}), abs(seq_un), OutX, OutY);
            end
            CNR(l,:) = tempCNR;
            CNR2(l,:) = tempCNR2;
            PSNR3(l,:) = tempPSNR3;
            SNR4(l,:) = tempSNR;
            fprintf('SVF tau=%g CNR=%g\n',taus(l),max(tempCNR));
        end
        PSNR_SVF(j,k,:,:) = PSNR3(:,:);
        CNR_SVF(j,k,:,:) = CNR(:,:);
        CNR2_SVF(j,k,:,:) = CNR2(:,:);
        SNR_SVF(j,k,:,:) = SNR4(:,:);
        [l,p] = find(CNR == max(CNR(:)));
        fprintf('SVF tau=%g alpha=%g CNR=%g\n',taus(l),alphas(p),max(CNR(:)));
        
%         continue;
        
        % Run MCA
        % use default values for: P, OL_P, DICT_SIZE (2:1), 
        params = [];
        params.N = Nvalues(k)-1;
        params.M = Mvalue-1;
        params.OL_M = params.M;
%         params.OL_M = params.M/2+1;
        params.parallel = true;
        params.partition = 0.9;
        params.exact = 0;
        params.ATOMS_KSVD = round(0.1*(params.N+1)*(params.M+1));
%         params.DICT_SIZE = 4*(params.N+1)*(params.M+1);
        params.DICT_SIZE = 2*(params.N+1)*(params.M+1);
        
        % Learn from data
        tic
        D = KSVDCDoubleFieldRealSequence(seq, params);
        DtD = D'*D;
        toc

        % Remove clutter
%         CNR = zeros(numel(OMP_THRs),numel(ATOM_THRs));
%         CNR2 = zeros(numel(OMP_THRs),numel(ATOM_THRs));
        params.ATOM_THRs = ATOM_THRs;
%         params.OMP_MAX_ATOMS = round(0.2*(params.N+1)*(params.M+1));
        params.OMP_MAX_ATOMS = min(40, round(0.2*(params.N+1)*(params.M+1)));
        params.OL_M = params.M;
        
        % Compute NSTD
%         patches = extractPatches(seq, startFrame,  params.N, 0, params.M, 0, params.OL_M);
%         G = comp2(D'*patches, sum(conj(patches).*patches), DtD, sqrt(2*(params.M+1)*(params.N+1)), 'maxatoms', params.OMP_MAX_ATOMS); %params.OMP_MAX_ATOMS);
%         pp = D*G;
%         NSTD = mean(sqrt(sum((abs(pp-patches)).^2)))/sqrt(2*(params.M+1)*(params.N+1))/2
        
%         NSTD = 450;
        
%         params.OL_M = params.M/2+1;
        params.OL_M = params.M;
%         l=1;
%         parfor l = 1:numel(OMP_THRs)
            tempCNR =zeros(numel(ATOM_THRs),1);
            tempCNR2 =zeros(numel(ATOM_THRs),1);
            tempPSNR3=zeros(numel(ATOM_THRs),1); 
            tempSNR4=zeros(numel(ATOM_THRs),1);
            locparams = params;
            locparams.OMP_ERROR = sqrt(2*(locparams.M+1)*(locparams.N+1)); %*NSTD*OMP_THRs(l);
            seq_mca = KSVDCField2SupressClutterSeq(D, DtD, seq, locparams, startFrame, stopFrame, Mask);
            for p = 1:numel(ATOM_THRs)
                [tempCNR(p), tempCNR2(p), tempSNR4(p)] = computeSNRandCNRsequence(abs(seq_mca{p}), InX, InY, OutX, OutY);
                tempPSNR3(p) = computePSNRsequence(abs(seq_mca{p}), abs(seq_un), OutX, OutY);
            end
%             fprintf('MCA errThr=%g CNR=%g\n',OMP_THRs(l),max(tempCNR));
%             CNR(l,:) = tempCNR;
%             CNR2(l,:) = tempCNR2;
%         end
%         CNR_MCA(j,k,:,:) = CNR(:,:);
%         CNR2_MCA(j,k,:,:) = CNR2(:,:);
        CNR = tempCNR;
        CNR2 = tempCNR2;
        PSNR3 = tempPSNR3;
        SNR4 = tempSNR4;
        [~,p] = max(CNR);
        fprintf('MCA THR=%g CNR=%g\n',ATOM_THRs(p),max(CNR(:)));
        CNR_MCA(j,k,:) = CNR(:);
        CNR2_MCA(j,k,:) = CNR2(:);
        PSNR_MCA(j,k,:) = PSNR3(:);
        SNR_MCA(j,k,:) = SNR4(:);
    end  
    save search_temp_snr.mat
end

end