%% load data
load simulatedData40
artIQ = artifactIQdata;
IQdata = IQdata(:,:,1:20);
clear artifactIQdata;


%% send to server
Mvalue = 33;
Nvalue = 9;
taus = 0.35:0.025:0.95;
ATOM_THRs = 0.35:0.025:0.95;
OMP_ERROR = 2.3;

SNRs = 0:5:50;

% Mvalue = 5;
% Nvalue = 5;
% taus = 0.75:0.05:0.8;
% ATOM_THRs = 0.75:0.05:0.8;
% SNRs = [0 5];
% procs = 1;
% max_exper = 4;

max_exper = 100;
procs = 16;
poolprocs = 4-1;
iterperproc = max_exper/procs;
j = cell(procs,1);
c = parcluster('MINA');

last=0;
for i = 1:procs
    first = last + 1;
    last = min(round(i*iterperproc),max_exper);
    fprintf('Sending proc %2d  f %3d l %3d\n',[i first last]);
    j{i} = batch('figureNoiseBatch',nargout(@figureNoiseBatch),{IQdata, artIQ, first, last, SNRs, Mvalue, Nvalue,taus, OMP_ERROR, ATOM_THRs},'Pool',poolprocs,'AttachedFiles','ompbox/');
%     j{i} = batch('figureNoiseBatch',nargout(@figureNoiseBatch),{IQdata, artIQ, first, last, SNRs, Mvalue, Nvalue,taus, OMP_ERROR, ATOM_THRs},'AttachedFiles','ompbox/');
end

%% save the jobs info
jobsID = zeros(procs,1);
for i = 1:procs
    jobsID(i) = j{i}.ID;
end
save('results\Noise-jobs.mat','jobsID','max_exper','procs','Mvalue','Nvalue','taus','OMP_ERROR','ATOM_THRs','SNRs');
diary(j{1});


%% merge the results

clear
load('results\Noise-jobs.mat');

nSNRs = numel(SNRs);
CNR_FIR  = zeros(max_exper,nSNRs);
PSNR_FIR = zeros(max_exper,nSNRs);

CNR_BMODE  = zeros(max_exper,nSNRs);
PSNR_BMODE = zeros(max_exper,nSNRs);

CNR_PF  = zeros(max_exper,nSNRs);
PSNR_PF = zeros(max_exper,nSNRs);

cluster = parcluster('MINA');

for i = 1:procs
    job = cluster.findJob('ID',jobsID(i));
    r = fetchOutputs(job);
    first_exper = r{1};
    last_exper = r{2};
    if i==1
        [~,a,b] = size(r{3});
        CNR_KSVDC  = zeros(max_exper,a,b);
        CNR2_KSVDC = zeros(max_exper,a,b);
        PSNR_KSVDC = zeros(max_exper,a,b);
        MSE_KSVDC  = zeros(max_exper,a,b);
        PSNR2_KSVDC = zeros(max_exper,a,b);
        MSE2_KSVDC  = zeros(max_exper,a,b);
        
        [~,a,b] = size(r{9});
        CNR_SVF  = zeros(max_exper,a,b);
        CNR2_SVF = zeros(max_exper,a,b);
        PSNR_SVF = zeros(max_exper,a,b);
        MSE_SVF  = zeros(max_exper,a,b);
        PSNR2_SVF = zeros(max_exper,a,b);
        MSE2_SVF  = zeros(max_exper,a,b);
    end
    CNR_KSVDC(first_exper:last_exper,:,:) = r{3}(first_exper:last_exper,:,:);
    CNR2_KSVDC(first_exper:last_exper,:,:) = r{4}(first_exper:last_exper,:,:);
    PSNR_KSVDC(first_exper:last_exper,:,:) = r{5}(first_exper:last_exper,:,:);
    MSE_KSVDC(first_exper:last_exper,:,:) = r{6}(first_exper:last_exper,:,:);
    PSNR2_KSVDC(first_exper:last_exper,:,:) = r{7}(first_exper:last_exper,:,:);
    MSE2_KSVDC(first_exper:last_exper,:,:) = r{8}(first_exper:last_exper,:,:);
    CNR_SVF(first_exper:last_exper,:,:) = r{9}(first_exper:last_exper,:,:);
    CNR2_SVF(first_exper:last_exper,:,:) = r{10}(first_exper:last_exper,:,:);
    PSNR_SVF(first_exper:last_exper,:,:) = r{11}(first_exper:last_exper,:,:);
    MSE_SVF(first_exper:last_exper,:,:) = r{12}(first_exper:last_exper,:,:);
    PSNR2_SVF(first_exper:last_exper,:,:) = r{13}(first_exper:last_exper,:,:);
    MSE2_SVF(first_exper:last_exper,:,:) = r{14}(first_exper:last_exper,:,:);
    CNR_FIR(first_exper:last_exper,:) = r{15}(first_exper:last_exper,:);
    CNR2_FIR(first_exper:last_exper,:) = r{16}(first_exper:last_exper,:);
    PSNR_FIR(first_exper:last_exper,:) = r{17}(first_exper:last_exper,:);
    MSE_FIR(first_exper:last_exper,:) = r{18}(first_exper:last_exper,:);
    PSNR2_FIR(first_exper:last_exper,:) = r{19}(first_exper:last_exper,:);
    MSE2_FIR(first_exper:last_exper,:) = r{20}(first_exper:last_exper,:);
    CNR_BMODE(first_exper:last_exper,:) = r{21}(first_exper:last_exper,:);
    CNR2_BMODE(first_exper:last_exper,:) = r{22}(first_exper:last_exper,:);
    PSNR_BMODE(first_exper:last_exper,:) = r{23}(first_exper:last_exper,:);
    MSE_BMODE(first_exper:last_exper,:) = r{24}(first_exper:last_exper,:);
%     PSNR2_BMODE(first_exper:last_exper,:) = r{25}(first_exper:last_exper,:);
%     MSE2_BMODE(first_exper:last_exper,:) = r{26}(first_exper:last_exper,:);
    CNR_PF(first_exper:last_exper,:) = r{27}(first_exper:last_exper,:);
    CNR2_PF(first_exper:last_exper,:) = r{28}(first_exper:last_exper,:);
    PSNR_PF(first_exper:last_exper,:) = r{29}(first_exper:last_exper,:);
    MSE_PF(first_exper:last_exper,:) = r{30}(first_exper:last_exper,:);
end

%save the results
clear cluster job r 
save('results\noise.mat');

%% plot figures:
close all

load('results\noise.mat');
folder = 'results';
fontsz = 20;

CNR_PF = CNR_PF(:,2:end);
CNR_BMODE = CNR_BMODE(:,2:end);
CNR_FIR = CNR_FIR(:,2:end,:);
CNR_KSVDC = CNR_KSVDC(:,2:end,:,:);
CNR_SVF = CNR_SVF(:,2:end,:,:);
SNRs = SNRs(2:end);

PSNR_BMODE = PSNR_BMODE(:,2:end);
PSNR_FIR = PSNR_FIR(:,2:end,:);
PSNR_KSVDC = PSNR_KSVDC(:,2:end,:,:);
PSNR_SVF = PSNR_SVF(:,2:end,:,:);

% CNR 
figure;
% errorbar(SNRs, mean(CNR_PF), std(CNR_PF), '--k','linewidth',2);
plot(SNRs, mean(CNR_PF), '--k','linewidth',2);
hold on;
% errorbar(SNRs, mean(CNR_BMODE), std(CNR_BMODE),':k','linewidth',2,'markersize',15);
plot(SNRs, mean(CNR_BMODE),':k','linewidth',2,'markersize',15);
errorbar(SNRs, mean(CNR_FIR), std(CNR_FIR), 'o-k','linewidth',2,'markersize',15);

[CNR, IDX1] = max(squeeze(mean(CNR_KSVDC)),[],2);
STD1 = squeeze(std(CNR_KSVDC));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(SNRs, CNR, STD1, 's-k','linewidth',2,'markersize',15);

[CNR, IDX1] = max(squeeze(mean(CNR_SVF)),[],2);
STD1 = squeeze(std(CNR_SVF));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(SNRs, CNR, STD1,'d-k','linewidth',2,'markersize',15);

xlim([2 54]);
ylim([0 6]);
xlabel('Electronic SNR','FontSize',fontsz);
ylabel('CNR (dB)','FontSize',fontsz);
legend({'Perfect Filtering','Unfiltered','FIR','MCA','SVF'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
print([folder '\figureSNRCNR.eps'],'-depsc2','-r600');
print([folder '\figureSNRCNR.png'],'-dpng');

%PSNR
figure;
% errorbar(SNRs, mean(PSNR_BMODE), std(PSNR_BMODE),':k','linewidth',2,'markersize',15);
plot(SNRs, mean(PSNR_BMODE),':k','linewidth',2,'markersize',15);
hold on;
errorbar(SNRs, mean(PSNR_FIR), std(PSNR_FIR), 'o-k','linewidth',2,'markersize',15);

[PSNR, IDX1] = max(squeeze(mean(PSNR_KSVDC)),[],2);
STD1 = squeeze(std(PSNR_KSVDC));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(SNRs, PSNR, STD1, 's-k','linewidth',2,'markersize',15);

[PSNR, IDX1] = max(squeeze(mean(PSNR_SVF)),[],2);
STD1 = squeeze(std(PSNR_SVF));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(SNRs, PSNR, STD1,'d-k','linewidth',2,'markersize',15);

xlim([2 54]);
ylim([20 55]);
xlabel('Electronic SNR','FontSize',fontsz);
ylabel('PSNR (dB)','FontSize',fontsz);
legend({'Unfiltered','FIR','MCA','SVF'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
print([folder '\figureSNRPSNR.eps'],'-depsc2','-r600');
print([folder '\figureSNRPSNR.png'],'-dpng');
