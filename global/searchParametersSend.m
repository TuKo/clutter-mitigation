%% load data
load simulatedData40
artIQ = artifactIQdata;
IQdata = IQdata(:,:,1:20);
clear artifactIQdata;


%% send to server
Mvalues = 13:4:33;
Nvalues = [5,7,9];
taus = 0.35:0.025:0.95;
OMP_ERRORs = 0.8:0.1:2.3;
ATOM_THRs = 0.35:0.025:0.95;

c = parcluster('MINA');

max_exper = 100;
procs = 16;
poolprocs = 4-1;
iterperproc = max_exper/procs;
j = cell(procs,1);

last=0;
for i = 1:procs
    first = last + 1;
    last = min(round(i*iterperproc),max_exper);
    fprintf('Sending proc %2d  f %3d l %3d\n',[i first last]);
    j{i} = batch('searchParametersBatch',nargout(@searchParametersBatch),{IQdata, artIQ, first, last, Mvalues, Nvalues,taus, OMP_ERRORs, ATOM_THRs},'Pool',poolprocs,'AttachedFiles','ompbox/');
end

%% save the jobs info
jobsID = zeros(procs,1);
for i = 1:procs
    jobsID(i) = j{i}.ID;
end
save('results\Parameters-jobs.mat','jobsID','max_exper','procs','Mvalues','Nvalues','taus','OMP_ERRORs','ATOM_THRs');
diary(j{1});


%% merge the results

clear
load('results\Parameters-jobs.mat');

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

cluster = parcluster('MINA');

for i = 1:procs
    job = cluster.findJob('ID',jobsID(i));
    r = fetchOutputs(job);
    first_exper = r{1};
    last_exper = r{2};
    if i==1
        [~,a,b,c,d] = size(r{3});
        CNR_KSVDC  = zeros(max_exper,a,b,c,d);
        CNR2_KSVDC = zeros(max_exper,a,b,c,d);
        PSNR_KSVDC = zeros(max_exper,a,b,c,d);
        MSE_KSVDC  = zeros(max_exper,a,b,c,d);
        PSNR2_KSVDC = zeros(max_exper,a,b,c,d);
        MSE2_KSVDC  = zeros(max_exper,a,b,c,d);
        
        [~,a,b,c] = size(r{9});
        CNR_SVF  = zeros(max_exper,a,b,c);
        CNR2_SVF = zeros(max_exper,a,b,c);
        PSNR_SVF = zeros(max_exper,a,b,c);
        MSE_SVF  = zeros(max_exper,a,b,c);
        PSNR2_SVF = zeros(max_exper,a,b,c);
        MSE2_SVF  = zeros(max_exper,a,b,c);
    end
    CNR_KSVDC(first_exper:last_exper,:,:,:,:) = r{3}(first_exper:last_exper,:,:,:,:);
    CNR2_KSVDC(first_exper:last_exper,:,:,:,:) = r{4}(first_exper:last_exper,:,:,:,:);
    PSNR_KSVDC(first_exper:last_exper,:,:,:,:) = r{5}(first_exper:last_exper,:,:,:,:);
    MSE_KSVDC(first_exper:last_exper,:,:,:,:) = r{6}(first_exper:last_exper,:,:,:,:);
    PSNR2_KSVDC(first_exper:last_exper,:,:,:,:) = r{7}(first_exper:last_exper,:,:,:,:);
    MSE2_KSVDC(first_exper:last_exper,:,:,:,:) = r{8}(first_exper:last_exper,:,:,:,:);
    CNR_SVF(first_exper:last_exper,:,:,:) = r{9}(first_exper:last_exper,:,:,:);
    CNR2_SVF(first_exper:last_exper,:,:,:) = r{10}(first_exper:last_exper,:,:,:);
    PSNR_SVF(first_exper:last_exper,:,:,:) = r{11}(first_exper:last_exper,:,:,:);
    MSE_SVF(first_exper:last_exper,:,:,:) = r{12}(first_exper:last_exper,:,:,:);
    PSNR2_SVF(first_exper:last_exper,:,:,:) = r{13}(first_exper:last_exper,:,:,:);
    MSE2_SVF(first_exper:last_exper,:,:,:) = r{14}(first_exper:last_exper,:,:,:);
    CNR_FIR(first_exper:last_exper) = r{15}(first_exper:last_exper);
    CNR2_FIR(first_exper:last_exper) = r{16}(first_exper:last_exper);
    PSNR_FIR(first_exper:last_exper) = r{17}(first_exper:last_exper);
    MSE_FIR(first_exper:last_exper) = r{18}(first_exper:last_exper);
    PSNR2_FIR(first_exper:last_exper) = r{19}(first_exper:last_exper);
    MSE2_FIR(first_exper:last_exper) = r{20}(first_exper:last_exper);
    CNR_BMODE(first_exper:last_exper) = r{21}(first_exper:last_exper);
    CNR2_BMODE(first_exper:last_exper) = r{22}(first_exper:last_exper);
    PSNR_BMODE(first_exper:last_exper) = r{23}(first_exper:last_exper);
    MSE_BMODE(first_exper:last_exper) = r{24}(first_exper:last_exper);
%     PSNR2_BMODE(first_exper:last_exper) = r{25}(first_exper:last_exper);
%     MSE2_BMODE(first_exper:last_exper) = r{26}(first_exper:last_exper);
    CNR_PF(first_exper:last_exper) = r{27}(first_exper:last_exper);
    CNR2_PF(first_exper:last_exper) = r{28}(first_exper:last_exper);
    PSNR_PF(first_exper:last_exper) = r{29}(first_exper:last_exper);
    MSE_PF(first_exper:last_exper) = r{30}(first_exper:last_exper);
end

%save the results
clear cluster job r 
save('results\parameters.mat');

%% find the best parameters

mean_KSVDC = squeeze(mean(CNR_KSVDC));
mean_SVF   = squeeze(mean(CNR_SVF));
mean_PF    = mean(CNR_PF);
mean_BMODE = mean(CNR_BMODE);
mean_FIR   = mean(CNR_FIR);

best_KSVD = find(abs(mean_KSVDC - max(mean_KSVDC(:))) < 1e-10);
[n,m,x,y] = ind2sub(size(mean_KSVDC),best_KSVD);
fprintf('Best parameters K-SVDC: M=%d, N=%d, ERR=%g, BETA=%g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(x),ATOM_THRs(y));

best_SVF = find(abs(mean_SVF - max(mean_SVF(:))) < 1e-10);
[rr,ss,tt] = ind2sub(size(mean_SVF),best_SVF);
fprintf('Best parameters SVF: M=%d, N=%d, TAU=%g\n',Mvalues(ss),Nvalues(rr),taus(tt));

%% plot figures:
clear
close all
load('results\parameters.mat');
folder = 'results';
fontsz = 20;

best_m = 6;
best_n = 3;

% CNR 
figure;
plot(ATOM_THRs,mean(CNR_PF)*ones(numel(ATOM_THRs),1), '--k','linewidth',2);
hold on;
plot(ATOM_THRs,mean(CNR_BMODE)*ones(numel(ATOM_THRs),1),':k','linewidth',2,'markersize',15);

[CNR, IDX1] = max(squeeze(mean(CNR_KSVDC(:,best_n,best_m,:,:))),[],1);
STD1 = squeeze(std(CNR_KSVDC(:,best_n,best_m,:,:)));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(ATOM_THRs,CNR, STD1, 's-k','linewidth',2,'markersize',15);
errorbar(ATOM_THRs,CNR, STD1, 'o-g','linewidth',2,'markersize',15);

xlim([min(ATOM_THRs)-0.025 max(ATOM_THRs)+0.025]);
ylim([-2.5 6]);
xlabel('Cut-off value \beta','FontSize',fontsz);
ylabel('CNR (dB)','FontSize',fontsz);
legend({'Perfect Filtering','Unfiltered','MCA'},'FontSize',fontsz,'Location','SouthEast');
legend({'Perfect Filt.','Unfiltered','ON-MCA'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
print([folder '\figureThresholdCNR.eps'],'-depsc2','-r600');
print([folder '\figureThresholdCNR.png'],'-dpng');

%PSNR
figure;
% plot(ATOM_THRs,mean(PSNR_PF)*ones(numel(ATOM_THRs),1), '--k','linewidth',2);
plot(ATOM_THRs,mean(PSNR_BMODE)*ones(numel(ATOM_THRs),1),':k','linewidth',2,'markersize',15);
hold on;

[PSNR, IDX1] = max(squeeze(mean(PSNR_KSVDC(:,best_n,best_m,:,:))),[],1);
STD1 = squeeze(std(PSNR_KSVDC(:,best_n,best_m,:,:)));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(ATOM_THRs,PSNR, STD1, 's-k','linewidth',2,'markersize',15);
errorbar(ATOM_THRs,PSNR, STD1, 'o-g','linewidth',2,'markersize',15);

xlim([min(ATOM_THRs)-0.025 max(ATOM_THRs)+0.025]);
xlabel('Cut-off value \beta','FontSize',fontsz);
ylabel('PSNR (dB)','FontSize',fontsz);
legend({'Unfiltered','MCA'},'FontSize',fontsz,'Location','SouthEast');
legend({'Unfiltered','ON-MCA'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
print([folder '\figureThresholdPSNR.eps'],'-depsc2','-r600');
print([folder '\figureThresholdPSNR.png'],'-dpng');

