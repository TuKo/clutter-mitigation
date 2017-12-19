%% load data
load ..\datasets\Simulated\simulatedData40.mat
artIQ = artifactIQdata;
IQdata = IQdata(:,:,1:20);
clear artifactIQdata;


%% send to server
Mvalues = 13:4:33;
Nvalues = [5,7,9];
OMP_ERRORs = 0.8:0.1:2.3;
cutoffs = 0.05:0.025:0.65;
ATOM_THRs = 0.35:0.025:0.95;
taus = 0.35:0.025:0.95;

c = parcluster('GIPCluster');

max_exper = 50;
procs = 24;
poolprocs = 2-1;
iterperproc = max_exper/procs;

j = cell(procs,1);

last=0;
for i = 1:procs
    first = last + 1;
    last = min(round(i*iterperproc),max_exper);
    fprintf('Sending proc %2d  f %3d l %3d\n',[i first last]);
    j{i} = batch(c,'figPatchSize',nargout(@figPatchSize),{IQdata, artIQ, first, last, Mvalues, Nvalues, OMP_ERRORs, cutoffs, ATOM_THRs, taus},'Pool',poolprocs,'AttachedFiles','ompbox/');
end

%% save the jobs info
jobsID = zeros(procs,1);
for i = 1:procs
    jobsID(i) = j{i}.ID;
end
save('results\Parameters-jobs.mat','jobsID','max_exper','procs','Mvalues','Nvalues','OMP_ERRORs','taus','ATOM_THRs','cutoffs');
diary(j{1});


%% merge the results

clear
load('results\Parameters-jobs.mat');

CNR_BMODE  = zeros(max_exper,1);
PSNR_BMODE = zeros(max_exper,1);

CNR_PF     = zeros(max_exper,1);
PSNR_PF    = zeros(max_exper,1);

CNR_FIR    = zeros(max_exper,1);
PSNR_FIR   = zeros(max_exper,1);

cluster = parcluster('GIPCluster');

for i = 1:procs
    job = cluster.findJob('ID',jobsID(i));
    r = fetchOutputs(job);
    first_exper = r{1};
    last_exper = r{2};
    fprintf('f %2d l %2d\n',first_exper,last_exper);
    if i==1
        [~,a,b,c] = size(r{3});
        CNR_GMCA  = zeros(max_exper,a,b,c);
        PSNR_GMCA = zeros(max_exper,a,b,c);
        [~,a,b,c,d] = size(r{9});
        CNR_AMCA  = zeros(max_exper,a,b,c,d);
        PSNR_AMCA = zeros(max_exper,a,b,c,d);
        [~,a,b,c,d] = size(r{11});
        CNR_FOMCA  = zeros(max_exper,a,b,c,d);
        PSNR_FOMCA = zeros(max_exper,a,b,c,d);
        [~,a,b,c] = size(r{14});
        CNR_SVF  = zeros(max_exper,a,b,c);
        PSNR_SVF = zeros(max_exper,a,b,c);
    end
    CNR_GMCA(first_exper:last_exper,:,:,:) = r{3}(first_exper:last_exper,:,:,:);
    PSNR_GMCA(first_exper:last_exper,:,:,:) = r{6}(first_exper:last_exper,:,:,:);
    CNR_AMCA(first_exper:last_exper,:,:,:,:) = r{9}(first_exper:last_exper,:,:,:,:);
    PSNR_AMCA(first_exper:last_exper,:,:,:,:) = r{10}(first_exper:last_exper,:,:,:,:);
    CNR_FOMCA(first_exper:last_exper,:,:,:,:) = r{11}(first_exper:last_exper,:,:,:,:);
    PSNR_FOMCA(first_exper:last_exper,:,:,:,:) = r{12}(first_exper:last_exper,:,:,:,:);
    CNR_SVF(first_exper:last_exper,:,:,:) = r{13}(first_exper:last_exper,:,:,:);
    PSNR_SVF(first_exper:last_exper,:,:,:) = r{14}(first_exper:last_exper,:,:,:);
    CNR_BMODE(first_exper:last_exper) = r{5}(first_exper:last_exper);
    PSNR_BMODE(first_exper:last_exper) = r{8}(first_exper:last_exper);
    CNR_PF(first_exper:last_exper) = r{4}(first_exper:last_exper);
    PSNR_PF(first_exper:last_exper) = r{7}(first_exper:last_exper);
    CNR_FIR(first_exper:last_exper) = r{15}(first_exper:last_exper);
    PSNR_FIR(first_exper:last_exper) = r{16}(first_exper:last_exper);
end

%save the results
clear cluster job r 
save('results\parameters.mat');

%% find the best parameters by CNR
load('results\parameters.mat');
fprintf('Best Parameters for CNR:\n');

mean_GMCA  = squeeze(mean(CNR_GMCA));
mean_AMCA  = squeeze(mean(CNR_AMCA));
mean_FOMCA = squeeze(mean(CNR_FOMCA));
mean_SVF   = squeeze(mean(CNR_SVF));
mean_PF = mean(CNR_PF);
mean_UN = mean(CNR_BMODE);

fprintf('PF: \tCNR=%.4g\n',mean_PF);
fprintf('UN: \tCNR=%.4g\n',mean_UN);

best_GMCA = find(abs(mean_GMCA - max(mean_GMCA(:))) < 1e-10);
[n,m,e] = ind2sub(size(mean_GMCA),best_GMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) mean_GMCA(best_GMCA)]
fprintf('GMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tCNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),mean_GMCA(n,m,e));

best_AMCA = find(abs(mean_AMCA - max(mean_AMCA(:))) < 1e-10);
[n,m,e,c] = ind2sub(size(mean_AMCA),best_AMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) cutoffs(c) mean_AMCA(best_AMCA)]
fprintf('AMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tCUT=%g, \tCNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),cutoffs(c),mean_AMCA(n,m,e,c));

best_FOMCA = find(abs(mean_FOMCA - max(mean_FOMCA(:))) < 1e-10);
[n,m,e,b] = ind2sub(size(mean_FOMCA),best_FOMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) ATOM_THRs(b) mean_FOMCA(best_FOMCA)]
fprintf('FOMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tBET=%g, \tCNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),ATOM_THRs(b), mean_FOMCA(n,m,e,b));

best_SVF = find(abs(mean_SVF - max(mean_SVF(:))) < 1e-10);
[n,m,t] = ind2sub(size(mean_SVF),best_SVF);
% [Nvalues(n) Mvalues(m) taus(t) mean_SVF(best_SVF)]
fprintf('SVF: \tM=%d, \tN=%d, \tTAU=%4.2g, \tCNR=%.4g\n',Mvalues(m),Nvalues(n),taus(t), mean_SVF(n,m,t));

%% find the best parameters by PSNR

fprintf('Best Parameters for PSNR:\n');

mean_GMCA  = squeeze(mean(PSNR_GMCA));
mean_AMCA  = squeeze(mean(PSNR_AMCA));
mean_FOMCA = squeeze(mean(PSNR_FOMCA));
mean_SVF   = squeeze(mean(PSNR_SVF));
mean_PF = mean(PSNR_PF);
mean_UN = mean(PSNR_BMODE);

fprintf('PF: \tPSNR=%.4g\n',mean_PF);
fprintf('UN: \tPSNR=%.4g\n',mean_UN);

best_GMCA = find(abs(mean_GMCA - max(mean_GMCA(:))) < 1e-10);
[n,m,e] = ind2sub(size(mean_GMCA),best_GMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) mean_GMCA(best_GMCA)]
fprintf('GMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),mean_GMCA(n,m,e));

best_AMCA = find(abs(mean_AMCA - max(mean_AMCA(:))) < 1e-10);
[n,m,e,c] = ind2sub(size(mean_AMCA),best_AMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) cutoffs(c) mean_AMCA(best_AMCA)]
fprintf('AMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tCUT=%g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),cutoffs(c),mean_AMCA(n,m,e,c));

best_FOMCA = find(abs(mean_FOMCA - max(mean_FOMCA(:))) < 1e-10);
[n,m,e,b] = ind2sub(size(mean_FOMCA),best_FOMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) ATOM_THRs(b) mean_FOMCA(best_FOMCA)]
fprintf('FOMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tBET=%g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),ATOM_THRs(b), mean_FOMCA(n,m,e,b));

best_SVF = find(abs(mean_SVF - max(mean_SVF(:))) < 1e-10);
[n,m,t] = ind2sub(size(mean_SVF),best_SVF);
% [Nvalues(n) Mvalues(m) taus(t) mean_SVF(best_SVF)]
fprintf('SVF: \tM=%d, \tN=%d, \tTAU=%4.2g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),taus(t), mean_SVF(n,m,t));

%% temp
fprintf('Comparison:\n');
mean_GMCA  = squeeze(mean(PSNR_GMCA));
mean_AMCA  = squeeze(mean(PSNR_AMCA));
mean_FOMCA = squeeze(mean(PSNR_FOMCA));
mean_SVF   = squeeze(mean(PSNR_SVF));

% mean_GMCA  = squeeze(mean(CNR_GMCA));
% mean_AMCA  = squeeze(mean(CNR_AMCA));
% mean_FOMCA = squeeze(mean(CNR_FOMCA));
% mean_SVF   = squeeze(mean(CNR_SVF));

n = 3; %N = 9
m = 2; %M = 17
m = 6; %M = 33
e = 16; %N = 2.3
% e = 3;

mean_GMCA(n,m,e)
max(squeeze(mean_AMCA(n,m,e,:)))
max(squeeze(mean_FOMCA(n,m,e,:)))
max(squeeze(mean_SVF(n,m,:)))

% max(squeeze(mean_GMCA(n,m,:)))
% figure, mesh(squeeze(mean_FOMCA(n,m,:,:)))
% figure, mesh( cutoffs, OMP_ERRORs, squeeze(mean_AMCA(n,m,:,:)))

% get best cutoff for these n,m values
mm_AMCA = squeeze(mean_AMCA(n,m,:,:));
best_AMCA = find(abs(mm_AMCA - max(mm_AMCA(:))) < 1e-10);
[e,c] = ind2sub(size(mm_AMCA),best_AMCA);
fprintf('AMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tCUT=%g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),cutoffs(c),mm_AMCA(e,c));

best_GMCA = find(abs(mean_GMCA - max(mean_GMCA(:))) < 1e-10);
[n,m,e] = ind2sub(size(mean_GMCA),best_GMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) mean_GMCA(best_GMCA)]
fprintf('GMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),mean_GMCA(n,m,e));

best_AMCA = find(abs(mean_AMCA - max(mean_AMCA(:))) < 1e-10);
[n,m,e,c] = ind2sub(size(mean_AMCA),best_AMCA);
% [Nvalues(n) Mvalues(m) OMP_ERRORs(e) cutoffs(c) mean_AMCA(best_AMCA)]
fprintf('AMCA: \tM=%d, \tN=%d, \tERR=%4.2g, \tCUT=%g, \tPSNR=%.4g\n',Mvalues(m),Nvalues(n),OMP_ERRORs(e),cutoffs(c),mean_AMCA(n,m,e,c));

%% plot figures:
% load('results\parameters.mat');
% folder = 'results';
% 
% best_m = 6;
% best_n = 3;
% 
% % CNR 
% figure;
% plot(ATOM_THRs,mean(CNR_PF)*ones(numel(ATOM_THRs),1), '--k','linewidth',2);
% hold on;
% plot(ATOM_THRs,mean(CNR_BMODE)*ones(numel(ATOM_THRs),1),':k','linewidth',2,'markersize',15);
% 
% [CNR, IDX1] = max(squeeze(mean(CNR_GMCA(:,best_n,best_m,:,:))),[],1);
% STD1 = squeeze(std(CNR_GMCA(:,best_n,best_m,:,:)));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(ATOM_THRs,CNR, STD1, 's-k','linewidth',2,'markersize',15);
% 
% xlim([min(ATOM_THRs)-0.025 max(ATOM_THRs)+0.025]);
% ylim([-2.5 6]);
% xlabel('Cut-off value \beta','FontSize',14);
% ylabel('CNR (dB)','FontSize',14);
% legend({'Perfect Filtering','Unfiltered','MCA'},'FontSize',14,'Location','SouthEast');
% set(gca,'FontSize',14);
% print([folder '\figureThresholdCNR.eps'],'-depsc2','-r600');
% print([folder '\figureThresholdCNR.png'],'-dpng');
% 
% %PSNR
% figure;
% % plot(ATOM_THRs,mean(PSNR_PF)*ones(numel(ATOM_THRs),1), '--k','linewidth',2);
% plot(ATOM_THRs,mean(PSNR_BMODE)*ones(numel(ATOM_THRs),1),':k','linewidth',2,'markersize',15);
% hold on;
% 
% [PSNR, IDX1] = max(squeeze(mean(PSNR_GMCA(:,best_n,best_m,:,:))),[],1);
% STD1 = squeeze(std(PSNR_GMCA(:,best_n,best_m,:,:)));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(ATOM_THRs,PSNR, STD1, 's-k','linewidth',2,'markersize',15);
% 
% xlim([min(ATOM_THRs)-0.025 max(ATOM_THRs)+0.025]);
% xlabel('Cut-off value \beta','FontSize',14);
% ylabel('PSNR (dB)','FontSize',14);
% legend({'Unfiltered','MCA'},'FontSize',14,'Location','SouthEast');
% set(gca,'FontSize',14);
% print([folder '\figureThresholdPSNR.eps'],'-depsc2','-r600');
% print([folder '\figureThresholdPSNR.png'],'-dpng');

