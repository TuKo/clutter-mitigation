%% load data
load ..\datasets\Simulated\simulatedData40.mat
artIQ = artifactIQdata;
IQdata = IQdata(:,:,1:20);
clear artifactIQdata;


%% send to server
Mvalue = 33;
Nvalue = 9;
OMP_ERROR = [2.3, 2.3, 2.3]; % GMCA, AMCA, FOMCA
cutoffs = 0.05:0.025:0.65;
ATOM_THRs = 0.35:0.025:0.95;
taus = 0.35:0.025:0.95;

shifts = (0:0.25:2)*80;

max_exper = 50;
procs = 9;
poolprocs = 5-1;
iterperproc = max_exper/procs;

j = cell(procs,1);
c = parcluster('GIPCluster');

last=0;
for i = 1:procs
    first = last + 1;
    last = min(round(i*iterperproc),max_exper);
    fprintf('Sending proc %2d  f %3d l %3d\n',[i first last]);
    j{i} = batch(c,'figMotion',nargout(@figMotion),{IQdata, artIQ, first, last, shifts, Mvalue, Nvalue, OMP_ERROR, cutoffs, taus, ATOM_THRs},'Pool',poolprocs,'AttachedFiles','ompbox/');
end

%% save the jobs info
jobsID = zeros(procs,1);
for i = 1:procs
    jobsID(i) = j{i}.ID;
end
save('results\Motion-jobs.mat','jobsID','max_exper','procs','Mvalue','Nvalue','OMP_ERROR','shifts','cutoffs', 'taus', 'ATOM_THRs');
diary(j{1});


%% merge the results

clear
load('results\Motion-jobs.mat');

nvalues = numel(shifts);

CNR_BMODE  = zeros(max_exper,nvalues);
PSNR_BMODE = zeros(max_exper,nvalues);

CNR_PF  = zeros(max_exper,nvalues);
PSNR_PF = zeros(max_exper,nvalues);

CNR_FIR   = zeros(max_exper,nvalues);
PSNR_FIR  = zeros(max_exper,nvalues);

CNR_GMCA  = zeros(max_exper,nvalues);
PSNR_GMCA = zeros(max_exper,nvalues);

CNR_GMCA2 = zeros(max_exper,nvalues);
PSNR_GMCA2= zeros(max_exper,nvalues);

cluster = parcluster('GIPCluster');

for i = 1:procs
    job = cluster.findJob('ID',jobsID(i));
    r = fetchOutputs(job);
    first_exper = r{1};
    last_exper = r{2};
    if i==1
        [~,a,b] = size(r{9});
        CNR_AMCA  = zeros(max_exper,a,b);
        PSNR_AMCA = zeros(max_exper,a,b);
        [~,a,b] = size(r{11});
        CNR_FOMCA  = zeros(max_exper,a,b);
        PSNR_FOMCA = zeros(max_exper,a,b);
        [~,a,b] = size(r{14});
        CNR_SVF  = zeros(max_exper,a,b);
        PSNR_SVF = zeros(max_exper,a,b);
    end
    CNR_GMCA(first_exper:last_exper,:)     = r{3}(first_exper:last_exper,:);
    PSNR_GMCA(first_exper:last_exper,:)    = r{6}(first_exper:last_exper,:);
    CNR_GMCA2(first_exper:last_exper,:)    = r{17}(first_exper:last_exper,:);
    PSNR_GMCA2(first_exper:last_exper,:)   = r{18}(first_exper:last_exper,:);
    CNR_AMCA(first_exper:last_exper,:,:)   = r{9}(first_exper:last_exper,:,:);
    PSNR_AMCA(first_exper:last_exper,:,:)  = r{10}(first_exper:last_exper,:,:);
    CNR_FOMCA(first_exper:last_exper,:,:)  = r{11}(first_exper:last_exper,:,:);
    PSNR_FOMCA(first_exper:last_exper,:,:) = r{12}(first_exper:last_exper,:,:);
    CNR_SVF(first_exper:last_exper,:,:)  = r{13}(first_exper:last_exper,:,:);
    PSNR_SVF(first_exper:last_exper,:,:) = r{14}(first_exper:last_exper,:,:);
    CNR_BMODE(first_exper:last_exper,:)  = r{5}(first_exper:last_exper,:);
    PSNR_BMODE(first_exper:last_exper,:) = r{8}(first_exper:last_exper,:);
    CNR_PF(first_exper:last_exper,:)     = r{4}(first_exper:last_exper,:);
    PSNR_PF(first_exper:last_exper,:)    = r{7}(first_exper:last_exper,:);
    CNR_FIR(first_exper:last_exper,:)    = r{15}(first_exper:last_exper,:);
    PSNR_FIR(first_exper:last_exper,:)   = r{16}(first_exper:last_exper,:);
end

%save the results
clear cluster job r 
save('results\motion.mat');

%% plot figure

close all;
clear
load('results\motion-corrected.mat');

mean_GMCA  = squeeze(mean(PSNR_GMCA));
mean_GMCA2 = squeeze(mean(PSNR_GMCA2));
mean_AMCA  = squeeze(mean(PSNR_AMCA));
mean_FOMCA = squeeze(mean(PSNR_FOMCA));
mean_SVF   = squeeze(mean(PSNR_SVF));
mean_PF    = squeeze(mean(PSNR_PF));
mean_UN    = squeeze(mean(PSNR_BMODE));
mean_FIR    = squeeze(mean(PSNR_FIR));

shifts = shifts/80;
folder = 'results';
fontsz = 20;

%PSNR
figure;
plot(shifts, mean_UN,':k','linewidth',2,'markersize',15);
hold on;
plot(shifts, mean_PF,'--k','linewidth',2,'markersize',15);

errorbar(shifts, mean_GMCA, std(PSNR_GMCA), 's-r','linewidth',2,'markersize',15);
errorbar(shifts, mean_GMCA2, std(PSNR_GMCA2), '^-r','linewidth',2,'markersize',15);

[PSNR, IDX1] = max(mean_AMCA,[],2);
STD1 = squeeze(std(PSNR_AMCA));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(shifts, PSNR, STD1, 'v-b','linewidth',2,'markersize',15);

[PSNR, IDX1] = max(mean_FOMCA,[],2);
STD1 = squeeze(std(PSNR_FOMCA));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(shifts, PSNR, STD1, 'o-g','linewidth',2,'markersize',15);


[PSNR, IDX1] = max(mean_SVF,[],2);
STD1 = squeeze(std(PSNR_SVF));
STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(shifts, PSNR, STD1,'x-c','linewidth',2,'markersize',15);

% errorbar(shifts, mean_FIR, std(PSNR_FIR), 'd-m','linewidth',2,'markersize',15);

xlim([min(shifts)-0.2 max(shifts)+0.2]);
ylim([30,55]);
xlabel('Tissue Axial Motion (Elements/Frame)','FontSize',fontsz);
ylabel('PSNR (dB)','FontSize',fontsz);
legend({'Unfiltered','Perfect Filt.','OFF-MCA','OFF-MCA impr.','TA-MCA','ON-MCA','SVF'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
print([folder '\figureMotion.eps'],'-depsc2','-r600');
print([folder '\figureMotion.png'],'-dpng');


