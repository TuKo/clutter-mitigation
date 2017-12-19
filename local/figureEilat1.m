clear
close all
load('results\correlation.mat');

mean_GMCA  = squeeze(mean(PSNR_GMCA));
mean_GMCA2 = squeeze(mean(PSNR_GMCA2));
mean_AMCA  = squeeze(mean(PSNR_AMCA));
mean_FOMCA = squeeze(mean(PSNR_FOMCA));
mean_SVF   = squeeze(mean(PSNR_SVF));
mean_PF    = squeeze(mean(PSNR_PF));
mean_UN    = squeeze(mean(PSNR_BMODE));
mean_FIR    = squeeze(mean(PSNR_FIR));

folder = 'results';
fontsz = 20;
idx = 7;

%PSNR
figure;
plot(cutoffs,mean_UN(idx)*ones(numel(cutoffs),1),':k','linewidth',2,'markersize',15);
hold on;
plot(cutoffs,mean_PF(idx)*ones(numel(cutoffs),1),'--k','linewidth',2,'markersize',15);

% xs = 1:numel(rhos);
% errorbar(xs,mean_GMCA, std(PSNR_GMCA), 's-r','linewidth',2,'markersize',15);
% errorbar(xs,mean_GMCA2, std(PSNR_GMCA2), '^-r','linewidth',2,'markersize',15);
plot(cutoffs,mean_GMCA(idx)*ones(numel(cutoffs),1),'-r','linewidth',2,'markersize',15);
% plot(cutoffs,mean_GMCA2(idx)*ones(numel(cutoffs),1),'--r','linewidth',2,'markersize',15);

[PSNR] = mean_AMCA(idx,:);
STD1 = squeeze(std(squeeze(PSNR_AMCA(:,idx,:))));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
errorbar(cutoffs,PSNR, STD1, 'v-b','linewidth',2,'markersize',15);

% [PSNR, IDX1] = max(mean_FOMCA,[],2);
% STD1 = squeeze(std(PSNR_FOMCA));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(xs, PSNR, STD1, 'o-g','linewidth',2,'markersize',15);

% [PSNR, IDX1] = max(mean_SVF,[],2);
% STD1 = squeeze(std(PSNR_SVF));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(xs,PSNR, STD1,'x-c','linewidth',2,'markersize',15);

xlim([min(cutoffs) max(cutoffs)]);
ylim([30,55]);
xlabel('Threshold value \tau','FontSize',fontsz);
ylabel('PSNR (dB)','FontSize',fontsz);
legend({'Unfiltered','Perfect Filt.','OFF-MCA','TA-MCA'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
% set(gca,'XTick',cutoffs)
% set(gca,'XTickLabel',{cutoffs})
print([folder '\figureTauParameter.eps'],'-depsc2','-r600');
print([folder '\figureTauParameter.png'],'-dpng');

%% CNR
mean_GMCA  = squeeze(mean(CNR_GMCA));
mean_GMCA2 = squeeze(mean(CNR_GMCA2));
mean_AMCA  = squeeze(mean(CNR_AMCA));
mean_PF    = squeeze(mean(CNR_PF));
mean_UN    = squeeze(mean(CNR_BMODE));


figure;
plot(cutoffs,mean_UN(idx)*ones(numel(cutoffs),1),':k','linewidth',2,'markersize',15);
hold on;
plot(cutoffs,mean_PF(idx)*ones(numel(cutoffs),1),'--k','linewidth',2,'markersize',15);

% xs = 1:numel(rhos);
% errorbar(xs,mean_GMCA, std(PSNR_GMCA), 's-r','linewidth',2,'markersize',15);
% errorbar(xs,mean_GMCA2, std(PSNR_GMCA2), '^-r','linewidth',2,'markersize',15);
plot(cutoffs,mean_GMCA(idx)*ones(numel(cutoffs),1),'-r','linewidth',2,'markersize',15);
% plot(cutoffs,mean_GMCA2(idx)*ones(numel(cutoffs),1),'--r','linewidth',2,'markersize',15);

[CNR] = mean_AMCA(idx,:);
STD1 = squeeze(std(squeeze(CNR_AMCA(:,idx,:))));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(cutoffs,CNR, STD1, 'v-b','linewidth',2,'markersize',15);
plot(cutoffs,CNR, 'v-b','linewidth',2,'markersize',15);

% [PSNR, IDX1] = max(mean_FOMCA,[],2);
% STD1 = squeeze(std(PSNR_FOMCA));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(xs, PSNR, STD1, 'o-g','linewidth',2,'markersize',15);

% [PSNR, IDX1] = max(mean_SVF,[],2);
% STD1 = squeeze(std(PSNR_SVF));
% STD1 = STD1(sub2ind(size(STD1),ones(size(IDX1)),IDX1));
% errorbar(xs,PSNR, STD1,'x-c','linewidth',2,'markersize',15);

xlim([min(cutoffs) max(cutoffs)]);
ylim([4.74,4.78]);
xlabel('Threshold value \tau','FontSize',fontsz);
ylabel('CNR (dB)','FontSize',fontsz);
legend({'Unfiltered','Perfect Filt.','OFF-MCA','TA-MCA'},'FontSize',fontsz,'Location','SouthEast');
set(gca,'FontSize',fontsz);
% set(gca,'XTick',cutoffs)
% set(gca,'XTickLabel',{cutoffs})
print([folder '\figureTauParameterCNR.eps'],'-depsc2','-r600');
print([folder '\figureTauParameterCNR.png'],'-dpng');