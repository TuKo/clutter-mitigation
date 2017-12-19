addpath('helpers/');
addpath('ompbox/');

% Normalizar las imagenes:

datasets = {'datasets_norm\set1\','datasets_norm\set2\','datasets_norm\set3\',...
            'datasets_norm\set4\','datasets_norm\set5\','datasets_norm\set6\',...
            'datasets_norm\set7\','datasets_norm\set8\','datasets_norm\set9\',...
            'datasets_norm\set10\','datasets_norm\set11\','datasets_norm\set12\',...
            'datasets_norm\set13\','datasets_norm\set14\','datasets_norm\set15\',...
};
folder = 'datasets_norm\';

seq_base_num = 2;
[seq_base, nMLAs] = readDataset(datasets{seq_base_num});
DR = 30;
GAIN = -40;
height = 400;
width = 400;

for j = 1:numel(datasets)
    j
    if j == seq_base_num
        writeVideo([folder 'set' num2str(j) '.avi'], convertMovie(scanConvertMovie(tissueProcessing(stb(seq_base,nMLAs,nMLAs/2),DR,GAIN),height,width)), 20);
        continue;
    end
    [seq, nMLAs] = readDataset(datasets{j});
    display(max(abs(seq_base(:))) / max(abs(seq(:))))
    seq = (max(abs(seq_base(:)))/ max(abs(seq(:))) )*seq;
    writeVideo([folder 'set' num2str(j) '-a.avi'], convertMovie(tissueProcessing(stb(seq,nMLAs,nMLAs/2),DR,GAIN)),20)
%     continue;
%     writeDataset([datasets{j}],seq,nMLAs);
    writeVideo([folder 'set' num2str(j) '.avi'], convertMovie(scanConvertMovie(tissueProcessing(stb(seq,nMLAs,nMLAs/2),DR,GAIN),height,width)), 20);
end



%%
%parpool('local');

% Nvalues = 5:2:19;
% Mvalue = 15;
% datasets = {'datasets\set1\','datasets\set2\','datasets\set3\',...
%             'datasets\set4\','datasets\set5\','datasets\set6\',...
%             'datasets\set7\','datasets\set8\','datasets\set9\'};
datasets = {'datasets_norm\set1\','datasets_norm\set2\','datasets_norm\set3\',...
            'datasets_norm\set4\','datasets_norm\set5\','datasets_norm\set6\',...
            'datasets_norm\set7\','datasets_norm\set8\','datasets_norm\set9\',...
            'datasets_norm\set10\','datasets_norm\set11\','datasets_norm\set12\',...
            'datasets_norm\set13\','datasets_norm\set14\','datasets_norm\set15\',...
            };

taus = 0.1:0.025:0.8;
alphas = 20:5:30;

OMPerrs = 1;
ATOMthrs = 0.05:0.025:0.95;

% datasets = {'datasets\set8\'};
% ATOMthrs = [0.35, 0.40, 0.45, 0.5];
% alphas = 10:5:50;
% alphas = [25 1000];
% taus = 0.4;
% taus = 0.35;
% OMPerrs = 0.05:0.05:1.5;
% ATOMthrs = 0.5;
% ATOMthrs = 0.45;


Nvalues = 15;
Mvalue = 15;
% alphas = 25; % Parametro SVF

[CNR_MCA, CNR2_MCA, CNR_SVF, CNR2_SVF, CNR_FIR, CNR2_FIR, CNR_UNF, CNR2_UNF,PSNR_MCA, PSNR_SVF, PSNR_FIR, SNR_UNF, SNR_MCA, SNR_SVF, SNR_FIR] = searchParametersRealdataBatch(datasets, Mvalue, Nvalues, taus, alphas, ATOMthrs, OMPerrs);


save results\realdata_results_norm_snr
return;

%% Figure to compare best results CNR from MCA and SVF
load results\realdata_results_norm

CNR_SVF = CNR_SVF(1:end-2,:,:,:);
CNR_MCA = CNR_MCA(1:end-2,:,:,:);
CNR_UNF = CNR_UNF(1:end-2,:,:,:);
PSNR_SVF = squeeze(PSNR_SVF(1:end-2,:,:,:));
PSNR_MCA = squeeze(PSNR_MCA(1:end-2,:,:,:));

CNR_SVF_ORIG = CNR_SVF;
CNR_SVF = squeeze(CNR_SVF);
best_alpha = 1;

CNR_MCA_ORIG = CNR_MCA;
CNR_MCA = squeeze(CNR_MCA);

% plot the CNR improvement graph as a function of the parameter
figure
plot(ATOMthrs, mean(squeeze(CNR_MCA) - repmat(CNR_UNF,[1 numel(ATOMthrs)])),'-b')
hold on , plot(taus, mean(squeeze(CNR_SVF(:,:,best_alpha))- repmat(CNR_UNF,[1 numel(taus)])),'-r');
% hold on , plot(taus, mean(max(CNR_SVF,[],3)- repmat(CNR_UNF,[1 numel(taus)])),'-g');
xlabel('Parameter value (\alpha and \beta)');
ylabel('CNR improvement over unfiltered image');

%%

PSNR_SVF_ORIG = PSNR_SVF;
PSNR_SVF = squeeze(PSNR_SVF);
best_alpha = 1;

PSNR_MCA_ORIG = PSNR_MCA;
PSNR_MCA = squeeze(PSNR_MCA);

% plot the CNR improvement graph as a function of the MSE (or "psnr")
%PSNR = 10*log10(MAX^2/MSE) with MAX = 10;
figure
plot(sqrt(10^2./10.^(mean(PSNR_MCA)/10)), mean(squeeze(CNR_MCA) - repmat(CNR_UNF,[1 numel(ATOMthrs)])),'-b')
hold on , plot(sqrt(10^2./10.^(mean(squeeze(PSNR_SVF(:,:,best_alpha)))/10)), mean(squeeze(CNR_SVF(:,:,best_alpha))- repmat(CNR_UNF,[1 numel(taus)])),'-r');
xlabel('MSE');
ylabel('CNR improvement over unfiltered image');

%%
best_tau = find(taus >= 0.325,1);
best_beta = find(ATOMthrs >= 0.45,1);
figure
% plot(sqrt(10^2./10.^((PSNR_MCA(:,best_beta))/10))-sqrt(10^2./10.^(squeeze(PSNR_SVF(:,best_tau,best_alpha))/10)), CNR_MCA(:,best_beta) -  CNR_SVF(:,best_tau,best_alpha),'+b')
plot(sqrt(10^2./10.^((PSNR_MCA(:,best_beta))/10)), CNR_MCA(:,best_beta) - CNR_UNF,'+b')
hold on , plot(sqrt(10^2./10.^(squeeze(PSNR_SVF(:,best_tau,best_alpha))/10)), CNR_SVF(:,best_tau,best_alpha)- CNR_UNF,'+r');
plot([0,max([ sqrt(10^2./10.^((PSNR_MCA(:,best_beta))/10));sqrt(10^2./10.^(squeeze(PSNR_SVF(:,best_tau,best_alpha))/10))])],[mean(CNR_MCA(:,best_beta) - CNR_UNF) mean(CNR_MCA(:,best_beta) - CNR_UNF)],'--b');
plot([0,max([ sqrt(10^2./10.^((PSNR_MCA(:,best_beta))/10));sqrt(10^2./10.^(squeeze(PSNR_SVF(:,best_tau,best_alpha))/10))])],[mean(CNR_SVF(:,best_tau,best_alpha)- CNR_UNF) mean(CNR_SVF(:,best_tau,best_alpha)- CNR_UNF)],'--r');
plot([mean(sqrt(10^2./10.^((PSNR_MCA(:,best_beta))/10))) mean(sqrt(10^2./10.^((PSNR_MCA(:,best_beta))/10)))],[0 max([(CNR_MCA(:,best_beta) - CNR_UNF);(CNR_SVF(:,best_tau,best_alpha)- CNR_UNF)])],'--b');
plot([mean(sqrt(10^2./10.^((PSNR_SVF(:,best_tau,best_alpha))/10))) mean(sqrt(10^2./10.^((PSNR_SVF(:,best_tau,best_alpha))/10)))],[0 max([(CNR_MCA(:,best_beta) - CNR_UNF);(CNR_SVF(:,best_tau,best_alpha)- CNR_UNF)])],'--r');
grid on
xlabel('MSE');
ylabel('CNR improvement over unfiltered image');
% xlabel('better MCA <--- MSE diff --> better SVF');
% ylabel('better SVF <--- CNR diff --> better MCA');

%% Find the correspondent value of tau based on PSNR/MSE for the current value of Beta
clear val1 val2
for i = 1:numel(ATOMthrs)
    % best_beta = find(ATOMthrs >= ATOMthrs(i),1);
    best_beta = i;
    % mean(PSNR_MCA(:,best_beta))
    best_tau = find(mean(PSNR_SVF(:,:,best_alpha)) >=mean(PSNR_MCA(:,best_beta)),1);
    if isempty(best_tau)
        best_tau = numel(taus);
    end
    val1(i) = mean(squeeze(CNR_SVF(:,best_tau,best_alpha))- CNR_UNF);
    val2(i) = mean(squeeze(CNR_MCA(:,best_beta))- CNR_UNF);
end
figure
plot(val1,val2, '+b')
% hold on, plot([0, max()
xlabel('SVF improvement');
ylabel('MCA improvement');

% figure
% mesh(ATOMthrs,Nvalues,squeeze(mean(CNR_MCA-repmat(CNR_UNF,[1,size(CNR_MCA,2),size(CNR_MCA,3)]))))
% xlabel('Cut-off value \beta','FontSize',14);
% ylabel('Number of frames','FontSize',14);
% zlabel('CNR (dB) improvement','FontSize',14);
% set(gca,'FontSize',14);
% % need to rotate the graph manually!
% print(['results\figureParametersRealData.eps'],'-depsc2','-r600');
% print(['results\figureParametersRealData.png'],'-dpng');


% figure
% mesh(taus,Nvalues,squeeze(mean(max(CNR_SVF,[],4)-repmat(CNR_UNF,[1,size(CNR_SVF,2),size(CNR_SVF,3)]))))

% figure
% mesh(ATOMthrs,Nvalues,squeeze(mean(CNR2_MCA-repmat(CNR2_UNF,[1,size(CNR2_MCA,2),size(CNR2_MCA,3)]))))

% figure
% mesh(taus,Nvalues,squeeze(mean(CNR2_SVF-repmat(CNR2_UNF,[1,size(CNR2_SVF,2),size(CNR2_SVF,3)]))))

%% Crear grafico con los resultados
load results\realdata_results_norm

imp_FIR = CNR_FIR - CNR_UNF;

% imp_MCA = squeeze(max(max(CNR_MCA,[],3),[],2)) - CNR_UNF;
% imp_SVF = squeeze(max(max(max(CNR_SVF,[],3),[],2),[],4)) - CNR_UNF;

best_thr = find(ATOMthrs>=0.44,1);
best_alpha = 1;
best_tau = find(taus>=0.35,1);
imp_MCA = CNR_MCA(:,:,best_thr) - CNR_UNF;
imp_SVF = squeeze(CNR_SVF(:,:,best_tau,best_alpha)) - CNR_UNF;
imp_MCA = squeeze(max(max(CNR_MCA(:,:,ATOMthrs>=0.3),[],3),[],2)) - CNR_UNF;
imp_SVF = squeeze(max(max(max(CNR_SVF(:,:,taus>=0.2,:),[],3),[],2),[],4)) - CNR_UNF;

% imp_MCA(14) = [];
% imp_SVF(14) = [];
% imp_FIR(14) = [];
imp_MCA = imp_MCA(1:end-2);
imp_SVF = imp_SVF(1:end-2);
imp_FIR = imp_FIR(1:end-2);

% imp_MCA = squeeze(squeeze(CNR_MCA(:,:,ATOMthrs==0.45))) - CNR_UNF;
% imp_SVF = squeeze(squeeze(max(CNR_SVF(:,:,taus==0.35,:),[],4))) - CNR_UNF;

folder = 'results\';

figure;
errorbar([mean(imp_MCA), mean(imp_FIR), mean(imp_SVF)],[std(imp_MCA), std(imp_FIR), std(imp_SVF)],'+k','linewidth',2,'markersize',15)
hold on
bar([mean(imp_MCA), mean(imp_FIR), mean(imp_SVF)],'k')
ylabel('CNR improvement (dB)','FontSize',20);
ylim([0 3]);
set(gca,'FontSize',20);
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'MCA','FIR','SVF'});
print([folder '\figureHeartComparison.eps'],'-depsc2','-r600');
print([folder '\figureHeartComparison.png'],'-dpng');


%% ejecutar en cada dataset con los parametros promedio y guardar cada video
%del output, el dataset limpio, y una imagen del antes y el despues de cada
%uno (imagen del medio?)

% datasets = {'datasets_norm\set1\','datasets_norm\set2\','datasets_norm\set3\',...
%             'datasets_norm\set4\','datasets_norm\set5\','datasets_norm\set6\',...
%             'datasets_norm\set7\','datasets_norm\set8\','datasets_norm\set9\',...
%             'datasets_norm\set10\','datasets_norm\set11\','datasets_norm\set12\',...
%             'datasets_norm\set13\'};

% datasets = {'datasets_norm\set14\','datasets_norm\set15\'};
% datasets = {'datasets_norm\set1\','datasets_norm\set2\'};
% folder = 'resultstemp\';

datasets = {'datasets_norm\set8\','datasets_norm\set8\','datasets_norm\set8\','datasets_norm\set9\','datasets_norm\set9\','datasets_norm\set15\'};
folder = 'results\hindawi\';

% folder = 'C:\Temp\results\';

%Best visual values for K-SVD (set1): 0.45 to 0.5
%Best visual values for SVF (set1): 0.4

%Best visual values for K-SVD (set8): 0.45
%Best visual values for SVF (set8): 0.35

% datasets = {'datasets\p6s2\'};
% folder = ['resultstemp\test6s2\'];
GAIN = [-45,-50,-55,-45,-50,-40];


thr = 0.45;
% folder = ['resultstemp\' num2str(thr) '\'];

% Parametros K-SVD
% ATOMthrs = 0.45;
ATOMthrs = thr;
N = 15;
M = 15;

% Parametros SVF
tau = 0.35;
% tau = thr;
alpha = 20;

% GAIN = [-55, -40, -40, -40, -20, -55, -55, -55];
% GAIN = [-40];

visualFrame = 1:10;
% visualFrame =[1, 3, 7, 10];
obtainRealdataOutput(datasets, folder, ATOMthrs, N, M, tau, alpha, visualFrame, GAIN)

date
