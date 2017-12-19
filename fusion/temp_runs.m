c = parcluster('GIPCluster');
test = 'dual1';
% test = 'dual2';

[I1clutter, nMLAs] = readDataset('data\clutter1-1\');
I1clutter = stb(I1clutter, nMLAs, nMLAs/2);
[I2clutter, nMLAs] = readDataset('data\clutter2-2\');
I2clutter = stb(I2clutter, nMLAs, nMLAs/2);
% [I2tissue, nMLAs] = readDataset('data\tissue1-2\');
[I1, I2, nMLAs] = readDatasetDualHarm(['data\' test '\']);
I1 = stb(I1, nMLAs, nMLAs/2);
I2 = stb(I2, nMLAs, nMLAs/2);


D1C = [];
D2C = [];
D = [];
DT = [];


%%

% j1 = batch(c,@experiment1,nargout(@experiment1),{I1,I2,I1clutter,I2clutter},'Pool',7,'AttachedFiles',{'ompbox\'});
j2 = batch(c,@experiment2,nargout(@experiment2),{I1,I2,I1clutter,I2clutter, D1C, D2C, D},'Pool',7,'AttachedFiles',{'ompbox\'});
j3 = batch(c,@experiment3,nargout(@experiment3),{I1,I2,I1clutter,I2clutter},'Pool',7,'AttachedFiles',{'ompbox\'});
% j4 = batch(c,@experiment4,nargout(@experiment4),{I1,I2,I1clutter,I2clutter},'Pool',7,'AttachedFiles',{'ompbox\'});
% j5 = batch(c,@experiment5,nargout(@experiment5),{I1,I2,I1clutter,I2clutter, D1C, D2C, D, DT},'Pool',7,'AttachedFiles',{'ompbox\'});


%%
THRs = [0.05:0.05:1];
j22 = cell(numel(THRs),1);
i=1;
j22{i} = batch(c,@experiment2thr,nargout(@experiment2thr),{I1,I2,I1clutter,I2clutter, D1C, D2C, D, THRs(i)},'Pool',7,'AttachedFiles',{'ompbox\'});
wait(j22{i});
rr = fetchOutputs(j22{i});
D1C = rr{4};
D2C = rr{5};
D = rr{6};
for i = 2:numel(THRs)
    j22{i} = batch(c,@experiment2thr,nargout(@experiment2thr),{I1,I2,I1clutter,I2clutter, D1C, D2C, D, THRs(i)},'Pool',7,'AttachedFiles',{'ompbox\'});
end
%%

% r = fetchOutputs(j1);
% r = fetchOutputs(j2);
% r = fetchOutputs(j3);
% r = fetchOutputs(j4);
r = fetchOutputs(j5);

%%
I1N = (max(abs(I2(:)))/ max(abs(I1(:))) )*I1;
GAIN1 = -40;
GAIN2 = GAIN1+20;
M = 15;
N = 15;
H = 400;
W = 400;
[height,width,nFrames] = size(I1);
startFrame = (N-1)/2+1;
stopFrame = nFrames - (N-1)/2;
I1N = I1N(:,:,startFrame:stopFrame);
I2N = I2(:,:,startFrame:stopFrame);

%%
k = 2;

k = 3;

k = 4;


k=5;


%%
out1 = r{1};
out2 = r{2};
D1C = r{3};
D2C = r{4};
D = r{5};
% DC = r{3};
% D = r{4};
% AT = r{5};
tissue1 = r{6};
tissue2 = r{7};

% out1 = out1(:,:,startFrame:stopFrame);
% out2 = out2(:,:,startFrame:stopFrame);
%% with stb
% mixBefore = scanConvertMovie(tissueProcessing(0.5*stb(I1N,nMLAs,2)+0.5*stb(I2N,nMLAs,2),30,GAIN2),H,W);
% mixAfter  = scanConvertMovie(tissueProcessing(0.5*stb(out1,nMLAs,2)+0.5*stb(out2,nMLAs,2),30,GAIN2),H,W);
% originals = [scanConvertMovie(tissueProcessing(stb(I1N,nMLAs,2),30,GAIN2),H,W), scanConvertMovie(tissueProcessing(stb(I2N,nMLAs,2),30,GAIN2),H,W)];
% % tissueonly = [scanConvertMovie(tissueProcessing(stb(tissue1,nMLAs,2),30,GAIN2),H,W), scanConvertMovie(tissueProcessing(stb(tissue2,nMLAs,2),30,GAIN2),H,W)];
% results   = [scanConvertMovie(tissueProcessing(stb(out1,nMLAs,2),30,GAIN2),H,W), scanConvertMovie(tissueProcessing(stb(out2,nMLAs,2),30,GAIN2),H,W)];
% allmix    = [originals, mixBefore];
% allresult = [results,   mixAfter];
% 
% writeVideo(['results\exper' num2str(k) '-mix.avi'], convertMovie([mixBefore, mixAfter]), 20);
% % writeVideo(['results\exper' num2str(k) '-tissueOnly.avi'], convertMovie([originals; tissueonly]), 20);
% writeVideo(['results\exper' num2str(k) '-results.avi'], convertMovie([originals; results]), 20);
% writeVideo(['results\exper' num2str(k) '-all.avi'], convertMovie([allmix; allresult]), 20);
% % writeVideo(['results\exper' num2str(k) '-miximage-mixpatches.avi'], convertMovie([mixBefore mixAfter scanConvertMovie(tissueProcessing(stb(tissue2,nMLAs,2),30,GAIN2),H,W)]), 20);

%% without stb
mixBefore = scanConvertMovie(tissueProcessing(0.5*I1N+0.5*I2N,30,GAIN2),H,W);
mixAfter  = scanConvertMovie(tissueProcessing(0.5*out1+0.5*out2,30,GAIN2),H,W);
% mixBefore = scanConvertMovie(tissueProcessing(0.5*abs(I1N)+0.5*abs(I2N),30,GAIN2),H,W);
% mixAfter  = scanConvertMovie(tissueProcessing(0.5*abs(out1)+0.5*abs(out2),30,GAIN2),H,W);
originals = [scanConvertMovie(tissueProcessing(I1N,30,GAIN2),H,W), scanConvertMovie(tissueProcessing(I2N,30,GAIN2),H,W)];
% tissueonly = [scanConvertMovie(tissueProcessing((tissue1),30,GAIN2),H,W), scanConvertMovie(tissueProcessing((tissue2),30,GAIN2),H,W)];
results   = [scanConvertMovie(tissueProcessing(out1,30,GAIN2),H,W), scanConvertMovie(tissueProcessing(out2,30,GAIN2),H,W)];
allmix    = [originals, mixBefore];
allresult = [results,   mixAfter];

writeVideo(['results\' test '-exper' num2str(k) '-mix-nostb.avi'], convertMovie([mixBefore, mixAfter]), 20);
% writeVideo(['results\' test '-exper' num2str(k) '-tissueOnly-nostb.avi'], convertMovie([originals; tissueonly]), 20);
writeVideo(['results\' test '-exper' num2str(k) '-results-nostb.avi'], convertMovie([originals; results]), 20);
writeVideo(['results\' test '-exper' num2str(k) '-all-nostb.avi'], convertMovie([allmix; allresult]), 20);
% writeVideo(['results\' test '-exper' num2str(k) '-miximage-mixpatches-nostb.avi'], convertMovie([mixBefore mixAfter scanConvertMovie(tissueProcessing((tissue2),30,GAIN2),H,W)]), 20);
% writeVideo(['results\' test '-exper' num2str(k) '-onlymix-nostb.avi'], convertMovie([ mixBefore results]), 20);


%%

out1 = cell(numel(THRs),1);
out2 = cell(numel(THRs),1);
for i = 1:numel(THRs)
    r = fetchOutputs(j22{i});
    out1{i} = r{1};
    out2{i} = r{2};
    if i==1
        D1C = r{4};
        D2C = r{5};
        D = r{6};
    end
end
k = 22;
%%
% D = DT;

mm = size(DC,2);
D1C = DC(1:(N*M),1:mm/2);
D2C = DC((N*M+1):end,(mm/2+1):end);

m = size(D,2);
D1T = D(1:(N*M),1:m/2);
D2T = D((N*M+1):end,(m/2+1):end);



AT1 = comp(D1C'*D1T, D1C'*D1C, M);
AT2 = comp(D2C'*D2T, D2C'*D2C, M);

AT = comp(DC'*D, DC'*DC, M);
ATnn = comp(DC'*D, DC'*DC, 2*M);

DDT1 = D(:,~(sum(abs(D - DC*AT).^2)<0.2));
DDC1 = D(:,(sum(abs(D - DC*AT).^2)<0.2));

DDT1 = D(:,~(sum(abs(D - DC*ATnn).^2)<0.1));
DDC1 = D(:,(sum(abs(D - DC*ATnn).^2)<0.1));

DDT1 = D(:,~(sum(abs(D1T - D1C*AT1).^2)<0.2 & sum(abs(D2T - D2C*AT2).^2)<0.2));
DDC1 = D(:,(sum(abs(D1T - D1C*AT1).^2)<0.2 & sum(abs(D2T - D2C*AT2).^2)<0.2));


% queda como tissue
figure;
showdict(abs(DDT1(1:M*N,:)),[M,N],floor(sqrt(size(DDT1,2))), ceil(sqrt(size(DDT1,2))),'highcontrast','whitelines');
figure;
showdict(abs(DDT1(M*N+(1:M*N),:)),[M,N],floor(sqrt(size(DDT1,2))), ceil(sqrt(size(DDT1,2))),'highcontrast','whitelines');

% eliminado x parecido a clutter
figure;
showdict(abs(DDC1(1:M*N,:)),[M,N],floor(sqrt(size(DDC1,2))), ceil(sqrt(size(DDC1,2))),'highcontrast','whitelines');
figure;
showdict(abs(DDC1(M*N+(1:M*N),:)),[M,N],floor(sqrt(size(DDC1,2))), ceil(sqrt(size(DDC1,2))),'highcontrast','whitelines');





AT = comp(D2C'*D, D2C'*D2C, M);
DDT2 = D(:,~(sum(abs(D - D2C*AT).^2)<0.2));
DDC2 = D(:,(sum(abs(D - D2C*AT).^2)<0.2));

m = size(DC,2);
figure, showdict(abs(DC(1:M*N,1:m/2)),[M,N],floor(sqrt(m/2)), floor(sqrt(m/2)),'highcontrast','whitelines');
figure, showdict(abs(DC(M*N+(1:M*N),m/2+(1:m/2))),[M,N],floor(sqrt(m/2)), floor(sqrt(m/2)),'highcontrast','whitelines');


figure;
showdict(abs(DDT2),[M,N],floor(sqrt(size(DDT2,2))), ceil(sqrt(size(DDT2,2))),'highcontrast','whitelines');
figure;
showdict(abs(DDC2),[M,N],floor(sqrt(size(DDC2,2))), ceil(sqrt(size(DDC2,2))),'highcontrast','whitelines');


%%



