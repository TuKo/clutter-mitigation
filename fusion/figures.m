%% This script creates the tables and figures for the ICASSP paper
clc

table = zeros(9,4);
table2 = zeros(8,4);
tableCNR = zeros(11,2);
%% Load TAMCA results for comparison

load('results\mat-nosync\TAMCA-dual2.mat');
out1TAMCA = out1;
out2TAMCA = out2;

%% Work on Dual 2
load('results\mat-nosync\dual2-exper2-results.mat')


%% Setup of input signals + other output parameters
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

%% Region 1 - ~ventriculo grande
fprintf('\nRegion 1\n');
% RegX = [1 25]; % ~10% improvement
% RegY = [290 290 + 45];
RegX = [1 25]; % ~12% improvement 
RegY = [290 290 + 30]+10;
% RegX = [5 5+20]; % ~6% improvement
% RegY = [437 437+50];

% FUSION 
[snr,ssnr, st, sd] = SpeckleSNR(abs(I1N),RegX,RegY);
table(1,1) = st; table2(1,1) = sd;
fprintf('I1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(I2N),RegX,RegY);
table(2,1) = st; table2(2,1) = sd;
fprintf('I2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY);
fprintf('Im: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1),RegX,RegY);
table(3,1) = st; table2(3,1) = sd;
fprintf('O1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2),RegX,RegY);
table(4,1) = st; table2(4,1) = sd;
fprintf('O2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1+out2)/2),RegX,RegY);
table(5,1) = st; table2(5,1) = sd;
fprintf('Om: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
% TAMCA
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1TAMCA),RegX,RegY);
table(6,1) = st; table2(6,1) = sd;
fprintf('O1 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2TAMCA),RegX,RegY);
table(7,1) = st; table2(7,1) = sd;
fprintf('O2 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1TAMCA+out2TAMCA)/2),RegX,RegY);
table(8,1) = st; table2(8,1) = sd;
fprintf('Om TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);

% medidas sobre el promedio de los envelopes
% [snr,ssnr] = SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY);
% [snr,ssnr] = SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY);
fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs((out1+out2)/2),RegX,RegY)/SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY));
% fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY)/SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY));
fprintf('#Pixels: %d\n', (RegX(2)-RegX(1))*(RegY(2)-RegY(1)));
table(end,1) = (RegX(2)-RegX(1))*(RegY(2)-RegY(1));

%% Region 2 - ventriculo chico
fprintf('\nRegion 2\n');
% RegX = [201 220]-124; % ~8% improvement
% RegY = [313 361];
% RegX = [200 210+24]-124; % ~23% impr with high std!
% RegY = [335 335+45];
RegX = [72 72+25]; % ~7% impr. 
RegY = [330 330+45];

[snr,ssnr, st, sd] = SpeckleSNR(abs(I1N),RegX,RegY);
table(1,2) = st; table2(1,2) = sd;
fprintf('I1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(I2N),RegX,RegY);
table(2,2) = st; table2(2,2) = sd;
fprintf('I2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY);
fprintf('Im: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1),RegX,RegY);
table(3,2) = st; table2(3,2) = sd;
fprintf('O1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2),RegX,RegY);
table(4,2) = st; table2(4,2) = sd;
fprintf('O2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1+out2)/2),RegX,RegY);
table(5,2) = st; table2(5,2) = sd;
fprintf('Om: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
% TAMCA
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1TAMCA),RegX,RegY);
table(6,2) = st; table2(6,2) = sd;
fprintf('O1 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2TAMCA),RegX,RegY);
table(7,2) = st; table2(7,2) = sd;
fprintf('O2 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1TAMCA+out2TAMCA)/2),RegX,RegY);
table(8,2) = st; table2(8,2) = sd;
fprintf('Om TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);

% medidas sobre el promedio de los envelopes
% [snr,ssnr] = SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY);
% [snr,ssnr] = SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY);
fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs((out1+out2)/2),RegX,RegY)/SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY));
% fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY)/SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY));
fprintf('#Pixels: %d\n', (RegX(2)-RegX(1))*(RegY(2)-RegY(1)));
% return;
table(end,2) = (RegX(2)-RegX(1))*(RegY(2)-RegY(1));

%% CNR Measurements
fprintf('\nCNR\n');
% InX = [100 124];
% InY = [297 297+44];

InX =[229 229+16]-124;
InY = [248 248+29];

OutX = [139 139+16]-124;
OutY = [211 211+29];

[zoneInY, zoneInX] = meshgrid(InY(1):InY(2), InX(1):InX(2));
[zoneOutY, zoneOutX] = meshgrid(OutY(1):OutY(2), OutX(1):OutX(2));

[cnr,stdCNR] = CNR(abs(I1N),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(1,1) = cnr;
fprintf('I1: \t%.3g \t%.3g \n', cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(I2N),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(2,1) = cnr;
fprintf('I2: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(3,1) = cnr;
fprintf('Im: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out1),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(4,1) = cnr;
fprintf('O1: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(5,1) = cnr;
fprintf('O2: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2/2+out1/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(6,1) = cnr;
fprintf('Om: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out1TAMCA),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(7,1) = cnr;
fprintf('O1 TA: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2TAMCA),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(8,1) = cnr;
fprintf('O2 TA: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2TAMCA/2+out1TAMCA/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(9,1) = cnr;
fprintf('Om TA: \t%.3g \t%.3g\n',cnr,stdCNR);

fprintf('Improvement Om/Im: %g\n',CNR(abs(out2/2+out1/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY));
tableCNR(end-1,1) = CNR(abs(out2/2+out1/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
fprintf('Improvement Om/Im TA: \t%g\n',CNR(abs(out2TAMCA/2+out1TAMCA/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY));
tableCNR(end,1) = CNR(abs(out2TAMCA/2+out1TAMCA/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY);

%% Load TAMCA results for comparison

load('results\mat-nosync\TAMCA-dual1.mat');
out1TAMCA = out1;
out2TAMCA = out2;
%% Work on Dual 1
load('results\mat-nosync\dual1-exper2-results.mat')

%% Setup of input signals + other output parameters
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

%% Region 1 - ventriculo grande
fprintf('\nRegion 1\n');
RegX = [5 5+28]; % ~11% improvement
RegY = [437 437+32];

[snr,ssnr, st, sd] = SpeckleSNR(abs(I1N),RegX,RegY);
table(1,3) = st; table2(1,3) = sd;
fprintf('I1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(I2N),RegX,RegY);
table(2,3) = st; table2(2,3) = sd;
fprintf('I2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY);
fprintf('Im: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1),RegX,RegY);
table(3,3) = st; table2(3,3) = sd;
fprintf('O1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2),RegX,RegY);
table(4,3) = st; table2(4,3) = sd;
fprintf('O2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1+out2)/2),RegX,RegY);
table(5,3) = st; table2(5,3) = sd;
fprintf('Om: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
% TAMCA
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1TAMCA),RegX,RegY);
table(6,3) = st; table2(6,3) = sd;
fprintf('O1 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2TAMCA),RegX,RegY);
table(7,3) = st; table2(7,3) = sd;
fprintf('O2 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1TAMCA+out2TAMCA)/2),RegX,RegY);
table(8,3) = st; table2(8,3) = sd;
fprintf('Om TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);

% medidas sobre el promedio de los envelopes
% [snr,ssnr] = SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY);
% [snr,ssnr] = SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY);
fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs((out1+out2)/2),RegX,RegY)/SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY));
% fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY)/SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY));
fprintf('#Pixels: %d\n', (RegX(2)-RegX(1))*(RegY(2)-RegY(1)));
table(end,3) = (RegX(2)-RegX(1))*(RegY(2)-RegY(1));

%% Region 2 - ventriculo chico
fprintf('\nRegion 2\n');
RegX = [82 82+26]; % ~10% improvement
RegY = [338 338+32];

[snr,ssnr, st, sd] = SpeckleSNR(abs(I1N),RegX,RegY);
table(1,4) = st; table2(1,4) = sd;
fprintf('I1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(I2N),RegX,RegY);
table(2,4) = st; table2(2,4) = sd;
fprintf('I2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY);
fprintf('Im: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1),RegX,RegY);
table(3,4) = st; table2(3,4) = sd;
fprintf('O1: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2),RegX,RegY);
table(4,4) = st; table2(4,4) = sd;
fprintf('O2: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1+out2)/2),RegX,RegY);
table(5,4) = st; table2(5,4) = sd;
fprintf('Om: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
% TAMCA
[snr,ssnr, st, sd] = SpeckleSNR(abs(out1TAMCA),RegX,RegY);
table(6,4) = st; table2(6,4) = sd;
fprintf('O1 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs(out2TAMCA),RegX,RegY);
table(7,4) = st; table2(7,4) = sd;
fprintf('O2 TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);
[snr,ssnr, st, sd] = SpeckleSNR(abs((out1TAMCA+out2TAMCA)/2),RegX,RegY);
table(8,4) = st; table2(8,4) = sd;
fprintf('Om TA: \t%.3g \t%.3g\t%.3g\n',snr,ssnr,st);

% medidas sobre el promedio de los envelopes
% [snr,ssnr] = SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY);
% [snr,ssnr] = SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY);
fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs((out1+out2)/2),RegX,RegY)/SpeckleSNR(abs((I1N+I2N)/2),RegX,RegY));
% fprintf('Improvement Om/Im: %g\n', SpeckleSNR(abs(out1)/2+abs(out2)/2,RegX,RegY)/SpeckleSNR(abs(I1N)/2+abs(I2N)/2,RegX,RegY));
fprintf('#Pixels: %d\n', (RegX(2)-RegX(1))*(RegY(2)-RegY(1)));
table(end,4) = (RegX(2)-RegX(1))*(RegY(2)-RegY(1));

%% CNR Measurements
fprintf('\nCNR\n');
InX =[101 101+21];
InY = [231 231+44];

% OutX = [8 8+30];
% OutY = [220 220+25];
% 
% OutX = [11  11+18];
% OutY = [598 598+38];

OutX = [6  6+48];
OutY = [211 211+28];

[zoneInY, zoneInX] = meshgrid(InY(1):InY(2), InX(1):InX(2));
[zoneOutY, zoneOutX] = meshgrid(OutY(1):OutY(2), OutX(1):OutX(2));

[cnr,stdCNR] = CNR(abs(I1N),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(1,2) = cnr;
fprintf('I1: \t%.3g \t%.3g \n', cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(I2N),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(2,2) = cnr;
fprintf('I2: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(3,2) = cnr;
fprintf('Im: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out1),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(4,2) = cnr;
fprintf('O1: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(5,2) = cnr;
fprintf('O2: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2/2+out1/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(6,2) = cnr;
fprintf('Om: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out1TAMCA),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(7,2) = cnr;
fprintf('O1 TA: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2TAMCA),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(8,2) = cnr;
fprintf('O2 TA: \t%.3g \t%.3g\n',cnr,stdCNR);
[cnr,stdCNR] = CNR(abs(out2TAMCA/2+out1TAMCA/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
tableCNR(9,2) = cnr;
fprintf('Om TA: \t%.3g \t%.3g\n',cnr,stdCNR);

fprintf('Improvement Om/Im: \t\t%g\n',CNR(abs(out2/2+out1/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY));
tableCNR(end-1,2) = CNR(abs(out2/2+out1/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY);
fprintf('Improvement Om/Im TA: \t%g\n',CNR(abs(out2TAMCA/2+out1TAMCA/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY));
tableCNR(end,2) = CNR(abs(out2TAMCA/2+out1TAMCA/2),zoneInX, zoneInY, zoneOutX, zoneOutY)-CNR(abs(I2N/2+I1N/2),zoneInX, zoneInY, zoneOutX, zoneOutY);

%% Save the table as latex
table
matrix2latex(table, 'results/specklenoiseSpace.tex', 'rowLabels', {'Orig. 1st Harm.','Orig. 2nd Harm.','Filt. 1st Harm.','Filt. 2nd Harm.','Compounded','TA-MCA Filt. 1st Harm.','TA-MCA Filt. 2nd Harm.','TA-MCA Compounded','# Pixels'},'columnLabels',{'Dataset 1 - 1','Dataset 1 - 2','Dataset 2 - 1','Dataset 2 - 2'});

table2
matrix2latex(table2, 'results/specklenoiseTime.tex', 'rowLabels', {'Orig. 1st Harm.','Orig. 2nd Harm.','Filt. 1st Harm.','Filt. 2nd Harm.','Compounded','TA-MCA Filt. 1st Harm.','TA-MCA Filt. 2nd Harm.','TA-MCA Compounded'},'columnLabels',{'Dataset 1 - 1','Dataset 1 - 2','Dataset 2 - 1','Dataset 2 - 2'});

tableCNR
matrix2latex(tableCNR, 'results/CNR.tex', 'rowLabels', {'Orig. 1st Harm.','Orig. 2nd Harm.','Compounded','Filt. 1st Harm.','Filt. 2nd Harm.','Filt. Compounded','TA-MCA Filt. 1st Harm.','TA-MCA Filt. 2nd Harm.','TA-MCA Compounded','Improvement','Improvement TA-MCA'},'columnLabels',{'Dataset 1 - 1','Dataset 1 - 2','Dataset 2 - 1','Dataset 2 - 2'});

%% load dual 2

load('results\mat-nosync\dual2-exper2-results.mat')
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

%% write figures

frame = 15;
color = [255 0 0];

% mark clutter in input images
Im1st = tissueProcessing(I1N(:,:,frame),30,GAIN2);
Im2nd = tissueProcessing(I2N(:,:,frame),30,GAIN2);

InX =[229 229+16]-124;
InY = [248 248+29];
OutX = [139 139+16]-124;
OutY = [211 211+29];
Im1st(:,:,1:3) = repmat(Im1st,[1,1,3]);
Im2nd(:,:,1:3) = repmat(Im2nd,[1,1,3]);
Im1st = drawBox(Im1st, [InX(1) InY(1)], [InX(2) InY(2)], [1 0 1]);
Im1st = drawBox(Im1st, [OutX(1) OutY(1)], [OutX(2) OutY(2)], [0 0 1]);
Im2nd = drawBox(Im2nd, [InX(1) InY(1)], [InX(2) InY(2)], [1 0 1]);
Im2nd = drawBox(Im2nd, [OutX(1) OutY(1)], [OutX(2) OutY(2)],[0 0 1] );

Im1st = scanConvertMovie(Im1st,H,W);
Im2nd = scanConvertMovie(Im2nd,H,W);

% p1 = [249, 179];
% p2 = [261, 154];
% delta = [-10 -30];
% Im1st = drawArrow(Im1st, p1, p2, color);
% Im1st = drawArrow(Im1st, p1+delta, p2+delta, color);
% Im1st = drawArrow(Im1st, p1-delta, p2-delta, color);
% Im2nd = drawArrow(Im2nd, p1, p2, color);
% Im2nd = drawArrow(Im2nd, p1+delta, p2+delta, color);
% Im2nd = drawArrow(Im2nd, p1-delta, p2-delta, color);


p1 = [139 299];
p2 = [125, 274];
color = [0, 1, 0];
outCompound = scanConvertMovie(tissueProcessing(0.5*out1(:,:,frame)+0.5*out2(:,:,frame),30,GAIN2),H,W);
outCompound(:,:,1:3) = repmat(outCompound,[1,1,3]);
outCompound = drawArrow(outCompound, p1, p2, color);

for i = 1:3
    Im1st(:,:,i) = fliplr(Im1st(:,:,i));
    Im2nd(:,:,i) = fliplr(Im2nd(:,:,i));
    outCompound(:,:,i) = fliplr(outCompound(:,:,i));
end

% imshow(outCompound);
%% save files
imwrite(outCompound,'results/figSet1Mixed.png');
imwrite(fliplr(scanConvertMovie(tissueProcessing(out1(:,:,frame),30,GAIN2),H,W)),'results/figSet1Filtered1st.png');
imwrite(fliplr(scanConvertMovie(tissueProcessing(out2(:,:,frame),30,GAIN2),H,W)),'results/figSet1Filtered2nd.png');
imwrite(Im1st,'results/figSet1Orig1st.png');
imwrite(Im2nd,'results/figSet1Orig2nd.png');

% same from Eilat
load('results\mat-nosync\TAMCA-dual2.mat')
imwrite(fliplr(scanConvertMovie(tissueProcessing(out1(:,:,frame),30,GAIN2),H,W)),'results/figSet1Filtered1stTAMCA.png');
imwrite(fliplr(scanConvertMovie(tissueProcessing(out2(:,:,frame),30,GAIN2),H,W)),'results/figSet1Filtered2ndTAMCA.png');
imwrite(fliplr(scanConvertMovie(tissueProcessing(0.5*out1(:,:,frame)+0.5*out2(:,:,frame),30,GAIN2),H,W)),'results/figSet1MixedTAMCA.png');

%% load dual 1
load('results\mat-nosync\dual1-exper2-results.mat')
I1N = (max(abs(I2(:)))/ max(abs(I1(:))) )*I1;
GAIN1 = -40;
GAIN2 = GAIN1+20;
H = 400;
W = 400;
[height,width,nFrames] = size(I1);
startFrame = (N-1)/2+1;
stopFrame = nFrames - (N-1)/2;
I1N = I1N(:,:,startFrame:stopFrame);
I2N = I2(:,:,startFrame:stopFrame);

%% write figures

frame = 28;
color = [255 0 0];

InX =[101 101+21];
InY = [231 231+44];
OutX = [6  6+48];
OutY = [211 211+28];

% mark clutter in input images
Im1st = tissueProcessing(I1N(:,:,frame),30,GAIN2);
Im2nd = tissueProcessing(I2N(:,:,frame),30,GAIN2);

Im1st(:,:,1:3) = repmat(Im1st,[1,1,3]);
Im2nd(:,:,1:3) = repmat(Im2nd,[1,1,3]);
Im1st = drawBox(Im1st, [InX(1) InY(1)], [InX(2) InY(2)], [0 0 0]);
Im1st = drawBox(Im1st, [OutX(1) OutY(1)], [OutX(2) OutY(2)], [0 0 1]);
Im2nd = drawBox(Im2nd, [InX(1) InY(1)], [InX(2) InY(2)], [0 0 0]);
Im2nd = drawBox(Im2nd, [OutX(1) OutY(1)], [OutX(2) OutY(2)],[0 0 1] );

Im1st = scanConvertMovie(Im1st,H,W);
Im2nd = scanConvertMovie(Im2nd,H,W);

p1 = [249, 179];
p2 = [261, 154];
delta = [-10 -30];
Im1st = drawArrow(Im1st, p1, p2, color);
Im1st = drawArrow(Im1st, p1+delta, p2+delta, color);
Im1st = drawArrow(Im1st, p1-delta, p2-delta, color);
Im2nd = drawArrow(Im2nd, p1, p2, color);
Im2nd = drawArrow(Im2nd, p1+delta, p2+delta, color);
Im2nd = drawArrow(Im2nd, p1-delta, p2-delta, color);


for i = 1:3
    Im1st(:,:,i) = fliplr(Im1st(:,:,i));
    Im2nd(:,:,i) = fliplr(Im2nd(:,:,i));
end
%% write files
imwrite(fliplr(scanConvertMovie(tissueProcessing(0.5*out1(:,:,frame)+0.5*out2(:,:,frame),30,GAIN2),H,W)),'results/figSet2Mixed.png');
% imwrite(scanConvertMovie(tissueProcessing(out1(:,:,frame),30,GAIN2),H,W),'results/figSet1Filtered1st.png');
% imwrite(scanConvertMovie(tissueProcessing(out2(:,:,frame),30,GAIN2),H,W),'results/figSet1Filtered2nd.png');
imwrite(Im1st,'results/figSet2Orig1st.png');
imwrite(Im2nd,'results/figSet2Orig2nd.png');


%% write image clutter 1st and 2nd
frame = 20;
I1C = (max(abs(I2clutter(:)))/ max(abs(I1clutter(:))) )*I1clutter;
GAIN1 = -50;
GAIN2 = GAIN1+20;
H = 400;
W = 400;
I1C = I1C(:,:,frame);
I2C = I2clutter(:,:,frame);

imwrite(fliplr(scanConvertMovie(tissueProcessing(I1C,30,GAIN2),H,W)),'results/figClutter1st.png');
imwrite(fliplr(scanConvertMovie(tissueProcessing(I2C,30,GAIN2),H,W)),'results/figClutter2nd.png');


%% Get some atoms from each dictionary
load('results\mat-nosync\dual2-exper2-results.mat');

N = 15;
M = 15;

% clean DT first
cutoff = 0.2;
A2 = comp(D2C'*D, D2C'*D2C, M);
DDT = D(:,~(sum(abs(D - D2C*A2).^2)< cutoff));

natoms = 5;

rng(1);
atoms = randperm(size(D1C,2));
atoms = atoms(1:natoms);
D1 = tissueProcessing(D1C(:,atoms),1.2,0);
D2 = tissueProcessing(D2C(:,atoms),1.2,0);
DT = tissueProcessing(D(:,atoms),1.5,0);

for i = 1:natoms
    imwrite(reshape(D1(:,i),[N M]),['results/DC1-atom-' num2str(i) '.png']);
end

for i = 1:natoms
    imwrite(reshape(D2(:,i),[N M]),['results/DC2-atom-' num2str(i) '.png']);
end

for i = 1:natoms
    imwrite(reshape(DT(:,i),[N M]),['results/DT-atom-' num2str(i) '.png']);
end

D1 = tissueProcessing(D1C,1.2,0);
dN = ceil(sqrt(size(D1,2)));
dM = ceil(size(D1,2)/dN);
imD1C = showdict(D1,[15 15],dN, dM,'highcontrast','whitelines');

D2 = tissueProcessing(D2C,1.2,0);
dN = ceil(sqrt(size(D2,2)));
dM = ceil(size(D2,2)/dN);
imD2C = showdict(D2,[15 15],dN, dM,'highcontrast','whitelines');

DT = tissueProcessing(DDT,1.5,0);
dN = ceil(sqrt(size(DT,2)));
dM = ceil(size(DT,2)/dN);
imDT = showdict(DT,[15 15],dN, dM,'highcontrast','whitelines');

imwrite(imDT, 'results/Dictionary-DT.png');
imwrite(imD1C,'results/Dictionary-D1C.png');
imwrite(imD2C,'results/Dictionary-D2C.png');

% Get the phase of the atoms 
D1 = getPhase(D1C, N,M);
imD1C = showdict(D1,[15 15],dN, dM,'highcontrast','whitelines');
% imshow(imD1C); 

D2 = getPhase(D1C, N,M);
imD2C = showdict(D2,[15 15],dN, dM,'highcontrast','whitelines');
% imshow(imD2C); 

DT = getPhase(DDT, N,M);
imDT = showdict(DT,[15 15],dN, dM,'highcontrast','whitelines');
% imshow(imDT); 


%% example figures for presentation (signal model slides)

load('results\mat-nosync\dual2-exper2-results.mat');

imwrite(scanConvertMovie(tissueProcessing(mean(I1N - out1,3),30,GAIN2),H,W),'results/clutter1st.png');
imwrite(scanConvertMovie(tissueProcessing(mean(I2N - out2,3),30,GAIN2),H,W),'results/clutter2nd.png');


[II1, mla] = readDataset('data/exampleonly1st/');
II1 = stb(II1,mla,mla/2);
[II2, mla] = readDataset('data/exampleonly2nd/');
II2 = stb(II2,mla,mla/2);

imwrite(scanConvertMovie(tissueProcessing(II1(:,:,5),30,-30),H,W),'results/tissue1st.png');
imwrite(scanConvertMovie(tissueProcessing(II2(:,:,5),30,-30),H,W),'results/tissue2nd.png');

