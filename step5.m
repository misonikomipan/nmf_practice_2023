clear; close all; clc;

% 行列サイズとランク
nRow = 100;
nCol = 100;
rank = 10;

% 欠損値割合
nanRate = 0.4;
nNan = nRow * nCol * nanRate;

% 階数=rankの非負観測行列
x1 = rand(nRow, rank);
x2 = rand(rank, nCol);
X = x1 * x2;

% 欠損値を上書き
losX = X;
losX(randperm(nRow * nCol, nNan)) = nan;

[W_EU, H_EU] = NMFforMissData(losX,rank,"typeCostFunction","EU");
[W_KL, H_KL] = NMFforMissData(losX,rank,"typeCostFunction","KL");
[W_IS, H_IS] = NMFforMissData(losX,rank,"typeCostFunction","IS");

% 近似された観測行列の表示
Xhat_EU = W_EU * H_EU;
Xhat_KL = W_KL * H_KL;
Xhat_IS = W_IS * H_IS;

figure; imagesc(X);
title("完全な観測行列");

figure; imagesc(losX);
title("欠損値あり観測行列");

figure; imagesc(Xhat_EU);
t = 'EUNMF $\hat{X}_{EU}$';
title(t,'interpreter','latex');

figure; imagesc(abs(Xhat_EU-X)./X*100);
t = 'EUNMF $|\hat{X}_{EU}|\otimes X^{-1} \times 100$';
title(t,'interpreter','latex');

figure; imagesc(Xhat_KL);
t = 'KLNMF $\hat{X}_{KL}$';
title(t,'interpreter','latex');

figure; imagesc(abs(Xhat_KL-X)./X*100);
t = 'KLNMF $|\hat{X}_{KL}|\otimes X^{-1} \times 100$';
title(t,'interpreter','latex');

figure; imagesc(Xhat_IS);
t = 'ISNMF $\hat{X}_{IS}$';
title(t,'interpreter','latex');

figure; imagesc(abs(Xhat_IS-X)./X*100);
t = 'ISNMF $|\hat{X}_{IS}|\otimes X^{-1} \times 100$';
title(t,'interpreter','latex');

% ---精度評価---
% 誤差平均(絶対値)
errAVG_EU = sum(abs(Xhat_EU-X),"all")/(nRow*nCol);
errAVG_KL = sum(abs(Xhat_KL-X),"all")/(nRow*nCol);
errAVG_IS = sum(abs(Xhat_IS-X),"all")/(nRow*nCol);

% 誤差(絶対値)の標準偏差
errDev_EU = (sum((abs(Xhat_EU-X)-errAVG_EU)^2,'all')/(nRow*nCol))^(0.5);
errDev_KL = (sum((abs(Xhat_KL-X)-errAVG_KL)^2,'all')/(nRow*nCol))^(0.5);
errDev_IS = (sum((abs(Xhat_IS-X)-errAVG_IS)^2,'all')/(nRow*nCol))^(0.5);

errAVG = [errAVG_EU,errAVG_KL,errAVG_IS];
errDev = [errDev_EU,errDev_KL,errDev_IS];

figure;bar(errAVG);
title("誤差平均:1-EU,2-KL,3-IS");
figure;bar(errDev);
title("誤差の標準偏差:1-EU,2-KL,3-IS");