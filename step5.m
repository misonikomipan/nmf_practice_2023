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

[W, H] = NMFforMissData(losX,rank,"convergenceCurve",true,"typeCostFunction","KL");

% 近似された観測行列の表示
Xhat = W * H;
figure; imagesc(X);
title("完全な観測行列");
figure; imagesc(Xhat);
title("近似行列");
figure; imagesc(abs(X-Xhat));
title("誤差[abs(X - Xhat)]");