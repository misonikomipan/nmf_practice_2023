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
title("近似行列 EUNMF");

figure; imagesc(abs(Xhat_EU-X)./X*100);
title("誤差率 EUNMF");

figure; imagesc(Xhat_KL);
title("近似行列 KLNMF");

figure; imagesc(abs(Xhat_KL-X)./X*100);
title("誤差率 KLNMF");

figure; imagesc(Xhat_IS);
title("近似行列 ISNMF");

figure; imagesc(abs(Xhat_IS-X)./X*100);
title("誤差率 ISNMF");