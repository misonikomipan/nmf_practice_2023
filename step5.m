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

% マスク行列の生成
M = not(isnan(losX));

% NMF変数の乱数初期化
W = rand(nRow, rank);
H = rand(rank, nCol);
