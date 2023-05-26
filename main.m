clear; close all; clc;

% 行列のサイズ
I = 100; % Xの行数
J = 100; % Xの列数
K = 10; % 基底数

% パラメータ
nItr = 1000; % 更新式の反復回数

% 非負観測行列の生成
trueW = rand(I, K); % 非負乱数（開区間(0, 1)）
trueH = rand(K, J); % 非負乱数（開区間(0, 1)）
X = trueW * trueH; % ランクKの非負観測行列

% Xの表示
figure; imagesc(X);

[W_EU, H_EU] = NMF(X, K, "convergenceCurve",true,"nItr",nItr,"typeCostFunction","KL");

% 近似された観測行列の表示
Xhat_EU = W_EU*H_EU;
figure; imagesc(Xhat_EU);