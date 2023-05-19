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

% NMF変数の乱数初期化
W = rand(I, K);
H = rand(K, J);

% --- IS_NMF ---
W_IS = W;
H_IS = H;
% コスト関数値格納行列定義
cost_ISNMF = zeros(nItr + 1,1);

% コスト関数値の初期値格納(フロベニウスノルムの二乗値計算)
err = X - W_IS*H_IS;
traceErr = sum(err.*err, 'all');
cost_ISNMF(1) = sqrt(traceErr);
