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

% 初期値の観測行列
Xhat = W*H;

figure; imagesc(Xhat);

% Eu-NMFの更新式
for iItr = 1:nItr
    W = W .* ( (X*H.') ./ (W*(H*H.')) );
    H = H .* ( (W.'*X) ./ ((W.'*W)*H) );
end

% 近似された観測行列
Xhat = W*H;

figure; imagesc(Xhat);

% 宿題（ステップ1）
% Eu-NMFのコスト関数（フロベニウスノルムの二乗値）を
% 反復ごとに計算して「横軸反復回数vs縦軸コスト関数値」
% のグラフを描いてみてください