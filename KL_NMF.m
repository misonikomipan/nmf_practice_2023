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

W_KL = W;
H_KL = H;
% コスト関数値格納行列定義
cost_KLNMF = zeros(nItr + 1,1);

% コスト関数値の初期値格納(フロベニウスノルムの二乗値計算)
err = X - W_KL*H_KL;
traceErr = sum(err.*err, 'all');
cost_KLNMF(1) = sqrt(traceErr);

% onesの事前準備
ONE = ones(I,J);

% KL-NMFの更新式
for iItr = 1:nItr
    W_KL = W_KL .* (((X./(W_KL*H_KL))*H_KL.')./ (ONE*H_KL.'));
    H_KL = H_KL .* ((W_KL.'*(X./(W_KL*H_KL))) ./ (W_KL.'*ONE));

    % Wの列毎に正規化係数を計算（列の総和）
    % todo
    nomalC_KL = sum(W_KL,1);

    % Wに正規化係数を適用
    % todo
    W_KL = W_KL ./ nomalC_KL;
       
    % Hに正規化係数を適用
    % todo
    for i = 1:K
        for j = 1:J
            H_KL(i,j) = nomalC_KL(i) * H_KL(i,j);
        end
    end

    nomalWH_KL = W_KL * H_KL;
    % フロベニウスノルムの二乗値計算
    err_KL = X - nomalWH_KL;
    traceErr_KL = sum(err_KL.*err_KL, 'all');
    cost_KLNMF(iItr+1) = sqrt(traceErr_KL);
end

% 近似された観測行列の表示
Xhat_KL = nomalWH_KL;
figure; imagesc(Xhat_KL);

% コスト関数値のグラフ描画(線形)
figure; plot(cost_KLNMF);
xlabel("反復回数", "FontSize", 14);
ylabel("コスト関数値(線形軸)", "FontSize", 14);

% コスト関数値のグラフ描画(対数軸)
figure; semilogy(cost_KLNMF);
xlabel("反復回数", "FontSize", 14);
ylabel("コスト関数値(対数軸)", "FontSize", 14);