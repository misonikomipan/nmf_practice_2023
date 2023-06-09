% 宿題（ステップ1）
% Eu-NMFのコスト関数（フロベニウスノルムの二乗値）を
% 反復ごとに計算して「横軸反復回数vs縦軸コスト関数値」
% のグラフを描いてみてください


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

% コスト関数値格納行列定義
cost = zeros(nItr + 1,1);

% コスト関数値の初期値格納(フロベニウスノルムの二乗値計算)
err = X - W*H;
traceErr = sum(err.*err, 'all');
cost(1) = sqrt(traceErr);

% Eu-NMFの更新式
for iItr = 1:nItr
    W = W .* ( (X*H.') ./ (W*(H*H.')) );
    H = H .* ( (W.'*X) ./ ((W.'*W)*H) );    
    
    % Wの列毎に正規化係数を計算（列の総和）
    % todo
    
    % Wに正規化係数を適用
    % todo

    % Hに正規化係数を適用
    % todo

    % フロベニウスノルムの二乗値計算
    err = X - W*H;
    traceErr = sum(err.*err, 'all');
    cost(iItr+1) = sqrt(traceErr);
end

% 近似された観測行列の表示
Xhat = W*H;
figure; imagesc(Xhat);

% コスト関数値のグラフ描画(線形)
figure; plot(cost);
xlabel("反復回数", "FontSize", 14);
ylabel("コスト関数値(線形軸)", "FontSize", 14);

% コスト関数値のグラフ描画(対数軸)
figure; semilogy(cost);
xlabel("反復回数", "FontSize", 14);
ylabel("コスト関数値(対数軸)", "FontSize", 14);