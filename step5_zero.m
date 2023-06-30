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

% 以降，欠損値(nan)は0に変換
losX = fillmissing(losX, "constant", 0);

% NMF変数の乱数初期化
W = rand(nRow, rank);
H = rand(rank, nCol);

% 更新回数
nItr = 1000;

% EuNMF 
% コスト関数値格納行列
cost = zeros(nItr + 1, 1);

% コスト関数値の初期値格納(フロベニウスノルムの二乗値計算)
err = M .* (losX - W*H);
cost(1) = sqrt(sum(err.^2, "all"));

for iItr = 1:nItr
    % EuNMFの更新式
    W = W .* ((losX * H.') ./ (W*H.*M*H.'));
    H = H .* ((W.' * losX) ./ (W.' * (W*H.*M)) );


    % 正規化
    normalC = sum(W, 1);
    W = W ./ normalC;
    H = normalC.' .* H;

    % コスト計算
    err = M .* (losX - W*H);
    cost(iItr+1) = sqrt(sum(err.^2, "all"));
end

% 近似された観測行列の表示
Xhat = W * H;
figure; plot(cost);
title("コスト関数値推移");
xlabel("反復回数");
ylabel("コスト");
figure; imagesc(X);
title("完全な観測行列");
figure; imagesc(Xhat);
title("近似行列");
figure; imagesc(abs(X-Xhat));
title("誤差[abs(X - Xhat)]");