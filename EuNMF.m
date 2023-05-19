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

% ------ Eu-NMF ------
% コスト関数値格納行列定義
cost_EuNMF = zeros(nItr + 1,1);
W_Eu = W;
H_Eu = H;

% コスト関数値の初期値格納(フロベニウスノルムの二乗値計算)
err = X - W_Eu*H_Eu;
traceErr = sum(err.*err, 'all');
cost_EuNMF(1) = sqrt(traceErr);

% 正規化によるWHが変化しないことの確認
cntErr = 0;

% Eu-NMFの更新式
for iItr = 1:nItr
    W_Eu = W_Eu .* ( (X*H_Eu.') ./ (W_Eu*(H_Eu*H_Eu.')) );
    H_Eu = H_Eu .* ( (W_Eu.'*X) ./ ((W_Eu.'*W_Eu)*H_Eu) );    
    
    % WHが変化しないことの確認
    prevWH = W_Eu*H_Eu;

    % Wの列毎に正規化係数を計算（列の総和）
    % todo
    nomalC_Eu = sum(W_Eu,1);

    % Wに正規化係数を適用
    % todo
    W_Eu = W_Eu ./ nomalC_Eu;
       
    % Hに正規化係数を適用
    % todo
    H_Eu = nomalC_Eu.' .*H_Eu;

    % WHが変化しないことの確認
    nomalWH_Eu = W_Eu*H_Eu;
    if(nomalWH_Eu - prevWH > 10^-15)
        cntErr = cntErr + 1;
    end

    % フロベニウスノルムの二乗値計算
    err = X - nomalWH_Eu;
    traceErr = sum(err.*err, 'all');
    cost_EuNMF(iItr+1) = sqrt(traceErr);
end

% 近似された観測行列の表示
Xhat_Eu = nomalWH_Eu;
figure; imagesc(Xhat_Eu);

% コスト関数値のグラフ描画(線形)
figure; plot(cost_EuNMF);
xlabel("反復回数", "FontSize", 14);
ylabel("コスト関数値(線形軸)", "FontSize", 14);

% コスト関数値のグラフ描画(対数軸)
figure; semilogy(cost_EuNMF);
xlabel("反復回数", "FontSize", 14);
ylabel("コスト関数値(対数軸)", "FontSize", 14);
