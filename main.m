clear; close all; clc;

% 基底
K = 10;

% 音声ファイル
audiofile = "kitamuravoice.wav";

% 音声ファイルの振幅スペクトログラム変換(非負観測行列)
X = STFT(audiofile);
Nyquist = fix(size(X,1)/2) + 1;
X = X(1:Nyquist,:);

% 行列サイズの取得
[I, J] = size(X);

% パラメータ
nItr = 1000; % 更新式の反復回数

% Xの表示
figure; imagesc(X);

[W, H] = NMF(X, K, "convergenceCurve",true,"nItr",nItr,"typeCostFunction","EU");

%imagesc(W_EU,H_EU);
% 近似された観測行列の表示
Xhat = W*H;
figure; imagesc(Xhat);

