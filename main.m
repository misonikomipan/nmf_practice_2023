clear; close all; clc;

% 基底
K = 2;

% 音声ファイル
audiofile = "sound.wav";

% 音声ファイルの振幅スペクトログラム変換(非負観測行列)
[X, fs] = STFT(audiofile);
Nyquist = fix(size(X,1)/2) + 1;
X = X(1:Nyquist,:);

% 行列サイズの取得
[I, J] = size(X);

% パラメータ
nItr = 1000; % 更新式の反復回数

% Xの表示
yGrid = linspace(0, fs/2, Nyquist);
xGrid = linspace(0, ((I+1)*J)/(fs), J);
colorlim = [0, 80];
ylimit = [0, 5000];

figure; imagesc(xGrid, yGrid, X);
axis xy;
fontsize(gca, 12, "pixels");
ylabel("周波数 [Hz]");
xlabel("時間 [s]");
caxis(colorlim);
ylim(ylimit);
colorbar;

[W, H] = NMF(X, K, "convergenceCurve",true,"nItr",nItr,"typeCostFunction","EU");


% 分解音の表示
Y=W(:,1)*H(1,:);
figure; imagesc(xGrid, yGrid, W(:,1)*H(1,:));
axis xy;
fontsize(gca, 12, "pixels");
ylabel("周波数 [Hz]");
xlabel("時間 [s]");
caxis(colorlim);
ylim(ylimit);
colorbar;

figure; imagesc(xGrid, yGrid, W(:,2)*H(2,:));
axis xy;
fontsize(gca, 12, "pixels");
ylabel("周波数 [Hz]");
xlabel("時間 [s]");
caxis(colorlim);
ylim(ylimit);
colorbar;

