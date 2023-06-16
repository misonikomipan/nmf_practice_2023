clear; close all; clc;

% 基底
K = 2;

% 音声ファイル
audiofile = "sound.wav";

% 音声ファイルの振幅スペクトログラム変換(非負観測行列)
[X, fs] = STFT(audiofile, "fftSize", 2048);
Nyquist = fix(size(X,1)/2) + 1; % ナイキスト周波数にあたるデータ位置
X = X(1:Nyquist,:); % ナイキスト周波数までを抽出

% 行列サイズの取得
[I, J] = size(X);

% パラメータ
nItr = 1000; % 更新式の反復回数

[W, H] = NMF(X, K, "convergenceCurve",true,"nItr",nItr,"typeCostFunction","KL");

% 分解音の表示
% 縦軸を0からナイキスト周波数までI個に分割
yGrid = linspace(0, fs/2, I);
% numRow=(signalLength - fftSize)/shiftSize+1
% numRow=J,fftSize=2*shiftSize=I
% 再生時間=signalLength/fs, 以上から再生時間算出
% 横軸として再生時間をJ個に分割
xGrid = linspace(0, ((I+1)*J)/(fs), J);
colorlim = [0, 80];
ylimit = [0, 10000];
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

