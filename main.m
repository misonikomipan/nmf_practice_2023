clear; close all; clc;

% 基底
K = 2;

% 音声ファイルの読み込み
audiofile = "sound.wav";
[y, fs] = audioread(audiofile);
if(size(y,2) == 2)
    y = ( y(:, 1) + y(:, 2) ) / 2;
end

% DGTtoolのインスタンス生成
F = DGTtool("windowShift", 1024, "windowLength", 2048, "FFTnum", 2048, "windowName", "blackman");

% 音声ファイルのSTFT
X = F(y);
Y = STFT(audiofile, "fftSize", 2048);

% 行列サイズの取得
[I, J] = size(X);

Nyquist = fix(I/2) + 1; % ナイキスト周波数にあたるデータ位置
X = X(1:Nyquist,:); % ナイキスト周波数までを抽出

% 振幅スペクトルに変換
ampX = abs(X);

% 位相スペクトルに変換
phaseX = angle(X);

% パラメータ
nItr = 1000; % 更新式の反復回数

[W, H] = NMF(ampX, K, "convergenceCurve",true,"nItr",nItr,"typeCostFunction","KL");

% 分離振幅スペクトル
Xhat1 = W(:,1)*H(1,:);
Xhat2 = W(:,2)*H(2,:);

% 分離複素スペクトル
cpxXhat1 = Xhat1 .* exp(1i*phaseX);
cpxXhat2 = Xhat2 .* exp(1i*phaseX);

% 分離音声の復元
sigXhat1 = F.pinv(cpxXhat1);
sigXhat2 = F.pinv(cpxXhat2);

% 分解音の表示
figure; plot(y);
figure; plot(sigXhat1);
figure; plot(sigXhat2);