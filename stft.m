%myvoiceの読み込み
[y, fs] = audioread("myvoice.wav");
y = ( y(:, 1) + y(:, 2) ) / 2;

%信号長取得
signalLength = size(y, 1);

%fftSize定義
fftSize = 1024;

%shiftSize定義
shiftSize = fftSize / 2;

%ハン窓作成
window = hann(fftSize);

%行サイズ計算
numRow = ceil((signalLength - fftSize) / shiftSize) + 1;

%paddingSize計算
paddingSize = fftSize - 1;

%zeros生成
padding = zeros(paddingSize, 1);

%padding結合
yPadding = [y;padding];

%spec定義
spec = zeros(fftSize, numRow);
 
for n = 1:numRow
    %yから抽出
    vec = yPadding(1 + (n - 1)*shiftSize:fftSize + (n - 1)*shiftSize,1);

    %ハン窓乗算
    vecWindow = vec .* window;

    %fft
    yDft = fft(vecWindow);

    %結果格納
    spec(:, n) = yDft;
end

%---パワースペクトログラムの表示部---

%絶対値の計算
ampSpec = sqrt(real(spec).^2 + imag(spec).^2);
%specAbs = abs(spec);

%db変換
powerSpec = 20 * log10(ampSpec);

%縦軸生成(周波数[Hz])
%fsのfftsize分割
yGrid = linspace(0,fs,fftSize);

%横軸生成(時間[s])
%1秒のデータ数=fs個
%全体の時間=signalLength/fs
xGrid = linspace(0, signalLength / fs, numRow);

%スペクトログラム描画
imagesc(xGrid, yGrid, powerSpec);
axis xy;
fontsize(gca, 15, "pixels")
ylabel("周波数 [Hz]", 'FontSize', 18);
xlabel("時間 [s]", 'FontSize', 18);
caxis([-60 30]);
ylim([0 fs/2]);
colorbar;

