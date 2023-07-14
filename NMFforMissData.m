function [W, H] = NMFforMissData(X,K,args)
% EuNMF: 入力の2次元実数非負行列と基底数から,指定されたNMFで分解したW,Hを出力する関数
%   [Input]
%       X: 2次元実数非負観測行列(3次元配列や複素行列の場合はエラー)
%       K: 基底ベクトル数(正の自然数スカラーのみ)
%       nItr: 更新式の反復回数(正の自然数スカラーのみ)
%       typeCostFunction: コスト関数の距離の種類("EU", "KL", "IS")
%       convergenceCurve: コスト関数値の収束カーブの表示のon/off(true, false)
%   [Output]
%       W: 分解された行列(XのサイズがI*JのときI*K)
%       H: 分解された行列(XのサイズがI*JのときK*J)
%
%

arguments
    X (:,:) double {mustBeReal}
    K (1,1) double {mustBePositive, mustBeInteger}
    args.nItr (1,1) double {mustBePositive, mustBeInteger} = 1000
    args.typeCostFunction string {mustBeMember(args.typeCostFunction, ["EU","KL","IS"])} = "EU"
    args.convergenceCurve (1,1) logical = false
end
nItr = args.nItr;
typeCostFunction = args.typeCostFunction;
convergenceCurve = args.convergenceCurve;

% ---入力行列の検証---
% Xが負値を含むときエラー出力(欠損値を0に置換した上で負値を含むとき)
mustBeNonnegative(fillmissing(X, "constant", 0));

% 配列サイズの取得
[I, J] = size(X);

% NMF変数の乱数初期化
W = rand(I, K);
H = rand(K, J);

if(typeCostFunction == "EU")
    [W,H,cost] = local_EUNMF(X,W,H,args);
end
if(typeCostFunction == "KL")
    [W,H,cost] = local_KLNMF(X,W,H,args);
end
if(typeCostFunction == "IS")
    [W,H,cost] = local_ISNMF(X,W,H,args);
end

if(convergenceCurve)
    % コスト関数値のグラフ描画(線形)
    figure; plot(cost);
    xlabel("反復回数", "FontSize", 14);
    ylabel("コスト関数値(線形軸)", "FontSize", 14);

    % コスト関数値のグラフ描画(対数軸)
    figure; semilogy(cost);
    xlabel("反復回数", "FontSize", 14);
    ylabel("コスト関数値(対数軸)", "FontSize", 14);
end
end

%% local function
function cost = local_FrobeniusNorm(X, Xhat, M)
err = M.*(X - Xhat);
cost = sqrt(sum(err.^2, "all"));
end

function [nomalW, nomalH] = local_normalization(W,H)
nomalC = sum(W, 1);
nomalW = W ./ nomalC;
nomalH = nomalC.'.*H;
end

function [W, H, cost] = local_EUNMF(X,W,H,args)
% マスク行列の生成
M = not(isnan(X));
% 入力行列の欠損値を0に置換
X = fillmissing(X, "constant", 0);
% コスト関数値の格納行列定義，初期値格納
cost = 0;
if(args.convergenceCurve)
    cost = zeros(args.nItr + 1,1);
    cost(1) = local_FrobeniusNorm(X,W*H,M);
end
for iItr = 1:args.nItr
    W = W .* ((X * H.') ./ (W*H.*M*H.'));
    H = H .* ((W.' * X) ./ (W.' * (W*H.*M)) );        

    % 正規化
    [W, H] = local_normalization(W,H);

    % コスト関数値格納
    if(args.convergenceCurve)
        cost(iItr+1) = local_FrobeniusNorm(X,W*H,M);
    end
end    
end

function [W, H, cost] = local_KLNMF(X,W,H,args)
% マスク行列の生成
M = not(isnan(X));
% 入力行列の欠損値を0に置換
X = fillmissing(X, "constant", 0);
% コスト関数値の格納行列定義，初期値格納
cost = 0;
if(args.convergenceCurve)
    cost = zeros(args.nItr + 1,1);
    cost(1) = local_FrobeniusNorm(X,W*H,M);
end
ONE = ones(size(X));
for iItr = 1:args.nItr
    W = W .* (((X./(W*H))*H.')./((M.*ONE)*H.'));
    H = H .* ((W.'*(X./(W*H)))./(W.'*(ONE.*M)));

    % 正規化
    [W, H] = local_normalization(W,H);

    % コスト関数値格納
    if(args.convergenceCurve)
        cost(iItr+1) = local_FrobeniusNorm(X,W*H,M);
    end
end
end

function [W, H, cost] = local_ISNMF(X,W,H,args)
% マスク行列の生成
M = not(isnan(X));
% 入力行列の欠損値を0に置換
X = fillmissing(X, "constant", 0);
% コスト関数値の格納行列定義，初期値格納
cost = 0;
if(args.convergenceCurve)
    cost = zeros(args.nItr + 1,1);
    cost(1) = local_FrobeniusNorm(X,W*H,M);
end
ONE = ones(size(X));
for iItr = 1:args.nItr
    W = W.*((X./(W*H).^2 )*H.' ./ (M.*ONE./(W*H)*H.')).^(0.5);
    H = H.*(W.'*(X./((W*H).^2)) ./ (W.'*((ONE./(W*H)).*M))).^(0.5);

    % 正規化
    [W, H] = local_normalization(W,H);

    % コスト関数値格納
    if(args.convergenceCurve)
        cost(iItr+1) = local_FrobeniusNorm(X,W*H,M);
    end
end
end