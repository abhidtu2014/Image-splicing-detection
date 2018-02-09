function [h] = ILBP81ri(X)
%Im=rgb2gray(X);
%ILBP81ri texture features
%
%INPUT:
%X: grey-scale image
%
%OUTPUT:
%h:     featureVector
%

    % Image size
    [L M] = size(X);

    codesILBP81 = codeILBP81(X);
    histCodesILBP81 = sum(hist(codesILBP81,0:511)');

    lut81 = lutRotInv('CCR81ri');

    nbins = length(unique(lut81));

    histCodesILBP81ri = zeros(1,nbins);

    for p = 1:length(histCodesILBP81);
        r = lut81(p);
        histCodesILBP81ri(r+1) =  histCodesILBP81ri(r+1) + histCodesILBP81(p);
    end

    h  =  histCodesILBP81ri / ((L-2)*(M-2));

    %Remove the first bin (is always 0)
    h(1) = [];

end

function [Xcod] = codeILBP81(X)
%
% [Xcod] = Features_CCR81(X,T)
% 
% Computes CCR81 codes of a single channel texture image
%
% Inputs:
%   X    -  Single channel texture image (at least 3x3 pixels)
%   T    -  Binarization threshold
%
% Outputs:
%   Xcod  -  Image with CCR81 codes


    % Coeficients used in bilinear interpolation
    sqrt_2 = 1.4142; 
    center  =  (1-1/sqrt_2)^2;
    corner  =  (1/sqrt_2)^2;
    diagon  =  (1-1/sqrt_2)*(1/sqrt_2);

    % Conversion to avoid errors when using sort, unique...
    X = double(X);

    % Image size
    [L M] = size(X);

    % Displacement directions
    north = 1:L;   % N: North
    south = 3:L+2; % S: South
    equad = 2:L+1; % equator
    east  = 3:M+2; % E: East
    west  = 1:M;   % W: West
    meri  = 2:M+1; % meridian
    % 0: No displacement

    [X0, XN, XNE, XE, XSE, XS, XSW, XNW, XW] = deal(zeros(L+2,M+2));    
    [pN, pNE, pE, pSE, pS, pSW, pW, pNW, p0] = deal(zeros(L+2,M+2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pN (north,meri) = diagon*X ;
    pS (south,meri) = diagon*X ;
    pE (equad,east) = diagon*X ;
    pW (equad,west) = diagon*X ;

    pNE(north,east) = X ;
    pNW(north,west) = X ;
    pSE(south,east) = X ;
    pSW(south,west) = X ;

    p0 (equad,meri) = center*X ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    XN(north,meri)  = X;
    XS(south,meri)  = X;
    XE(equad,east)  = X;
    XW(equad,west)  = X;

    X0(equad,meri)  = X;

    XSE = round(corner*pSE + pS + pE + p0);
    XSW = round(corner*pSW + pS + pW + p0);
    XNW = round(corner*pNW + pN + pW + p0);
    XNE = round(corner*pNE + pN + pE + p0);

    T = round((XN + XS + XE + XW + XNE + XNW + XSE + XSW + X0)/9);

    X0 = 256*(XSE>=T) + 128*(XS>=T) + 64*(XSW>=T) + 32*(XE>=T) + 16*(X0>=T) + 8*(XW>=T) + 4*(XNE>=T) + 2*(XN>=T) + (XNW>=T);

    Xcod = X0(3:L,3:M);

end