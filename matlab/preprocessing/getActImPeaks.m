function thetaf = getActImPeaks(actIM, peakth, bgIM, peakFuncOpt, exclusionMask)
if nargin < 4
    peakFuncOpt = 1;
end
if nargin < 5
    exclusionMask = false(size(actIM));
end

switch peakFuncOpt
    case 1
        peakFunc = @gaussianPeaks;
    case 2
        peakFunc = @gaussianPeaksIntegrated;
end

mu_bg = median(actIM,'all','omitmissing');
sigma_bg = mad(actIM(~isnan(actIM)),1,'all') ./ 0.6741891400433162;
peak_thresh = mu_bg + peakth * sigma_bg;

explored = actIM .* ~exclusionMask; pTmp = ordfilt2(explored,8,ones(3)) > peak_thresh & explored == ordfilt2(explored, 9, ones(3));
thetaf = zeros(0,4);
pLocs = zeros(0,2);

if sum(pTmp(:))
    nPeaks = sum(pTmp(:));

    [pY, pX] = ind2sub(size(pTmp),find(pTmp));
    switch peakFuncOpt
        case 1
            amp = actIM(pTmp);
        case 2
            amp = actIM(pTmp) ./ 0.75;
    end
    widths = 0.35 * ones(nPeaks,1);
    
    actSelPix = imdilate(pTmp, ones(9)) & ~isnan(actIM);
    [actSelY, actSelX] = ind2sub(size(pTmp),find(actSelPix));
    
    theta0 = [thetaf; [amp, pY, pX, widths]];
    pLocs = [pY, pX];

    ub = [Inf*ones(size(theta0,1),1),min(size(actIM,1)+0.5,pLocs(:,1)+1.5),min(size(actIM,2)+0.5,pLocs(:,2)+1.5),5*ones(size(theta0,1),1)];
    lb = [zeros(size(theta0,1),1),max(0.5,pLocs(:,1)-1.5),max(0.5,pLocs(:,2)-1.5),zeros(size(theta0,1),1)];
    opts = optimset('MaxFunEvals',5000,'Display','off');

    thetaf = lsqcurvefit(peakFunc,theta0,[actSelY,actSelX],actIM(actSelPix)-mu_bg,lb,ub,opts);

    pIM = false(size(actIM));
    pIM(sub2ind(size(actIM),round(thetaf(:,2)),round(thetaf(:,3)))) = true;
    bufferMask = pIM; %imdilate(pIM,ones(3,3));

    fitIM = zeros(size(actIM));
    fitIM(actSelPix) = peakFunc(thetaf,[actSelY,actSelX]);
    resIM = actIM - fitIM - mu_bg;
    resIM = resIM ./ ((fitIM+bgIM+mu_bg) ./ bgIM);

    explored = resIM .* ~bufferMask .* ~exclusionMask; % .* (resIM > actIM ./ 3);
    pTmp = zeros(size(explored));
    if max(explored(fitIM > 1e-3)) > (peak_thresh-mu_bg)
        pTmp = explored == max(explored(fitIM > 1e-3)); % & explored == ordfilt2(explored, 9, ones(3)) & (fitIM > 1e-3);
    end

    while sum(pTmp(:))
        idxNew = find(pTmp,1,'first');
        [pY, pX] = ind2sub(size(pTmp),idxNew);
        switch peakFuncOpt
            case 1
                amp = actIM(pTmp);
            case 2
                amp = actIM(pTmp) ./ 0.75;
        end

        thetaf = [thetaf; [amp, pY, pX, 0.35]];
        pLocs = [pLocs; [pY, pX]];

        CC = bwconncomp(actSelPix);
        for i = 1:CC.NumObjects
            if ismember(idxNew, CC.PixelIdxList{i})
                actSelPix = false(size(actSelPix));
                actSelPix(CC.PixelIdxList{i}) = true;
                [actSelY, actSelX] = ind2sub(size(pTmp),find(actSelPix));
                
                % Find which theta indices correspond to peaks in this connected component
                peakIndices = sub2ind(size(actIM), round(thetaf(:,2)), round(thetaf(:,3)));
                thetaIdxsToFit = actSelPix(peakIndices);

                break;
            end
        end

        ub = [Inf*ones(sum(thetaIdxsToFit),1),min(size(actIM,1)+0.5,pLocs(thetaIdxsToFit,1)+1.5),min(size(actIM,2)+0.5,pLocs(thetaIdxsToFit,2)+1.5),5*ones(sum(thetaIdxsToFit),1)];
        lb = [zeros(sum(thetaIdxsToFit),1),max(0.5,pLocs(thetaIdxsToFit,1)-1.5),max(0.5,pLocs(thetaIdxsToFit,2)-1.5),zeros(sum(thetaIdxsToFit),1)];
        opts = optimset('MaxFunEvals',5000,'Display','off');
        thetaf(thetaIdxsToFit,:) = lsqcurvefit(peakFunc,thetaf(thetaIdxsToFit,:),[actSelY,actSelX],actIM(actSelPix)-mu_bg,lb,ub,opts);
        
        pIM = false(size(actIM));
        pIM(sub2ind(size(actIM),round(thetaf(:,2)),round(thetaf(:,3)))) = true;
        bufferMask = pIM; %imdilate(pIM,ones(3,3));

        fitIM(actSelPix) = peakFunc(thetaf,[actSelY,actSelX]);
        resIM = actIM - fitIM - mu_bg;
        resIM = resIM ./ ((fitIM+bgIM+mu_bg) ./ bgIM);
    
        explored = resIM .* ~bufferMask .* ~exclusionMask; % .* (resIM > actIM ./ 3);
        pTmp = zeros(size(explored));
        if max(explored(fitIM > 1e-3)) > (peak_thresh-mu_bg)
            pTmp = explored == max(explored(fitIM > 1e-3)); % & explored == ordfilt2(explored, 9, ones(3)) & (fitIM > 1e-3);
        end
    end

    switch peakFuncOpt
        case 1
            adj_thresh = peak_thresh ./ exp(-0.25./thetaf(:,4).^2);
        case 2
            adj_thresh = peak_thresh ./ (pi/2.*thetaf(:,4).^2.*erf(1./(sqrt(2).*thetaf(:,4))).^2);
    end
    small_peaks = thetaf(:,1) < adj_thresh;
    
    thetaf(small_peaks,:) = [];
end

end

function val = gaussianPeaks(theta, yxdata)
% theta: [N x 4] with columns [amp, mu_x, mu_y, sigma]
% yxdata: [M x 2] with columns [x, y]
%
% val: [M x 1]

y = yxdata(:,1);          % Mx1
x = yxdata(:,2);          % Mx1

A  = theta(:,1).';        % 1xN
my = theta(:,2).';        % 1xN
mx = theta(:,3).';        % 1xN
s  = theta(:,4).';        % 1xN

inv2s2 = 1 ./ (2 * s.^2); % 1xN

dx = x - mx;              % MxN (implicit expansion)
dy = y - my;              % MxN

E = exp( -(dx.^2 + dy.^2) .* inv2s2 );  % MxN
val = E * A.';                           % (MxN)*(Nx1) -> Mx1
end

function val = gaussianPeaksIntegrated(theta, yxdata)
% Same interface as original:
% theta  : N×4  [A, mux, muy, sigma]
% yxdata : M×2  pixel centers [x, y]
% val    : M×1  integrated Gaussian over each pixel

x = yxdata(:,1);
y = yxdata(:,2);

% --- infer pixel size from grid ---
xu = unique(x);
yu = unique(y);

dx = median(diff(xu));
dy = median(diff(yu));

% build edges
xEdges = [xu(1)-dx/2; (xu(1:end-1)+xu(2:end))/2; xu(end)+dx/2];
yEdges = [yu(1)-dy/2; (yu(1:end-1)+yu(2:end))/2; yu(end)+dy/2];

W = numel(xu);
H = numel(yu);

% reshape index map
[~, xIdx] = ismember(x, xu);
[~, yIdx] = ismember(y, yu);

% --- parameters ---
A  = theta(:,1).';   % 1×N
mx = theta(:,2).';
my = theta(:,3).';
s  = max(theta(:,4).', eps);

c   = sqrt(pi/2);
rt2 = sqrt(2);

% --- integrate in x (W×N) ---
xL = xEdges(1:end-1);
xR = xEdges(2:end);
Ix = c .* s .* ( ...
    erf((xR - mx)./(rt2*s)) - ...
    erf((xL - mx)./(rt2*s)) );

% --- integrate in y (H×N) ---
yB = yEdges(1:end-1);
yT = yEdges(2:end);
Iy = c .* s .* ( ...
    erf((yT - my)./(rt2*s)) - ...
    erf((yB - my)./(rt2*s)) );

% --- combine to full image ---
img = (Iy .* A) * Ix.';   % H×W

% --- return in original point ordering ---
val = img(sub2ind([H W], yIdx, xIdx));
end
