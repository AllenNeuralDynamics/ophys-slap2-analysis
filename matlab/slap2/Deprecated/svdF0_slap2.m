function F0 = svdF0_slap2(IM, nPCs, denoiseWindow, baselineWindow)
%compute F0 for all pixels in IM, ignoring nans, using an SVD-based
%approach. This is good if individual pixels are noisy
%time is the first dimension
if nargin<2
    nPCs = 8;
end
if nargin<3
    denoiseWindow = 33;
end
if nargin<4
    baselineWindow = 101;
end

denoiseWindow = 2*floor(denoiseWindow/2)+1; %ensure denoiseWindow is odd
baselineWindow = 2*floor(baselineWindow/2)+1; %ensure baselineWindow is odd

origsz = size(IM);
IM = reshape(IM, origsz(1), []);
nanframes = all(isnan(IM),2);

IMsmooth = smoothdata(IM, 1, 'movmean',denoiseWindow, 'omitnan');
selPix = all(~isnan(IMsmooth(~nanframes,:)),1);

IMsel = IMsmooth;
%IMsel = IM; IMsel(isnan(IMsel)) = IMsmooth(isnan(IMsel));
IMsel = IMsel(~nanframes, selPix);

% IMsel = smoothdata(IMsel, 'movmean', denoiseWindow, 'omitnan');
% IMsel = smoothExp(IMsel, 'movmedian', baselineWindow);

IMsel = IMsel-mean(IMsel,1);

[UU,~,VV] = svds(IMsel',nPCs);
spaceFac = mean(abs(UU),1)./sqrt(mean(UU.^2,1));
selU = (spaceFac>0.6) | (spaceFac>mean(spaceFac)); %select components that are spatially distributed

warning('off','stats:statrobustfit:IterationLimit');
Y = IM(~nanframes,:);
for px = 1:size(IM,2)
    y = Y(:,px);
    selT = ~isnan(y);
    if sum(selT)>2*sum(selU) && sum(selT)>(2*sum(selU))
        %regular regression
        % b = regress(y(selT), [VV(selT,selU) ones(sum(selT),1)]);
        % Y(selT, px) = min(IMsmooth(selT,px), [VV(selT,selU) ones(sum(selT),1)]*b);

        %iteratively reweighted least squares
        b2 = robustfit([VV(selT,selU) ones(sum(selT),1)],y(selT),[],2,'off') ;
        Y(selT, px) = min(IMsmooth(selT,px), [VV(selT,selU) ones(sum(selT),1)]*b2);
    else
        Y(:,px) = nan; %very few timepoints valid
    end
end

F0 = nan(size(IM));
F0(~nanframes,:) = max(Y,0);
F0(isnan(IM)) = nan;
F0 = reshape(F0, origsz);
end