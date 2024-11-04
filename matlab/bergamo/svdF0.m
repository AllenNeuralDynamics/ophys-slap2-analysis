function F0 = svdF0(IM, nPCs, baselineWindow, algo, denoiseWindow)
%compute F0 for all pixels in IM, ignoring nans, using an SVD-based
%approach. This is good if individual pixels are noisy
%time is the first dimension
if nargin<2
    nPCs = 4;
end
if nargin<3
    baselineWindow = 101;
end
if nargin<4
    algo = 2;
end
if nargin<5
    denoiseWindow = 35;
end
origsz = size(IM);
IM = reshape(IM, origsz(1), []);
nanframes = all(isnan(IM),2);
selPix = mean(~isnan(IM(~nanframes,:)),1)>0.9;

IMsel = IM(~nanframes, selPix);

IMsmooth = smoothdata(IMsel,1, 'movmean', 2*round(size(IMsel,1)/8)+1, 'omitnan');
IMsel(isnan(IMsel)) = IMsmooth(isnan(IMsel));
IMsel(isnan(IMsel)) = 0; %anything left over

IMsel = smoothExp(IMsel, 'movmedian', baselineWindow);
% IMsel = IMsel-mean(IMsel,1);

%%
if algo == 1
    rng(1);
    k = nPCs;
    opts1 = statset('MaxIter', 5, 'Display','final');%, 'UseParallel', true);
    [VV, UU] = nnmf2(IMsel,k,'Algorithm','mult','Options',opts1);
    VV_filt = max(0,computeF0(VV,denoiseWindow, baselineWindow));

    for iters = 1:10
        [VV,UU] = nnmf2(IMsel,k,'w0',VV_filt,'h0',UU,'Options',opts1);

        keepU = mean(abs(UU),2)./sqrt(mean(UU.^2,2));
        keepU = keepU > 0.5;
        k = sum(keepU);

        UU = UU(keepU,:);
        VV_filt = max(0,computeF0(VV(:,keepU),denoiseWindow, baselineWindow));
    end

    warning("off", 'stats:statrobustfit:IterationLimit')
    Y = IM(~nanframes,:);
    for px = 1:size(IM,2)
        y = Y(:,px);
        selT = ~isnan(y);
        if sum(selT)>2*size(VV_filt,2) && sum(selT)>5
            %use iteratively reweighted least squares to fit
            b = robustfit(VV_filt(selT,:), y(selT), 'bisquare', 2.5, true);
            Y(selT, px) = [ones(sum(selT),1) VV_filt(selT,:)  ]*b;
        else
            Y(selT,px) = nan; %very few timepoints valid
        end
    end
    nnF0 = nan(size(IM));
    nnF0(~nanframes,:) = Y;
    nnF0 = reshape(nnF0, origsz);

    resIM = IM - nnF0;

    resIMsel = resIM(~nanframes, selPix);

    IMsmooth = smoothdata(resIMsel,1, 'movmean', size(resIMsel,1)/4, 'omitnan');
    resIMsel(isnan(resIMsel)) = IMsmooth(isnan(resIMsel));
    resIMsel(isnan(resIMsel)) = 0; %anything left over

    resIMsel = smoothExp(resIMsel, 'movmedian', baselineWindow);
    resIMsel = resIMsel-mean(resIMsel,1);

    [UU,~,VV] = svds(resIMsel',nPCs);
    selU = mean(abs(UU),1)./sqrt(mean(UU.^2,1));
    selU = selU > 0.5;

    warning("off", 'stats:statrobustfit:IterationLimit')
    Y = resIM(~nanframes,:);
    for px = 1:size(resIM,2)
        y = Y(:,px);
        selT = ~isnan(y);
        if sum(selT)>2*sum(selU) && sum(selT)>5
            %use iteratively reweighted least squares to fit
            b = robustfit(VV(selT,selU), y(selT), 'bisquare', 2.5, true);
            Y(selT, px) = [ones(sum(selT),1) VV(selT,selU)  ]*b;
        else
            Y(selT,px) = nan; %very few timepoints valid
        end
    end
    pcaF0 = nan(size(resIM));
    pcaF0(~nanframes,:) = Y;
    %F0 = computeF0(F0, denoiseWindow,baselineWindow,2);
    pcaF0 = reshape(pcaF0, origsz);

    F0 = pcaF0 + nnF0;

elseif algo == 2
    IMsel = IMsel-mean(IMsel,1);
    [UU,~,VV] = svds(IMsel',nPCs);
    selU = mean(abs(UU),1)./sqrt(mean(UU.^2,1));
    selU = selU>mean(selU); %0.67; %select only components that are spatially distributed

    VV2 = smoothExp(VV, 'movmedian', baselineWindow);

    warning("off", 'stats:statrobustfit:IterationLimit')
    Y = IM(~nanframes,:);
    for px = 1:size(IM,2)
        y = Y(:,px);
        selT = ~isnan(y);
        if sum(selT)>2*sum(selU) && sum(selT)>5
            %use iteratively reweighted least squares to fit
            b = robustfit(VV2(selT,selU), y(selT), 'bisquare', 2.5, true);
            Y(selT, px) = [ones(sum(selT),1) VV2(selT,selU)  ]*b;
        else
            Y(selT,px) = nan; %very few timepoints valid
        end
    end
    F0 = nan(size(IM));
    F0(~nanframes,:) = Y;
    %F0 = computeF0(F0, denoiseWindow,baselineWindow,2);
    F0 = reshape(F0, origsz);
end

end

