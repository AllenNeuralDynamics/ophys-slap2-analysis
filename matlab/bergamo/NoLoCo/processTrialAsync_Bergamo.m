function exptSummary = processTrialAsync_Bergamo(dr, fnRaw, ~, ~, W0, F0selDS, selPix, discardFrames, ~, ~, motOutput, ~, params)
    fn = fnRaw;
    activityChannel = params.activityChannel;
    numChannels = params.numChannels;

    disp('Loading high-res data for file:')
    disp([dr filesep fn])

    %load the high time resolution tiff
    fn = fnRaw;
    [IM, desc, meta] = networkScanImageTiffReader([dr filesep fn]);
    IM = double(IM);

    selPx2D = any(selPix,3);
    sz = size(selPx2D);

    %rearrange IM into correct dimensions
    IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []);

    if activityChannel>1
        IM = IM(:,:,[activityChannel:end, 1:activityChannel-1],:);
        disp('Reordering channels for analysis!')
    end

    IM1 = squeeze(IM(:,:,1,:));
    if numChannels==2
        IM2 =  squeeze(IM(:,:,2,:));
        clear IM;
        IM2sel = interpArray(IM2, any(selPix,3), motOutput); %interpolate the movie at the shifted coordinates
        clear IM2;
        IM2sel = IM2sel - min(0, min(mean(IM2sel,2, 'omitnan')));%ensure that the baseline is not overestimated
        discard = reshape(repmat(discardFrames{trialIx}(:), 1,params.dsFac)', 1,[]); %upsample the discard frames
        IM2sel(:,discard) = nan;     %throw away movement frames as above
    else %1 channel
        clear IM;
    end 

    IMsel = interpArray(IM1, any(selPix,3), motOutput); %interpolate the movie at the shifted coordinates
    
    % %compute global ROI activity
    % mIM = meanIM(:,:,cix);
    % yLabeled = double(Y(labeled)); nans= isnan(yLabeled);
    % FF = (sum(yLabeled(~nans))./sum(mLabeled{cix}(~nans))).*sum(mLabeled{cix});
    % exptSummary.global.F(fix,cix) = FF;
    %IMsel = IMsel - min(0, min(mean(IMsel,2, 'omitnan'))); %ensure that the baseline is not overestimated
    
    discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,2)))>0; %1,params.dsFac)', 1,[]); %upsample the discard frames
    IMsel(:,discard) = nan;     %throw away movement frames as above
    %exptSummary.global.F(discard,:) = nan;

    [IMselFilt, IMselRaw, F0sel,  W, selNans] = prepareNMFproblem(IMsel, W0, F0selDS, params);
    
    %least squares solve the match filtered movie
    H0 = W\IMselFilt;
    
    %fit F0
    f0 = W\F0sel;
    F0_H = computeF0(H0', ceil(params.denoiseWindow_s*params.analyzeHz), ceil(params.baselineWindow_Glu_s*params.analyzeHz), 1)';
    H = H0-F0_H;
    f0 = f0+F0_H;

    Hsub = 3*std(H,0,2, 'omitmissing');
    H = H+Hsub;

    Wfull = nan([sz(1)*sz(2) size(W,2)]);
    Wfull(any(selPix,3),:) = W;
    Wfull = reshape(Wfull, sz(1),sz(2), []);

    %Frames where more than 25% of pixels in a source are nan
    nanFramesH = false(size(H));
    for sourceIx = 1:size(H,1)
        support = W(:,sourceIx)>0;
        nanFramesH(sourceIx,:) = mean(selNans(support,:)) > 0.25;
    end
    nanFramesH(:,end+(-floor(params.tau_full):0)) = 1;
    deconvWeights = 1-double(nanFramesH);

    %deconvolve out matched filter we applied earlier
    kernel = [zeros(1,ceil(8*params.tau_full)) exp(-(0:ceil(8*params.tau_full))/params.tau_full)];
    kernel = kernel./sum(kernel);
    fKernel = fliplr(kernel);
    doubleKernel = conv(kernel, fKernel, 'same');
    doubleKernel = doubleKernel./sum(doubleKernel);
    
    %perform deconvolution, filling in NaNs with reconstructed values every
    %few iterations:
    J = deconvlucy({H},doubleKernel, 20,[], deconvWeights);
    Js = deconvlucy({H},fKernel, 10,[],deconvWeights);
    for iter = 1:3
            recon = convn(J{2}, doubleKernel, 'same');
            J{1}(nanFramesH) = recon(nanFramesH);
            J = deconvlucy(J,doubleKernel, 20, [], deconvWeights);

            recon2 =  convn(Js{2}, fKernel, 'same');
            Js{1}(nanFramesH) = recon2(nanFramesH);
            Js = deconvlucy(Js,fKernel, 10, [], deconvWeights);
    end
    H2 = J{2} - Hsub; %an estimate of the underlying spikes; deconvolving out the double kernel from the match-filtetred data
    H3 = convn(H2, kernel, 'same'); %the denoised trace, composed of the events convolved witht he forward kernel
    H4 = Js{2}-Hsub; %an estimate of raw fluorescence (no denoising), deconvolving out the matched kernel that was applied prior to NMF
    H5 = W \ IMselRaw;
    %assess residuals
    %assessResiduals(W*H3, IMsel, selPx2D)
    %assess F0

    %compute F0
    setNan = imdilate(nanFramesH, ones(1,2*floor(params.tau_full)+1));

    % F0mean = repmat(mean(F0sel,2, 'omitnan'),1,size(F0sel,2));
    % F0sel(~isfinite(F0sel)) = max(0, F0mean(~isfinite(F0sel)));
    % F0= (W./max(W,[],1))' *F0sel + sum(W,1)'.*F0_H;%(W./max(W,[],1))' *F0sel;
    F0 = max(f0,0) + median(f0(f0>0))/4; %some regularization and ensure F0 is positive
    F0(setNan) = nan;

    %NaN out invalid data
    H(setNan) = nan; %The match filtered signal
    H2(setNan) = nan; %The detected events
    H3(setNan) = nan; %The denoised activity
    H4(setNan) = nan; %raw Fluorescence estimate
    H5(setNan) = nan; %raw Fluorescence estimate, no matched filter

    exptSummary.matchFilt(:,:,1) = sum(W,1)'.*H; %[source#, time, channel]
    exptSummary.events(:,:,1) = sum(W,1)'.*H2; %[source#, time, channel]
    exptSummary.dF(:,:,1) = sum(W,1)'.*H3; %[source#, time, channel]
    exptSummary.dF2(:,:,1) = sum(W,1)'.*H4;
    exptSummary.dFls(:,:,1) = sum(W,1)'.*H5;
    exptSummary.F0(:,:,1) = F0;
    exptSummary.footprints(:,:,1:size(W,2)) = Wfull;
    
    % %compute channel 2 signals
    % if numChannels==2
    %     IM2sel(isnan(IM2sel)) = 0;
    %     F_2 = (W./max(W,[],1))' * IM2sel;
    %     F_2(nanFramesH) = nan;
    % 
    %     F0_2 = nan(size(IM2sel));
    %     for rix = 1:size(F0selDS{trialIx}, 1)
    %         F0_2(rix,:) = interp(F02selDS(rix,:), params.dsFac)./params.dsFac;
    %     end
    % 
    %     F02mean = repmat(mean(F0_2,2, 'omitnan'),1,size(F0_2,2));
    %     F0_2(~isfinite(F0_2)) = max(0, F02mean(~isfinite(F0_2)));
    %     F0_2 = (W./max(W,[],1))' *F0_2;
    %     F0_2(nanFramesH) = nan;
    % 
    %     exptSummary.dF(:,:,2) = F_2 - F0_2;
    %     exptSummary.dF2(:,:,2) = F_2 - F0_2;
    %     exptSummary.dFls(:,:,2) = F_2 - F0_2;
    %     exptSummary.F0(:,:,2) = F0_2;
    % end

    exptSummary.dFF = exptSummary.dF./exptSummary.F0;
    exptSummary.dFF2 = exptSummary.dF2./exptSummary.F0;
    exptSummary.dFFls = exptSummary.dFls./exptSummary.F0;
end

function [IMsel, IMselRaw, F0,  W, selNans] = prepareNMFproblem(IMsel, W0, F0selDS, params)
    F0 = nan(size(IMsel));
    for rix = 1:size(F0selDS, 1)
        F0(rix,:) = interp(F0selDS(rix,:), params.dsFac)./params.dsFac;
    end
    IMsel = IMsel - F0;
    selNans = isnan(IMsel);
    IMsel(selNans) = 0;
    IMselRaw = IMsel;
    IMsel = matchedExpFilter(IMsel, params.tau_full);
    %IMsel = IMsel- computeF0(IMsel', params.denoiseWindow_samps*params.dsFac+1,baselineWindow, 1)';
    W = W0;
end