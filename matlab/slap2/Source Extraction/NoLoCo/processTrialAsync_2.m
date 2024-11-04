function exptSummary = processTrialAsync_2(dr, fnRaw, startLine, endLine, W0, F0selDS, selPix, discardFrames, alignData, meanIM, motOutput, roiData, params)
    fn = fnRaw;
    disp('Loading high-res data for file:')
    disp([dr filesep fn])

    %load the high time resolution data
    S2data = slap2.Slap2DataFile([dr filesep fn]);
    meta = loadMetadata([dr filesep fn]);
    numChannels = S2data.numChannels;
    linerateHz = 1/meta.linePeriod_s;
    dt = linerateHz/params.analyzeHz;
    frameLines = ceil(startLine:dt:endLine);
    nFrames= length(frameLines);
    selPx2D = any(selPix,3);
    sz = size(selPx2D);

    %upsample motion
    motionC = interp1(alignData.DSframes, alignData.motionDSc, frameLines, 'pchip', 'extrap') + motOutput(2);
    motionR = interp1(alignData.DSframes, alignData.motionDSr, frameLines, 'pchip', 'extrap') + motOutput(1);
    
    labeled = medfilt2(meanIM(:,:,1), [3 3]);
    labeled = ~isnan(meanIM(:,:,1)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
    for cix = numChannels:-1:1
        tmp = double(meanIM(:,:,cix));
        mLabeled{cix} = tmp(labeled);
    end

    Ysz = [length(alignData.trimRows) length(alignData.trimCols)];
    selPx2D = selPx2D(1:size(alignData.viewC,1), 1:size(alignData.viewC,2));

    %interpolate the raw data at the shifted coordinates
    for fix = nFrames:-1:1
        for cix = 1:numChannels
            Y = S2data.getImage(cix, frameLines(fix), ceil(dt), 1);
            Y = Y(alignData.trimRows, alignData.trimCols);
            Y = interp2(1:Ysz(2), 1:Ysz(1), Y,alignData.viewC+motionC(fix), alignData.viewR+motionR(fix), 'linear', nan);
            if cix==1
                IMsel(:, fix) = Y(selPx2D);
            elseif cix==2
                IM2sel(:, fix) = Y(selPx2D);
            end

            %compute global ROI activity
            mIM = meanIM(:,:,cix);
            yLabeled = double(Y(labeled)); nans= isnan(yLabeled);
            FF = (sum(yLabeled(~nans))./sum(mLabeled{cix}(~nans))).*sum(mLabeled{cix});
            exptSummary.global.F(fix,cix) = FF;

            %compute user ROI activity
            for rix = 1:length(roiData)
                mask = roiData{rix}.mask;
                tmp1 = Y(mask); tmp2 = mIM(mask);   tmp2(isnan(tmp2)) = 0;
                nans= isnan(tmp1);
                Fpx{rix}(:,fix,cix) = tmp1;
                FF = (sum(tmp1(~nans))./sum(tmp2(~nans))).*sum(tmp2);
                exptSummary.ROIs.F(rix, fix,cix) = FF;
            end

            % if fix<300 && cix==1
            %     sanitycheckTMP(:,:,301-fix) = Y;
            % end
        end
    end

    %perform SVD on user ROIs to denoise
    exptSummary.ROIs.Fsvd = nan(length(roiData), nFrames, numChannels);
    for rix = 1:length(roiData)
        for cix = 1:numChannels 
            Dtmp = double(Fpx{rix}(:,:,cix));
            nanFrames = all(isnan(Dtmp),1);
            nanPx = any(isnan(Dtmp(:, ~nanFrames)),2);
            Dtmp = Dtmp(~nanPx, ~nanFrames);
            [UU,SS,VV] = svds(Dtmp,1);
            exptSummary.ROIs.Fsvd(rix,~nanFrames,cix) = sum(UU)*SS*VV';
        end
    end
    clear Fpx;
    
    discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,2)))>0; %1,params.dsFac)', 1,[]); %upsample the discard frames
    IMsel(:,discard) = nan;     %throw away movement frames as above
    exptSummary.global.F(discard,:) = nan;
    exptSummary.ROIs.F(:, discard,:) = nan;
    exptSummary.ROIs.Fsvd(:,discard,:) = nan;

    [IMselFilt, IMselRaw, F0sel, W, selNans] = prepareNMFproblem(IMsel, W0, F0selDS, params);
    
    %least squares solve the match filtered movie
    H0 = W\IMselFilt;

    %W = W1;

    
    %fit F0
    f0 = W\F0sel;
    F0_H = computeF0(H0', ceil(params.denoiseWindow_s*params.analyzeHz), ceil(params.baselineWindow_Glu_s*params.analyzeHz), 1)';
    H = H0-F0_H;
    f0 = f0+F0_H;

    Hsub = 3*std(H,0,2, 'omitmissing');
    H = H+Hsub;
    % %NMF
    % H0 = ones(size(W1,2), size(IMselFilt,2));
    % for iter = 1:3 %perform nonnegative matrix division to initialize H0 with W0 constant
    %     H0 = H0.*(max(0,W1'*IMselFilt))./((W1'*W1)*(H0 + mean(H0(:)/100))); %confirm this is right; per Lee and Seung NIPS 2000 'Algorithms for nonnegative matrix factorization'
    % end
    % %'h0', H0
    % H0 = max(max(H0(:)/100), H0 + max(H0(:)/100));
    % opts = statset('MaxIter', 3,  'Display', 'final');
    % [W,H] = nnmf2(IMselFilt, size(W1,2),'algorithm', 'mult', 'w0', W1, 'h0', H0, 'options', opts); %!!nnmf2 has been modified to keep the ordering of the provided factors
    Wfull = nan([sz(1)*sz(2) size(W,2)]);
    Wfull(any(selPix,3),:) = W;
    Wfull = reshape(Wfull, sz(1),sz(2), []);

    %Visualize the difference between the initialized footprints and fit
    %footprints
    % W1full = nan([sz(1)*sz(2) size(W1,2)]);
    % W1full(any(selPix,3),:) = W1;
    % W1full = reshape(W1full, sz(1),sz(2), []);

    % figure('name', 'original'), imshow3D(Wfull./sum(Wfull,[1 2], 'omitnan'));
    % figure('name', 'refit'), imshow3D(W1full./sum(W1full,[1 2], 'omitnan'));
    % figure('name', 'difference'), imshow3D(1000*((Wfull./sum(Wfull,[1 2], 'omitnan'))-(W1full./sum(W1full,[1 2], 'omitnan'))));

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
    
    %Compute a 'least squares' dF
    H5 = W \ IMselRaw;
    F0_H5 = computeF0(H5', ceil(params.denoiseWindow_s*params.analyzeHz), ceil(params.baselineWindow_Glu_s*params.analyzeHz), 1)';
    H5 = H5-F0_H5;

    %compute F0
    setNan = imdilate(nanFramesH, ones(1,2*floor(params.tau_full)+1));

    % F0mean = repmat(mean(F0sel,2, 'omitnan'),1,size(F0sel,2));
    % F0sel(~isfinite(F0sel)) = max(0, F0mean(~isfinite(F0sel)));
    % F0= (W./max(W,[],1))' *F0sel + sum(W,1)'.*F0_H;%(W./max(W,[],1))' *F0sel;
    F0 = max(f0,0) + median(f0(f0>0))/4; %some regularization and ensure F0 is positive
    F0(setNan) = nan;
    
    % %ensure F0 is positive
    % Fmin = min(F0(:));
    % desMin = prctile(F0(:),1) - Fmin;
    % if Fmin<desMin
    %     F0 = F0 + desMin - Fmin;
    % end

    %NaN out invalid data
    H(setNan) = nan; %The match filtered signal
    H2(setNan) = nan; %The detected events
    H3(setNan) = nan; %The denoised activity
    H4(setNan) = nan; %raw Fluorescence estimate
    H5(setNan) = nan; %raw Fluorescence estimate

    %exptSummary.dFerr = sum(W,1)'.*errH;
    exptSummary.matchFilt(:,:,1) = sum(W,1)'.*H; %[source#, time, channel]
    exptSummary.events(:,:,1) = sum(W,1)'.*H2; %[source#, time, channel]
    exptSummary.denoised(:,:,1) = sum(W,1)'.*H3; %[source#, time, channel]
    exptSummary.dFraw(:,:,1) = sum(W,1)'.*H4; %[source#, time, channel]
    exptSummary.dFls(:,:,1) = sum(W,1)'.*H5; %[source#, time, channel]
    exptSummary.F0(:,:,1) = F0;
    exptSummary.footprints = Wfull;
    exptSummary.discardFrames = discard;
    %exptSummary.sanityCheck = mean(sanitycheckTMP,3, 'omitnan');

    %compute channel 2 signals
    if numChannels==2
        for rix = size(W,2):-1:1
            sel = W(:,rix)>0;
            Wr = W(sel,rix)./max(W(sel,rix));
            F_2(rix,:) = Wr'*IM2sel(sel,:);
        end
        F0_2 = computeF0(F_2', ceil(params.denoiseWindow_s*params.analyzeHz), ceil(params.baselineWindow_Ca_s*params.analyzeHz),1)';
        exptSummary.dFraw(:,:,2) = F_2;
        exptSummary.F0(:,:,2) = F0_2;
    end

    exptSummary.dFF = exptSummary.dFraw./exptSummary.F0;
    exptSummary.dFFls = exptSummary.dFls./exptSummary.F0;    
end