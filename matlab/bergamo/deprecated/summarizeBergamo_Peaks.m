function summarizeBergamo_Peaks(dr, fns, activityChannel)
    %TO DO:
    %PARAMETER SENSITIVITY ANALYSIS

    %add some global NMF components?
        %initialize with random values

%cd('\\allen\aind\scratch\ophys\Adrian\iGluSnFR testing\datasets\699943\10-11-23\scans\scan_00002_20231011_115018')
%cd('C:\temp\SYNAPSES\Test');
%cd('C:\Users\kaspar.podgorski\OneDrive - Allen Institute\Documents\GitHub\ophys-slap2-analysis\matlab\rasterDataProcessing\Bergamo\simulations\data\SIMULATIONS')

params.tau_s = 0.027; % time constant in seconds for glutamate channel; from Aggarwal et al 2023 Fig 5
params.sigma_px = 1.33;   % space constant in pixels
params.eventRateThresh_hz = 1/10; % minimum event rate in Hz
params.sparseFac = 0.1; %sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser
params.nmfIter = 5; %number of iterations of NMF refinement
params.dXY = 3; %how large sources can be (radius), pixels
params.upsample = 3; %how many times to upsample the imaging resolution for finding local maxima to identify sources; affects maximum source density
params.nmfBackgroundComps = 0; % <=4, max number of background components to use for NMF. If 0, we compute F0 instead of fitting background
params.denoiseWindow_samps = 35; %number of samples to average together for denoising
params.baselineWindow_Glu_s = 2; %timescale for calculating F0 in glutamate channel, seconds
params.baselineWindow_Ca_s = 2; %timescale for calculating F0 in calcium channel, seconds

if nargin == 0
    %select a set of aligned downsampled recordings (trials)
    [fns, dr] = uigetfile('*REGISTERED*.tif', 'multiselect', 'on');
end
if nargin<3
    activityChannel = 2;
end

if ~iscell(fns)
    fns = {fns};
end

fnStemEnd = strfind(fns{1}, '_REGISTERED') -1;
fnStem = fns{1}(1:fnStemEnd);
savedr = dr;
% savedr = [dr filesep 'ExperimentSummary'];
% if ~exist(savedr, 'dir')
%     mkdir(savedr);
% end
if nargin == 0
    [fnsave, drsave] = uiputfile([savedr filesep 'Summary.mat'], 'Set filename for saving summary for this condition');
else
    fnsave = [fnStem '_EXPTSUMMARY.mat'];
    drsave = savedr;
end
%load some metadata
load([dr filesep fnStem '_ALIGNMENTDATA.mat'], 'aData');
numChannels = aData.numChannels;

%generate a concensus alignment across trials for further analysis
meanIM = nan(1,1,numChannels,1);
actIM = nan(1,1,numChannels,1);
disp('Loading data and performing localizations...')
for trialIx = length(fns):-1:1
    fn = fns{trialIx};

    %ensure the high res file exists
    ind =strfind(fn, '_DOWNSAMPLED');
    fnRaw{trialIx} = [fn(1:ind) 'RAW.tif'];
    assert(exist([dr filesep fnRaw{trialIx}], 'file'), ['No corresponding RAW tiff recording found for file:' fn])

    %load the tiff
    [IM, desc, meta] = networkScanImageTiffReader([dr filesep fn]);
    IM = double(IM);
    % A = ScanImageTiffReader([dr filesep fn]);
    % IM = double(A.data);
    if size(IM,3)<100
        error(['The file:' fn 'is very short. You should probably not include it?']);
    end
    IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;
    
    %Reorder channels so the activity channel is first
    if activityChannel>1
        IM = IM(:,:,[activityChannel:end, 1:activityChannel-1],:);
        disp('Reordering channels for analysis!')
    end

    meanIM(end:size(IM,1),:,:,:) = nan;
    meanIM(:, end:size(IM,2),:,:) = nan;
    meanIM(:,:,:,trialIx) = nan;
    meanIM(1:size(IM,1),1:size(IM,2),:,trialIx) = mean(IM,4, 'omitnan');

    %load alignment data
    fnStemEnd = strfind(fn, '_REGISTERED') -1;
    load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
    aData.dsFac = round(length(aData.motionC)./length(aData.motionDSc));
    params.dsFac = aData.dsFac;
    params.frametime = aData.frametime;

    %discard motion frames
    % tmp = aData.aRankCorr(:)-smoothExp(aData.aRankCorr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorr,2, 'movmedian', ceil(2/aData.frametime));
    % discardFrames{trialIx} = imdilate(tmp<-(4*std(tmp)), ones(1,5));
    tmp = aData.recNegErr(:)- median(aData.recNegErr(:)); %smoothExp(aData.recNegErr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
    tmp = tmp./(median(aData.recNegErr(:)) - prctile(aData.recNegErr(:),1)); 
    thresh = 3; %decrease thresh to be more stringent on motion correction
    window = 2*ceil(0.1/params.frametime)+1;% a window in time to censor aronud motion events, ~0.2 seconds;
    discardFrames{trialIx} = imclose(imdilate(tmp>thresh, ones(window,1)) | (tmp>(thresh/2) & imdilate(tmp>thresh, ones(2*window+1,1))), ones(window,1));

    % tmp = aData.aRankCorr(:)-smoothExp(aData.aRankCorr(:),'movmedian', ceil(10/(aData.frametime*aData.dsFac)));
    % filtTmp = smoothExp(tmp, 'movmean',ceil(.2/(aData.frametime*aData.dsFac)));
    % discardFrames{trialIx} = imdilate(filtTmp<-(4*std(filtTmp)), ones(1,5));
    rawIMs{trialIx} = squeeze(IM(:,:,1,:));
    rawIMs{trialIx}(:,:,discardFrames{trialIx}) = nan;
    if numChannels==2
        rawIM2s{trialIx} = squeeze(IM(:,:,2,:));
        rawIM2s{trialIx}(:,:,discardFrames{trialIx}) = nan;
    end

    [IMc, peaks(trialIx), params] = localizeFlashesBergamo(rawIMs{trialIx}, aData, params);

    %calculate correlation image
    actIM(end:size(IMc,1),:,:,:) = nan;
    actIM(:,end:size(IMc,2),:,:) = nan;
    actIM(:,:,:,trialIx) = nan;
    actIM(1:size(IMc,1),1:size(IMc,2),1,trialIx) = IMc;
end
params.sz = size(meanIM, [1 2]);

clear aData

%Make template
disp('Making template for aligning across trials...')
maxshift = 12;
M = squeeze(sum(meanIM, 3));
template = makeTemplateMultiRoi(M, maxshift);

%align all mean images to template
disp('Aligning across trials...')
Mpad = nan([size(template) size(M,3)]);
Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
for trialIx = length(fns):-1:1
    disp(['trial: ' int2str(trialIx)])
    mot1 = xcorr2_nans(Mpad(:,:,trialIx), template, [0 ; 0], maxshift);
    [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
    [rr,cc] = ndgrid(1:size(meanIM,1), 1:size(meanIM,2));
    for chIx = 1:size(meanIM,3)
        meanAligned(:,:,chIx,trialIx) = interp2(meanIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
        actAligned(:,:,chIx,trialIx) = interp2(actIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
    end
    if trialIx==length(fns)
        peaksCat.row = peaks(trialIx).row - motOutput(1,trialIx);
        peaksCat.col = peaks(trialIx).col - motOutput(2,trialIx);
        peaksCat.val = peaks(trialIx).val;
        peaksCat.t = peaks(trialIx).t;
         figure, imagesc(template)
         hold on, scatter( maxshift+ peaks(trialIx).col - motOutput(2,trialIx), maxshift+ peaks(trialIx).row - motOutput(1,trialIx));
    else
        peaksCat.row = cat(2, peaks(trialIx).row - motOutput(1,trialIx),peaksCat.row);
        peaksCat.col = cat(2, peaks(trialIx).col - motOutput(2,trialIx), peaksCat.col);
        peaksCat.val = cat(1, peaks(trialIx).val, peaksCat.val);
        peaksCat.t = peaksCat.t + size(rawIMs{trialIx},3);
        peaksCat.t = cat(1, peaks(trialIx).t,peaksCat.t);
         hold on, scatter( maxshift+ peaks(trialIx).col - motOutput(2,trialIx), maxshift+ peaks(trialIx).row - motOutput(1,trialIx));
    end
end

%identify outliers in alignment quality
cc = corrCoeff;
corrThresh = min(0.96, median(cc)-2*std(cc));

validTrials= find(cc>corrThresh);
exptSummary.meanIM = mean(meanAligned,4, 'omitnan');
exptSummary.actIM = mean(actAligned, 4, 'omitnan');

%cluster localizations
totalFrames = sum(cellfun(@(x)(size(x,3)), rawIMs));
params.minEvents = totalFrames*params.frametime*params.dsFac*params.eventRateThresh_hz;
[~, P, sources.R, sources.C, params] = clusterLocalizations(peaksCat, params);
k = length(sources.R); %number of sources
sz = params.sz;

%Generate IMsel; the data only in the selected region, aligned across movies
selPix = false(sz(1:2));
for sourceIx = k:-1:1
    rr = max(1, round(sources.R(sourceIx)-params.dXY)):min(sz(1),  round(sources.R(sourceIx))+params.dXY);
    cc = max(1, round(sources.C(sourceIx)-params.dXY)):min(sz(2),  round(sources.C(sourceIx))+params.dXY);
    selPix(rr,cc,sourceIx) = true;
end
nSelPix = sum(any(selPix,3), 'all'); %number of selected pixels
dFsel = nan(nSelPix, totalFrames); %extize
if numChannels == 2
    dF2sel = nan(nSelPix, totalFrames); %extize
end
baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));
frameInd = 0;
for trialIx = validTrials
    szTmp = size(rawIMs{trialIx});
    IMrawSel = interpArray(rawIMs{trialIx}, any(selPix,3), motOutput(:,trialIx)); %interpolates the movie at the shifted coordinates
    F0selDS{trialIx} = svdF0(IMrawSel', 5, baselineWindow, 1, params.denoiseWindow_samps)'; %#ok<AGROW>
    dFsel(:,frameInd+(1:szTmp(3))) = IMrawSel - F0selDS{trialIx};
    if numChannels == 2
        IM2rawSel = interpArray(rawIM2s{trialIx}, any(selPix,3), motOutput(:,trialIx)); %interpolates the movie at the shifted coordinates
        F02selDS{trialIx} = svdF0(IM2rawSel', 5, baselineWindow, 1, params.denoiseWindow_samps)';
        dF2sel(:,frameInd+(1:szTmp(3))) = IM2rawSel - F02selDS{trialIx};
    end
    frameInd = frameInd+szTmp(3);
end
clear rawIMs

%extract sources from the downsampled movies
[W0,~] = extractSourcesLoRes(dFsel, P, sources, selPix, params);

%for each file, load high res data and refine
params.tau_full=params.tau_frames*params.dsFac;
exptSummary.params = params;
for trialIx = validTrials
    fn = fnRaw{trialIx};
    
    %load the high time resolution tiff
    [IM, desc, meta] = networkScanImageTiffReader([dr filesep fn]);
    IM = double(IM);
    % A = ScanImageTiffReader([dr filesep fn]);
    % IM = double(A.data);

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
        IM2sel = interpArray(IM2, any(selPix,3), motOutput(:,trialIx)); %interpolate the movie at the shifted coordinates
        clear IM2;
        IM2sel = IM2sel - min(0, min(mean(IM2sel,2, 'omitnan')));%ensure that the baseline is not overestimated
        discard = reshape(repmat(discardFrames{trialIx}(:), 1,params.dsFac)', 1,[]); %upsample the discard frames
        IM2sel(:,discard) = nan;     %throw away movement frames as above
    else %1 channel
        clear IM;
    end 

    IMsel = interpArray(IM1, any(selPix,3), motOutput(:,trialIx)); %interpolate the movie at the shifted coordinates
    %IMsel = IMsel - min(0, min(mean(IMsel,2, 'omitnan'))); %ensure that the baseline is not overestimated
    
    discard = reshape(repmat(discardFrames{trialIx}(:), 1,params.dsFac)', 1,[]); %upsample the discard frames
    IMsel(:,discard) = nan;     %throw away movement frames as above

    [IMsel, IMselRaw, F0sel, W1, selNans] = prepareNMFproblem(IMsel, W0, F0selDS{trialIx}, params);

    %NMF
    H0 = ones(size(W1,2), size(IMsel,2));
    for iter = 1:3 %perform nonnegative matrix division to initialize H0 with W0 constant
        H0 = H0.*(max(0,W1'*IMsel))./((W1'*W1)*(H0 + mean(H0(:)/100))); %confirm this is right; per Lee and Seung NIPS 2000 'Algorithms for nonnegative matrix factorization'
    end
    opts = statset('MaxIter', 10,  'Display', 'final');
    [W,H] = nnmf2(IMsel, size(W1,2),'algorithm', 'mult', 'w0', W1, 'h0', H0, 'options', opts); %!!nnmf2 has been modified to keep the ordering of the provided factors
    Wfull = nan([sz(1)*sz(2) size(W,2)]);
    Wfull(any(selPix,3),:) = W;
    Wfull = reshape(Wfull, sz(1),sz(2), []);
    
    %Frames where more than 25% of pixels in a source are nan
    nanFramesH = false(size(H));
    for sourceIx = 1:size(H,1)
        support = W(:,sourceIx)>0;
        nanFramesH(sourceIx,:) = mean(selNans(support,:)) > 0.25;
    end

    %deconvolve out matched filter we applied earlier
    kernel = [zeros(1,ceil(8*params.tau_full)) exp(-(0:ceil(8*params.tau_full))/params.tau_full)];
    kernel = kernel./sum(kernel);
    doubleKernel = conv(kernel, fliplr(kernel), 'same');
    doubleKernel = doubleKernel./sum(doubleKernel);

    flipKernel = fliplr(kernel);
    
    %perform deconvolution, filling in NaNs with reconstructed values every
    %few iterations:
    J = deconvlucy({H},doubleKernel, 20);
    for iter = 1:5
            recon = convn(J{2}, doubleKernel, 'same');
            J{1}(nanFramesH) = recon(nanFramesH);
            J = deconvlucy(J,doubleKernel, 25);
    end
    H2 = J{2};
    H3 = convn(H2, kernel, 'same');


    J = deconvlucy({H},flipKernel, 20);
    for iter = 1:5
            recon = convn(J{2}, flipKernel, 'same');
            J{1}(nanFramesH) = recon(nanFramesH);
            J = deconvlucy(J,flipKernel, 25);
    end
    H4 = J{2};

    H5 = W \ IMselRaw;

    errH = (H - convn(H2, doubleKernel, 'same')).^2;
    errH = sqrt(convn(errH, doubleKernel, 'same')); % uncertainty at each point

    %compute F0
    F0mean = repmat(mean(F0sel,2, 'omitnan'),1,size(F0sel,2));
    F0sel(~isfinite(F0sel)) = max(0, F0mean(~isfinite(F0sel)));
    F0= (W./max(W,[],1))' *F0sel;
    F0(nanFramesH) = nan;
    
    %ensure F0 is positive
    Fmin = min(F0(:));
    desMin = prctile(F0(:),1) - Fmin;
    if Fmin<desMin
        F0 = F0 + desMin - Fmin;
    end

    %NaN out invalid data
    H(nanFramesH) = nan; %The match filtered signal
    H2(nanFramesH) = nan; %The detected events
    H3(nanFramesH) = nan; %The denoised activity
    H4(nanFramesH) = nan;
    H5(nanFramesH) = nan;

    exptSummary.dFerr{trialIx} = sum(W,1)'.*errH;
    exptSummary.matchFilt{trialIx}(:,:,1) = sum(W,1)'.*H; %[source#, time, channel]
    exptSummary.events{trialIx}(:,:,1) = sum(W,1)'.*H2; %[source#, time, channel]
    exptSummary.dF{trialIx}(:,:,1) = sum(W,1)'.*H3; %[source#, time, channel]
    exptSummary.dF2{trialIx}(:,:,1) = sum(W,1)'.*H4;
    exptSummary.dFls{trialIx}(:,:,1) = sum(W,1)'.*H5;
    exptSummary.F0{trialIx}(:,:,1) = F0;
    exptSummary.footprints(:,:,1:size(W0,2),trialIx) = Wfull;
    
    %compute channel 2 signals
    if numChannels==2
        IM2sel(isnan(IM2sel)) = 0;
        F_2 = (W./max(W,[],1))' * IM2sel;
        F_2(nanFramesH) = nan;

        F0_2 = nan(size(IM2sel));
        for rix = 1:size(F0selDS{trialIx}, 1)
            F0_2(rix,:) = interp(F02selDS{trialIx}(rix,:), params.dsFac)./params.dsFac;
        end

        F02mean = repmat(mean(F0_2,2, 'omitnan'),1,size(F0_2,2));
        F0_2(~isfinite(F0_2)) = max(0, F02mean(~isfinite(F0_2)));
        F0_2 = (W./max(W,[],1))' *F0_2;
        F0_2(nanFramesH) = nan;

        exptSummary.dF{trialIx}(:,:,2) = F_2 - F0_2;
        exptSummary.dF2{trialIx}(:,:,2) = F_2 - F0_2;
        exptSummary.dFls{trialIx}(:,:,2) = F_2 - F0_2;
        exptSummary.F0{trialIx}(:,:,2) = F0_2;
    end

    exptSummary.dFF{trialIx} = exptSummary.dF{trialIx}./exptSummary.F0{trialIx};
    exptSummary.dFF2{trialIx} = exptSummary.dF2{trialIx}./exptSummary.F0{trialIx};
    exptSummary.dFFls{trialIx} = exptSummary.dFls{trialIx}./exptSummary.F0{trialIx};
end

%prepare file for saving
exptSummary.fns = fns;
exptSummary.dr = dr;

%per-trial images
exptSummary.originalActivityChannel = activityChannel;
exptSummary.peaks = peaks;
exptSummary.perTrialMeanIMs = meanIM;
exptSummary.perTrialActIms = actIM;
exptSummary.perTrialAlignmentOffsets = motOutput; %the alignment vector for each trial

%save
save([drsave filesep fnsave], 'exptSummary','-v7.3');

disp('Done summarizeBergamo_Peaks')
end


function [IMsel, IMselRaw, F0,  W, selNans] = prepareNMFproblem(IMsel, W0, F0selDS, params)
if ~params.nmfBackgroundComps
    F0 = nan(size(IMsel));
    for rix = 1:size(F0selDS, 1)
        F0(rix,:) = interp(F0selDS(rix,:), params.dsFac)./params.dsFac;
    end
    %baselineWindow = ceil(params.baselineWindow_Glu_s/params.frametime);
    %F0 = computeF0(IMsel', params.denoiseWindow_samps*params.dsFac+1, baselineWindow, 1)';
    IMsel = IMsel - F0;
    selNans = isnan(IMsel);
    IMsel(selNans) = 0;
    IMselRaw = IMsel;
    IMsel = matchedExpFilter(IMsel, params.tau_full);
    %IMsel = IMsel- computeF0(IMsel', params.denoiseWindow_samps*params.dsFac+1,baselineWindow, 1)';
    W = W0;
else
    error('fitting global background components not implemented')
    % NOT IMPLEMENTED! 
    %GOALS
    % %compute a background (F0) per pixel
    % %NMF the background, and add those components to the NMF spatial components
    % %temporal filter the residual (dF)
    % %IMsel is the temporal filtered residual + background
    % F0sel = computeF0(IMsel', params.denoiseWindow*params.dsFac+1, params.baselineWindow*params.dsFac+1, 2)'; %algo2 doesn't underestimate F0 as badly in noisy conditions
    % dFsel = IMsel - F0sel;    
    % dFselTf = matchedExpFilter(dFsel, params.tau_full);
    % IMsel = dFselTf+F0sel;
    % 
    % selNans = isnan(IMsel);
    % nValid = sum(~selNans,2); %the number of values in each row
    % selRows = nValid>median(nValid)*0.9; %rows that we will use to estimate the background
    % 
    % IMBG = IMsel;
    % IMBG = IMBG - mean(IMBG,2, 'omitnan');
    % IMBG(selNans) = 0;
    % [~,~,VV] = svds(IMBG(selRows,:), params.nmfBackgroundComps);
    % pred = cat(2, VV, ones(size(VV,1),1));
    % IMBG2 = smoothdata(IMBG, 2, 'movmedian', params.baselineWindow, 'omitnan');
    % %for each row, fill nans with best fit
    % for rix = 1:size(IMBG,1)
    %     b = regress(IMBG2(rix, ~selNans(rix,:))', pred(~selNans(rix,:),:));
    %     fitVal = pred*b;
    %     IMBG(rix, selNans(rix,:)) = fitVal(selNans(rix,:));
    % end
    % 
    % %initialize components using PCA
    % tmp = IMsel - mean(IMsel,2, 'omitnan');
    % tmp(isnan(tmp)) = 0;
    % [U,~,VV] = svds(tmp,2);
    % clear tmp
    % UU = cat(2, U(:,1), -U(:,1), U(:,2), -U(:,2));
    % UU(UU<0) = 0;
    % UU = UU./sqrt(mean(UU.^2,1));
    % selU = mean(UU,1) > 0.66*mean(UU,'all');%select only components that are more spatially distributed
    % UU= UU(:,selU)+0.1;
    % 
    % 
    % IMBG =  smoothdata(IMBG,2,'movmean', ceil(params.baselineWindow*params.dsFac/2), 'omitnan');
    % opts = statset('MaxIter', 30,  'Display', 'final');
    % [Wbg,~] = nnmf(max(0,IMBG),size(UU,2), 'algorithm', 'mult', 'w0', UU, 'options', opts);
    % 
    % 
    % W = cat(2, W0, Wbg);
end
end

function IMsel = interpArray (IM, sel, shiftRC)
%linearly interpolate the 3D matrix IM in each 2D plane at the selected
%pixels sel, shifted by shiftRC
%returns a 2D array: [sum(sel) x size(IM,3)]
sz = size(IM);
IMsel = nan(sum(sel(:)), sz(3));
IM = reshape(IM, sz(1)*sz(2), []);

inds = zeros(size(sel));
inds(sel) = 1:sum(sel(:));

sel = sel(1:sz(1), 1:sz(2));
inds = inds(1:sz(1), 1:sz(2));

intShift = floor(shiftRC);
shiftRC = shiftRC-intShift;
sel = imtranslate(sel, [intShift(2) intShift(1)]);
inds = imtranslate(inds, [intShift(2) intShift(1)]);

%ensure that all the masks have the same number of values:
if shiftRC(1)>0.05
    sel(end,:) = false; inds(end,:) = false;
end
if shiftRC(2)>0.05
    sel(:,end) = false; inds(:,end) = false;
end
inds = inds(sel);

mask00 = sel;                    %unshifted
mask10 = imtranslate(sel,[0 1]); %shifted 1 row
mask01 = imtranslate(sel,[1 0]); %shifted 1 col
mask11 = imtranslate(sel,[1 1]); %shifted 1 row and 1 col

if shiftRC(1)>0.05 %the subpixel shift is nonnegligible, so interpolate
    R0 = (1-shiftRC(1)).*IM(mask00(:),:) + shiftRC(1).*IM(mask10(:),:);
    R1 = (1-shiftRC(1)).*IM(mask01(:),:) + shiftRC(1).*IM(mask11(:),:);
else %subpixel shift is negligible, use the unshifted data (this prevents NaNing out good data at edges)
    R0 = IM(mask00(:),:);
    R1 = IM(mask01(:),:);
end
if shiftRC(2)>0.05
    IMsel(inds,:) = (1-shiftRC(2)).*R0 + shiftRC(2).*R1;
else
    IMsel(inds,:) = R0;
end
end


function [density, peaks, sourceR, sourceC, params] = clusterLocalizations(peaks, params)
sz = params.sz; %XY size of the summary image, should be at least as large as the largest individual session; ideally all sessions are same shape 
upsample = params.upsample;
ampThresh = min(peaks.val);
nEventsTotal = length(peaks.val);
kernelSigma = 1/upsample;
rGrid = linspace(1,sz(1), upsample*sz(1)+1);
cGrid = linspace(1,sz(2), upsample*sz(2)+1);
[rr,cc] = ndgrid(rGrid,cGrid);
density = zeros(size(rr));

for eIx = 1:nEventsTotal
    rMin = peaks.row(eIx)-(8*kernelSigma);
    rMax = peaks.row(eIx)+(8*kernelSigma);
    cMin = peaks.col(eIx)-(8*kernelSigma);
    cMax = peaks.col(eIx)+(8*kernelSigma);
    rSel = rGrid>rMin & rGrid<rMax;
    cSel = cGrid>cMin & cGrid<cMax;

    rLoc = rr(rSel,cSel); cLoc = cc(rSel,cSel);
    A = mvnpdf([rLoc(:) cLoc(:)], [peaks.row(eIx) peaks.col(eIx)], kernelSigma.*eye(2));
    A = (peaks.val(eIx)./max(A)).*A;

    density(rSel,cSel) = density(rSel,cSel)+reshape(A, size(rLoc));
end

%figure, imagesc(cGrid,rGrid,density); hold on, scatter(peaks.col, peaks.row, 'r');

deconvSigma =sqrt(2)*upsample*params.sigma_px./sqrt(ampThresh); % sigmaXY/sqrt(amp) is the localization precision;
filtSize = 2*ceil(3*deconvSigma)+1;
selRows = imdilate(any(density,2), ones(filtSize,1));
selCols = imdilate(any(density,1), ones(1,filtSize));
PSF = fspecial('gaussian',filtSize,deconvSigma); %prior on loclaization accuracy
IMest = density;
IMest(selRows,selCols) = deconvlucy(density(selRows,selCols),PSF, 20); %should replace with our own algorithm

BW = imregionalmax(IMest);
BW = BW & IMest>(mean(IMest(IMest>0)));
[maxR, maxC] = find(BW);
[V, sortorder] = sort(IMest(BW), 'descend');
maxR = maxR(sortorder);
maxC = maxC(sortorder);

%compute pairwise distances to cull spurious maxima
if length(maxR)>1
    keep = true(1,length(V));
    dMaxima = squareform(pdist([maxR maxC]));
    dMaxima(eye(size(dMaxima), 'logical')) = inf;
    for vIx = 1:length(V)
        if ~isnan(dMaxima(vIx,vIx))
            sel = dMaxima(vIx,:)<(upsample);
            dMaxima(sel,:) = nan;
            dMaxima(:,sel) = nan;
            keep(sel) = false;
        end
    end
    maxR = maxR(keep);
    maxC = maxC(keep);
    V = V(keep);
end
k= length(V);
sourceR = rGrid(maxR);
sourceC = cGrid(maxC);

%Perform assignment using k-means-like approach
weights = ones(1,k);
assignments = zeros(1, nEventsTotal);
done = false;
ksigma = (upsample*params.sigma_px./sqrt(2*ampThresh));
while ~done
    zScore = sqrt((sourceR - peaks.row').^2 + (sourceC - peaks.col').^2)/ksigma; %squared distance from events to centers
    zScore(zScore>2) = inf;
    likelihoods = normpdf(zScore).*weights;
    [maxVal, maxInds] = max(likelihoods,[],2);
    maxInds(maxVal==0) = 0; %unassigned points
    
    % figure,
    % for ix = 1:length(sourceR)
    %     color = rand(1,3);
    %     sel = assignments==ix;
    %     scatter(sourceR(ix),sourceC(ix), 100,'k','x');
    %     hold on,
    %     scatter(peaks.row(sel), peaks.col(sel), 10,color, 'o')
    % end

    if all(maxInds(:)==assignments(:))
        done = true;
        keepEvents = min(zScore,[],2)<2;
        %remove extra sources
        for ix = k:-1:1
            keepSources(ix) = sum(assignments==ix)>=params.minEvents;
            keepEvents(assignments==ix) = keepSources(ix);
        end
        peaks.row = peaks.row(keepEvents);
        peaks.col = peaks.col(keepEvents);
        peaks.val = peaks.val(keepEvents);
        peaks.t = peaks.t(keepEvents);

        sourceR = sourceR(keepSources);
        sourceC =sourceC(keepSources);
        weights = weights(keepSources);
        k = sum(keepSources);

        %recalculate with sources removed
        zScore = sqrt((sourceR - peaks.row').^2 + (sourceC - peaks.col').^2)/ksigma; %squared distance from events to centers
        zScore(zScore>2) = inf;
        likelihoods = normpdf(zScore).*weights;
        [~, assignments] = max(likelihoods,[],2);

        assignProbs = likelihoods./sum(likelihoods,2);
    else
        assignments = maxInds;
        for ii = 1:k
            sel = maxInds==ii;
            sourceR(ii) = mean(peaks.row(sel));
            sourceC(ii) = mean(peaks.col(sel));
            weights(ii) = sqrt(sum(sel));
        end
    end
end
peaks.assignments = assignments;
peaks.assignProbs = assignProbs;

%Plot event assignments to sources
figure('name', 'Event assignments to sources'), imagesc(cGrid,rGrid,density); hold on;
colors = hsv(k);
colors = colors(randperm(k),:);
for sourceIx = 1:k
    sel = assignProbs(:,sourceIx)>0.5;
    scatter(sourceC(sourceIx), sourceR(sourceIx), 300, 'marker', 'x', 'markeredgecolor',colors(sourceIx,:), 'linewidth', 2); hold on;
    scatter(peaks.col(sel), peaks.row(sel),'markeredgecolor', colors(sourceIx,:));
end
end

function D = matchedExpFilter(D, tau_frames)
gamma = exp(-1/tau_frames);
selNans = isnan(D);
D(selNans) = 0;
for t = size(D,2)-1:-1:1
    D(:,t) = max(0,gamma*D(:,t+1)) + (1-gamma)*D(:,t);
end
D(selNans) = nan;
end

function [W0,H0] = extractSourcesLoRes(dFsel, peaks, sources, selPix, params)
k = length(sources.R);
sz = size(selPix, [1 2]);
anySel = any(selPix,3);
nSelPix = sum(anySel(:));

%Temporal filter the movie
dFselTf = matchedExpFilter(dFsel, params.tau_frames);
%baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));
%dFselTf = dFselTf - computeF0(dFselTf', params.denoiseWindow_samps, baselineWindow, 1)'; %subtracting F0 again to emphasize large events; we could uniformly subtract a quantile instead
selNans = isnan(dFselTf);
dFselTf(selNans) = 0;

nComp = k; %we could use extra components for background if desired
W0 = nan(nSelPix, nComp);
spatialScore = nan(1,nComp);
for sourceIx = 1:k
    thresh = min(prctile(peaks.assignProbs(peaks.assignments==sourceIx,sourceIx), 33), 0.5);
    selPeaks = find(peaks.assignProbs(:,sourceIx)>thresh);

    %create a selection vector for indexing into the smaller dFsel matrix
    selSel = selPix(:,:,sourceIx);
    [rr,cc] = find(selSel);
    spatialTemplate = mvnpdf(cat(2,rr,cc), [mean(rr) mean(cc)], params.sigma_px*eye(2)); % we will score each source based on how much the average event looks like a centered gaussian
    selSel = selSel(anySel);

    allEventsThisSite = nan(sum(selSel(:)), length(selPeaks));
    for eIx = 1:length(selPeaks)
        T = peaks.t(selPeaks(eIx));
        allEventsThisSite(:,eIx) = dFselTf(selSel, T);
    end
    avgEvent = mean(allEventsThisSite,2, 'omitnan');
    spatialScore(sourceIx) = corr(avgEvent,spatialTemplate).*sum(avgEvent);
    thresh = max(avgEvent(:))*1e-5;
    avgEvent(avgEvent<thresh) = thresh;
    W0(selSel,sourceIx) = avgEvent;
end

%discard sources with poor spatial profiles/response amplitudes
spaceThresh = median(spatialScore)/3;
W0 = W0(:, spatialScore>spaceThresh);
W0(isnan(W0)) = 0;
nComp = size(W0,2);

W0full = zeros(sz(1)*sz(2),nComp);
W0full(anySel(:),:) = W0;
W0full = reshape(W0full, sz(1),sz(2),[]);
W0full = min(W0full, imgaussfilt(W0full, params.sigma_px/2));
W0full = reshape(W0full, sz(1)*sz(2),[]);
W0 =  W0full(anySel,:);

%Use multiplicative updates NMF, which makes it easy to zero out pixels
opts1 = statset('MaxIter', 6,  'Display', 'final');%, 'UseParallel', true);
[W0,H0] = nnmf(dFselTf, nComp,'algorithm', 'mult', 'w0',W0, 'options', opts1); %!!nnmf has been modified to allow it to take more than rank(Y) inputs
for bigIter = 1:(params.nmfIter+3)
    disp(['outer loop ' int2str(bigIter) ' of ' int2str(params.nmfIter)]);

    %apply sparsity
    W0 = max(0, W0-params.sparseFac.*max(W0,[],1));
    
    W0full = zeros(sz(1)*sz(2),nComp);
    W0full(anySel(:),:) = W0;
    W0full = reshape(W0full, sz(1),sz(2),[]);
    W0full = min(W0full, imgaussfilt(W0full, params.sigma_px/2)); %sculpt the spatial profiles
    W0full = reshape(W0full, sz(1)*sz(2),[]);
    [W0,H0] = nnmf(dFselTf, nComp,'algorithm', 'mult', 'w0', W0full(anySel,:), 'h0', H0, 'options', opts1);

    if bigIter == params.nmfIter
        disp('Merging sources...')
        [W0,H0] = mergeSources (W0,H0,dFselTf,anySel);
        disp(['Kept ' int2str(size(W0,2)) ' of ' int2str(nComp) 'sources']);
        nComp = size(W0,2);
    end
end

%evaluate how well fit the data are
%figure, imshow3D(cat(3,dFselTf, W0*H0, dFselTf-W0*H0));
end

% function [W,H] = mergeSources (W,H)
% %sort by variance; nnmf usually does htis automatically but we disabled it
% %in nnmf2
% [~, sortorder] = sort(sum(W.^2,1), 'descend');
% W = W(:,sortorder);
% H = H(sortorder,:);
% 
% % compute correlation in activity
% C = corr(H'); C(logical(eye(size(C)))) = nan; 
% overlap = ((W>0)'*(W>0))>2;
% k = size(W,2);
% keep = true(1,k);
% merged = false(1,k); 
% for sourceIx = k:-1:2 % assumes inputs are already sorted by variance explained
%     if merged(sourceIx)
%         continue %this source already has had others merged into it; don't merge recursively
%     end
%     %if a source is better predicted by an overlapping higher-variance source than by any
%     %non-overlapping source, and vice versa, merge them
%     [maxC, maxInd] = max(C(sourceIx, 1:sourceIx-1).^2);
%     if overlap(sourceIx, maxInd) && maxC==max(C(maxInd, ~overlap(maxInd,:)).^2)
%        %merge 
%        keep(sourceIx) = false;
%        merged(maxInd) = true;
%        W(:,maxInd) = W(:,maxInd) + W(:,sourceIx);
%        H(maxInd,:) = H(maxInd,:) +  H(sourceIx,:);
%     end
% end
% W = W(:,keep);
% H(merged,:) = H(merged,:)./sum(H(merged,:).^2, 2); %mormalize the activities that we merged
% H = H(keep,:); 
% end
