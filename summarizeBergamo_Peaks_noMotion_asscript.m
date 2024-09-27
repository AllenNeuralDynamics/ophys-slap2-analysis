addpath(genpath('matlab'));
%dr = '/local/data/scan_00003_20240816_100000'
%fns = 'scan_00003_20240816_100000_REGISTERED_DOWNSAMPLED-2x.tif'

%dr = '/local/data/scan_00003_20240816_001000'
%fns = 'scan_00003_20240816_001000_REGISTERED_DOWNSAMPLED-2x.tif'

%dr = '/local/data/scan_00001_20240720_100000/'
%fns = 'scan_00001_20240720_100000_REGISTERED_DOWNSAMPLED-2x.tif'
% summarizeBergamo_Peaks(dr, fn)

dr = '/local/data/iGluSnFR-simulation/113'
fns = 'SIMULATION_scan_00001_113_Trial1.tif'

Tmax = 50000;


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

% begin ADDED
params.frametime = 0.0023;
params.dsFac = 1;
activityChannel = 1;
numChannels = 1
% end ADDED

if ~iscell(fns)
    fns = {fns};
end

% fnStemEnd = strfind(fns{1}, '_REGISTERED') -1;
% fnStem = fns{1}(1:fnStemEnd);
% savedr = dr;
% % savedr = [dr filesep 'ExperimentSummary'];
% % if ~exist(savedr, 'dir')
% %     mkdir(savedr);
% % end
% fnsave = [fnStem '_EXPTSUMMARY.mat'];
% drsave = savedr;
% %load some metadata
% load([dr filesep fnStem '_ALIGNMENTDATA.mat'], 'aData');
% numChannels = aData.numChannels;


%generate a concensus alignment across trials for further analysis
meanIM = nan(1,1,numChannels,1);
actIM = nan(1,1,numChannels,1);
disp('Loading data and performing localizations...')
for trialIx = length(fns):-1:1
    fn = fns{trialIx};

    % %ensure the high res file exists
    % ind =strfind(fn, '_DOWNSAMPLED');
    % fnRaw{trialIx} = [fn(1:ind) 'RAW.tif'];
    % assert(exist([dr filesep fnRaw{trialIx}], 'file'), ['No corresponding RAW tiff recording found for file:' fn])

    %load the tiff
    [IM, desc, meta] = networkScanImageTiffReader([dr filesep fn]);
    IM = IM(:,:,1:min(Tmax, size(IM, 3)));
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

    % %load alignment data
    % fnStemEnd = strfind(fn, '_REGISTERED') -1;
    % load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
    % aData.dsFac = round(length(aData.motionC)./length(aData.motionDSc));
    % params.dsFac = aData.dsFac;
    % params.frametime = aData.frametime;

    % %discard motion frames
    % % tmp = aData.aRankCorr(:)-smoothExp(aData.aRankCorr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorr,2, 'movmedian', ceil(2/aData.frametime));
    % % discardFrames{trialIx} = imdilate(tmp<-(4*std(tmp)), ones(1,5));
    % tmp = aData.aRankCorr(:)-smoothExp(aData.aRankCorr(:),'movmedian', ceil(10/(aData.frametime*aData.dsFac)));
    % filtTmp = smoothExp(tmp, 'movmean',ceil(.2/(aData.frametime*aData.dsFac)));
    % discardFrames{trialIx} = imdilate(filtTmp<-(4*std(filtTmp)), ones(1,5));
    rawIMs{trialIx} = squeeze(IM(:,:,1,:));
    % rawIMs{trialIx}(:,:,discardFrames{trialIx}) = nan;
    % if numChannels==2
    %     rawIM2s{trialIx} = squeeze(IM(:,:,2,:));
    %     rawIM2s{trialIx}(:,:,discardFrames{trialIx}) = nan;
    % end

    % begin ADDED
    aData = {};
    rawIMs{trialIx}(1,:,:) = nan;
    rawIMs{trialIx}(end,:,:) = nan;
    rawIMs{trialIx}(:,1,:) = nan;
    rawIMs{trialIx}(:,end,:) = nan;
    % end ADDED

    disp('Localizing Flashes')
    tic
    [IMc, peaks(trialIx), params] = localizeFlashesBergamo(rawIMs{trialIx}, aData, params);
    toc

    %calculate correlation image
    actIM(end:size(IMc,1),:,:,:) = nan;
    actIM(:,end:size(IMc,2),:,:) = nan;
    actIM(:,:,:,trialIx) = nan;
    actIM(1:size(IMc,1),1:size(IMc,2),1,trialIx) = IMc;
end
params.sz = size(meanIM, [1 2]);

% clear aData

%Make template
disp('Making template for aligning across trials...')
maxshift = 12;
M = squeeze(sum(meanIM, 3));
template = makeTemplateMultiRoi(M, maxshift);

% %align all mean images to template
% disp('Aligning across trials...')
% Mpad = nan([size(template) size(M,3)]);
% Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
% for trialIx = length(fns):-1:1
%     disp(['trial: ' int2str(trialIx)])
%     mot1 = xcorr2_nans(Mpad(:,:,trialIx), template, [0 ; 0], maxshift);
%     [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
%     [rr,cc] = ndgrid(1:size(meanIM,1), 1:size(meanIM,2));
%     for chIx = 1:size(meanIM,3)
%         meanAligned(:,:,chIx,trialIx) = interp2(meanIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
%         actAligned(:,:,chIx,trialIx) = interp2(actIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
%     end
%     if trialIx==length(fns)
%         peaksCat.row = peaks(trialIx).row - motOutput(1,trialIx);
%         peaksCat.col = peaks(trialIx).col - motOutput(2,trialIx);
%         peaksCat.val = peaks(trialIx).val;
%         peaksCat.t = peaks(trialIx).t;
%          figure, imagesc(template)
%          hold on, scatter( maxshift+ peaks(trialIx).col - motOutput(2,trialIx), maxshift+ peaks(trialIx).row - motOutput(1,trialIx));
%     else
%         peaksCat.row = cat(2, peaks(trialIx).row - motOutput(1,trialIx),peaksCat.row);
%         peaksCat.col = cat(2, peaks(trialIx).col - motOutput(2,trialIx), peaksCat.col);
%         peaksCat.val = cat(1, peaks(trialIx).val, peaksCat.val);
%         peaksCat.t = peaksCat.t + size(rawIMs{trialIx},3);
%         peaksCat.t = cat(1, peaks(trialIx).t,peaksCat.t);
%          hold on, scatter( maxshift+ peaks(trialIx).col - motOutput(2,trialIx), maxshift+ peaks(trialIx).row - motOutput(1,trialIx));
%     end
% end

% %identify outliers in alignment quality
% cc = corrCoeff;
% corrThresh = min(0.96, median(cc)-2*std(cc));
% 
% validTrials= find(cc>corrThresh);
% exptSummary.meanIM = mean(meanAligned,4, 'omitnan');
% exptSummary.actIM = mean(actAligned, 4, 'omitnan');


%cluster localizations
totalFrames = sum(cellfun(@(x)(size(x,3)), rawIMs));
params.minEvents = totalFrames*params.frametime*params.dsFac*params.eventRateThresh_hz;
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

P = peaks;
sources.R = sourceR;
sources.C = sourceC;

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

% function IMsel = interpArray (IM, sel, shiftRC)
% %linearly interpolate the 3D matrix IM in each 2D plane at the selected
% %pixels sel, shifted by shiftRC
% %returns a 2D array: [sum(sel) x size(IM,3)]
% sz = size(IM);
% IMsel = nan(sum(sel(:)), sz(3));
% IM = reshape(IM, sz(1)*sz(2), []);
% 
% inds = zeros(size(sel));
% inds(sel) = 1:sum(sel(:));
% 
% sel = sel(1:sz(1), 1:sz(2));
% inds = inds(1:sz(1), 1:sz(2));
% 
% intShift = floor(shiftRC);
% shiftRC = shiftRC-intShift;
% sel = imtranslate(sel, [intShift(2) intShift(1)]);
% inds = imtranslate(inds, [intShift(2) intShift(1)]);
% size(sel)
% %ensure that all the masks have the same number of values:
% if shiftRC(1)>0.05
%     sel(end,:) = false; inds(end,:) = false;
% end
% if shiftRC(2)>0.05
%     sel(:,end) = false; inds(:,end) = false;
% end
% inds = inds(sel);
% 
% mask00 = sel;                    %unshifted
% mask10 = imtranslate(sel,[0 1]); %shifted 1 row
% mask01 = imtranslate(sel,[1 0]); %shifted 1 col
% mask11 = imtranslate(sel,[1 1]); %shifted 1 row and 1 col
% 
% if shiftRC(1)>0.05 %the subpixel shift is nonnegligible, so interpolate
%     R0 = (1-shiftRC(1)).*IM(mask00(:),:) + shiftRC(1).*IM(mask10(:),:);
%     R1 = (1-shiftRC(1)).*IM(mask01(:),:) + shiftRC(1).*IM(mask11(:),:);
% else %subpixel shift is negligible, use the unshifted data (this prevents NaNing out good data at edges)
%     R0 = IM(mask00(:),:);
%     R1 = IM(mask01(:),:);
% end
% if shiftRC(2)>0.05
%     IMsel(inds,:) = (1-shiftRC(2)).*R0 + shiftRC(2).*R1;
% else
%     IMsel(inds,:) = R0;
% end
% end

% begin ADDED 
validTrials = 1:length(fns);
% end ADDED 
for trialIx = validTrials
    szTmp = size(rawIMs{trialIx});
    % begin REPLACED
    % IMrawSel = interpArray(rawIMs{trialIx}, any(selPix,3), motOutput(:,trialIx)); %interpolates the movie at the shifted coordinates
    IMrawSel = reshape(rawIMs{trialIx}(repmat(any(selPix,3), 1, 1, szTmp(3))), [], szTmp(3));
    % end REPLACED 
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

