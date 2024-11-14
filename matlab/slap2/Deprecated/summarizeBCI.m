function summarizeBCI(dr)
    %TO DO:
    %PARAMETER SENSITIVITY ANALYSIS

    %add some global NMF components?
        %initialize with random values
%cd('\\allen\aind\scratch\ophys\BCI\inactive_mice\709390_SBCI12\SLAP2\slap2_709390_2024-01-10_12-17-31')
%cd('C:\temp\SYNAPSES\Test');
%cd('C:\Users\kaspar.podgorski\OneDrive - Allen Institute\Documents\GitHub\ophys-slap2-analysis\matlab\rasterDataProcessing\Bergamo\simulations\data\SIMULATIONS')
if isempty(gcp('nocreate'))
    nWorkers = 10;
    if nWorkers<15
        warning('You are using few parallel workers! adjust this in summarizeBCI.m');
    end
    disp(['Parallel workers:' int2str(nWorkers)])
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end
disp(['## SUMMARIZEBCI ##' newline 'Folder:'])
disp(dr)
nDMDs = 2;

%general params
params.analyzeHz = 200; %frame rate used for analysis
params.discardInitial_s = 0.1; %discard the first short period of each trial as the beam stabilization locks on and the imaging system warms up

%filtering params
params.sigma_px = 1.5; %1.33;   % space constant in pixels
done = false;
while ~done
    params.tau_s = str2double(inputdlg('Time constant of glutamate indicator (s)?', 'tau', 1, {'0.05'}));
    done = params.tau_s>0 && params.tau_s<1;
end
%params.tau_s = 0.027; % time constant in seconds for glutamate channel; from Aggarwal et al 2023 Fig 5
%params.denoiseWindow_samps = 20; %number of samples to average together for denoising; REFACTOR THIS OUT AND USE SECONDS!
params.denoiseWindow_s = 0.25; %number of samples to average together for denoising
params.baselineWindow_Glu_s = 2; %timescale for calculating F0 in glutamate channel, seconds
params.baselineWindow_Ca_s = 2; %timescale for calculating F0 in calcium channel, seconds

%localization params
params.eventRateThresh_hz = 0.01; % minimum event rate in Hz
params.tilesizeLoc = 64; %Tile size used for computing statistics for localization; Accounts for image statistics varying slowly across FOV.
params.threshSNRloc = 12; %6; % Z-score threshold for calling a localization
params.threshSKloc = 4; %threshold on skewness-based summary statistic for pixels to consider for localizations
params.upsample = 3; %how many times to upsample the imaging resolution for finding local maxima to identify sources; affects maximum source density

%NMF params
params.sparseFac = 0.08; %sparsity factor for shrinking sources in space, 0-1, higher value makes things sparser
params.nmfIter = 5; %number of iterations of NMF refinement
params.dXY = 4; %how large sources can be (radius), pixels
params.nmfBackgroundComps = 0; % <=4, max number of background components to use for NMF. If 0, we compute F0 instead of fitting background

%load trial table
load([dr filesep 'trialTable.mat'], 'trialTable');

savedr = [dr filesep 'ExperimentSummary'];
if ~exist(savedr, 'dir')
    mkdir(savedr);
end
fnsave = [savedr filesep 'Summary-' datestr(now, 'YYmmDD-HHMMSS') '.mat'];

%confirm that all files exist for both DMDs
nTrials = length(trialTable.trueTrialIx);
keepTrials = true(nDMDs, nTrials);
for trialIx = nTrials:-1:1
    for DMDix = 1:nDMDs
        trialStr = ['E' int2str(trialTable.epoch(trialIx)) 'T' int2str(trialIx) 'DMD' int2str(DMDix)];
        tiffFn = [trialStr '_REGISTERED_DOWNSAMPLED-*Hz.tif'];
        if isempty(dir([dr filesep tiffFn]))
            disp(['Missing tiff file:' tiffFn])
            keepTrials(DMDix,trialIx) = false;
        end
        alignFn =  [trialStr '_ALIGNMENTDATA.mat'];
        if ~exist([dr filesep alignFn], 'file')
            disp(['Missing alignData file:' alignFn])
            keepTrials(DMDix,trialIx) = false;
        end
        fnRaw{trialIx} = trialTable.(strcat('DMD', int2str(DMDix), 'filename')){trialIx};
        if ~exist([dr filesep fnRaw{trialIx}], 'file')
            disp(['Missing Dat file:' alignFn])
            keepTrials(DMDix, trialIx) = false;
        end
    end
end
if ~all(keepTrials)
    disp(['Files were missing for ' int2str(sum(~keepTrials(:))) ' [DMDs x trials]; likely failed alignments. Proceeding without them.']);
end
if ~any(keepTrials)
    error('All trials were rejected due to missing alignment files!');
end

%call up a GUI for the user to define Soma ROI and regions to exclude
fnAnn = [dr filesep 'ANNOTATIONS.mat'];
if exist(fnAnn, 'file')
    load(fnAnn, 'ROIs')
else
    for DMDix = 1:nDMDs
        %load image data
        trialStr = ['E' int2str(trialTable.epoch(1)) 'T1' 'DMD' int2str(DMDix)];
        fobj = dir([dr filesep trialStr '_REGISTERED_DOWNSAMPLED-*Hz.tif']);
        fn = fobj(1).name;
        IM = copyReadDeleteScanImageTiff([dr filesep fn]);
        IM = squeeze(mean(IM,[3 4], 'omitnan'));
        hROIs(DMDix) = drawROIs(sqrt(max(0,IM)), dr, fn);
        ROIs(DMDix).dr = dr;
        ROIs(DMDix).fn = fn;
    end
    for DMDix = 1:nDMDs
        waitfor(hROIs(DMDix).hF);
        ROIs(DMDix).roiData = hROIs(DMDix).roiData;
    end
    save(fnAnn, 'ROIs'); clear hROIs;
end

%load some metadata
firstValidTrial = find(keepTrials(DMDix,:),1,'first');
trialStr = ['E' int2str(trialTable.epoch(firstValidTrial)) 'T' int2str(firstValidTrial) 'DMD' int2str(DMDix)];
load([dr filesep trialStr '_ALIGNMENTDATA.mat'], 'aData');
numChannels = aData.numChannels;
params.numChannels = numChannels;
params.dsFac = 1; %aData.dsFac;
params.frametime = aData.frametime;

%generate a concensus alignment across trials for further analysis
for DMDix = nDMDs:-1:1
%meanIM = nan(1,1,numChannels,1);
%actIM = nan(1,1,numChannels,1);
disp('Loading data and performing localizations...')

%Perform Localizations
mIM = cell(1, nTrials); aIM = cell(1,nTrials); alignData = cell(1, nTrials); peaks = []; discardFrames = cell(1,nTrials); %rawIMs = cell(1,nTrials);
parfor trialIx = 1:nTrials
    if keepTrials(DMDix,trialIx)
        trialStr = ['E' int2str(trialTable.epoch(trialIx)) 'T' int2str(trialIx) 'DMD' int2str(DMDix)];
        fobj = dir([dr filesep trialStr '_REGISTERED_DOWNSAMPLED-*Hz.tif']);
        fn = fobj(1).name;
        [~, mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks(trialIx), discardFrames{trialIx}]= loadAndLocalizeTrialAsync(dr, fn, numChannels, params); %rawIMs{trialIx}
    end
end
%Assemble same-sized mean images from different-sized trial means
szm1 = max(cellfun(@(x)size(x,1),mIM)); szm2 = max(cellfun(@(x)size(x,2), aIM));
meanIM = nan(szm1,szm2,numChannels, nTrials); actIM = nan(szm1,szm2,1, nTrials);
for trialIx = 1:nTrials
    tmp =  mIM{trialIx};
    meanIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp; 
    tmp =  aIM{trialIx};
    actIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp; 
end
params.sz = size(meanIM, [1 2]);

clear aData

%Make template
disp('Making template for aligning across trials...')
maxshift = 5;
M = squeeze(sum(meanIM, 3));
template = makeTemplateMultiRoi(M(:,:,keepTrials(DMDix,:)), maxshift);

%align all mean images to template
disp('Aligning across trials...')
meanAligned = []; actAligned = []; 
corrCoeff = nan(1,nTrials);
motOutput = nan(2,nTrials);
peaksCat = struct;
Mpad = nan([size(template) size(M,3)]);
Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
for trialIx = nTrials:-1:1
    if ~keepTrials(DMDix,trialIx)
        continue %skip
    end
    disp(['trial: ' int2str(trialIx)])
    mot1 = xcorr2_nans(Mpad(:,:,trialIx), template, [0 ; 0], maxshift);
    [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
    [rr,cc] = ndgrid(1:size(meanIM,1), 1:size(meanIM,2));
    
    for chIx = 1:size(meanIM,3)
        meanAligned(:,:,chIx,trialIx) = interp2(meanIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
    end
    actAligned(:,:,1,trialIx) = interp2(actIM(:,:,1,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
end

%identify outliers in alignment quality to determine valid trials
ccf = corrCoeff;
if nargin>1 && forceCorrThresh>0
    corrThresh = forceCorrThresh;
else
    corrThresh = min(0.90, median(ccf, 'omitnan')-2*std(ccf, 'omitmissing'));
end
validTrials= find(ccf>corrThresh);
exptSummary.meanIM{DMDix} = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');
exptSummary.actIM{DMDix} = mean(actAligned(:,:,:,validTrials), 4, 'omitnan');

%accumulate peaks, only from valid trials
nValidFramesDS = 0;
for vTrialIx = 1:length(validTrials)
    trialIx = validTrials(vTrialIx);
    if ~isfield(peaksCat, 'row')
        peaksCat.row = peaks(trialIx).row - motOutput(1,trialIx);
        peaksCat.col = peaks(trialIx).col - motOutput(2,trialIx);
        peaksCat.val = peaks(trialIx).val;
        peaksCat.t = peaks(trialIx).t;
         % figure, imagesc(template)
         % hold on, scatter( maxshift+ peaks(trialIx).col - motOutput(2,trialIx), maxshift+ peaks(trialIx).row - motOutput(1,trialIx));
    else
        peaksCat.row = cat(2, peaksCat.row, peaks(trialIx).row - motOutput(1,trialIx));
        peaksCat.col = cat(2, peaksCat.col, peaks(trialIx).col - motOutput(2,trialIx));
        peaksCat.val = cat(1, peaksCat.val, peaks(trialIx).val);
        peaksCat.t = peaksCat.t + size(discardFrames{trialIx},1);
        peaksCat.t = cat(1,peaksCat.t, peaks(trialIx).t+nValidFramesDS);
         % hold on, scatter( maxshift+ peaks(trialIx).col - motOutput(2,trialIx), maxshift+ peaks(trialIx).row - motOutput(1,trialIx));
    end
    nValidFramesDS=  nValidFramesDS + size(discardFrames{trialIx},1);
end

%cluster localizations
totalFrames = sum(cellfun(@(x)(sum(~x)), discardFrames));
params.minEvents = totalFrames*params.frametime*params.dsFac*params.eventRateThresh_hz;
[~, P, sources.R, sources.C, params] = clusterLocalizations(peaksCat, params);
k = length(sources.R); %number of sources
sz = params.sz;

%Generate IMsel; the data only in the selected region, aligned across movies
selPix = false([sz(1:2) k]);
for sourceIx = k:-1:1
    rr = max(1, round(sources.R(sourceIx)-params.dXY)):min(sz(1),  round(sources.R(sourceIx))+params.dXY);
    cc = max(1, round(sources.C(sourceIx)-params.dXY)):min(sz(2),  round(sources.C(sourceIx))+params.dXY);
    selPix(rr,cc,sourceIx) = true;
end
pxAlwaysValid = mean(isnan(meanAligned(:,:,1,validTrials)),4)<0.1;
selPix = selPix & repmat(pxAlwaysValid, 1, 1, k); %ADJUST SELECTED PIXELS NOT TO INCLUDE POORLY MEASURED PIXELS

%prune any sources that got clipped by pixel selection process
keepSources = sum(selPix, [1 2])>5;
k = sum(keepSources, 'all');
sources.R = sources.R(keepSources);
sources.C = sources.C(keepSources);
selPix = selPix(:,:,keepSources);

%accumulate dF of selected pixels for initial NMF of downsampled movies
F0selDS = cell(nTrials,1); dFsel = cell(1,nTrials);
parfor trialIx = 1:nTrials
    if any(validTrials==trialIx)
        trialStr = ['E' int2str(trialTable.epoch(trialIx)) 'T' int2str(trialIx) 'DMD' int2str(DMDix)];
        fobj = dir([dr filesep trialStr '_REGISTERED_DOWNSAMPLED-*Hz.tif']);
        fn = fobj(1).name;
        [dFsel{trialIx}, F0selDS{trialIx}] = DFselAsync([dr filesep fn], discardFrames{trialIx}, selPix, motOutput(:,trialIx), params); %path, discardFrames, selPix, motOutput, params
    end
end
dFsel = cell2mat(dFsel); %collapse into an array

%error checking
assert(nValidFramesDS==size(dFsel,2));
assert(max(P.t)<nValidFramesDS);

%extract sources from downsampled movie dF
try
    [W0,~] = extractSourcesLoRes(dFsel, P, sources, selPix, params);
catch ME
    disp('fatal Error extracting sources')
    rethrow(ME)
end

%for each file, load high res data and refine
params.tau_full=params.tau_s*params.analyzeHz;
roiData = ROIs(DMDix).roiData;
E = cell(nTrials,1);
parfor trialIx = 1:nTrials
    if any(validTrials==trialIx)
        fnRaw = trialTable.(strcat('DMD', int2str(DMDix), 'filename')){trialIx};
        startLine=trialTable.(strcat('DMD', int2str(DMDix), 'firstLine'))(trialIx);
        endLine =trialTable.(strcat('DMD', int2str(DMDix), 'lastLine'))(trialIx);
        E{trialIx} = processTrialAsync(dr, fnRaw, startLine, endLine, W0, F0selDS{trialIx}, selPix, discardFrames{trialIx}, alignData{trialIx}, mIM{trialIx}, motOutput(:,trialIx), roiData, params);
    end
end

%per-trial images
exptSummary.E(:,DMDix) = E; %experiment data
exptSummary.aData(:,DMDix) = alignData;
exptSummary.userROIs{DMDix} = roiData;
exptSummary.peaks{DMDix}= peaks;
exptSummary.perTrialMeanIMs{DMDix} = meanIM;
exptSummary.perTrialMeanIMsAligned{DMDix} = meanAligned;
exptSummary.perTrialActIms{DMDix} = actIM;
exptSummary.perTrialActIMsAligned{DMDix} = actAligned;
exptSummary.perTrialAlignmentOffsets{DMDix} = motOutput; %the alignment vector for each trial
end

%prepare file for saving
exptSummary.params = params;
exptSummary.trialTable = trialTable;
exptSummary.dr = dr;


% -- Lucas' Contribution ----
% behavioralData = extract_behavior_from_harp(dr);
% exptSummary.behavior = behavioralData;
% ---------------------------


%save
save(fnsave, 'exptSummary', "-v7.3");
disp('Done summarizeMultiROI_Peaks')
end

function IMsel = interpArray (IM, sel, shiftRC)
%linearly interpolate the 3D matrix IM in each 2D plane at the selected
%pixels sel, shifted by shiftRC
%returns a 2D array: [sum(sel) x size(IM,3)]
sz = size(IM);
assert(all(size(sel)==sz(1:2)));
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

function [W0,H0] = extractSourcesLoRes(dFsel, peaks, sources, selPix, params)
%dFsel: delta fluorescence over selected pixels;  has dimensions [pixels time]

k = length(sources.R);
sz = size(selPix, [1 2]);
anySel = any(selPix,3);
nSelPix = sum(anySel(:));

%Temporal filter the movie
tau = params.tau_s./(params.frametime*params.dsFac); %time constant in frames
dFselTf = matchedExpFilter(dFsel, tau);
%baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));
%dFselTf = dFselTf - computeF0(dFselTf', params.denoiseWindow_samps, baselineWindow, 1)'; %subtracting F0 again to emphasize large events; we could uniformly subtract a quantile instead

selNans = imclose(isnan(dFselTf), ones(1, 2*ceil(params.denoiseWindow_s/params.frametime)+1));
meanDF = repmat(mean(dFselTf,2, 'omitnan'), 1, size(dFsel,2));
meanDF(isnan(meanDF)) = 0;
dFselTf(selNans) = meanDF(selNans);

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
    avgEvent = mean(allEventsThisSite,2, 'omitmissing');
    spatialScore(sourceIx) = corr(avgEvent,spatialTemplate).*sum(avgEvent);
    thresh = max(avgEvent(:))*1e-5;
    avgEvent(avgEvent<thresh) = thresh;
    W0(selSel,sourceIx) = avgEvent;
end

%discard sources with poor spatial profiles/response amplitudes
spaceThresh = median(spatialScore, 'omitmissing')/3;
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
opts1 = statset('MaxIter', 6,  'Display', 'final');
[W0,H0] = nnmf2(dFselTf, nComp,'algorithm', 'mult', 'w0',W0, 'options', opts1); %!!nnmf has been modified to allow it to take more than rank(Y) inputs
for bigIter = 1:(params.nmfIter+3)
    H0 = H0+sqrt(0.1/size(H0,2)); %get rid of zeroing out effect
    disp(['outer loop ' int2str(bigIter) ' of ' int2str(params.nmfIter)]);

    %apply sparsity
    W0 = max(0, W0-params.sparseFac.*max(W0,[],1));
    
    W0full = zeros(sz(1)*sz(2),nComp);
    W0full(anySel(:),:) = W0;
    W0full = reshape(W0full, sz(1),sz(2),[]);
    W0full = min(W0full, imgaussfilt(W0full, params.sigma_px/2)); %sculpt the spatial profiles
    W0full = reshape(W0full, sz(1)*sz(2),[]);
    [W0,H0] = nnmf2(dFselTf, nComp,'algorithm', 'mult', 'w0', W0full(anySel,:), 'h0', H0, 'options', opts1);

    if bigIter == params.nmfIter
        disp('Merging sources...')
        [W0,H0] = mergeSources(W0,H0, dFselTf, anySel);
        disp(['Kept ' int2str(size(W0,2)) ' of ' int2str(nComp) 'sources']);
        nComp = size(W0,2);
    end
end

%evaluate how well fit the data are
%figure, imshow3D(cat(3,dFselTf, W0*H0, dFselTf-W0*H0));
end

function [dFsel, F0selDS] = DFselAsync(path, discardFrames, selPix, motOutput, params)
    baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));
    denoiseWindow = ceil(params.denoiseWindow_s/(params.frametime*params.dsFac)); %params.denoiseWindow_samps;
    
    rawIM = copyReadDeleteScanImageTiff(path);
    rawIM = reshape(rawIM, size(rawIM,1), size(rawIM,2), params.numChannels, []); %deinterleave;
    rawIM = squeeze(rawIM(:,:,1,:));
    rawIM(:,:,discardFrames) = nan;
    
    szTmp = size(rawIM);
    sel2D = any(selPix,3); 
    sel2D = sel2D(1:szTmp(1), 1:szTmp(2));
    
    IMrawSel = interpArray(rawIM, sel2D, motOutput); %interpolates the movie at the shifted coordinates
    F0selDS = svdF0_slap2(IMrawSel', 8, denoiseWindow, baselineWindow)';
    dFsel = IMrawSel - F0selDS;
end