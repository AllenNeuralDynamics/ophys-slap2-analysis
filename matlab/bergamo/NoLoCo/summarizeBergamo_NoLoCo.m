function summarizeBergamo_NoLoCo(fns, dr, paramsIn)

%PARAMETER SETTING
if nargin>2
    params = setParams('summarizeBergamo_NoLoCo', paramsIn);
else
    params = setParams('summarizeBergamo_NoLoCo');
end

%SELECT FILES IF NOT PASSED IN
if ~nargin || isempty(fns)
    [fns, dr] = uigetfile('*REGISTERED*.tif', 'multiselect', 'on');
end
if ~iscell(fns)
    fns = {fns};
end
drsave = dr;
fnStemEnd = strfind(fns{1}, '_REGISTERED') -1;
fnStem = fns{1}(1:fnStemEnd);
fnsave = [fnStem '_EXPTSUMMARY_NOLOCO.mat'];

%LOAD METADATA
load([dr filesep fnStem '_ALIGNMENTDATA.mat'], 'aData');
numChannels = aData.numChannels;
params.numChannels = numChannels;
nTrials = length(fns);

%SUMMARIZE ACTIVITY FOR EACH TRIAL
meanIM = nan(1,1,numChannels,1);
actIM = nan(1,1,numChannels,1);
disp('Loading data and performing localizations...')
for trialIx = nTrials:-1:1
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
    if params.activityChannel>1
        IM = IM(:,:,[params.activityChannel:end, 1:params.activityChannel-1],:);
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
    params.analyzeHz = 1/params.frametime;
    params.alignHz = params.analyzeHz./params.dsFac;

    %discard motion frames
    tmp = aData.aRankCorr(:)-smoothExp(aData.aRankCorr(:),'movmedian', ceil(10/(aData.frametime*aData.dsFac)));
    filtTmp = smoothExp(tmp, 'movmean',ceil(.2/(aData.frametime*aData.dsFac)));
    discardFrames{trialIx} = imdilate(filtTmp<-(4*std(filtTmp)), ones(1,5));
    rawIMs{trialIx} = squeeze(IM(:,:,1,:));
    rawIMs{trialIx}(:,:,discardFrames{trialIx}) = nan;
    if numChannels==2
        rawIM2s{trialIx} = squeeze(IM(:,:,2,:));
        rawIM2s{trialIx}(:,:,discardFrames{trialIx}) = nan;
    end

    [IMc, peaks{trialIx}] = localizeSourcesSLAP2(rawIMs{trialIx}, aData, params);

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

%align trial images to template
disp('Aligning across trials...')
meanAligned = []; actAligned = []; 
corrCoeff = nan(1,nTrials);
motOutput = nan(2,nTrials);
Mpad = nan([size(template) size(M,3)]);
Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
for trialIx = nTrials:-1:1
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
corrThresh = min(0.90, median(ccf, 'omitnan')-2*std(ccf, 'omitmissing'));
validTrials= find(ccf>corrThresh);
exptSummary.meanIM = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');
actIM = mean(actAligned(:,:,:,validTrials), 4, 'includenan');
exptSummary.actIM = actIM;

%accumulate peaks, only from valid trials
peaksCat = struct;
for vTrialIx = 1:length(validTrials)
    trialIx = validTrials(vTrialIx);
    if ~isfield(peaksCat, 'row')
        peaksCat.row = peaks{trialIx}.row - motOutput(1,trialIx);
        peaksCat.col = peaks{trialIx}.col - motOutput(2,trialIx);
        peaksCat.val = peaks{trialIx}.val;
    else
        peaksCat.row = cat(1, peaksCat.row, peaks{trialIx}.row - motOutput(1,trialIx));
        peaksCat.col = cat(1, peaksCat.col, peaks{trialIx}.col - motOutput(2,trialIx));
        peaksCat.val = cat(1, peaksCat.val, peaks{trialIx}.val);
    end
end

%select sources
%strategy 1: find peaks directly on aligned activity image
actIM = mean(actAligned(:,:,:,validTrials), 4, 'includenan');
pIM = actIM == ordfilt2(actIM, 9, ones(3)); 
p = actIM(pIM);
sortedP = sort(p, 'descend');
totalPix = sum(~isnan(actIM(:)));
threshP = 2*sortedP(ceil(totalPix/100)); %maximum synapse density
pp = actIM; pp(~pIM) = 0; pp(pp<threshP) = 0;
[sources.R,sources.C,sources.V] = find(pp);
sz = size(pp);
k = length(sources.R);

%Generate dFsel and F0selDS; the data only in the selected region, aligned across movies
selPix = false(sz(1:2));
for sourceIx = k:-1:1
    rr = max(1, round(sources.R(sourceIx)-params.dXY)):min(sz(1),  round(sources.R(sourceIx))+params.dXY);
    cc = max(1, round(sources.C(sourceIx)-params.dXY)):min(sz(2),  round(sources.C(sourceIx))+params.dXY);
    selPix(rr,cc,sourceIx) = true;
end
F0selDS = cell(nTrials,1); dFsel = cell(1,nTrials);
for trialIx = 1:nTrials %parfor
    if any(validTrials==trialIx)
        [dFsel{trialIx}, F0selDS{trialIx}] = DFselAsync(rawIMs{trialIx}, discardFrames{trialIx}, selPix, motOutput(:,trialIx), params); %path, discardFrames, selPix, motOutput, params
    end
end
dFsel = cell2mat(dFsel); %collapse into an array


% for trialIx = validTrials
%     szTmp = size(rawIMs{trialIx});
%     IMrawSel = interpArray(rawIMs{trialIx}, any(selPix,3), motOutput(:,trialIx)); %interpolates the movie at the shifted coordinates
%     F0selDS{trialIx} = svdF0(IMrawSel', 5, baselineWindow, 1, params.denoiseWindow_samps)'; %#ok<AGROW>
%     dFsel(:,frameInd+(1:szTmp(3))) = IMrawSel - F0selDS{trialIx};
%     if numChannels == 2
%         IM2rawSel = interpArray(rawIM2s{trialIx}, any(selPix,3), motOutput(:,trialIx)); %interpolates the movie at the shifted coordinates
%         F02selDS{trialIx} = svdF0(IM2rawSel', 5, baselineWindow, 1, params.denoiseWindow_samps)';
%         dF2sel(:,frameInd+(1:szTmp(3))) = IM2rawSel - F02selDS{trialIx};
%     end
%     frameInd = frameInd+szTmp(3);
% end
clear rawIMs

%extract sources from the downsampled movies
[W0,~] = extractSourcesLoRes(dFsel, sources, selPix, params);

%for each file, process the high res data
params.tau_full=params.tau_s*params.analyzeHz;
E = cell(nTrials,1);
parfor trialIx = 1:nTrials
    if any(validTrials==trialIx)
        E{trialIx} = processTrialAsync_Bergamo(dr, fnRaw{trialIx}, [], [], W0, F0selDS{trialIx}, selPix, discardFrames{trialIx}, [], [], motOutput(:,trialIx), [], params);
    end
end

%prepare file for saving
exptSummary.E = E;
exptSummary.fns = fns;
exptSummary.dr = dr;

%per-trial images
exptSummary.originalActivityChannel = params.activityChannel;
exptSummary.peaks = peaks;
exptSummary.perTrialMeanIMs = meanIM;
exptSummary.perTrialActIms = actIM;
exptSummary.perTrialAlignmentOffsets = motOutput; %the alignment vector for each trial

%save
save([drsave filesep fnsave], 'exptSummary','-v7.3');

disp('Done summarizeBergamo_NoLoCo')
end

function [W0,H0] = extractSourcesLoRes(dFsel, sources, selPix, params)
%dFsel: delta fluorescence over selected pixels;  has dimensions [pixels time]

k = length(sources.R);
sz = size(selPix, [1 2]);
anySel = any(selPix,3);
nSelPix = sum(anySel(:));

%Temporal filter the movie
tau = params.tau_s./(params.frametime*params.dsFac); %time constant in frames
dFselTf = matchedExpFilter(dFsel, tau);
selNans = imclose(isnan(dFselTf), ones(1, 2*ceil(params.denoiseWindow_s/params.frametime)+1));
meanDF = repmat(mean(dFselTf,2, 'omitnan'), 1, size(dFsel,2));
meanDF(isnan(meanDF)) = 0;
dFselTf(selNans) = meanDF(selNans);

nComp = k; %we could use extra components for background if desired
W0 = zeros(nSelPix, nComp);
for sourceIx = 1:k
    %create a selection vector for indexing into the smaller dFsel matrix
    selSel = selPix(:,:,sourceIx);
    [rr,cc] = find(selSel);
    spatialTemplate = mvnpdf(cat(2,rr,cc), [sources.R(sourceIx) sources.C(sourceIx)], params.sigma_px*eye(2)); % we will score each source based on how much the average event looks like a centered gaussian
    selSel = selSel(anySel);
    W0(selSel,sourceIx) = spatialTemplate./max(spatialTemplate(:));
end

%discard sources with poor spatial profiles/response amplitudes
nComp = size(W0,2);
W0full = zeros(sz(1)*sz(2),nComp);
W0full(anySel(:),:) = W0;
W0full = reshape(W0full, sz(1),sz(2),[]);
W00full = imgaussfilt(double(W0full>0), 1).^3;
W00full(W0full==0) = 0;

%Use multiplicative updates NMF, which makes it easy to zero out pixels
opts1 = statset('MaxIter', 4,  'Display', 'final');
[W0,H0] = nnmf2(dFselTf, nComp,'algorithm', 'mult', 'w0',W0, 'options', opts1); %!!nnmf has been modified to allow it to take more than rank(Y) inputs
for bigIter = 1:params.nmfIter
    H0 = H0+sqrt(0.1/size(H0,2)); %get rid of any zeroing out effect on activity
    disp(['outer loop ' int2str(bigIter) ' of ' int2str(params.nmfIter)]);   
    W0 = max(0, W0-params.sparseFac.*max(W0,[],1)); %sparsify
    [W0,H0] = nnmf2(dFselTf, nComp,'algorithm', 'mult', 'w0', W0, 'h0', H0, 'options', opts1); %'h0', H0
    
    % W0full = zeros(sz(1)*sz(2),nComp);
    % W0full(anySel(:),:) = W0;
    %W0full = reshape(W0full, sz(1),sz(2),[]);
    %W0full = min(W0full, imgaussfilt(W0full, params.sigma_px/2)); %sculpt the spatial profiles
    % W0full = reshape(W0full, sz(1)*sz(2),[]);
    %[W0,H0] = nnmf2(dFselTf, nComp,'algorithm', 'mult', 'w0', W0full(anySel,:), 'h0', H0, 'options', opts1);

    if bigIter == 2
        disp('Merging sources...')
        [W0,H0] = mergeSources_2(W0,H0,dFselTf, anySel);
        disp(['Kept ' int2str(size(W0,2)) ' of ' int2str(nComp) ' sources']);
        nComp = size(W0,2);

        W0full = zeros(sz(1)*sz(2),nComp);
        W0full(anySel(:),:) = W0;
        W0full = reshape(W0full, sz(1),sz(2),[]);
        W00full = imgaussfilt(double(W0full>0), 1).^3;
        W00full(W0full==0) = 0;
    else
        W0full = zeros(sz(1)*sz(2),nComp);
        W0full(anySel(:),:) = W0;
        W0full = reshape(W0full, sz(1),sz(2),[]);
        W0full = W0full.*W00full; %shrink edges
        W0full = reshape(W0full, sz(1)*sz(2),[]);
        W0 =  W0full(anySel,:);
    end
end

% W0full = zeros(sz(1)*sz(2),nComp);
% W0full(anySel(:),:) = W0;
% W0full = reshape(W0full, sz(1),sz(2),[]);
% % W00 = imgaussfilt(double(W0full>0), 3).^3;
% % W00(W0full==0) = 0;
% figure, imagesc(sum(W0full,3));
%evaluate how well fit the data are
%figure, imshow3D(cat(3,dFselTf, W0*H0, dFselTf-W0*H0));
end

function [dFsel, F0selDS] = DFselAsync(rawIM, discardFrames, selPix, motOutput, params)
    baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));

    rawIM(:,:,discardFrames) = nan;
    
    szTmp = size(rawIM);
    sel2D = any(selPix,3); 
    sel2D = sel2D(1:szTmp(1), 1:szTmp(2));
    
    IMrawSel = interpArray(rawIM, sel2D, motOutput); %interpolates the movie at the shifted coordinates
    F0selDS =   smoothdata(IMrawSel,2, 'movmean',baselineWindow,'omitmissing');%svdF0_slap2(IMrawSel', 8, denoiseWindow, baselineWindow)';
    dFsel = IMrawSel - F0selDS;
end
