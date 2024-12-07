function summarize_NoLoCo(dr_or_pathToTrialTable, paramsIn)
%PARAMETER SETTING
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('summarize_NoLoCo', paramsIn);
else
    params = setParams('summarize_NoLoCo');
end
if ~nargin
    [trialTablefn, dr] = uigetfile('*trialTable*.mat');
else
    %parse dr
    %_or_pathToTrialTable
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        trialTablefn = 'trialTable.mat';
    else
        [dr trialTablefn ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end

copyReadDeleteScanImageTiff([]); %make sure we can use the function in parallel loops 

%confirm that all files exist
[trialTable, keepTrials] = verifyFiles(trialTablefn, dr, params);
nDMDs = size(trialTable.filename,1); %the trial table has size #DMDs x # trials; Bergamo is treated as '1 DMD'
nTrials = size(trialTable.filename,2);
firstValidTrial = find(all(keepTrials,1),1,'first');

%parameters that depend only on the microscope, hidden from GUI
switch params.microscope
    case 'SLAP2'
        params.discardInitial_s = 0.1;
        trialTable.fnRaw = trialTable.filename;
    case 'bergamo'
        params.discardInitial_s = 0;
        params.analyzeHz = nan;
end

disp(['## SUMMARIZING' newline 'Folder:'])
disp(dr)

savedr = [dr filesep 'ExperimentSummary'];
% on CodeOcean /data is read-only and we save to /results
is_CodeOcean = ~(getenv("CO_CPUS") == "");
if is_CodeOcean
    savedr = strrep(savedr, '/data', '/results');
end
if ~exist(savedr, 'dir')
    mkdir(savedr);
end
fnsave = [savedr filesep 'Summary-' datestr(now, 'YYmmDD-HHMMSS') '.mat'];

%call up a GUI for the user to define Soma ROI and regions to exclude
if params.drawUserRois
    fnAnn = [dr filesep 'ANNOTATIONS.mat'];
    if exist(fnAnn, 'file')
        load(fnAnn, 'ROIs')
    else
        for DMDix = 1:nDMDs
            %load image data
            [~, fn, ext] = fileparts(trialTable.fnRegDS{DMDix,firstValidTrial});
            IM = copyReadDeleteScanImageTiff([dr filesep fn ext]);
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
else
    ROIs = [];
end

%load some metadata
[~, fn, ext] = fileparts(trialTable.fnAdata{1,firstValidTrial});
load([dr filesep fn ext], 'aData');
numChannels = aData.numChannels;
params.numChannels = numChannels;
params.alignHz = aData.alignHz;
if ~strcmpi(params.microscope, 'SLAP2')
    params.analyzeHz = 1/aData.frametime; %analyze conventional recordings at the acquisitoin framerate
end

%PROCESS DATA
for DMDix = nDMDs:-1:1
    %set up parallelization
    if params.nParallelWorkers>1
        p = gcp('nocreate');
        if isempty(p)
            poolsize = 0;
        else
            poolsize = p.NumWorkers;
        end
        nWorkers = min(params.nParallelWorkers, size(trialTable.filename,2));
        if poolsize~=nWorkers
            delete(gcp('nocreate'));
            parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
        end
    else
        delete(gcp('nocreate'));
    end

    %Perform Localizations
    disp('Loading data and performing localizations...')
    mIM = cell(1, nTrials); aIM = cell(1,nTrials); alignData = cell(1, nTrials); peaks = cell(1, nTrials); discardFrames = cell(1,nTrials); %rawIMs = cell(1,nTrials)
    fns = trialTable.fnRegDS(DMDix, :);
    parfor trialIx = 1:nTrials
        if keepTrials(DMDix,trialIx)
            [~,fn, ext] = fileparts(fns{trialIx}); fn = [fn ext];
            [~, mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks{trialIx}, discardFrames{trialIx}]= loadAndProcessTrialAsync(dr, fn, numChannels, params); %rawIMs{trialIx}
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
    meanAligned = []; actAligned = []; pAligned = [];
    corrCoeff = nan(1,nTrials);
    motOutput = nan(2,nTrials);
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
    corrThresh = min(0.90, median(ccf, 'omitnan')-2*std(ccf, 'omitmissing'));
    actValidPix = squeeze(mean(~isnan(actAligned(:,:,1,:)), [1 2]));
    validTrials= find(ccf(:)>corrThresh & actValidPix(:)>mean(actValidPix)/2);
    exptSummary.meanIM{DMDix} = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');
    actIM = mean(actAligned(:,:,:,validTrials), 4, 'omitnan');
    nanFrac = mean(isnan(actAligned(:,:,:,validTrials)), 4);
    actIM(nanFrac>0.6) = nan;
    exptSummary.actIM{DMDix} = actIM;

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
    %Mask out somata from activity image
    somaMask = false(size(actIM));
    if ~isempty(ROIs)
        for rix = 1:numel(ROIs(DMDix).roiData)
            if contains(upper(ROIs(DMDix).roiData{rix}.Label), 'SOMA')
                tmp = ROIs(DMDix).roiData{rix}.mask;
                somaMask(1:size(tmp,1), 1:size(tmp,2)) = somaMask(1:size(tmp,1), 1:size(tmp,2)) | tmp;
            end
        end
    end
    pIM(somaMask) = 0;
    p = actIM(pIM);
    sortedP = sort(p, 'descend');
    totalPix = sum(~isnan(actIM(:)) & ~somaMask(:));
    if totalPix == 0
        k = 0;
    else
        threshP = 2*sortedP(min(end, ceil(totalPix*params.maxSynapseDensity))); %maximum synapse density
        pp = actIM; pp(~pIM) = 0; pp(pp<threshP) = 0;
        [sources.R,sources.C,sources.V] = find(pp);
        sz = size(pp);
        k = length(sources.R);
    end
    %strategy 2: use peaks
    %not implemented

    if k>0
    %Generate IMsel; the data only in the selected region, aligned across movies
    selPix = false([sz(1:2) k]);
    for sourceIx = k:-1:1
        rr = max(1, round(sources.R(sourceIx)-params.dXY)):min(sz(1),  round(sources.R(sourceIx))+params.dXY);
        cc = max(1, round(sources.C(sourceIx)-params.dXY)):min(sz(2),  round(sources.C(sourceIx))+params.dXY);
        selPix(rr,cc,sourceIx) = true;
    end
    pxAlwaysValid = mean(isnan(meanAligned(:,:,1,validTrials)),4)<params.nanThresh;
    selPix = selPix & repmat(pxAlwaysValid, 1, 1, k); %ADJUST SELECTED PIXELS NOT TO INCLUDE POORLY MEASURED PIXELS

    %prune any sources that got clipped by pixel selection process
    keepSources = sum(selPix, [1 2])>5;
    k = sum(keepSources, 'all');
    sources.R = sources.R(keepSources);
    sources.C = sources.C(keepSources);
    selPix = selPix(:,:,keepSources);

    %accumulate dF of selected pixels for initial NMF of downsampled movies
    F0selDS = cell(nTrials,1); dFsel = cell(1,nTrials);
    fns = trialTable.fnRegDS(DMDix, :);
    parfor trialIx = 1:nTrials
        if any(validTrials==trialIx)
            [~,fn, ext] = fileparts(fns{trialIx}); fn = [fn ext];
            [dFsel{trialIx}, F0selDS{trialIx}] = DFselAsync([dr filesep fn], discardFrames{trialIx}, selPix, motOutput(:,trialIx), params); %path, discardFrames, selPix, motOutput, params
        end
    end
    dFsel = cell2mat(dFsel); %collapse into an array

    %to do:
    %option to use only top 10% most active frames
    %extract sources from downsampled movie dF
    try
        [W0,~] = extractSourcesLoRes(dFsel, sources, selPix, params);
    catch ME
        disp('fatal Error extracting sources')
        rethrow(ME)
    end

    %for each file, load high res data and refine
    params.tau_full=params.tau_s*params.analyzeHz;
    if isempty(ROIs) || isempty(ROIs(DMDix))
        roiData =[];
    else
        roiData = ROIs(DMDix).roiData;
    end
    E = cell(nTrials,1);
    fns = trialTable.fnRaw(DMDix,:);
    if size(W0, 2) > 0
    if strcmpi(params.microscope, 'SLAP2')
        if params.nParallelWorkers>1
            newN = min(min(24,nWorkers*5), size(trialTable.filename,2)); %this step is not very memory-demanding for SLAP2; increase parallel workers
            if nWorkers~=newN
                delete(gcp('nocreate'));
                parpool('processes',newN); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
            end
        end

        fls = trialTable.firstLine(DMDix,:);
        els = trialTable.lastLine(DMDix,:);
        parfor trialIx = 1:nTrials
            if any(validTrials==trialIx)
                fnRaw = fns{trialIx};
                E{trialIx} = processTrialAsync(dr, fnRaw, fls(trialIx), els(trialIx), W0, F0selDS{trialIx}, selPix, discardFrames{trialIx}, alignData{trialIx}, mIM{trialIx}, motOutput(:,trialIx), roiData, params);
            end
        end
    else %BERGAMO
        parfor trialIx = 1:nTrials
            if any(validTrials==trialIx)
                fnRaw = fns{trialIx};
                E{trialIx} = processTrialAsync(dr, fnRaw, [], [], W0, F0selDS{trialIx}, selPix, discardFrames{trialIx}, alignData{trialIx}, mIM{trialIx}, motOutput(:,trialIx), roiData, params);
            end
            end
        end
    end

    %per-trial images
    exptSummary.E(:,DMDix) = E; %experiment data
    else
    roiData =[];
    end
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

%save
save(fnsave, 'exptSummary', "-v7.3");
disp('Done summarize_NoLoCo')
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

function [W0,H0] = extractSourcesLoRes(dFsel, sources, selPix, params)
%dFsel: delta fluorescence over selected pixels;  has dimensions [pixels time]

k = length(sources.R);
sz = size(selPix, [1 2]);
anySel = any(selPix,3);
nSelPix = sum(anySel(:));

%Temporal filter the movie
tau = params.tau_s.*params.alignHz; %time constant in frames
dFselTf = matchedExpFilter(dFsel, tau);

selNans = imclose(isnan(dFselTf), ones(1, 2*ceil(params.denoiseWindow_s.*params.alignHz)+1));
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

%W0full = min(W0full, imgaussfilt(W0full, params.sigma_px/2));
%W0full = reshape(W0full, sz(1)*sz(2),[]);
%W0 =  W0full(anySel,:);

% H0 = W0\dFselTf;
% nActive = sum(H0>prctile(H0,99,2),1);
% selTime = nActive>

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

function [dFsel, F0selDS] = DFselAsync(path, discardFrames, selPix, motOutput, params)
baselineWindow = ceil(params.baselineWindow_Glu_s .* params.alignHz);
denoiseWindow = ceil(params.denoiseWindow_s .* params.alignHz); %params.denoiseWindow_samps;

rawIM = copyReadDeleteScanImageTiff(path);
rawIM = reshape(rawIM, size(rawIM,1), size(rawIM,2), params.numChannels, []); %deinterleave;
rawIM = squeeze(rawIM(:,:,1,:));
rawIM(:,:,discardFrames) = nan;

szTmp = size(rawIM);
sel2D = any(selPix,3);
sel2D = sel2D(1:szTmp(1), 1:szTmp(2));

IMrawSel = interpArray(rawIM, sel2D, motOutput); %interpolates the movie at the shifted coordinates
F0selDS =   smoothdata(IMrawSel,2, 'movmean',baselineWindow,'omitmissing');%svdF0_slap2(IMrawSel', 8, denoiseWindow, baselineWindow)';
dFsel = IMrawSel - F0selDS;
end