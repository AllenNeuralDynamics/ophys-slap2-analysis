function summarize_LoCo(dr_or_pathToTrialTable, paramsIn)
%PARAMETER SETTING
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('summarize_LoCo', paramsIn);
else
    params = setParams('summarize_LoCo');
end
if ~nargin
    [trialTablefn, dr] =  uigetfile('*.mat', 'Select a trialTable file', '*trialTable*.mat' );
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

if params.makeJSON
    pythonenv_dir = uigetdir(getenv("USERPROFILE"),'Select Python Environment Directory');
    disp(['python env: ' pythonenv_dir])
end

params.startTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

copyReadDeleteScanImageTiff([]); %make sure we can use the function in parallel loops

%confirm that all files exist
[trialTable, keepTrials] = verifyFiles(trialTablefn, dr, params);
% for dmdIx = 1:numel(trialTable.refStack)
%     trialTable.refStack{dmdIx}.IM = []; %this uses a lot of memory and we won't need it
%end
nDMDs = size(trialTable.filename,1); %the trial table has size #DMDs x # trials; Bergamo is treated as '1 DMD'
nTrials = size(trialTable.filename,2);

%parameters that depend only on the microscope, hidden from GUI
switch params.microscope
    case 'SLAP2'
        trialTable.fnRaw = trialTable.filename;
    case 'bergamo'
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
fnsave = [savedr filesep 'SummaryLoCo-' datestr(now, 'YYmmDD-HHMMSS') '.mat'];

% Persist authoritative SLAP2 acquisition metadata before any trial-level
% processing. This mirrors the voltage extractor so downstream Python can
% reconcile source acquisition epochs to HARP without inferring duration from
% processed traces or rescaling the imaging timebase.
if strcmpi(params.microscope, 'SLAP2')
    slap2AcqInfo = collectSlap2AcquisitionInfo(trialTable, dr);
    exptSummary.trialEpoch = slap2AcqInfo.trialEpoch;
    exptSummary.trialFilePrefix = slap2AcqInfo.trialFilePrefix;
    exptSummary.epochTable = slap2AcqInfo.epochTable;
    exptSummary.nEpochs = slap2AcqInfo.nEpochs;
    exptSummary.multiEpochAcquisition = slap2AcqInfo.nEpochs > 1;
    exptSummary.dmd = slap2AcqInfo.dmd;
    exptSummary.acquisitionMetadataSchemaVersion = '1.0';
else
    exptSummary.trialEpoch = ones(1, nTrials);
    exptSummary.trialFilePrefix = repmat({''}, 1, nTrials);
    exptSummary.epochTable = table(1, {''}, 1, nTrials, nTrials, ...
        'VariableNames', {'epochIdx', 'filePrefix', 'firstTrial', 'lastTrial', 'nTrials'});
    exptSummary.nEpochs = 1;
    exptSummary.multiEpochAcquisition = false;
    exptSummary.dmd = repmat(struct(), 1, nDMDs);
    exptSummary.acquisitionMetadataSchemaVersion = 'not_slap2';
end

%call up a GUI for the user to define Soma ROI and regions to exclude
if params.drawUserRois
    fnPostAcqAnn = [dr filesep 'postAcqANNOTATIONS.mat'];
    fnAnn = [dr filesep 'ANNOTATIONS.mat'];
    ROIs = [];
    if exist(fnPostAcqAnn, 'file')
        annData = load(fnPostAcqAnn);
        if isfield(annData, 'ROIs')
            ROIs = annData.ROIs;
        else
            warning('postAcqANNOTATIONS.mat exists but does not contain ROIs. Opening ROI drawing GUI.');
        end
    end

    if ~isempty(ROIs)
        updatedAnnotations = false;
        for DMDix = 1:nDMDs
            %load image data
            firstValidTrial = find(keepTrials(DMDix,:),1,"first");
            [~, fn, ext] = fileparts(trialTable.fnRegDS{DMDix,firstValidTrial});
            ROIs(DMDix).dr = dr;
            ROIs(DMDix).fn = fn;
            targetKey = trialTable.fnRegDS{DMDix,firstValidTrial};

            IMtargetRaw = copyReadDeleteScanImageTiff([dr filesep fn ext]);
            IMtarget = collapseImageForRegistration(IMtargetRaw);
            [IMann, hasAnnotationImage] = extractAnnotationImage(ROIs(DMDix));
            alreadyAligned = isfield(ROIs(DMDix), 'alignedToFnRegDS') && strcmp(ROIs(DMDix).alignedToFnRegDS, targetKey);
            hasRoiData = isfield(ROIs(DMDix), 'roiData') && ~isempty(ROIs(DMDix).roiData);
            if hasAnnotationImage && hasRoiData && ~alreadyAligned
                [shiftRC, ok] = estimateRegistrationShift(IMtarget, IMann);
                if ok
                    newRoiData = shiftRoiData(ROIs(DMDix).roiData, shiftRC, size(IMtarget,[1 2]));
                    ROIs(DMDix).roiData = newRoiData;
                    ROIs(DMDix).alignedToFnRegDS = targetKey;
                    ROIs(DMDix).alignmentShiftRC = shiftRC;
                    updatedAnnotations = true;
                else
                    warning(['Could not register ANNOTATIONS image for DMD' int2str(DMDix) '. Using stored ROI coordinates as-is.']);
                end
            end
        end
        save(fnAnn, 'ROIs');
        if ~updatedAnnotations
            warning('No post-acquisition annotation transforms were applied.');
        end
    else
        if ~exist(fnPostAcqAnn, 'file')
            warning('postAcqANNOTATIONS.mat not found. Opening ROI drawing GUI.');
        end

        for DMDix = 1:nDMDs
            %load image data
            firstValidTrial = find(keepTrials(DMDix,:),1,"first");
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

%PROCESS DATA
for DMDix = nDMDs:-1:1
    %load some metadata
    firstValidTrial = find(keepTrials(DMDix,:),1,"first");
    fn = trialTable.fnAdata{DMDix,firstValidTrial};
    load([dr filesep fn], 'aData');
    numChannels = aData.numChannels;
    params.numChannels = numChannels;
    params.alignHz = aData.alignHz;
    if ~strcmpi(params.microscope, 'SLAP2')
        params.analyzeHz = 1/aData.frametime; %analyze conventional recordings at the acquisitoin framerate
    end
    if isfield(aData, 'Z')
        exptSummary.Z(DMDix) = aData.Z;
    else
        exptSummary.Z(DMDix) = nan;
        warning('Alignment data missing Z-plane, likely out of date!!')
    end
    clear aData

    %set up parallelization
    if params.nParallelWorkers>1
        p = gcp('nocreate');
        if isempty(p)
            poolsize = 0;
        else
            poolsize = p.NumWorkers;
        end
        dd = dir([dr filesep trialTable.fnRegDS{DMDix, firstValidTrial}]);
        try
            fileSize = dd.bytes;
        catch
            error(['Error loading registered tiff:' trialTable.fnRegDS{DMDix, firstValidTrial} '\n' 'Are paths in your trial table valid?']);
        end
        if ispc
            userMemInfo = memory;
            memAvailable = userMemInfo.MemAvailableAllArrays;
        else
            [~, result] = unix('grep MemAvailable /proc/meminfo | awk ''{print $2}''');
            memAvailable = str2double(result) * 1024;  % Convert KB to bytes
        end
        maxWorkers = max(1,min(size(trialTable.filename,2), floor(0.13*memAvailable/fileSize)));
        nWorkers = min(params.nParallelWorkers, maxWorkers);
        
        if poolsize~=nWorkers ||  ~strcmpi(class(p), 'parallel.ProcessPool')
            delete(gcp('nocreate'));
            parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM 
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
            [~, mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks{trialIx}, discardFrames{trialIx}]= loadAndProcessTrialAsync(dr, fns{trialIx}, numChannels, params); %rawIMs{trialIx}
        end
    end
    %Assemble same-sized mean images from different-sized trial means
    szm1 = max(cellfun(@(x)size(x,1),mIM)); szm2 = max(cellfun(@(x)size(x,2), aIM));
    meanIM = nan(szm1,szm2,numChannels, nTrials); activIM = nan(szm1,szm2,1, nTrials);
    for trialIx = 1:nTrials
        tmp =  mIM{trialIx};
        meanIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
        tmp =  aIM{trialIx};
        activIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
    end
    params.sz = size(meanIM, [1 2]);

    %Make template
    disp('Making template for aligning across trials...')
    maxshift = 5;
    M = squeeze(sum(meanIM, 3));
    samples = find(keepTrials(DMDix,:)); samples = samples(unique(round(linspace(1,length(samples),20))));
    template = makeTemplateMultiRoi(M(:,:,samples), maxshift);

    %align all mean images to template
    disp('Aligning across trials...')
    meanAligned = [];
    actAligned = nan(size(meanIM,1), size(meanIM,2),1,nTrials);
    corrCoeff = nan(1,nTrials);
    motOutput = nan(2,nTrials);
    Mpad = nan([size(template) size(M,3)]);
    Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
    %clear M

    fillval = min(template(:),[], 'omitnan')-1;
    tFFT = fft2(max(template, fillval));
    for trialIx = nTrials:-1:1
        if ~keepTrials(DMDix,trialIx) || all(isnan(activIM(:,:,1,trialIx)), 'all')
            disp(['skipping trial, dmd:' int2str(trialIx) ' ' int2str(DMDix)])
            continue %skip
        end
        disp(['trial: ' int2str(trialIx)])
        
        output1 = dftregistration_clipped(tFFT, fft2(max(Mpad(:,:,trialIx), fillval)),1,80);
        mot1 = [-output1(3) -output1(4)]; %xcorr2_nans(Mpad(:,:,trialIx), template, [-output1(3) ; -output1(4)], maxshift);
        [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
        [rr,cc] = ndgrid(1:size(meanIM,1), 1:size(meanIM,2));

        for chIx = 1:size(meanIM,3)
            meanAligned(:,:,chIx,trialIx) = interp2(meanIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
        end
        actAligned(:,:,1,trialIx) = interp2(activIM(:,:,1,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
    end
    %clear Mpad activIM

    %identify outliers in alignment quality to determine valid trials
    ccf = corrCoeff;
    corrThresh = min(0.90, median(ccf, 'omitnan')-2*std(ccf, 'omitmissing'));
    actValidPix = squeeze(mean(~isnan(actAligned(:,:,1,:)), [1 2]));
    validTrials= find(ccf(:)>corrThresh & actValidPix(:)>mean(actValidPix)/2);
    exptSummary.meanIM{DMDix} = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');
    actIM = prctile(actAligned(:,:,:,validTrials),80,4);  %mean(actAligned(:,:,:,validTrials), 4, 'omitnan');
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
    medIM = nanmedfilt2(actIM, (2*ceil(1.5*params.dXY)+1).*[1 1]);
    actIM = actIM-medIM; %subtract a local baseline
    explored = actIM; pTmp = explored>0 & explored == ordfilt2(explored, 9, ones(3));
    pIM = false(size(actIM));
    while any(pTmp(:))
        pIM = pIM | pTmp;
        explored(imdilate(pTmp, ones(5))) = 0;
        pTmp = explored>0 & explored == ordfilt2(explored, 9, ones(3));
    end

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
    if totalPix == 0 | isempty(p)
        k = 0;
    else
        threshP = 2*sortedP(min(end, ceil(totalPix*params.maxSynapseDensity))); %maximum synapse density
        pp = actIM; pp(~pIM) = 0; pp(pp<threshP) = 0;
        [sources.R,sources.C,sources.V] = find(pp);
        sz = size(pp);
        k = length(sources.R);
    end

    %select regions near synapses, aligned across movies
    selPix = false([sz(1:2) k]);
    params.selRadius = ceil(2*params.dXY);
    for sourceIx = k:-1:1
        rr = round(sources.R(sourceIx));
        cc = round(sources.C(sourceIx));
        selPix(rr,cc,sourceIx) = true;
        selPix(:,:,sourceIx) = imdilate(selPix(:,:,sourceIx), strel('disk',params.selRadius));
    end
    pxAlwaysValid = mean(isnan(meanAligned(:,:,1,validTrials)),4)<params.nanThresh;
    selPix = selPix & repmat(pxAlwaysValid, 1, 1, k); %ADJUST SELECTED PIXELS NOT TO INCLUDE POORLY MEASURED PIXELS

    %prune any sources that got clipped by pixel selection process
    keepSources = sum(selPix, [1 2])>5;
    sources.R = sources.R(keepSources);
    sources.C = sources.C(keepSources);
    selPix = selPix(:,:,keepSources);
    disp(['Number of sources: ' int2str(sum(keepSources))]);
        
    %for each file, load high res data and refine
    params.tau_full=params.tau_s*params.analyzeHz;
    params = setParamsExtractTrial(params);
    
    if isempty(ROIs) || isempty(ROIs(DMDix))
        roiData =[];
    else
        roiData = ROIs(DMDix).roiData;
    end

    if any(keepSources)
        fns = trialTable.fnRaw(DMDix,:);
            if strcmpi(params.microscope, 'SLAP2')
                fls = trialTable.firstLine(DMDix,:);
                els = trialTable.lastLine(DMDix,:);
                E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params);
            else %BERGAMO
                fls = cell(1,numel(fns)); %first frame; leave empty for most uses
                els = cell(1,numel(fns)); %last frame; leave empty for most uses
                E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params);
            end

        %per-trial images
        exptSummary.E(:,DMDix) = E; %experiment data
    end
    exptSummary.selPix{DMDix} = any(selPix,3);
    exptSummary.aData(:,DMDix) = alignData;
    exptSummary.userROIs{DMDix} = roiData;
    exptSummary.peaks{DMDix}= peaks;
    exptSummary.perTrialMeanIMs{DMDix} = meanIM;
    exptSummary.perTrialMeanIMsAligned{DMDix} = meanAligned;
    exptSummary.perTrialActIms{DMDix} = actIM;
    exptSummary.perTrialActIMsAligned{DMDix} = actAligned;
    exptSummary.perTrialAlignmentOffsets{DMDix} = motOutput; %the alignment vector for each trial

    clear meanAligned meanIM actAligned F0selDS E
end

params.endTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

%prepare file for saving
exptSummary.params = params;
exptSummary.trialTable = trialTable;
exptSummary.dr = dr;

%save
save(fnsave, 'exptSummary', "-v7.3");

if params.makeJSON
    try
        setenv("PYTHONHOME",pythonenv_dir)
        pyrunfile([fullfile(fileparts(mfilename('fullpath')), 'generate_processing_json_SLAP2_multiROI_raster.py') ' --mat_path "' fnsave '" --output_dir "' savedr '"']);
        disp('Saved processing.json')
    catch
        disp('Did not save processing.json')
    end
else
    disp('Did not save processing.json')
end

disp('Done summarize_LoCo')
end


function info = collectSlap2AcquisitionInfo(trialTable, dr)
%COLLECTSLAP2ACQUISITIONINFO Preserve raw SLAP2 acquisition/epoch metadata.
%
% The trial-wise SummaryLoCo traces are sampled at params.analyzeHz, but the
% authoritative acquisition duration lives in the raw SLAP2 files. For each DMD
% and acquisition epoch, open the corresponding MultiDataFiles object and retain
% actual nCycles/totalNumLines plus line-rate/parse-plan metadata. Downstream
% Python uses totalNumLines / lineRateHz for source-epoch duration QC, matching
% the voltage pipeline.

hasFilename = (isstruct(trialTable) && isfield(trialTable, 'filename')) || ...
    (istable(trialTable) && any(strcmp(trialTable.Properties.VariableNames, 'filename')));
if ~hasFilename
    error('summarize_LoCo:MissingFilename', ...
        'SLAP2 trialTable must contain filename to recover acquisition metadata.');
end
filename = trialTable.filename;
if isstring(filename)
    filename = cellstr(filename);
end
if ~iscell(filename)
    error('summarize_LoCo:InvalidFilename', ...
        'trialTable.filename must be a cell array or string array.');
end
[nDMDs, nTrials] = size(filename);
[trialEpoch, trialFilePrefix, epochTable] = inferSlap2AcquisitionEpochs(filename, trialTable);
nEpochs = height(epochTable);

info = struct();
info.trialEpoch = trialEpoch;
info.trialFilePrefix = trialFilePrefix;
info.epochTable = epochTable;
info.nEpochs = nEpochs;
info.dmd = repmat(struct(), 1, nDMDs);

for dmdIdx = 1:nDMDs
    epochRecords = repmat(struct( ...
        'epochIdx', [], 'filePrefix', '', 'firstTrial', [], 'lastTrial', [], ...
        'firstDatFile', '', 'totalNumLines', 0, 'nCycles', 0, ...
        'linesPerCycle', [], 'durationSec', nan, 'cycleRateHz', nan, ...
        'metadata', struct(), 'available', false), 1, nEpochs);

    canonicalMeta = struct();
    canonicalLinesPerCycle = [];
    firstCanonicalDatFile = '';

    for epochIdx = 1:nEpochs
        epochTrials = find(trialEpoch == epochIdx);
        validTrials = epochTrials(~cellfun(@isempty, filename(dmdIdx, epochTrials)));

        epochRecords(epochIdx).epochIdx = epochIdx;
        epochRecords(epochIdx).firstTrial = epochTrials(1);
        epochRecords(epochIdx).lastTrial = epochTrials(end);
        epochRecords(epochIdx).filePrefix = epochTable.filePrefix{epochIdx};

        if isempty(validTrials)
            warning('summarize_LoCo:MissingDmdEpoch', ...
                'DMD%d has no raw SLAP2 filename for acquisition epoch %d.', ...
                dmdIdx, epochIdx);
            continue
        end

        firstTrialThisEpoch = validTrials(1);
        firstDatFile = resolveSlap2DataFilePath(dr, filename{dmdIdx, firstTrialThisEpoch});
        fprintf('SLAP2 acquisition metadata DMD%d epoch%d: %s\n', ...
            dmdIdx, epochIdx, firstDatFile);
        hMDF = slap2.util.MultiDataFiles(firstDatFile);
        dmdMeta = getSlap2DmdMetadata(hMDF);

        if isempty(fieldnames(canonicalMeta))
            canonicalMeta = dmdMeta;
            canonicalLinesPerCycle = hMDF.header.linesPerCycle;
            firstCanonicalDatFile = firstDatFile;
        else
            warnSlap2EpochMetadataChange( ...
                dmdIdx, epochIdx, canonicalMeta, dmdMeta, ...
                canonicalLinesPerCycle, hMDF.header.linesPerCycle);
        end

        epochRecords(epochIdx).firstDatFile = firstDatFile;
        epochRecords(epochIdx).totalNumLines = double(hMDF.totalNumLines);
        epochRecords(epochIdx).nCycles = double(hMDF.numCycles);
        epochRecords(epochIdx).linesPerCycle = double(hMDF.header.linesPerCycle);
        epochRecords(epochIdx).metadata = dmdMeta;
        epochRecords(epochIdx).available = true;

        lineRateHz = scalarNumericField(dmdMeta, 'lineRateHz');
        if isfinite(lineRateHz) && lineRateHz > 0
            epochRecords(epochIdx).durationSec = double(hMDF.totalNumLines) / lineRateHz;
            if ~isempty(hMDF.header.linesPerCycle) && double(hMDF.header.linesPerCycle) > 0
                epochRecords(epochIdx).cycleRateHz = lineRateHz / double(hMDF.header.linesPerCycle);
            end
        end
    end

    info.dmd(dmdIdx).metadata = canonicalMeta;
    info.dmd(dmdIdx).epochs = epochRecords;
    info.dmd(dmdIdx).firstDatFile = firstCanonicalDatFile;
    info.dmd(dmdIdx).totalNumLines = sum([epochRecords.totalNumLines]);
    info.dmd(dmdIdx).nCycles = sum([epochRecords.nCycles]);
    info.dmd(dmdIdx).linesPerCycle = canonicalLinesPerCycle;
end
end


function [trialEpoch, trialFilePrefix, epochTable] = inferSlap2AcquisitionEpochs(filename, trialTable)
%INFERSLAP2ACQUISITIONEPOCHS Identify source acquisition blocks in trial order.
%
% Explicit trialEpoch/epoch labels are preferred when they are informative.
% Otherwise filename-prefix changes define acquisition restarts. Empty filenames
% are filled from neighboring trials so a missing/invalid trial does not create a
% spurious acquisition epoch.

[~, nTrials] = size(filename);
trialFilePrefix = repmat({''}, 1, nTrials);
for trialIdx = 1:nTrials
    fns = filename(:, trialIdx);
    fns = fns(~cellfun(@isempty, fns));
    if ~isempty(fns)
        trialFilePrefix{trialIdx} = slap2AcquisitionPrefixFromFilename(fns{1});
    end
end

% Fill missing prefixes from surrounding trials; raw-file absence should not by
% itself imply that SLAP2 acquisition restarted.
for trialIdx = 2:nTrials
    if isempty(trialFilePrefix{trialIdx})
        trialFilePrefix{trialIdx} = trialFilePrefix{trialIdx - 1};
    end
end
for trialIdx = nTrials-1:-1:1
    if isempty(trialFilePrefix{trialIdx})
        trialFilePrefix{trialIdx} = trialFilePrefix{trialIdx + 1};
    end
end

explicitEpoch = [];
try
    if istable(trialTable)
        vars = trialTable.Properties.VariableNames;
        if any(strcmp(vars, 'trialEpoch'))
            explicitEpoch = double(trialTable.trialEpoch(:))';
        elseif any(strcmp(vars, 'epoch'))
            explicitEpoch = double(trialTable.epoch(:))';
        end
    elseif isstruct(trialTable)
        if isfield(trialTable, 'trialEpoch')
            explicitEpoch = double(trialTable.trialEpoch(:))';
        elseif isfield(trialTable, 'epoch')
            explicitEpoch = double(trialTable.epoch(:))';
        end
    end
catch
    explicitEpoch = [];
end

finiteExplicit = explicitEpoch(isfinite(explicitEpoch));
if numel(explicitEpoch) == nTrials && numel(unique(finiteExplicit)) > 1
    trialEpoch = relabelSlap2EpochVectorStable(explicitEpoch);
else
    trialEpoch = ones(1, nTrials);
    currentEpoch = 1;
    prevPrefix = trialFilePrefix{1};
    for trialIdx = 2:nTrials
        thisPrefix = trialFilePrefix{trialIdx};
        if ~strcmp(thisPrefix, prevPrefix)
            currentEpoch = currentEpoch + 1;
            prevPrefix = thisPrefix;
        end
        trialEpoch(trialIdx) = currentEpoch;
    end
end

epochIds = unique(trialEpoch, 'stable');
epochIdxCol = epochIds(:);
filePrefixCol = cell(numel(epochIds), 1);
firstTrialCol = zeros(numel(epochIds), 1);
lastTrialCol = zeros(numel(epochIds), 1);
nTrialsCol = zeros(numel(epochIds), 1);
for k = 1:numel(epochIds)
    ep = epochIds(k);
    trials = find(trialEpoch == ep);
    firstTrialCol(k) = trials(1);
    lastTrialCol(k) = trials(end);
    nTrialsCol(k) = numel(trials);
    prefixes = trialFilePrefix(trials);
    prefixes = prefixes(~cellfun(@isempty, prefixes));
    if isempty(prefixes)
        filePrefixCol{k} = '';
    else
        filePrefixCol{k} = prefixes{1};
    end
end

epochTable = table(epochIdxCol, filePrefixCol, firstTrialCol, lastTrialCol, nTrialsCol, ...
    'VariableNames', {'epochIdx', 'filePrefix', 'firstTrial', 'lastTrial', 'nTrials'});
end


function trialEpoch = relabelSlap2EpochVectorStable(epochVals)
%RELABELSLAP2EPOCHVECTORSTABLE Convert arbitrary labels to consecutive 1..N.
trialEpoch = ones(size(epochVals));
seen = [];
current = 1;
for idx = 1:numel(epochVals)
    val = epochVals(idx);
    if ~isfinite(val)
        trialEpoch(idx) = current;
        continue
    end
    match = find(seen == val, 1, 'first');
    if isempty(match)
        seen(end + 1) = val; %#ok<AGROW>
        match = numel(seen);
    end
    current = match;
    trialEpoch(idx) = match;
end
end


function prefix = slap2AcquisitionPrefixFromFilename(fn)
%SLAP2ACQUISITIONPREFIXFROMFILENAME Strip DMD/CYCLE suffixes from raw filename.
[~, name, ~] = fileparts(char(fn));
prefix = regexprep(name, '[_-]DMD\d+.*$', '', 'ignorecase');
prefix = regexprep(prefix, '[_-]CYCLE[_-]?\d+.*$', '', 'ignorecase');
if isempty(prefix)
    prefix = name;
end
end


function datPath = resolveSlap2DataFilePath(dr, filename)
%RESOLVESLAP2DATAFILEPATH Resolve an absolute or session-relative raw .dat path.
filename = char(filename);
if exist(filename, 'file')
    datPath = filename;
elseif exist(fullfile(dr, filename), 'file')
    datPath = fullfile(dr, filename);
else
    error('summarize_LoCo:MissingDataFile', 'Could not find raw SLAP2 data file: %s', filename);
end
end


function dmdMeta = getSlap2DmdMetadata(hMDF)
%GETSLAP2DMDMETADATA Extract compact acquisition/parse-plan metadata.
dmdMeta = struct();
dmdMeta.dmdPixelsPerRow = getStructFieldOrEmpty(hMDF.metaData, 'dmdPixelsPerRow');
dmdMeta.dmdPixelsPerColumn = getStructFieldOrEmpty(hMDF.metaData, 'dmdPixelsPerColumn');
dmdMeta.linePeriod_s = getStructFieldOrEmpty(hMDF.metaData, 'linePeriod_s');
dmdMeta.samplesPerLine = getStructFieldOrEmpty(hMDF.metaData, 'samplesPerLine');
dmdMeta.channelsSave = getStructFieldOrEmpty(hMDF.metaData, 'channelsSave');
dmdMeta.acqDuration_s = getStructFieldOrEmpty(hMDF.metaData, 'acqDuration_s');
dmdMeta.acqDurationCycles = getStructFieldOrEmpty(hMDF.metaData, 'acqDurationCycles');
if isfield(hMDF.metaData, 'AcquisitionContainer') && ...
        isfield(hMDF.metaData.AcquisitionContainer, 'ParsePlan')
    pp = hMDF.metaData.AcquisitionContainer.ParsePlan;
    dmdMeta.parsePlanZs = getStructFieldOrEmpty(pp, 'zs');
    dmdMeta.lineRateHz = getStructFieldOrEmpty(pp, 'lineRateHz');
    dmdMeta.linesPerCycle = getStructFieldOrEmpty(pp, 'linesPerCycle');
    dmdMeta.linesPerFrame = getStructFieldOrEmpty(pp, 'linesPerFrame');
    dmdMeta.pixPerLine = getStructFieldOrEmpty(pp, 'pixPerLine');
end
if isempty(dmdMeta.linePeriod_s) && isfield(dmdMeta, 'lineRateHz') && ~isempty(dmdMeta.lineRateHz)
    dmdMeta.linePeriod_s = 1 ./ double(dmdMeta.lineRateHz);
end
end


function warnSlap2EpochMetadataChange(dmdIdx, epochIdx, canonicalMeta, epochMeta, canonicalLinesPerCycle, epochLinesPerCycle)
%WARNSLAP2EPOCHMETADATACHANGE Flag acquisition changes while preserving metadata.
if ~isequal(canonicalLinesPerCycle, epochLinesPerCycle)
    warning('summarize_LoCo:EpochLinesPerCycleChanged', ...
        'DMD%d epoch%d linesPerCycle changed from %s to %s.', ...
        dmdIdx, epochIdx, mat2str(canonicalLinesPerCycle), mat2str(epochLinesPerCycle));
end
for fieldName = {'lineRateHz', 'dmdPixelsPerRow', 'dmdPixelsPerColumn', 'samplesPerLine', 'channelsSave'}
    name = fieldName{1};
    a = getStructFieldOrEmpty(canonicalMeta, name);
    b = getStructFieldOrEmpty(epochMeta, name);
    if ~isempty(a) && ~isempty(b) && ~isequal(a, b)
        warning('summarize_LoCo:EpochMetadataChanged', ...
            'DMD%d epoch%d acquisition metadata field %s changed across epochs.', ...
            dmdIdx, epochIdx, name);
    end
end
end


function val = getStructFieldOrEmpty(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    val = s.(fieldName);
else
    val = [];
end
end


function value = scalarNumericField(s, fieldName)
value = nan;
if ~isstruct(s) || ~isfield(s, fieldName) || isempty(s.(fieldName))
    return
end
try
    arr = double(s.(fieldName));
    value = arr(1);
catch
    value = nan;
end
end


function [IMann, hasAnnotationImage] = extractAnnotationImage(roiEntry)
IMann = [];
hasAnnotationImage = false;
candidateFields = {'IM','im','image','annotationImage','annotationIM','refImage','referenceImage'};
for ix = 1:numel(candidateFields)
    fnm = candidateFields{ix};
    if isfield(roiEntry, fnm) && ~isempty(roiEntry.(fnm))
        IMann = collapseImageForRegistration(roiEntry.(fnm));
        hasAnnotationImage = ~isempty(IMann) && any(isfinite(IMann(:)));
        if hasAnnotationImage
            return
        end
    end
end
end

function IM = collapseImageForRegistration(IMraw)
IMraw = double(IMraw);
if isempty(IMraw)
    IM = [];
    return
end
if ndims(IMraw)<=2
    IM = IMraw;
else
    avgDims = 3:ndims(IMraw);
    IM = squeeze(mean(IMraw, avgDims, 'omitnan'));
end
end

function [shiftRC, ok] = estimateRegistrationShift(IMtarget, IMann)
shiftRC = [0 0];
ok = false;
if isempty(IMtarget) || isempty(IMann)
    return
end
targetSz = size(IMtarget,[1 2]);
annSz = size(IMann,[1 2]);
padSz = max(cat(1,targetSz,annSz),[],1);

targetPad = nan(padSz);
targetPad(1:targetSz(1),1:targetSz(2)) = IMtarget;
annPad = nan(padSz);
annPad(1:annSz(1),1:annSz(2)) = IMann;

fillval = 0;
targetPad(isnan(targetPad)) = fillval;
annPad(isnan(annPad)) = fillval;

clipShift = max(2, min([padSz(1)-100 padSz(2)-100]));

output1 = dftregistration_clipped( fft2(targetPad), fft2(annPad), 1, clipShift);
shiftRC = [-output1(3) -output1(4)];
corrCoeff = output1(1);

ok = corrCoeff>0.4;
end

function roiDataOut = shiftRoiData(roiDataIn, shiftRC, targetSz)
roiDataOut = roiDataIn;
rowShift = shiftRC(1);
colShift = shiftRC(2);
for rix = 1:numel(roiDataIn)
    roiEntry = roiDataIn{rix};
    if isfield(roiEntry, 'Center') && ~isempty(roiEntry.Center)
        roiEntry.Center = roiEntry.Center - [colShift rowShift];
    end
    if isfield(roiEntry, 'Position') && ~isempty(roiEntry.Position)
        roiEntry.Position = roiEntry.Position - [colShift rowShift];
    end
    if isfield(roiEntry, 'Vertices') && ~isempty(roiEntry.Vertices)
        roiEntry.Vertices = roiEntry.Vertices - [colShift rowShift];
    end
    if isfield(roiEntry, 'mask') && ~isempty(roiEntry.mask)
        roiEntry.mask = shiftMaskToTarget(roiEntry.mask, shiftRC, targetSz);
    end
    roiDataOut{rix} = roiEntry;
end
end

function maskOut = shiftMaskToTarget(maskIn, shiftRC, targetSz)
maskOut = false(targetSz);
if isempty(maskIn)
    return
end
[rr,cc] = find(maskIn);
rr2 = round(rr - shiftRC(1));
cc2 = round(cc - shiftRC(2));
valid = rr2>=1 & rr2<=targetSz(1) & cc2>=1 & cc2<=targetSz(2);
if any(valid)
    maskOut(sub2ind(targetSz, rr2(valid), cc2(valid))) = true;
end
end