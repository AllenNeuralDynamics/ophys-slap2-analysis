function summary = extractDendrites_new(dr_or_pathToTrialTable, paramsIn)
%EXTRACTDENDRITES_NEW Extract SLAP2 dendritic-voltage ROI traces.
%
%   summary = extractDendrites_new(dr_or_pathToTrialTable, paramsIn)
%
%   This pipeline reads a SLAP2 trial table, discovers DMD integration ROIs,
%   extracts voltage traces from the SLAP2 Trace backend, and writes large trace
%   arrays directly to HDF5 while saving lightweight session/ROI metadata to a
%   MATLAB summary file. It supports continuous CYCLE acquisitions and trial-file
%   acquisitions, with optional trial-sliced, continuous, or combined outputs.
%
%   Important conventions:
%       - ROI masks are stored in DMD/image pixel coordinates as
%         [dmdPixelsPerColumn x dmdPixelsPerRow x nROIs].
%       - When available, stored ROI mask images are preferred over masks
%         reconstructed from shapeData. This matches the summarize_LoCo.m access
%         pattern and avoids small coordinate offsets from shapeData conventions.
%       - Reference images are optional metadata/QC outputs; trace extraction only
%         requires masks and .dat/.meta files.
%
%   Typical usage:
%       params.outputMode = 'trial';      % 'trial', 'continuous', or 'both'
%       params.storageMode = 'h5';        % recommended
%       params.precision = 'single';
%       params.numWorkers = 4;
%       params.maxConcurrentROIs = 2;
%       summary = extractDendrites_new(sessionDir, params);
%
%   Outputs:
%       dendriticVoltageSummary-YYmmDD-HHMMSS.mat
%           Summary metadata, DMD ROI masks, reference images, trial/epoch
%           line ranges, ROI table, extraction status, and path to paired HDF5 traces.
%
%       dendriticVoltageTraces-YYmmDD-HHMMSS.h5
%           /traces/trial_####       [nTrialLines x nTotalROIs]
%           /traces/continuous/DMD#  [nDmdLines x nLocalROIs]
%
%   Epoch handling:
%       Multi-epoch sessions are detected from acquisition filename prefix changes
%       (for example acquisition_135515_DMD1 -> acquisition_140535_DMD1).
%       Epoch metadata are saved in summary.epochTable, summary.trialEpoch,
%       and summary.trialFilePrefix. For non-CYCLE multi-epoch acquisitions,
%       this function writes trial-sliced outputs; epoch-continuous HDF5 outputs
%       are intentionally deferred until downstream alignment can account for
%       behavior-continuous/imaging-discontinuous timebases.

import ScanImageTiffReader.ScanImageTiffReader

% -------------------------------------------------------------------------
% Resolve inputs, output locations, and run parameters
% -------------------------------------------------------------------------

if nargin < 1 || isempty(dr_or_pathToTrialTable)
    [trialTablefn, dr] = uigetfile({'*.mat', 'MAT-files (*.mat)'}, 'Select trialTable.mat');
    assert(~isnumeric(trialTablefn), 'No trialTable file selected.');
else
    [dr, trialTablefn] = resolveTrialTablePath(dr_or_pathToTrialTable);
end

if nargin < 2
    params = initializeParams();
else
    params = initializeParams(paramsIn);
end

timestamp = datestr(now, 'YYmmDD-HHMMSS');

if isempty(params.outputDir)
    if params.useOutputSubfolder
        if params.timestampOutputSubfolder
            params.outputDir = fullfile(dr, [params.outputSubfolderName '-' timestamp]);
        else
            params.outputDir = fullfile(dr, params.outputSubfolderName);
        end
    else
        params.outputDir = dr;
    end
end
if ~exist(params.outputDir, 'dir')
    mkdir(params.outputDir);
end

summaryPath = fullfile(params.outputDir, ['dendriticVoltageSummary-' timestamp '.mat']);
h5Path = fullfile(params.outputDir, ['dendriticVoltageTraces-' timestamp '.h5']);
if strcmpi(params.storageMode, 'h5') && exist(h5Path, 'file') && ~params.resume
    delete(h5Path);
end

% -------------------------------------------------------------------------
% Load trial table and initialize summary metadata
% -------------------------------------------------------------------------

fprintf('Loading trial table: %s\n', fullfile(dr, trialTablefn));
S = load(fullfile(dr, trialTablefn), 'trialTable');
trialTable = S.trialTable;
clear S

trialInfo = buildTrialInfo(trialTable);
nDMDs = size(trialInfo.filename, 1);
nTrials = size(trialInfo.filename, 2);
isContinuousAcq = any(contains(trialInfo.allFilenames, '-CYCLE-', 'IgnoreCase', true));

fprintf('Detected %d DMD(s), %d trial(s), continuous acquisition: %d\n', ...
    nDMDs, nTrials, isContinuousAcq);

if strcmpi(params.storageMode, 'memory')
    warning(['storageMode="memory" is provided only for compatibility and may ' ...
        'reintroduce memory pressure. Prefer storageMode="h5".']);
end

summary = struct();
summary.createdAt = char(datetime('now'));
summary.sourceTrialTable = fullfile(dr, trialTablefn);
summary.sessionDir = dr;
summary.params = params;
summary.outputMode = params.outputMode;
summary.storageMode = params.storageMode;
summary.outputDir = params.outputDir;
summary.outputH5 = h5Path;
summary.summaryPath = summaryPath;
summary.isContinuousAcquisition = isContinuousAcq;
summary.nDMDs = nDMDs;
summary.nTrials = nTrials;
summary.trialTable = trialInfo.minimal;
summary.trialLineRanges = trialInfo.lineRanges;
summary.trialEpoch = trialInfo.trialEpoch;
summary.trialFilePrefix = trialInfo.trialFilePrefix;
summary.epochTable = trialInfo.epochTable;
summary.nEpochs = trialInfo.nEpochs;
summary.multiEpochAcquisition = trialInfo.nEpochs > 1;
summary.outputPlan = planOutputs(params, isContinuousAcq, nTrials);
summary.masks = cell(1, nDMDs);
summary.maskImages = cell(1, nDMDs);
summary.refIM = cell(1, nDMDs);
summary.dmd = repmat(struct(), 1, nDMDs);

nAnalysisROIs = zeros(1, nDMDs);
roiRecords = struct('globalRoiIdx', {}, 'dmdIdx', {}, 'localRoiIdx', {}, ...
    'sourceRoiIdx', {}, 'nPixels', {}, 'isManual', {});

% -------------------------------------------------------------------------
% Discover DMD ROIs, reference images, and ROI bookkeeping
% -------------------------------------------------------------------------

fprintf('Discovering ROIs and reference images...\n');
for dmdIdx = 1:nDMDs
    firstValidTrial = findFirstValidTrial(trialInfo.filename, dmdIdx);
    firstDatFile = resolveDataFilePath(dr, trialInfo.filename{dmdIdx, firstValidTrial});

    fprintf('  DMD%d metadata: %s\n', dmdIdx, firstDatFile);
    hMDF = slap2.util.MultiDataFiles(firstDatFile);
    validateParsePlanCompatibility(hMDF, dmdIdx);

    dmdMeta = getDmdMetadata(hMDF);
    summary.dmd(dmdIdx).metadata = dmdMeta;
    summary.dmd(dmdIdx).firstDatFile = firstDatFile;
    summary.dmd(dmdIdx).totalNumLines = hMDF.totalNumLines;
    summary.dmd(dmdIdx).nCycles = hMDF.numCycles;
    summary.dmd(dmdIdx).linesPerCycle = hMDF.header.linesPerCycle;

    [masks, maskImage, sourceRoiIdx] = getIntegrationMasks(hMDF);
    [refIM, outlines] = loadReferenceImage(dr, dmdIdx, hMDF, masks, params);

    if params.manualROIs
        assert(~isempty(refIM), ['manualROIs=true requires a reference image. ' ...
            'Set saveRefImages/loadRefImages to true or provide a compatible reference TIFF.']);
        refFn = findReferenceTiff(dr, dmdIdx);
        hROIs = drawROIs(refIM, dr, refFn.name, outlines);
        waitfor(hROIs.hF);
        [masks, maskImage, sourceRoiIdx] = masksFromManualROIs(hROIs, hMDF, dmdIdx);
    end

    nAnalysisROIs(dmdIdx) = size(masks, 3);
    summary.masks{dmdIdx} = masks;
    summary.maskImages{dmdIdx} = maskImage;
    if params.saveRefImages
        summary.refIM{dmdIdx} = refIM;
    end

    if params.makePlots
        plotDmdPreview(dmdIdx, refIM, maskImage);
    end

    for localRoiIdx = 1:nAnalysisROIs(dmdIdx)
        globalRoiIdx = numel(roiRecords) + 1;
        roiRecords(globalRoiIdx).globalRoiIdx = globalRoiIdx;
        roiRecords(globalRoiIdx).dmdIdx = dmdIdx;
        roiRecords(globalRoiIdx).localRoiIdx = localRoiIdx;
        roiRecords(globalRoiIdx).sourceRoiIdx = sourceRoiIdx(localRoiIdx);
        roiRecords(globalRoiIdx).nPixels = nnz(masks(:, :, localRoiIdx));
        roiRecords(globalRoiIdx).isManual = logical(params.manualROIs);
    end

    delete(hMDF);
    clear hMDF
end

summary.nAnalysisROIs = nAnalysisROIs;
summary.nTotalROIs = sum(nAnalysisROIs);
summary.roiTable = struct2table(roiRecords);
summary.roiGlobalOffsets = [0, cumsum(nAnalysisROIs(1:end-1))];
summary.extractionStatus = initializeStatus(summary.roiTable);

fprintf('Detected %d imaging epoch(s) from trial-table filenames.', summary.nEpochs);
if summary.nEpochs > 1
    fprintf('  Epoch trial ranges: ');
    for epIdx = 1:height(summary.epochTable)
        fprintf('E%d=%d-%d ', summary.epochTable.epochIdx(epIdx), ...
            summary.epochTable.firstTrial(epIdx), summary.epochTable.lastTrial(epIdx));
    end
    fprintf('');
end
validateOutputPlan(summary.outputPlan, params, isContinuousAcq, nTrials);

% -------------------------------------------------------------------------
% Initialize output files, then extract traces
% -------------------------------------------------------------------------

if strcmpi(params.storageMode, 'h5')
    initializeH5Outputs(h5Path, summary, trialInfo, params);
end

save(summaryPath, 'summary', '-v7.3');
fprintf('Initialized summary: %s\n', summaryPath);
if strcmpi(params.storageMode, 'h5')
    fprintf('Initialized trace HDF5: %s\n', h5Path);
end

% Configure the parallel pool only after metadata/ROI discovery and HDF5
% initialization succeed. This avoids launching workers for runs that fail
% early because of path or metadata-parser problems.
params = configureParallelPool(params);
summary.params = params;
save(summaryPath, 'summary', '-v7.3');

if isContinuousAcq
    summary = extractContinuousAcquisition(summary, trialInfo, dr, params, summaryPath, h5Path);
else
    summary = extractTrialFiles(summary, trialInfo, dr, params, summaryPath, h5Path);
end

summary.completedAt = char(datetime('now'));
summary.complete = true;
save(summaryPath, 'summary', '-v7.3');
fprintf('Extraction complete. Summary saved to:\n  %s\n', summaryPath);
if strcmpi(params.storageMode, 'h5')
    fprintf('Trace data saved to:\n  %s\n', h5Path);
end
end

function outputPlan = planOutputs(params, isContinuousAcq, nTrials)
%PLANOUTPUTS Decide which HDF5 outputs are valid for this acquisition layout.
requestedTrial = wantsTrial(params.outputMode);
requestedContinuous = wantsContinuous(params.outputMode);

% Continuous output is safe for true -CYCLE acquisitions and for the edge case
% where a single non-CYCLE .dat file spans the whole acquisition. It is not yet
% safe for multi-trial/multi-epoch non-CYCLE sessions because behavior continues
% through imaging gaps and should be aligned downstream with explicit epochs.
canWriteContinuous = isContinuousAcq || nTrials == 1;

outputPlan = struct();
outputPlan.requestedMode = params.outputMode;
outputPlan.writeTrial = requestedTrial;
outputPlan.writeContinuous = requestedContinuous && canWriteContinuous;
outputPlan.canWriteContinuous = canWriteContinuous;
outputPlan.isContinuousAcquisition = isContinuousAcq;
outputPlan.nTrials = nTrials;
end

function validateOutputPlan(outputPlan, params, isContinuousAcq, nTrials)
%VALIDATEOUTPUTPLAN Fail early for unsupported output requests.
if wantsContinuous(params.outputMode) && ~outputPlan.writeContinuous
    if ~outputPlan.writeTrial
        error('extractDendrites:UnsupportedContinuousOutput', ...
            ['params.outputMode="continuous" is not supported for multi-trial ' ...
             'non-CYCLE acquisitions. Use outputMode="trial" to preserve trial ' ...
             'datasets plus summary.epochTable/trialEpoch, or outputMode="both" ' ...
             'if you want trial output while continuous writing is skipped.']);
    end
    warning(['Requested continuous output for a multi-trial non-CYCLE acquisition. ' ...
        'Continuous HDF5 datasets will not be written; trial datasets and epoch ' ...
        'metadata will be saved instead.']);
end

if isContinuousAcq && nTrials <= 0
    error('extractDendrites:InvalidTrialCount', ...
        'Continuous acquisition was detected, but no trials were found in trialInfo.');
end
end

% -------------------------------------------------------------------------
% Main extraction branches
% -------------------------------------------------------------------------

function summary = extractContinuousAcquisition(summary, trialInfo, dr, params, summaryPath, h5Path)
nDMDs = summary.nDMDs;

for dmdIdx = 1:nDMDs
    firstValidTrial = findFirstValidTrial(trialInfo.filename, dmdIdx);
    firstDatFile = resolveDataFilePath(dr, trialInfo.filename{dmdIdx, firstValidTrial});
    fprintf('\nExtracting continuous traces for DMD%d from %s\n', dmdIdx, firstDatFile);

    hMDF = slap2.util.MultiDataFiles(firstDatFile);
    localRois = 1:summary.nAnalysisROIs(dmdIdx);
    batches = makeBatches(localRois, params.maxConcurrentROIs);

    for batchIdx = 1:numel(batches)
        batch = batches{batchIdx};
        fprintf('  DMD%d batch %d/%d: ROI(s) %s\n', ...
            dmdIdx, batchIdx, numel(batches), mat2str(batch));

        [batchTraces, batchWeights, batchErrors] = extractRoiBatch( ...
            hMDF, summary.masks{dmdIdx}, batch, params);

        for j = 1:numel(batch)
            localRoiIdx = batch(j);
            globalRoiIdx = summary.roiGlobalOffsets(dmdIdx) + localRoiIdx;
            statusIdx = globalRoiIdx;
            summary.extractionStatus(statusIdx).startedAt = char(datetime('now'));

            if ~isempty(batchErrors{j})
                summary.extractionStatus(statusIdx).status = 'failed';
                summary.extractionStatus(statusIdx).errorMessage = batchErrors{j};
                if params.stopOnError
                    error('extractDendrites:ROIExtractionFailed', ...
                        'DMD%d ROI%d failed: %s', dmdIdx, localRoiIdx, batchErrors{j});
                end
                continue
            end

            trace = castTrace(batchTraces{j}, params.precision);
            summary.extractionStatus(statusIdx).nSamplesExtracted = numel(trace);
            summary.extractionStatus(statusIdx).weightClass = class(batchWeights{j});
            summary.extractionStatus(statusIdx).weightSize = sizeToString(size(batchWeights{j}));

            if summary.outputPlan.writeContinuous
                writeContinuousTrace(summary, h5Path, trace, dmdIdx, localRoiIdx, params);
            end
            if summary.outputPlan.writeTrial
                writeTrialSlices(summary, trialInfo, h5Path, trace, dmdIdx, globalRoiIdx, params);
            end

            summary.extractionStatus(statusIdx).status = 'complete';
            summary.extractionStatus(statusIdx).finishedAt = char(datetime('now'));
            summary.extractionStatus(statusIdx).errorMessage = '';

            clear trace
        end

        clear batchTraces batchWeights batchErrors
        if params.saveSummaryAfterEachBatch
            save(summaryPath, 'summary', '-v7.3');
        end
    end

    delete(hMDF);
    clear hMDF
end
end

function summary = extractTrialFiles(summary, trialInfo, dr, params, summaryPath, h5Path)
% Non-CYCLE acquisitions are either true trial-file outputs or a single file
% accidentally collected without a -CYCLE suffix. For the latter, continuous
% output is safe because the one .dat file spans the whole imaging epoch. For
% multi-file/multi-epoch acquisitions, keep trial outputs and preserve epoch
% metadata rather than stitching discontinuous imaging periods into one trace.
if wantsContinuous(params.outputMode) && ~summary.outputPlan.writeContinuous
    warning(['Continuous output for multi-trial non-CYCLE acquisitions is not implemented. ' ...
        'Writing trial-sliced output only; use summary.epochTable/trialEpoch downstream.']);
end

nDMDs = summary.nDMDs;
nTrials = summary.nTrials;

for trialIdx = 1:nTrials
    fprintf('\nExtracting trial %d/%d\n', trialIdx, nTrials);
    for dmdIdx = 1:nDMDs
        if isempty(trialInfo.filename{dmdIdx, trialIdx})
            continue
        end
        datFile = resolveDataFilePath(dr, trialInfo.filename{dmdIdx, trialIdx});
        fprintf('  DMD%d trial %d file: %s\n', dmdIdx, trialIdx, datFile);

        hMDF = slap2.util.MultiDataFiles(datFile);
        localRois = 1:summary.nAnalysisROIs(dmdIdx);
        batches = makeBatches(localRois, params.maxConcurrentROIs);

        for batchIdx = 1:numel(batches)
            batch = batches{batchIdx};
            fprintf('    batch %d/%d: ROI(s) %s\n', ...
                batchIdx, numel(batches), mat2str(batch));

            [batchTraces, batchWeights, batchErrors] = extractRoiBatch( ...
                hMDF, summary.masks{dmdIdx}, batch, params);

            for j = 1:numel(batch)
                localRoiIdx = batch(j);
                globalRoiIdx = summary.roiGlobalOffsets(dmdIdx) + localRoiIdx;
                statusIdx = globalRoiIdx;
                summary.extractionStatus(statusIdx).startedAt = char(datetime('now'));

                if ~isempty(batchErrors{j})
                    summary.extractionStatus(statusIdx).status = 'failed';
                    summary.extractionStatus(statusIdx).errorMessage = batchErrors{j};
                    if params.stopOnError
                        error('extractDendrites:ROIExtractionFailed', ...
                            'Trial %d DMD%d ROI%d failed: %s', ...
                            trialIdx, dmdIdx, localRoiIdx, batchErrors{j});
                    end
                    continue
                end

                trace = castTrace(batchTraces{j}, params.precision);
                summary.extractionStatus(statusIdx).nSamplesExtracted = ...
                    summary.extractionStatus(statusIdx).nSamplesExtracted + numel(trace);
                summary.extractionStatus(statusIdx).weightClass = class(batchWeights{j});
                summary.extractionStatus(statusIdx).weightSize = sizeToString(size(batchWeights{j}));

                if summary.outputPlan.writeContinuous
                    writeContinuousTrace(summary, h5Path, trace, dmdIdx, localRoiIdx, params);
                end
                if summary.outputPlan.writeTrial
                    writeOneTrialSlice(summary, trialInfo, h5Path, trace, ...
                        dmdIdx, trialIdx, globalRoiIdx, params);
                end

                summary.extractionStatus(statusIdx).status = 'in_progress';
                summary.extractionStatus(statusIdx).finishedAt = char(datetime('now'));
                summary.extractionStatus(statusIdx).errorMessage = '';
                clear trace
            end

            clear batchTraces batchWeights batchErrors
            if params.saveSummaryAfterEachBatch
                save(summaryPath, 'summary', '-v7.3');
            end
        end

        delete(hMDF);
        clear hMDF
    end
end

for idx = 1:numel(summary.extractionStatus)
    if strcmp(summary.extractionStatus(idx).status, 'in_progress')
        summary.extractionStatus(idx).status = 'complete';
    end
end
end

% -------------------------------------------------------------------------
% ROI extraction
% -------------------------------------------------------------------------

function [batchTraces, batchWeights, batchErrors] = extractRoiBatch(hMDF, masks, batch, params)
batchTraces = cell(1, numel(batch));
batchWeights = cell(1, numel(batch));
batchErrors = cell(1, numel(batch));

if params.useParallel
    futures = cell(1, numel(batch));
    hTraces = cell(1, numel(batch));

    for j = 1:numel(batch)
        try
            localRoiIdx = batch(j);
            hTrace = slap2.util.datafile.trace.Trace(hMDF, params.zIdx, params.chIdx);
            pixelMask = masks(:, :, localRoiIdx);
            hTrace.setPixelIdxs(pixelMask, pixelMask);
            futures{j} = hTrace.processAsync(params.windowWidth_lines, ...
                params.expectedWindowWidth_lines);
            hTraces{j} = hTrace;
        catch ME
            batchErrors{j} = getReport(ME, 'extended', 'hyperlinks', 'off');
        end
    end

    for j = 1:numel(batch)
        if ~isempty(batchErrors{j}) || isempty(futures{j})
            continue
        end
        try
            % Some installed SLAP2 Trace.processAsync variants return only
            % trace, while newer variants may also return weight. Fetch only
            % the first output for compatibility; weight is optional metadata.
            trace = fetchOutputs(futures{j});
            batchTraces{j} = trace;
            batchWeights{j} = [];
        catch ME
            batchErrors{j} = getReport(ME, 'extended', 'hyperlinks', 'off');
        end
    end

    clear futures hTraces
else
    for j = 1:numel(batch)
        try
            localRoiIdx = batch(j);
            hTrace = slap2.util.datafile.trace.Trace(hMDF, params.zIdx, params.chIdx);
            pixelMask = masks(:, :, localRoiIdx);
            hTrace.setPixelIdxs(pixelMask, pixelMask);
            % Some installed SLAP2 Trace.process variants return only trace,
            % while newer variants may also return weight. Request only one
            % output so this wrapper works with both APIs.
            trace = hTrace.process(params.windowWidth_lines, ...
                params.expectedWindowWidth_lines);
            batchTraces{j} = trace;
            batchWeights{j} = [];
            clear hTrace
        catch ME
            batchErrors{j} = getReport(ME, 'extended', 'hyperlinks', 'off');
        end
    end
end
end

% -------------------------------------------------------------------------
% HDF5 initialization/writes
% -------------------------------------------------------------------------

function initializeH5Outputs(h5Path, summary, trialInfo, params)
if summary.outputPlan.writeTrial
    for trialIdx = 1:summary.nTrials
        dset = trialDatasetName(trialIdx);

        % Size each trial dataset by the union of valid DMD line ranges.
        % This preserves acquisition-time alignment when DMD1 and DMD2 have
        % slightly different first/last line values for the same trial.
        nRows = trialInfo.trialGlobalNLines(trialIdx);
        nCols = summary.nTotalROIs;
        createH5DatasetIfNeeded(h5Path, dset, [nRows, nCols], params);
        h5writeatt(h5Path, dset, 'description', ...
            ['Time-aligned trial-sliced traces. Columns are global ROI indices ' ...
             'from summary.roiTable. Rows are relative to trialFirstLineGlobal.']);
        h5writeatt(h5Path, dset, 'trialIdx', trialIdx);
        h5writeatt(h5Path, dset, 'epochIdx', trialInfo.trialEpoch(trialIdx));
        h5writeatt(h5Path, dset, 'acquisitionPrefix', trialInfo.trialFilePrefix{trialIdx});
        h5writeatt(h5Path, dset, 'trialFirstLineGlobal', trialInfo.trialFirstLineGlobal(trialIdx));
        h5writeatt(h5Path, dset, 'trialLastLineGlobal', trialInfo.trialLastLineGlobal(trialIdx));
        h5writeatt(h5Path, dset, 'firstLineByDmd', trialInfo.firstLine(:, trialIdx));
        h5writeatt(h5Path, dset, 'lastLineByDmd', trialInfo.lastLine(:, trialIdx));
        h5writeatt(h5Path, dset, 'firstLineRoundedByDmd', trialInfo.firstLineRounded(:, trialIdx));
        h5writeatt(h5Path, dset, 'lastLineRoundedByDmd', trialInfo.lastLineRounded(:, trialIdx));
    end
end

if summary.outputPlan.writeContinuous
    for dmdIdx = 1:summary.nDMDs
        dset = continuousDatasetName(dmdIdx);
        nRows = summary.dmd(dmdIdx).totalNumLines;
        nCols = summary.nAnalysisROIs(dmdIdx);
        createH5DatasetIfNeeded(h5Path, dset, [nRows, nCols], params);
        h5writeatt(h5Path, dset, 'description', ...
            'Continuous traces for one DMD. Columns are local ROI indices for this DMD.');
        h5writeatt(h5Path, dset, 'dmdIdx', dmdIdx);
        h5writeatt(h5Path, dset, 'globalRoiIdx', ...
            summary.roiGlobalOffsets(dmdIdx) + (1:summary.nAnalysisROIs(dmdIdx)));
    end
end

h5writeatt(h5Path, '/', 'createdAt', char(datetime('now')));
h5writeatt(h5Path, '/', 'sourceTrialTable', summary.sourceTrialTable);
h5writeatt(h5Path, '/', 'outputMode', params.outputMode);
h5writeatt(h5Path, '/', 'writeTrial', summary.outputPlan.writeTrial);
h5writeatt(h5Path, '/', 'writeContinuous', summary.outputPlan.writeContinuous);
h5writeatt(h5Path, '/', 'nEpochs', summary.nEpochs);
h5writeatt(h5Path, '/', 'precision', params.precision);
end

function createH5DatasetIfNeeded(h5Path, dset, dims, params)
if h5DatasetExists(h5Path, dset)
    return
end

chunkLines = min(dims(1), params.h5ChunkLines);
chunkSize = [max(chunkLines, 1), 1];

args = {h5Path, dset, dims, 'Datatype', params.precision, ...
    'ChunkSize', chunkSize, 'FillValue', cast(NaN, params.precision)};
if params.h5Deflate > 0
    args = [args, {'Deflate', params.h5Deflate}]; %#ok<AGROW>
end
h5create(args{:});
end

function tf = h5DatasetExists(h5Path, dset)
tf = false;
if ~exist(h5Path, 'file')
    return
end
try
    h5info(h5Path, dset);
    tf = true;
catch
    tf = false;
end
end

function writeContinuousTrace(summary, h5Path, trace, dmdIdx, localRoiIdx, params)
if strcmpi(params.storageMode, 'h5')
    dset = continuousDatasetName(dmdIdx);
    nWrite = min(numel(trace), summary.dmd(dmdIdx).totalNumLines);
    h5write(h5Path, dset, trace(1:nWrite), [1, localRoiIdx], [nWrite, 1]);
else
    error('extractDendrites:MemoryStorageUnsupported', ...
        'Memory storage is not supported for continuous traces in this refactor. Use storageMode="h5".');
end
end

function writeTrialSlices(summary, trialInfo, h5Path, trace, dmdIdx, globalRoiIdx, params)
for trialIdx = 1:summary.nTrials
    writeOneTrialSlice(summary, trialInfo, h5Path, trace, ...
        dmdIdx, trialIdx, globalRoiIdx, params);
end
end

function writeOneTrialSlice(summary, trialInfo, h5Path, trace, dmdIdx, trialIdx, globalRoiIdx, params)
firstLine = trialInfo.firstLine(dmdIdx, trialIdx);
lastLine = trialInfo.lastLine(dmdIdx, trialIdx);
firstLineRounded = trialInfo.firstLineRounded(dmdIdx, trialIdx);
lastLineRounded = trialInfo.lastLineRounded(dmdIdx, trialIdx);
if isnan(firstLineRounded) || isnan(lastLineRounded) || ...
        firstLineRounded <= 0 || lastLineRounded < firstLineRounded
    return
end

startIdx = max(1, firstLineRounded);
stopIdx = min(numel(trace), lastLineRounded);
if stopIdx < startIdx
    return
end

seg = trace(startIdx:stopIdx);
nWrite = numel(seg);

if strcmpi(params.storageMode, 'h5')
    dset = trialDatasetName(trialIdx);

    % Preserve cross-DMD acquisition-time alignment within the trial dataset.
    % Each trial dataset starts at the earliest valid line across all DMDs.
    % This DMD's slice is written at an offset relative to that global start.
    trialFirstLine = trialInfo.trialFirstLineGlobal(trialIdx);
    if isnan(trialFirstLine) || trialFirstLine < 1
        return
    end
    rowStart = firstLineRounded - trialFirstLine + 1;

    info = h5info(h5Path, dset);
    nRows = info.Dataspace.Size(1);
    if rowStart < 1 || (rowStart + nWrite - 1) > nRows
        error('extractDendrites:H5TrialWriteOutOfBounds', ...
            ['Cannot write DMD%d ROI%d trial %d into %s: rowStart=%d, ' ...
             'nWrite=%d, datasetRows=%d, firstLine=%g, lastLine=%g, ' ...
             'trialFirstLineGlobal=%g, trialLastLineGlobal=%g, ' ...
             'firstLineRounded=%d, lastLineRounded=%d.'], ...
            dmdIdx, globalRoiIdx, trialIdx, dset, rowStart, nWrite, nRows, ...
            firstLine, lastLine, trialInfo.trialFirstLineGlobal(trialIdx), ...
            trialInfo.trialLastLineGlobal(trialIdx), firstLineRounded, lastLineRounded);
    end

    h5write(h5Path, dset, seg(:), [rowStart, globalRoiIdx], [nWrite, 1]);
else
    error('extractDendrites:MemoryStorageUnsupported', ...
        'Memory storage is not implemented for trial-sliced writes. Use storageMode="h5".');
end
end

function dset = trialDatasetName(trialIdx)
dset = sprintf('/traces/trial_%04d', trialIdx);
end

function dset = continuousDatasetName(dmdIdx)
dset = sprintf('/traces/continuous/DMD%d', dmdIdx);
end

% -------------------------------------------------------------------------
% ROI and reference image helpers
% -------------------------------------------------------------------------

function [masks, maskImage, sourceRoiIdx] = getIntegrationMasks(hMDF)
%GETINTEGRATIONMASKS Return DMD integration ROI masks in image coordinates.
%
% Prefer a stored ROI mask when present, matching summarize_LoCo.m and other
% ophys-slap2-analysis preprocessing code. Older SLAP2 metadata sometimes only
% stores shapeData; in that case reconstruct the mask using the legacy
% summarize_Voltage/extractDendrites convention.
meta = hMDF.metaData;
imagingRois = meta.AcquisitionContainer.ROIs.rois;
if ~iscell(imagingRois)
    imagingRois = num2cell(imagingRois);
end

isIntegration = false(1, numel(imagingRois));
for idx = 1:numel(imagingRois)
    roi = imagingRois{idx};
    mode = getRoiField(roi, 'imagingMode', '');
    isIntegration(idx) = strcmpi(char(mode), 'Integrate');
end

sourceRoiIdx = find(isIntegration);
integrationRois = imagingRois(isIntegration);
nRois = numel(integrationRois);

nRows = double(meta.dmdPixelsPerColumn);
nCols = double(meta.dmdPixelsPerRow);
masks = false(nRows, nCols, nRois);
maskImage = -1 .* ones(nRows, nCols, 'single');

for roiIdx = 1:nRois
    roi = integrationRois{roiIdx};

    storedMask = getRoiField(roi, 'mask', []);
    if ~isempty(storedMask)
        tmp = normalizeRoiMask(storedMask, nRows, nCols, roiIdx);
    else
        shape = getRoiField(roi, 'shapeData', []);
        if isempty(shape)
            error('extractDendrites:MissingRoiMask', ...
                'Integration ROI %d has neither mask nor shapeData.', sourceRoiIdx(roiIdx));
        end
        tmp = maskFromShapeData(shape, nRows, nCols, sourceRoiIdx(roiIdx));
    end

    masks(:, :, roiIdx) = tmp;
    maskImage(tmp) = roiIdx;
end
end

function tmp = normalizeRoiMask(mask, nRows, nCols, roiIdx)
%NORMALIZEROIMASK Coerce stored ROI mask to expected DMD image size.
tmp = logical(mask);
if ndims(tmp) > 2 && size(tmp, 3) == 1
    tmp = tmp(:, :, 1);
end
if isequal(size(tmp), [nRows, nCols])
    return
end
if isequal(size(tmp), [nCols, nRows])
    tmp = tmp.';
    return
end
error('extractDendrites:UnexpectedRoiMaskSize', ...
    'ROI %d stored mask has size %s, expected [%d %d].', ...
    roiIdx, mat2str(size(tmp)), nRows, nCols);
end

function tmp = maskFromShapeData(shape, nRows, nCols, sourceRoiIdx)
%MASKFROMSHAPEDATA Legacy fallback for metadata without stored ROI masks.
shape = double(shape);
if size(shape, 2) < 2
    error('extractDendrites:InvalidShapeData', ...
        'ROI %d shapeData must have at least two columns.', sourceRoiIdx);
end
shape = round(shape(:, 1:2));
valid = shape(:, 1) >= 1 & shape(:, 1) <= nRows & ...
        shape(:, 2) >= 1 & shape(:, 2) <= nCols;
if ~all(valid)
    warning('extractDendrites:ShapeDataOutOfBounds', ...
        'ROI %d shapeData contains %d out-of-bounds pixels; ignoring them.', ...
        sourceRoiIdx, nnz(~valid));
end
shape = shape(valid, :);
tmp = false(nRows, nCols);
if ~isempty(shape)
    tmp(sub2ind(size(tmp), shape(:, 1), shape(:, 2))) = true;
end
end

function val = getRoiField(roi, fieldName, defaultVal)
%GETROIFIELD Read a field/property from struct or object ROI metadata.
if nargin < 3
    defaultVal = [];
end
val = defaultVal;
if isstruct(roi) && isfield(roi, fieldName)
    val = roi.(fieldName);
elseif isobject(roi) && isprop(roi, fieldName)
    val = roi.(fieldName);
end
end

function [masks, maskImage, sourceRoiIdx] = masksFromManualROIs(hROIs, hMDF, dmdIdx) %#ok<INUSD>
nRows = double(hMDF.metaData.dmdPixelsPerColumn);
nCols = double(hMDF.metaData.dmdPixelsPerRow);
nRois = numel(hROIs.roiData);
masks = false(nRows, nCols, nRois);
maskImage = -1 .* ones(nRows, nCols, 'single');
sourceRoiIdx = nan(1, nRois);
for roiIdx = 1:nRois
    masks(:, :, roiIdx) = logical(hROIs.roiData{roiIdx}.mask);
    maskImage(masks(:, :, roiIdx)) = roiIdx;
end
end

function [refIM, outlines] = loadReferenceImage(dr, dmdIdx, hMDF, masks, params)
refIM = [];
outlines = {};

for roiIdx = 1:size(masks, 3)
    outlines = cat(1, outlines, bwboundaries(masks(:, :, roiIdx))); %#ok<AGROW>
end

if ~params.saveRefImages && ~params.manualROIs && ~params.makePlots
    return
end

try
    refFn = findReferenceTiff(dr, dmdIdx);
catch ME
    warning('Could not find reference TIFF for DMD%d: %s', dmdIdx, ME.message);
    return
end

try
    A = ScanImageTiffReader.ScanImageTiffReader(fullfile(refFn.folder, refFn.name));
    IDs = A.descriptions;
    z = zeros(numel(IDs), 1);
    ch = zeros(numel(IDs), 1);
    for imIdx = 1:numel(IDs)
        js = jsondecode(IDs{imIdx});
        z(imIdx) = double(js.z);
        ch(imIdx) = double(js.channel);
    end
    nChan = numel(unique(ch));
    Zs = unique(z);
    stack = A.data();
    stack = reshape(stack, size(stack, 1), size(stack, 2), nChan, []);

    if isfield(hMDF.metaData.AcquisitionContainer, 'ParsePlan') && ...
            isfield(hMDF.metaData.AcquisitionContainer.ParsePlan, 'zs')
        metaZ = double(hMDF.metaData.AcquisitionContainer.ParsePlan.zs);
    elseif ~isempty(hMDF.fastZs)
        metaZ = double(hMDF.fastZs(1));
    else
        metaZ = Zs(1);
    end

    [~, bestZix] = min(abs(Zs - metaZ));
    bestZix = bestZix(1);
    refIM = permute(stack(:, :, :, bestZix), [2, 1, 3]);
    clear stack A
catch ME
    warning('Failed to load reference image for DMD%d: %s', dmdIdx, ME.message);
    refIM = [];
end
end

function refFn = findReferenceTiff(dr, dmdIdx)
refFn = dir(fullfile(dr, '**', sprintf('*DMD%d*REFERENCE*.tif', dmdIdx)));
assert(numel(refFn) == 1, ...
    'Expected exactly one DMD%d reference TIFF; found %d.', dmdIdx, numel(refFn));
end

function plotDmdPreview(dmdIdx, refIM, maskImage)
%PLOTDMDPREVIEW Quick QC view of reference image and extracted ROI labels.
if ~isempty(refIM)
    figure('Name', sprintf('Reference + masks for DMD%d', dmdIdx));
    imshow(refIM(:, :, 1), []);
    hold on;
    B = bwboundaries(maskImage > 0);
    for idx = 1:numel(B)
        plot(B{idx}(:, 2), B{idx}(:, 1), 'c', 'LineWidth', 1.0);
    end
    title(sprintf('DMD%d reference image with ROI outlines', dmdIdx));
else
    figure('Name', sprintf('Masks for DMD%d', dmdIdx));
    imshow(maskImage, []);
    colormap('jet');
    title(sprintf('DMD%d ROI mask labels', dmdIdx));
end
drawnow;
end

% -------------------------------------------------------------------------
% Metadata and trial table helpers
% -------------------------------------------------------------------------

function [dr, trialTablefn] = resolveTrialTablePath(dr_or_pathToTrialTable)
if exist(dr_or_pathToTrialTable, 'dir')
    dr = char(dr_or_pathToTrialTable);
    trialTablefn = 'trialTable.mat';
else
    [dr, name, ext] = fileparts(char(dr_or_pathToTrialTable));
    if isempty(dr)
        dr = pwd;
    end
    trialTablefn = [name, ext];
end
end

function trialInfo = buildTrialInfo(trialTable)
filename = trialTable.filename;
if isstring(filename)
    filename = cellstr(filename);
end
if ~iscell(filename)
    error('extractDendrites:InvalidTrialTable', ...
        'trialTable.filename must be a cell array or string array.');
end

firstLine = double(trialTable.firstLine);
lastLine = double(trialTable.lastLine);

% The trial table line values can contain fractional/near-integer values, while
% Trace output is indexed with integer MATLAB subscripts. Use one consistent
% integer line-coordinate system everywhere HDF5 trial datasets are sized and
% written. This avoids off-by-one failures where raw ranges create a dataset
% that is one row shorter than the rounded write segment.
firstLineRounded = round(firstLine);
lastLineRounded = round(lastLine);

lineRanges = struct();
lineRanges.firstLine = firstLine;
lineRanges.lastLine = lastLine;
lineRanges.firstLineRounded = firstLineRounded;
lineRanges.lastLineRounded = lastLineRounded;
lineRanges.nLines = lastLineRounded - firstLineRounded + 1;
lineRanges.nLines(lineRanges.nLines < 0) = NaN;

valid = ~cellfun(@isempty, filename);
allFilenames = filename(valid);

% Per-trial union of DMD line ranges. These fields define the row coordinate
% system for /traces/trial_XXXX datasets. Rows are aligned to the earliest
% valid firstLine across DMDs for each trial, so DMDs with later starts are
% written with a positive row offset instead of being forced to row 1.
validLines = ~isnan(firstLineRounded) & ~isnan(lastLineRounded) & ...
    firstLineRounded > 0 & lastLineRounded >= firstLineRounded;
firstForMin = firstLineRounded;
firstForMin(~validLines) = NaN;
lastForMax = lastLineRounded;
lastForMax(~validLines) = NaN;

trialFirstLineGlobal = min(firstForMin, [], 1, 'omitnan');
trialLastLineGlobal = max(lastForMax, [], 1, 'omitnan');
invalidTrial = isnan(trialFirstLineGlobal) | isnan(trialLastLineGlobal) | ...
    trialLastLineGlobal < trialFirstLineGlobal;
trialFirstLineGlobal(invalidTrial) = 1;
trialLastLineGlobal(invalidTrial) = 1;
trialGlobalNLines = trialLastLineGlobal - trialFirstLineGlobal + 1;
trialGlobalNLines(trialGlobalNLines < 1 | isnan(trialGlobalNLines)) = 1;
trialGlobalNLines = round(trialGlobalNLines);

% Retain the old maximum per-DMD trial length for reference/backward
% compatibility, but do not use it to size HDF5 trial datasets.
maxTrialLines = max(lineRanges.nLines, [], 1, 'omitnan');
maxTrialLines(isnan(maxTrialLines) | maxTrialLines < 1) = 1;
maxTrialLines = round(maxTrialLines);

lineRanges.trialFirstLineGlobal = trialFirstLineGlobal;
lineRanges.trialLastLineGlobal = trialLastLineGlobal;
lineRanges.trialGlobalNLines = trialGlobalNLines;

[trialEpoch, trialFilePrefix, epochTable] = inferAcquisitionEpochs(filename, trialTable);

minimal = struct();
minimal.filename = filename;
minimal.firstLine = firstLine;
minimal.lastLine = lastLine;
minimal.firstLineRounded = firstLineRounded;
minimal.lastLineRounded = lastLineRounded;
minimal.trialEpoch = trialEpoch;
minimal.trialFilePrefix = trialFilePrefix;
% Preserve scalar/vector trial metadata but intentionally omit large image stacks.
trialFields = fieldnames(trialTable);
for tfIdx = 1:numel(trialFields)
    tfName = trialFields{tfIdx};
    if any(strcmp(tfName, {'refStack', 'filename', 'firstLine', 'lastLine'}))
        continue
    end
    minimal.(tfName) = trialTable.(tfName);
end

trialInfo = struct();
trialInfo.filename = filename;
trialInfo.firstLine = firstLine;
trialInfo.lastLine = lastLine;
trialInfo.firstLineRounded = firstLineRounded;
trialInfo.lastLineRounded = lastLineRounded;
trialInfo.lineRanges = lineRanges;
trialInfo.maxTrialLines = maxTrialLines;
trialInfo.trialFirstLineGlobal = trialFirstLineGlobal;
trialInfo.trialLastLineGlobal = trialLastLineGlobal;
trialInfo.trialGlobalNLines = trialGlobalNLines;
trialInfo.allFilenames = allFilenames(:);
trialInfo.trialEpoch = trialEpoch;
trialInfo.trialFilePrefix = trialFilePrefix;
trialInfo.epochTable = epochTable;
trialInfo.nEpochs = height(epochTable);
trialInfo.minimal = minimal;
end

function [trialEpoch, trialFilePrefix, epochTable] = inferAcquisitionEpochs(filename, trialTable)
%INFERACQUISITIONEPOCHS Infer imaging epochs from trial-table filenames.
%
% buildTrialTableSLAP2 sometimes leaves trialTable.epoch uninformative even when
% acquisition was paused/restarted. File prefixes are the most reliable local
% signal: acquisition_135515_DMD1-CYCLE-... and acquisition_135515_DMD2-CYCLE-...
% belong to the same epoch, while acquisition_140535_* starts a new epoch.
[~, nTrials] = size(filename);
trialFilePrefix = repmat({''}, 1, nTrials);
for trialIdx = 1:nTrials
    fns = filename(:, trialIdx);
    fns = fns(~cellfun(@isempty, fns));
    if isempty(fns)
        trialFilePrefix{trialIdx} = '';
    else
        trialFilePrefix{trialIdx} = acquisitionPrefixFromFilename(fns{1});
    end
end

% Use explicit epoch labels only if they are informative. Otherwise derive
% consecutive epochs from filename-prefix changes.
explicitEpoch = [];
try
    if istable(trialTable) && any(strcmp(trialTable.Properties.VariableNames, 'epoch'))
        explicitEpoch = double(trialTable.epoch(:))';
    elseif isstruct(trialTable) && isfield(trialTable, 'epoch')
        explicitEpoch = double(trialTable.epoch(:))';
    end
catch
    explicitEpoch = [];
end

if numel(explicitEpoch) == nTrials && numel(unique(explicitEpoch(~isnan(explicitEpoch)))) > 1
    trialEpoch = relabelEpochVectorStable(explicitEpoch);
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

function trialEpoch = relabelEpochVectorStable(epochVals)
%RELABELEPOCHVECTORSTABLE Convert arbitrary epoch labels to 1..N in stable order.
trialEpoch = nan(size(epochVals));
seen = [];
for idx = 1:numel(epochVals)
    val = epochVals(idx);
    if isnan(val)
        if isempty(seen)
            seen = val;
            trialEpoch(idx) = 1;
        else
            trialEpoch(idx) = numel(seen);
        end
        continue
    end
    match = find(seen == val, 1, 'first');
    if isempty(match)
        seen(end + 1) = val; %#ok<AGROW>
        match = numel(seen);
    end
    trialEpoch(idx) = match;
end
trialEpoch(isnan(trialEpoch)) = 1;
end

function prefix = acquisitionPrefixFromFilename(fn)
%ACQUISITIONPREFIXFROMFILENAME Strip DMD and CYCLE suffixes to get epoch key.
[~, name, ~] = fileparts(char(fn));
prefix = regexprep(name, '_DMD\d+.*$', '');
prefix = regexprep(prefix, '-DMD\d+.*$', '');
prefix = regexprep(prefix, '_CYCLE[_-]?\d+.*$', '');
prefix = regexprep(prefix, '-CYCLE[_-]?\d+.*$', '');
if isempty(prefix)
    prefix = name;
end
end

function firstValidTrial = findFirstValidTrial(filename, dmdIdx)
row = filename(dmdIdx, :);
firstValidTrial = find(~cellfun(@isempty, row), 1, 'first');
assert(~isempty(firstValidTrial), 'No valid filenames found for DMD%d.', dmdIdx);
end

function datPath = resolveDataFilePath(dr, filename)
filename = char(filename);
if exist(filename, 'file')
    datPath = filename;
elseif exist(fullfile(dr, filename), 'file')
    datPath = fullfile(dr, filename);
else
    error('extractDendrites:MissingDataFile', ...
        'Could not find data file: %s', filename);
end
end

function validateParsePlanCompatibility(hMDF, dmdIdx)
required = {'lineSuperPixelIDs', 'lineFastZIdxs', 'zPixelReplacementMaps'};
for idx = 1:numel(required)
    fieldName = required{idx};
    if ~isprop(hMDF, fieldName)
        error('extractDendrites:ParsePlanMissingField', ...
            ['DMD%d MultiDataFiles is missing %s. This usually means ' ...
             'slap2.util.DataFile/loadParsePlan is incompatible with the .meta file.'], ...
            dmdIdx, fieldName);
    end
end
if isempty(hMDF.zPixelReplacementMaps) || isempty(hMDF.lineSuperPixelIDs)
    error('extractDendrites:EmptyParsePlan', ...
        ['DMD%d has an empty parse plan after DataFile loading. Replace ' ...
         'loadParsePlan.m with the companion patched version.'], dmdIdx);
end
end

function dmdMeta = getDmdMetadata(hMDF)
dmdMeta = struct();
dmdMeta.dmdPixelsPerRow = getFieldOrEmpty(hMDF.metaData, 'dmdPixelsPerRow');
dmdMeta.dmdPixelsPerColumn = getFieldOrEmpty(hMDF.metaData, 'dmdPixelsPerColumn');
dmdMeta.linePeriod_s = getFieldOrEmpty(hMDF.metaData, 'linePeriod_s');
dmdMeta.samplesPerLine = getFieldOrEmpty(hMDF.metaData, 'samplesPerLine');
dmdMeta.channelsSave = getFieldOrEmpty(hMDF.metaData, 'channelsSave');
dmdMeta.acqDuration_s = getFieldOrEmpty(hMDF.metaData, 'acqDuration_s');
if isfield(hMDF.metaData, 'AcquisitionContainer') && ...
        isfield(hMDF.metaData.AcquisitionContainer, 'ParsePlan')
    pp = hMDF.metaData.AcquisitionContainer.ParsePlan;
    dmdMeta.parsePlanZs = getFieldOrEmpty(pp, 'zs');
    dmdMeta.lineRateHz = getFieldOrEmpty(pp, 'lineRateHz');
    dmdMeta.linesPerCycle = getFieldOrEmpty(pp, 'linesPerCycle');
    dmdMeta.linesPerFrame = getFieldOrEmpty(pp, 'linesPerFrame');
    dmdMeta.pixPerLine = getFieldOrEmpty(pp, 'pixPerLine');
end
if isempty(dmdMeta.linePeriod_s) && isfield(dmdMeta, 'lineRateHz') && ~isempty(dmdMeta.lineRateHz)
    dmdMeta.linePeriod_s = 1 ./ double(dmdMeta.lineRateHz);
end
end

function val = getFieldOrEmpty(s, fieldName)
if isstruct(s) && isfield(s, fieldName)
    val = s.(fieldName);
else
    val = [];
end
end

% -------------------------------------------------------------------------
% Status/progress helpers
% -------------------------------------------------------------------------

function status = initializeStatus(roiTable)
status = repmat(struct( ...
    'globalRoiIdx', [], ...
    'dmdIdx', [], ...
    'localRoiIdx', [], ...
    'status', 'pending', ...
    'startedAt', '', ...
    'finishedAt', '', ...
    'nSamplesExtracted', 0, ...
    'weightClass', '', ...
    'weightSize', '', ...
    'errorMessage', ''), height(roiTable), 1);

for idx = 1:height(roiTable)
    status(idx).globalRoiIdx = roiTable.globalRoiIdx(idx);
    status(idx).dmdIdx = roiTable.dmdIdx(idx);
    status(idx).localRoiIdx = roiTable.localRoiIdx(idx);
end
end

function batches = makeBatches(items, batchSize)
if isempty(items)
    batches = {};
    return
end
batchSize = max(1, round(batchSize));
nBatches = ceil(numel(items) / batchSize);
batches = cell(1, nBatches);
for batchIdx = 1:nBatches
    startIdx = (batchIdx - 1) * batchSize + 1;
    stopIdx = min(numel(items), batchIdx * batchSize);
    batches{batchIdx} = items(startIdx:stopIdx);
end
end

function out = sizeToString(sz)
out = sprintf('%dx', sz);
out = out(1:end-1);
end

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------

function params = initializeParams(paramsIn)
params = defaultParams();

if nargin < 1
    paramsIn = [];
end

if ischar(paramsIn) || isstring(paramsIn)
    paramsIn = jsondecode(char(paramsIn));
end

% Pull defaults from the local pipeline GUI using the extractDendrites_new
% case. If paramsIn is empty, setParams opens the GUI. If paramsIn is
% provided, setParams merges those values with GUI defaults and returns
% without opening the GUI.
try
    if exist('setParams', 'file') == 2
        if isempty(paramsIn)
            guiParams = setParams('extractDendrites_new');
        else
            guiParams = setParams('extractDendrites_new', paramsIn, false);
        end
        params = mergeStructs(params, guiParams);
    end
catch ME
    warning('Ignoring setParams(''extractDendrites_new'') defaults: %s', ME.message);
end

if ~isempty(paramsIn)
    params = mergeStructs(params, paramsIn);
end

% Accept a common plural typo from early test versions.
if isfield(params, 'outputMode') && strcmpi(params.outputMode, 'trials')
    params.outputMode = 'trial';
end

params.outputMode = validatestring(params.outputMode, {'trial', 'continuous', 'both'});
params.storageMode = validatestring(params.storageMode, {'h5', 'memory'});
params.precision = validatestring(params.precision, {'single', 'double'});
if ~(ischar(params.outputSubfolderName) || isstring(params.outputSubfolderName))
    error('extractDendrites:InvalidParams', ...
        'params.outputSubfolderName must be a char vector or string scalar.');
end
params.outputSubfolderName = char(params.outputSubfolderName);
params.maxConcurrentROIs = max(1, round(params.maxConcurrentROIs));
params.numWorkers = max(1, round(params.numWorkers));
params.h5ChunkLines = max(1, round(params.h5ChunkLines));
params.h5Deflate = max(0, min(9, round(params.h5Deflate)));
end

function params = defaultParams()
params = struct();
params.manualROIs = false;
params.chIdx = 1;
params.zIdx = 1;
params.windowWidth_lines = 16;
params.expectedWindowWidth_lines = 5000;
params.outputMode = 'trial';
params.storageMode = 'h5';
params.precision = 'single';
params.outputDir = '';
params.useOutputSubfolder = true;
params.outputSubfolderName = 'dendriticVoltageExtraction';
params.timestampOutputSubfolder = false;
params.useParallel = true;
params.numWorkers = 4;
params.restartPool = true;
params.maxConcurrentROIs = 2;
params.makePlots = false;
params.saveRefImages = true;
params.saveSummaryAfterEachBatch = true;
params.resume = false;
params.stopOnError = true;
params.h5ChunkLines = 100000;
params.h5Deflate = 0;
end

function out = mergeStructs(base, overrides)
out = base;
if isempty(overrides)
    return
end
if ~isstruct(overrides)
    error('extractDendrites:InvalidParams', 'paramsIn must be a struct or JSON string.');
end
fields = fieldnames(overrides);
for idx = 1:numel(fields)
    out.(fields{idx}) = overrides.(fields{idx});
end
end

function params = configureParallelPool(params)
if ~params.useParallel
    return
end

try
    pp = gcp('nocreate');
    if isempty(pp)
        fprintf('Starting parallel pool with %d workers...\n', params.numWorkers);
        parpool(params.numWorkers);
    elseif pp.NumWorkers ~= params.numWorkers
        if params.restartPool
            fprintf('Restarting parallel pool: existing=%d requested=%d\n', ...
                pp.NumWorkers, params.numWorkers);
            delete(pp);
            parpool(params.numWorkers);
        else
            warning(['Existing parallel pool has %d workers but params.numWorkers=%d. ' ...
                'Using existing pool.'], pp.NumWorkers, params.numWorkers);
        end
    else
        fprintf('Using existing parallel pool with %d workers.\n', pp.NumWorkers);
    end
catch ME
    warning('Could not configure parallel pool. Falling back to synchronous extraction: %s', ME.message);
    params.useParallel = false; %#ok<NASGU>
end
end

function tf = wantsTrial(outputMode)
tf = any(strcmpi(outputMode, {'trial', 'both'}));
end

function tf = wantsContinuous(outputMode)
tf = any(strcmpi(outputMode, {'continuous', 'both'}));
end

function trace = castTrace(trace, precision)
if strcmpi(precision, 'single')
    trace = single(trace(:));
else
    trace = double(trace(:));
end
end